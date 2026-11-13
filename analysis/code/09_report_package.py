#!/usr/bin/env python3
"""Pack the reports for reading off the server: results/DOLORES_report.zip.

The zip holds report_story.html (the short story) and documentation.html (the
project README: full documentation), with figures embedded so each page also
works on its own, plus every file they name, in the same layout as results/,
so the file links in the pages open the files.

Run from analysis/ with pandoc on PATH (it comes with the r_RNASEQ env):
    conda activate r_RNASEQ && python3 code/09_report_package.py
"""
import base64
import datetime
import html
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from urllib.parse import quote

ANALYSIS = Path(__file__).resolve().parent.parent
RES = ANALYSIS / "results"
PKG = "DOLORES_report"
OUT = RES / f"{PKG}.zip"

# Package path -> source. Folders are added whole, minus SKIP.
FILES = {
    "figs/dry_cycle_1.png": ANALYSIS.parent / "figs/dry_cycle_1.png",
    "figs/dry_cycle_2.png": ANALYSIS.parent / "figs/dry_cycle_2.png",
    "data/samples.tsv": ANALYSIS / "data/samples.tsv",
    "data/excluded_samples.tsv": ANALYSIS / "data/excluded_samples.tsv",
    "data/kegg_pathway_names.tsv": ANALYSIS / "data/kegg_pathway_names.tsv",
}
DIRS = {
    "plots": RES / "plots",
    "qc": RES / "qc",
    "de/group": RES / "de/group",
    "memory_genes": RES / "memory_genes",
    "enrichment": RES / "enrichment",
    "report_tables": RES / "report_tables",
    "code": ANALYSIS / "code",
    "h2_genotype_cycle": RES / "h2_genotype_cycle",
    "bud_batch": RES / "bud_batch",
    "h5_rank": RES / "h5_rank",
    "alleles": RES / "alleles",
    "wgcna": RES / "wgcna",
    "clock": RES / "clock",
    "lncrna_pilot": RES / "lncrna_pilot",
    "lncrna_genome": RES / "lncrna_genome",
}
SKIP = re.compile(r"\.rds$")          # R objects, 230 MB in de/group alone
SERVER_ONLY = ("de/gene",)             # allele-level tables, 531 MB

# Paths in the markdown that point outside results/ -> package paths.
# The README lives at the repository root, so its links carry the analysis/results/ prefix.
REWRITE = [("](analysis/results/", "]("), ("../../README.md", "documentation.html"),
           ("../../figs/", "figs/"), ("../code/", "code/"), ("../data/", "data/"),
           ("report_story.md", "report_story.html")]

# unlink=True: links to files that are not in the package (the rest of the repository) become plain text.
PAGES = [
    dict(md=RES / "report_story.md", out="report_story.html", unlink=False,
         kicker="DOLORES · drought memory in sugarcane · the short story",
         switch='Full documentation, with every analysis and figure: '
                '<a href="documentation.html">documentation.html</a>'),
    dict(md=ANALYSIS.parent / "README.md", out="documentation.html", unlink=True,
         kicker="DOLORES · drought memory in sugarcane · full documentation",
         switch='Short version, insights only: <a href="report_story.html">report_story.html</a>'),
]


def collect():
    files = {k: v for k, v in FILES.items()}
    for dest, src in DIRS.items():
        for p in sorted(src.rglob("*")):
            if p.is_file() and not SKIP.search(p.name) and "__pycache__" not in p.parts:
                files[f"{dest}/{p.relative_to(src).as_posix()}"] = p
    missing = [k for k, v in files.items() if not v.exists()]
    if missing:
        sys.exit(f"missing sources: {missing}")
    return files


def pandoc(md):
    exe = shutil.which("pandoc")
    if not exe:
        sys.exit("pandoc not found: run inside the r_RNASEQ conda env")
    fmt = "gfm-emoji"   # the reports are written GitHub-style (lists without a blank line before)
    return subprocess.run([exe, "-f", fmt, "-t", "html5", "--wrap=none"],
                          input=md, capture_output=True, text=True, check=True).stdout


def strip_tags(s):
    return html.unescape(re.sub(r"<[^>]+>", "", s)).strip()


def postprocess(body, files):
    dirs = {str(Path(f).parent) + "/" for f in files}
    dirs |= {"/".join(f.split("/")[:i]) + "/" for f in files for i in range(1, f.count("/"))}
    by_name = {}
    for f in files:
        by_name.setdefault(Path(f).name, []).append(f)

    def href(p):
        return quote(p, safe="/._-")

    # Figures: an image paragraph followed by its "Figure N." caption.
    def fig(m):
        src, alt, num, cap = m.group(1), m.group(2), m.group(3), m.group(4).strip()
        return (f'<figure id="fig-{num}">\n<img src="{src}" alt="{alt}" tabindex="0">\n'
                f'<figcaption><span class="label">Figure {num}.</span> {cap} '
                f'<a class="path" href="{href(html.unescape(src))}">{src}</a></figcaption>\n</figure>')
    body = re.sub(r'<p><img src="([^"]+)" alt="([^"]*)" /></p>\s*<p><em>Figure (\d+)\.(.*?)</em></p>',
                  fig, body, flags=re.S)
    body = re.sub(r'<p><img src="([^"]+)" alt="([^"]*)" /></p>',
                  r'<figure>\n<img src="\1" alt="\2" tabindex="0">\n</figure>', body)

    # Tables: let the browser size columns, scroll on narrow screens, style captions.
    body = re.sub(r"<colgroup>.*?</colgroup>\s*", "", body, flags=re.S)
    body = re.sub(r"<table>(.*?)</table>", r'<div class="table-wrap"><table>\1</table></div>',
                  body, flags=re.S)

    def tcap(m):
        cap = m.group(1).strip()
        t = re.match(r"Table (\d+)\.(.*)", cap, flags=re.S)
        if t:
            return (f'</table></div>\n<p class="caption" id="tab-{t.group(1)}">'
                    f'<span class="label">Table {t.group(1)}.</span>{t.group(2)}</p>')
        return f'</table></div>\n<p class="caption">{cap}</p>'
    body = re.sub(r"</table></div>\s*<p><em>(.*?)</em></p>", tcap, body, flags=re.S)

    # Concept boxes.
    def concept(m):
        rest = m.group(2).strip()
        out = (f'<blockquote class="concept">\n<p class="concept-title"><span class="tag">Concept</span>'
               f'{m.group(1).strip()}</p>')
        return out + (f"\n<p>{rest}" if rest else "\n<p>")
    body = re.sub(r"<blockquote>\s*<p><strong>Concept:(.*?)</strong>(.*?)(?=</p>)", concept, body,
                  flags=re.S)
    body = re.sub(r"<p>\s*</p>\n?", "", body)

    # Verdicts.
    body = re.sub(r"<p><strong>(H\d verdict)[.:]?</strong>(.*?)</p>(\s*<ul>.*?</ul>)?",
                  lambda m: (f'<div class="verdict"><p><span class="label">{m.group(1)}</span>'
                             f'{m.group(2)}</p>{m.group(3) or ""}</div>'), body, flags=re.S)

    # Lede: the first paragraph before the first section.
    first_h2 = body.find("<h2")
    m = re.search(r"<p>(.*?)</p>", body, flags=re.S)
    if m and (first_h2 < 0 or m.start() < first_h2):
        inner = m.group(1)
        e = re.fullmatch(r"<em>(.*)</em>", inner, flags=re.S)
        body = body[:m.start()] + f'<p class="lede">{e.group(1) if e else inner}</p>' + body[m.end():]

    # Section numbers in headings, and the number -> id map for cross-links.
    sections, toc = {}, []

    def head(m):
        lvl, hid, inner = m.group(1), m.group(2), m.group(3)
        n = re.match(r"(\d+(?:\.\d+)?)\.?\s+", inner)
        if n:
            sections[n.group(1)] = hid
            inner = f'<span class="num">{n.group(1)}</span>' + inner[n.end():]
        toc.append((int(lvl), hid, inner))
        return f'<h{lvl} id="{hid}">{inner}</h{lvl}>'
    body = re.sub(r'<h([23]) id="([^"]+)">(.*?)</h\1>', head, body, flags=re.S)
    fig_ids = set(re.findall(r'id="fig-(\d+)"', body))
    tab_ids = set(re.findall(r'id="tab-(\d+)"', body))

    def link_text(t):
        t = re.sub(r"\bFigure (\d+)\b",
                   lambda m: f'<a href="#fig-{m.group(1)}">{m.group(0)}</a>' if m.group(1) in fig_ids
                   else m.group(0), t)
        t = re.sub(r"\bTable (\d+)\b",
                   lambda m: f'<a href="#tab-{m.group(1)}">{m.group(0)}</a>' if m.group(1) in tab_ids
                   else m.group(0), t)

        def secs(m):
            if m.group(2) in sections:   # one number: link the whole "section N"
                return f'<a href="#{sections[m.group(2)]}">{m.group(0)}</a>'
            nums = re.sub(r"\d+(?:\.\d+)?",
                          lambda n: f'<a href="#{sections[n.group(0)]}">{n.group(0)}</a>'
                          if n.group(0) in sections else n.group(0), m.group(2))
            return m.group(1) + nums
        return re.sub(r"\b(sections? )(\d+(?:\.\d+)?(?:(?:, | and )\d+(?:\.\d+)?)*)", secs, t,
                      flags=re.I)

    # One pass over tags and text: cross-links in plain text, file links on code paths.
    toks = re.split(r"(<[^>]+>)", body)
    out, i = [], 0
    in_a = in_code = in_head = in_label = 0
    while i < len(toks):
        tok = toks[i]
        if tok.startswith("<"):
            name = re.match(r"</?([a-zA-Z0-9]+)", tok)
            name = name.group(1).lower() if name else ""
            close = tok.startswith("</")
            if name == "code" and not close and not in_a and i + 2 < len(toks) and toks[i + 2] == "</code>":
                p = html.unescape(toks[i + 1])
                target = None
                if p in files or p in dirs:
                    target = p
                elif "/" not in p and len(by_name.get(p, [])) == 1:
                    target = by_name[p][0]
                elif re.search(r"[<*…]", p) and "/" in p and p.rsplit("/", 1)[0] + "/" in dirs:
                    target = p.rsplit("/", 1)[0] + "/"   # a name pattern: link its folder
                if target:
                    out.append(f'<a class="path" href="{href(target)}"><code>{toks[i + 1]}</code></a>')
                    i += 3
                    continue
                if p.startswith(SERVER_ONLY):
                    out.append(f'<code class="server" title="Not in this package: on the server only">'
                               f'{toks[i + 1]}</code>')
                    i += 3
                    continue
            d = -1 if close else 1
            if name == "a":
                in_a += d
            elif name == "code":
                in_code += d
            elif name in ("h1", "h2", "h3", "h4"):
                in_head += d
            elif name == "span" and 'class="label"' in tok:
                in_label += 1
            elif name == "span" and close and in_label:
                in_label -= 1
            out.append(tok)
        else:
            out.append(tok if (in_a or in_code or in_head or in_label) else link_text(tok))
        i += 1
    return "".join(out), toc


def embed_images(body, files):
    def data_uri(m):
        src = html.unescape(m.group(1))
        if src not in files:
            sys.exit(f"image not in package: {src}")
        b64 = base64.b64encode(files[src].read_bytes()).decode()
        return f'<img src="data:image/png;base64,{b64}"'
    return re.sub(r'<img src="([^"]+)"', data_uri, body)


def toc_html(toc):
    items, open_sub = [], False
    for lvl, hid, inner in toc:
        if lvl == 2:
            if open_sub:
                items.append("</ol></li>")
                open_sub = False
            elif items:
                items.append("</li>")
            items.append(f'<li><a href="#{hid}">{inner}</a>')
        else:
            if not open_sub:
                items.append("<ol>")
                open_sub = True
            items.append(f'<li><a href="#{hid}">{inner}</a></li>')
    items.append("</ol></li>" if open_sub else "</li>")
    return "<ol>" + "".join(items) + "</ol>"


CSS = r"""
:root {
  --bg: #fbfaf7; --surface: #ffffff; --text: #1d2125; --muted: #5a626d; --line: #e2dfd7;
  --accent: #1d5c87; --accent-soft: #eaf1f6; --concept-bg: #f1f6f2; --concept-line: #2f6b4c;
  --tag-ink: #ffffff; --code-bg: #efece5; --verdict-line: #1d5c87;
  --serif: "Iowan Old Style", "Palatino Linotype", Palatino, "Book Antiqua", Georgia, serif;
  --sans: system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  --mono: ui-monospace, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #15171b; --surface: #1c1f24; --text: #e4e4e0; --muted: #a2a8b1; --line: #30343b;
    --accent: #8dbfe6; --accent-soft: #1d2731; --concept-bg: #18231d; --concept-line: #74b893;
    --tag-ink: #10140f; --code-bg: #262a30; --verdict-line: #8dbfe6;
  }
}
:root[data-theme="dark"] {
  --bg: #15171b; --surface: #1c1f24; --text: #e4e4e0; --muted: #a2a8b1; --line: #30343b;
  --accent: #8dbfe6; --accent-soft: #1d2731; --concept-bg: #18231d; --concept-line: #74b893;
  --tag-ink: #10140f; --code-bg: #262a30; --verdict-line: #8dbfe6;
}
* { box-sizing: border-box; }
html { scroll-behavior: smooth; }
body { margin: 0; background: var(--bg); color: var(--text); font-family: var(--serif);
  font-size: 1.125rem; line-height: 1.65; -webkit-text-size-adjust: 100%; }
.layout { max-width: 78rem; margin: 0 auto; padding: 0 16px; }
@media (min-width: 1100px) {
  .layout { display: grid; grid-template-columns: 16rem minmax(0, 46rem); gap: 4rem; justify-content: center; }
  .toc { position: sticky; top: 0; align-self: start; max-height: 100vh; overflow-y: auto; padding: 2.5rem 0 2rem; }
}
main { min-width: 0; padding: 2.5rem 0 6rem; max-width: 46rem; margin: 0 auto; }

/* contents */
.toc { font-family: var(--sans); font-size: .86rem; line-height: 1.4; }
.toc summary { cursor: pointer; font-weight: 650; letter-spacing: .06em; text-transform: uppercase;
  font-size: .75rem; color: var(--muted); padding: .4rem 0; }
.toc ol { list-style: none; margin: .3rem 0 0; padding: 0; }
.toc ol ol { margin: .15rem 0 .35rem .9rem; }
.toc a { display: block; padding: .22rem .6rem; color: var(--muted); text-decoration: none;
  border-left: 2px solid transparent; }
.toc a:hover { color: var(--text); }
.toc a.active { color: var(--accent); border-left-color: var(--accent); }
.toc .num { display: inline-block; min-width: 1.6em; font-variant-numeric: tabular-nums; }
@media (max-width: 1099px) {
  .toc { border-bottom: 1px solid var(--line); padding: 1rem 0 .5rem; }
}

/* header */
.kicker { font-family: var(--sans); font-size: .78rem; letter-spacing: .07em; text-transform: uppercase;
  color: var(--muted); margin: 0 0 .6rem; }
h1 { font-family: var(--sans); font-size: clamp(1.85rem, 4.5vw, 2.6rem); line-height: 1.15;
  letter-spacing: -.015em; margin: 0 0 1rem; text-wrap: balance; }
.switch { font-family: var(--sans); font-size: .88rem; color: var(--muted); margin: 0 0 1.5rem; }
.lede { font-size: 1.2rem; color: var(--muted); font-style: italic; margin: 0 0 1rem; }

/* text */
h2, h3, h4 { font-family: var(--sans); line-height: 1.25; scroll-margin-top: 1rem; text-wrap: balance; }
h2 { font-size: 1.55rem; margin: 3.2rem 0 1rem; letter-spacing: -.01em; }
h3 { font-size: 1.2rem; margin: 2.4rem 0 .8rem; }
.num { color: var(--accent); margin-right: .45em; font-variant-numeric: tabular-nums; }
p, ul, ol { margin: 0 0 1.05rem; }
li { margin: .25rem 0; }
li > ul, li > ol { margin: .3rem 0 .4rem; }
li > p { margin: .35rem 0; }
a { color: var(--accent); text-underline-offset: .16em; text-decoration-thickness: 1px; }
strong { font-weight: 650; }
hr { border: 0; border-top: 1px solid var(--line); margin: 3rem 0 0; }
code { font-family: var(--mono); font-size: .84em; background: var(--code-bg); padding: .08em .32em;
  border-radius: 4px; overflow-wrap: anywhere; }
a.path { text-decoration: none; }
a.path code { color: var(--accent); }
a.path:hover code { text-decoration: underline; }
code.server { opacity: .75; border: 1px dashed var(--line); }

/* concept boxes */
blockquote { margin: 1.8rem 0; padding: 1rem 1.25rem; border-left: 4px solid var(--line); }
blockquote.concept { background: var(--concept-bg); border-left-color: var(--concept-line);
  border-radius: 0 10px 10px 0; font-size: 1.02rem; }
blockquote > :last-child { margin-bottom: 0; }
.concept-title { font-family: var(--sans); font-weight: 650; font-size: 1rem; margin: 0 0 .55rem; }
.tag { display: inline-block; font-size: .68rem; font-weight: 700; letter-spacing: .08em;
  text-transform: uppercase; background: var(--concept-line); color: var(--tag-ink);
  padding: .15rem .45rem; border-radius: 4px; margin-right: .6rem; vertical-align: .12em; }

/* verdicts */
.verdict { background: var(--accent-soft); border-left: 4px solid var(--verdict-line);
  border-radius: 0 10px 10px 0; padding: .9rem 1.2rem; margin: 1.8rem 0; }
.verdict > :last-child { margin-bottom: 0; }
.verdict .label { font-family: var(--sans); font-weight: 700; margin-right: .5rem; }

/* figures */
figure { margin: 2.2rem 0; }
figure img { display: block; width: 100%; height: auto; background: #fff; border: 1px solid var(--line);
  border-radius: 8px; padding: 6px; cursor: zoom-in; }
figure img:focus-visible { outline: 3px solid var(--accent); outline-offset: 2px; }
figcaption, .caption { font-family: var(--sans); font-size: .88rem; line-height: 1.5; color: var(--muted);
  margin-top: .7rem; }
figcaption .label, .caption .label { color: var(--text); font-weight: 650; margin-right: .2rem; }
figcaption a.path { display: inline-block; font-family: var(--mono); font-size: .78rem; margin-left: .3rem; }
.caption { margin: .2rem 0 1.6rem; }

/* tables */
.table-wrap { overflow-x: auto; margin: 1.4rem 0 .4rem; -webkit-overflow-scrolling: touch; }
table { border-collapse: collapse; width: 100%; font-family: var(--sans); font-size: .9rem; line-height: 1.45;
  font-variant-numeric: tabular-nums; }
th, td { padding: .5rem .7rem; border-bottom: 1px solid var(--line); text-align: left; vertical-align: top; }
thead th { border-bottom: 2px solid var(--muted); font-weight: 650; white-space: nowrap; }
tbody tr:hover { background: var(--accent-soft); }
td code { font-size: .82em; }

/* lightbox */
.lightbox { position: fixed; inset: 0; z-index: 10; background: rgba(10, 12, 15, .88); display: flex;
  align-items: center; justify-content: center; padding: 2vh 2vw; cursor: zoom-out; }
.lightbox[hidden] { display: none; }
.lightbox img { max-width: 100%; max-height: 100%; background: #fff; border-radius: 6px; }
.lightbox p { position: fixed; bottom: .6rem; left: 0; right: 0; text-align: center; margin: 0;
  font-family: var(--sans); font-size: .8rem; color: #d6d6d6; }

@media print {
  .toc, .lightbox { display: none !important; }
  .layout { display: block; max-width: none; padding: 0; }
  main { max-width: none; padding: 0; }
  body { background: #fff; color: #000; font-size: 10.5pt; }
  a { color: inherit; }
  figure, table, blockquote, .verdict { break-inside: avoid; }
  figure img { border: 0; padding: 0; }
}
"""

JS = r"""
(function () {
  var details = document.querySelector('.toc details');
  if (details && window.matchMedia('(max-width: 1099px)').matches) details.open = false;
  var links = {};
  document.querySelectorAll('.toc a').forEach(function (a) { links[a.getAttribute('href').slice(1)] = a; });
  var heads = Array.prototype.filter.call(document.querySelectorAll('main h2[id], main h3[id]'),
    function (h) { return links[h.id]; });
  var current = null, ticking = false;
  function update() {
    ticking = false;
    var pick = heads[0];
    for (var i = 0; i < heads.length; i++) {
      if (heads[i].getBoundingClientRect().top < 140) pick = heads[i]; else break;
    }
    if (pick && links[pick.id] !== current) {
      if (current) current.classList.remove('active');
      current = links[pick.id];
      current.classList.add('active');
    }
  }
  window.addEventListener('scroll', function () {
    if (!ticking) { ticking = true; window.requestAnimationFrame(update); }
  }, { passive: true });
  update();

  var lb = document.getElementById('lightbox'), lbImg = lb.querySelector('img');
  function openImg(img) { lbImg.src = img.src; lbImg.alt = img.alt; lb.hidden = false; }
  function close() { lb.hidden = true; lbImg.removeAttribute('src'); }
  document.querySelectorAll('figure img').forEach(function (img) {
    img.addEventListener('click', function () { openImg(img); });
    img.addEventListener('keydown', function (e) {
      if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); openImg(img); }
    });
  });
  lb.addEventListener('click', close);
  document.addEventListener('keydown', function (e) { if (e.key === 'Escape' && !lb.hidden) close(); });
})();
"""

PAGE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>{css}</style>
</head>
<body>
<div class="layout">
<nav class="toc" aria-label="Contents"><details open><summary>Contents</summary>
{toc}
</details></nav>
<main>
<header>
<p class="kicker">{kicker}</p>
<h1>{h1}</h1>
<p class="switch">{switch}</p>
</header>
{body}
</main>
</div>
<div id="lightbox" class="lightbox" hidden><img alt=""><p>Click or press Esc to close. The full-size PNG is linked under each figure.</p></div>
<script>{js}</script>
</body>
</html>
"""

README = """DOLORES drought-memory report package (built {date})

Open report_story.html in a web browser: it is the short story, the main
insights only. documentation.html is the full documentation (the project
README), with every analysis, figure and method. Figures are embedded, so
both pages also work if you move them out of this folder.

The folders keep the layout of analysis/results/ on the server, so every file
the reports name is here, and the file links in the pages open it:
  plots/  qc/  de/group/  memory_genes/  enrichment/  report_tables/   results
  h2_genotype_cycle/  bud_batch/  h5_rank/  alleles/  wgcna/      results
  lncrna_pilot/    lncRNA pilot on chromosome 1
  lncrna_genome/   lncRNA search on the whole genome
  figs/            experiment schematics (Figures 1 and 2)
  data/            sample sheet and excluded libraries
  code/            the analysis scripts that made all of this

.tsv files are tab-separated tables (open with a spreadsheet program or R).
.tsv.gz files are the same, gzip-compressed: R reads them directly
(readr::read_tsv or data.table::fread); otherwise decompress them first.

Left out because of size; they stay on the server under
/Storage/data1/jorge.munoz/DOLORES/analysis/:
  results/de/gene/        allele-level DE tables (~530 MB)
  results/de/group/*.rds  DESeq2 objects (~230 MB)
  data/*.rds              imported counts and VST objects
"""


def main():
    files = collect()
    pages = {}
    for cfg in PAGES:
        md = cfg["md"].read_text(encoding="utf-8")
        first, _, rest = md.partition("\n")
        if not first.startswith("# "):
            sys.exit(f"{cfg['md'].name}: first line must be the '# ' title")
        h1 = first[2:].strip()
        for a, b in REWRITE:
            rest = rest.replace(a, b)
        body, toc = postprocess(pandoc(rest), files)
        def in_package(target):
            t = target.rstrip("/")
            return t in files or any(f.startswith(t + "/") for f in files) or \
                t in [p["out"] for p in PAGES]
        if cfg["unlink"]:
            body = re.sub(r'<a href="([^"#:][^"]*)">(.*?)</a>',
                          lambda m: m.group(0) if in_package(html.unescape(m.group(1))) else m.group(2),
                          body, flags=re.S)
        leftover = re.findall(r'(?:src|href)="\.\./[^"]*"', body)
        if leftover:
            sys.exit(f"{cfg['md'].name}: paths outside the package: {leftover}")
        for target in re.findall(r'href="([^"#:][^"]*)"', body):
            if not in_package(target):
                sys.exit(f"{cfg['md'].name}: link to a file not in the package: {target}")
        body = embed_images(body, files)
        h1_html = html.escape(h1)
        pages[cfg["out"]] = PAGE.format(title=h1_html, css=CSS, js=JS, toc=toc_html(toc),
                                        kicker=cfg["kicker"], h1=h1_html, switch=cfg["switch"],
                                        body=body)
        print(f"{cfg['out']}: {body.count('<figure')} figures, {body.count('<table')} tables, "
              f"{len(toc)} headings, {len(pages[cfg['out']]) / 1e6:.1f} MB")

    stored = (".png", ".gz", ".zip", ".pdf")
    with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for name, text in pages.items():
            z.writestr(f"{PKG}/{name}", text)
        z.writestr(f"{PKG}/README.txt", README.format(date=datetime.date.today().isoformat()))
        for arc, src in sorted(files.items()):
            z.write(src, f"{PKG}/{arc}",
                    compress_type=zipfile.ZIP_STORED if arc.endswith(stored) else zipfile.ZIP_DEFLATED)
    print(f"{OUT}: {len(files) + len(pages) + 1} files, {OUT.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
