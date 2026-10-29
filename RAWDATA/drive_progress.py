#!/usr/bin/env python3
"""Show the progress of download_drive_dolores.sh.

Reads only local files (no Google Drive calls):
  logs_drive/drive_manifest.tsv   expected files and sizes
  drive_provenance.tsv            files already finished (OK / failed)
  *.partial in RAWDATA/           file being downloaded right now
  logs_drive/drive_dolores.*.log  start time of the job

Usage:
  python3 drive_progress.py            # show once
  python3 drive_progress.py --watch 30 # refresh every 30 s (Ctrl+C to stop)
  python3 drive_progress.py --all      # also list every file
"""
import csv, glob, os, re, subprocess, sys, time
from datetime import datetime

RAW = '/Storage/data1/jorge.munoz/DOLORES/RAWDATA'
MANIFEST = f'{RAW}/logs_drive/drive_manifest.tsv'
PROVENANCE = f'{RAW}/drive_provenance.tsv'

GREEN, YELLOW, RED, DIM, BOLD, RESET = '\033[32m', '\033[33m', '\033[31m', '\033[2m', '\033[1m', '\033[0m'
if not sys.stdout.isatty():
    GREEN = YELLOW = RED = DIM = BOLD = RESET = ''


def gb(n): return f'{n / 1e9:.2f} GB'


def bar(frac, width=40):
    full = int(frac * width)
    return '█' * full + '░' * (width - full)


def hms(seconds):
    seconds = int(seconds)
    return f'{seconds // 3600}h{seconds % 3600 // 60:02d}m' if seconds >= 3600 else f'{seconds // 60}m{seconds % 60:02d}s'


def job_state():
    try:
        out = subprocess.run(['qstat', '-u', os.environ.get('USER', '')], capture_output=True, text=True, timeout=10).stdout
    except Exception:
        return None
    jobs = [l.split() for l in out.splitlines()[2:] if l.strip()]
    logs = sorted(glob.glob(f'{RAW}/logs_drive/drive_dolores.*.log'), key=os.path.getmtime)
    if logs:
        job_id = re.search(r'\.(\d+)\.log$', logs[-1]).group(1)
        for j in jobs:
            if j[0] == job_id:
                return f'job {job_id}: ' + {'r': 'running', 'qw': 'waiting in queue'}.get(j[4], j[4])
        return f'job {job_id}: not in queue (finished or killed)'
    return 'no job log found'


def start_time():
    logs = sorted(glob.glob(f'{RAW}/logs_drive/drive_dolores.*.log'), key=os.path.getmtime)
    if not logs:
        return None
    for line in open(logs[-1]):
        m = re.search(r'^DOWNLOAD .* (\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)$', line)
        if m:
            return datetime.strptime(m.group(1), '%Y-%m-%d %H:%M:%S')
    return None


def show(list_all):
    manifest = list(csv.DictReader(open(MANIFEST), delimiter='\t'))
    finished = {}
    if os.path.exists(PROVENANCE):
        for r in csv.DictReader(open(PROVENANCE), delimiter='\t'):
            finished[r['file']] = r['status']
    partials = {}
    for p in glob.glob(f'{RAW}/*.partial'):
        name = re.sub(r'\.[0-9a-f]+\.partial$', '', os.path.basename(p))
        partials[name] = os.path.getsize(p)

    total = sum(int(r['size_bytes']) for r in manifest)
    done_bytes = 0
    rows, n_ok, n_fail, n_skip = [], 0, 0, 0
    current = None
    for r in manifest:
        f, size = r['file'], int(r['size_bytes'])
        st = finished.get(f)
        if st == 'OK':
            n_ok += 1; done_bytes += size; label = f'{GREEN}✔ done{RESET}'
        elif st:
            n_fail += 1; label = f'{RED}✘ {st}{RESET}'
        elif f in partials:
            done_bytes += partials[f]; current = (f, partials[f], size)
            label = f'{YELLOW}⬇ {100 * partials[f] / size:5.1f}%{RESET}'
        elif os.path.exists(f'{RAW}/{f}'):
            n_skip += 1; label = f'{DIM}• already in RAWDATA{RESET}'
        else:
            label = f'{DIM}· pending{RESET}'
        rows.append((f, r['remote'], size, label))

    frac = done_bytes / total if total else 0
    print(f'{BOLD}Dolores Google Drive download{RESET}   {datetime.now():%Y-%m-%d %H:%M:%S}')
    print(f'{job_state()}\n')
    print(f'  {bar(frac)} {100 * frac:5.1f}%   {gb(done_bytes)} / {gb(total)}')
    print(f'  files: {GREEN}{n_ok} done{RESET}, {RED}{n_fail} failed{RESET}, {n_skip} skipped, '
          f'{len(manifest) - n_ok - n_fail - n_skip} to go (of {len(manifest)})')

    t0 = start_time()
    if t0 and done_bytes:
        elapsed = (datetime.now() - t0).total_seconds()
        speed = done_bytes / elapsed
        eta = (total - done_bytes) / speed if speed else 0
        print(f'  elapsed {hms(elapsed)}   average {speed / 1e6:.1f} MB/s   ETA {hms(eta) if n_ok + n_fail + n_skip < len(manifest) else "-"}')
    if current:
        f, got, size = current
        print(f'\n  now: {f}\n  {bar(got / size, 30)} {100 * got / size:5.1f}%   {gb(got)} / {gb(size)}')

    print()
    by_remote = {}
    for f, remote, size, label in rows:
        by_remote.setdefault(remote, []).append((f, size, label))
    for remote, items in by_remote.items():
        ok = sum('done' in l for _, _, l in items)
        print(f'  {remote}: {ok}/{len(items)} files')
        if list_all or any('✘' in l or '⬇' in l for _, _, l in items):
            for f, size, label in items:
                if list_all or '✘' in label or '⬇' in label:
                    print(f'      {f:<28} {gb(size):>9}  {label}')
    if n_fail:
        print(f'\n  {RED}Some files failed: see drive_provenance.tsv and logs_drive/. Resubmitting retries only the missing ones.{RESET}')


def main():
    args = sys.argv[1:]
    list_all = '--all' in args
    if '--watch' in args:
        i = args.index('--watch')
        every = int(args[i + 1]) if i + 1 < len(args) and args[i + 1].isdigit() else 30
        try:
            while True:
                os.system('clear')
                show(list_all)
                print(f'\n{DIM}refreshing every {every}s, Ctrl+C to stop{RESET}')
                time.sleep(every)
        except KeyboardInterrupt:
            print()
    else:
        show(list_all)


if __name__ == '__main__':
    main()
