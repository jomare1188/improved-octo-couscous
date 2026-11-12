# Colour-vision check for plot palettes (R port of the dataviz validator's core checks):
# adjacent-pair OKLab Delta E x100 under normal vision and simulated deutan/protan/tritan (colorspace, severity 1),
# and WCAG contrast against the white plot surface. Targets: CVD dE >= 8 (6-8 only with a second encoding),
# normal dE >= 15, contrast >= 3 for marks.
suppressMessages({ library(colorspace); library(farver) })
oklab <- function(hex) farver::convert_colour(t(col2rgb(hex)), "rgb", "oklab")
dE <- function(a, b) 100 * sqrt(sum((oklab(a) - oklab(b))^2))
lum <- function(hex) { x <- col2rgb(hex)[, 1] / 255; x <- ifelse(x <= 0.03928, x / 12.92, ((x + 0.055) / 1.055)^2.4); sum(c(0.2126, 0.7152, 0.0722) * x) }
contrast <- function(hex, bg = "#FFFFFF") { l <- sort(c(lum(hex), lum(bg)), TRUE); (l[1] + 0.05) / (l[2] + 0.05) }
check_palette <- function(pal, name = "", pairs = "adjacent") {
  sims <- list(normal = pal, deutan = deutan(pal), protan = protan(pal), tritan = tritan(pal))
  idx <- if (pairs == "all") t(combn(length(pal), 2)) else cbind(seq_len(length(pal) - 1), 2:length(pal))
  res <- sapply(sims, function(p) apply(idx, 1, function(i) dE(p[i[1]], p[i[2]])))
  res <- matrix(res, nrow = nrow(idx)); colnames(res) <- names(sims)
  worst_cvd <- min(res[, 2:4]); worst_norm <- min(res[, 1]); con <- sapply(pal, contrast)
  verdict <- if (worst_norm < 15 || worst_cvd < 6) "FAIL" else if (worst_cvd < 8 || min(con) < 3) "WARN" else "PASS"
  cat(sprintf("%-28s %s  min normal dE %.1f | min CVD dE %.1f | min contrast %.2f\n", name, verdict, worst_norm, worst_cvd, min(con)))
  invisible(verdict)
}
