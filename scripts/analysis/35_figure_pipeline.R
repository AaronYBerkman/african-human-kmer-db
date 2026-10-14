# Figure 1: construction of the African human k-mer database (scripts/db_build/01-11).
# Three input columns (microbial, African human, viral) are indexed and merged, then
# microbial and viral k-mers are subtracted from the African human set in turn.
# Drawn with grid in millimetres at Nature double-column width (183 mm), so every
# font size below is the true printed size (Helvetica, 7-8 pt).

library(grid)

repo_root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("--file=", "", args[grep("^--file=", args)])
  dirname(dirname(dirname(normalizePath(file_arg))))
})
outdir <- file.path(repo_root, "figures")

W <- 183; H <- 142                        # mm
pt_head <- 9; pt_title <- 8.5; pt_body <- 8
ink <- "grey10"; body_col <- "grey20"; detail <- "grey40"
arrow_col <- "grey45"

# Okabe-Ito orange and reddish purple for microbial and viral; the African column
# uses the "African DB only" cyan from Fig. 2 (#4DBBD5) so the figures link up.
hue <- list(
  microbial = c(line = "#B57B00", fill = "#FDF3E1", head = "#9A6A00"),
  african   = c(line = "#2A8CA3", fill = "#E6F5F9", head = "#1F7487"),
  viral     = c(line = "#A3527F", fill = "#F8EAF2", head = "#8C4570"),
  subtract  = c(line = "grey55",  fill = "grey96",  head = "grey30"),
  final     = c(line = "#1F7487", fill = "#CDEBF2", head = "#1F7487")
)

# Layout, y measured from the top in mm
xcol <- c(microbial = 30.5, african = W / 2, viral = W - 30.5); wcol <- 56
y_src <- 23.5; h_src <- 17
y_idx <- 43.5; h_idx <- 12
y_mrg <- 61; h_mrg <- 12
y_s1 <- 85.5; y_s2 <- 103; h_sub <- 12; w_sub <- 74
y_out <- 129.5; h_out <- 17; w_out <- 84

title <- function(s) list(s, "bold",  pt_title, ink)
body  <- function(s) list(s, "plain", pt_body,  body_col)
det   <- function(s) list(s, "plain", pt_body,  detail)

boxes <- list(
  list("microbial", xcol[["microbial"]], y_src, wcol, h_src, list(
    title("Microbial genomes"), body("60,294 bacterial + 707 archaeal"),
    det("NCBI RefSeq complete genomes"))),
  list("african", xcol[["african"]], y_src, wcol, h_src, list(
    title("African human assemblies"), body("3 diploid genomes, 6 haplotypes"),
    det("South Africa · Kenya · Burkina Faso"))),
  list("viral", xcol[["viral"]], y_src, wcol, h_src, list(
    title("Viral genomes"), body("14,626 genomes"),
    det("NCBI RefSeq complete genomes"))),
  list("microbial", xcol[["microbial"]], y_idx, wcol, h_idx, list(
    title("Dense 32-mer index"), body("Built per genome"))),
  list("african", xcol[["african"]], y_idx, wcol, h_idx, list(
    title("Dense 32-mer index"), body("4.5 billion k-mers"))),
  list("viral", xcol[["viral"]], y_idx, wcol, h_idx, list(
    title("Dense 32-mer index"), body("127.5 million k-mers"))),
  list("microbial", xcol[["microbial"]], y_mrg, wcol, h_mrg, list(
    title("Merge and deduplicate"), body("~16 billion unique k-mers"))),
  list("african", xcol[["african"]], y_mrg, wcol, h_mrg, list(
    title("Merge and deduplicate"), body("~1 billion unique k-mers"))),
  list("viral", xcol[["viral"]], y_mrg, wcol, h_mrg, list(
    title("Merge and deduplicate"), body("~94 million unique k-mers"))),
  list("subtract", W / 2, y_s1, w_sub, h_sub, list(
    title("Subtract microbial k-mers"),
    body("Remove k-mers shared with bacteria and archaea"))),
  list("subtract", W / 2, y_s2, w_sub, h_sub, list(
    title("Subtract viral k-mers"),
    body("Remove k-mers shared with viral genomes"))),
  list("final", W / 2, y_out, w_out, h_out, list(
    list("African human k-mer database", "bold", pt_head, ink),
    body("1.01 billion 32-mers · 8.1 GB"),
    det("Human k-mers absent from microbial and viral genomes")))
)

mm <- function(x) unit(x, "mm")
Y  <- function(y) mm(H - y)                    # top-down y to grid coordinates
arr <- arrow(length = mm(1.7), angle = 25, type = "closed")

seg <- function(x0, y0, x1, y1, col = arrow_col, lwd = 0.7, head = TRUE) {
  grid.segments(mm(x0), Y(y0), mm(x1), Y(y1), gp = gpar(col = col, fill = col, lwd = lwd),
                arrow = if (head) arr else NULL)
}

draw_box <- function(b) {
  h <- hue[[b[[1]]]]
  lwd <- if (b[[1]] == "final") 1.1 else 0.6
  grid.roundrect(mm(b[[2]]), Y(b[[3]]), mm(b[[4]]), mm(b[[5]]), r = mm(1.5),
                 gp = gpar(fill = h[["fill"]], col = h[["line"]], lwd = lwd))
  lines <- b[[6]]; n <- length(lines); gap <- 3.9
  ys <- b[[3]] - gap * ((n - 1) / 2 - (seq_len(n) - 1))
  for (i in seq_len(n)) {
    l <- lines[[i]]
    grid.text(l[[1]], mm(b[[2]]), Y(ys[i]),
              gp = gpar(fontface = l[[2]], fontsize = l[[3]], col = l[[4]]))
  }
}

# Section rule, broken where an arrow crosses it (gaps = x positions of crossings)
section <- function(label, y, gaps = numeric()) {
  grid.text(label, mm(4), Y(y), just = "left",
            gp = gpar(fontsize = pt_body, col = detail))
  w <- convertWidth(grobWidth(textGrob(label, gp = gpar(fontsize = pt_body))), "mm", TRUE)
  gaps <- sort(gaps)
  x0 <- c(4 + w + 2, gaps + 1.8); x1 <- c(gaps - 1.8, W - 3)
  keep <- x1 - x0 >= 3; x0 <- x0[keep]; x1 <- x1[keep]      # no stray stubs
  grid.segments(mm(x0), Y(y), mm(x1), Y(y), gp = gpar(col = "grey85", lwd = 0.5))
}

draw <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  section("Database construction", 4)
  section("k-mer subtraction", 73, c(xcol[["microbial"]], W / 2, xcol[["viral"]]))
  section("Output", 115, W / 2)

  heads <- c(microbial = "Microbial", african = "African human", viral = "Viral")
  for (k in names(heads))
    grid.text(heads[[k]], mm(xcol[[k]]), Y(11.5),
              gp = gpar(fontface = "bold", fontsize = pt_head, col = hue[[k]][["head"]]))

  # Column arrows
  for (k in names(xcol)) {
    seg(xcol[[k]], y_src + h_src / 2, xcol[[k]], y_idx - h_idx / 2 - 0.3)
    seg(xcol[[k]], y_idx + h_idx / 2, xcol[[k]], y_mrg - h_mrg / 2 - 0.3)
  }
  # Subtraction chain
  seg(W / 2, y_mrg + h_mrg / 2, W / 2, y_s1 - h_sub / 2 - 0.3)
  seg(W / 2, y_s1 + h_sub / 2,  W / 2, y_s2 - h_sub / 2 - 0.3)
  seg(W / 2, y_s2 + h_sub / 2,  W / 2, y_out - h_out / 2 - 0.3)

  # Elbow feeds from the microbial and viral sets
  mc <- hue$microbial[["line"]]; vc <- hue$viral[["line"]]
  seg(xcol[["microbial"]], y_mrg + h_mrg / 2, xcol[["microbial"]], y_s1, mc, 0.8, FALSE)
  seg(xcol[["microbial"]], y_s1, W / 2 - w_sub / 2 - 0.3, y_s1, mc, 0.8)
  seg(xcol[["viral"]], y_mrg + h_mrg / 2, xcol[["viral"]], y_s2, vc, 0.8, FALSE)
  seg(xcol[["viral"]], y_s2, W / 2 + w_sub / 2 + 0.3, y_s2, vc, 0.8)
  grid.text("subtract", mm((xcol[["microbial"]] + W / 2 - w_sub / 2) / 2), Y(y_s1 - 2.5),
            gp = gpar(fontsize = pt_body, col = hue$microbial[["head"]]))
  grid.text("subtract", mm((xcol[["viral"]] + W / 2 + w_sub / 2) / 2), Y(y_s2 - 2.5),
            gp = gpar(fontsize = pt_body, col = hue$viral[["head"]]))

  for (b in boxes) draw_box(b)
}

# Vector PDF with fonts embedded (cairo subsets the Helvetica-metric font), plus
# 1200 dpi PNG and LZW-compressed TIFF (line-art resolution for journal upload).
cairo_pdf(file.path(outdir, "pipeline_schematic.pdf"), width = W / 25.4, height = H / 25.4,
          family = "Helvetica")
draw(); invisible(dev.off())
png(file.path(outdir, "pipeline_schematic.png"), width = W / 25.4, height = H / 25.4,
    units = "in", res = 1200, type = "cairo", family = "Helvetica")
draw(); invisible(dev.off())
tiff(file.path(outdir, "pipeline_schematic.tiff"), width = W / 25.4, height = H / 25.4,
     units = "in", res = 1200, type = "cairo", family = "Helvetica", compression = "lzw")
draw(); invisible(dev.off())
cat("Wrote pipeline_schematic (", W, "x", H, "mm )\n")

# Check every text line fits its box with >= 2 mm padding per side (PDF metrics)
pdf(NULL, width = W / 25.4, height = H / 25.4, family = "Helvetica")
chk <- do.call(rbind, lapply(boxes, function(b) do.call(rbind, lapply(b[[6]], function(l) {
  w <- convertWidth(grobWidth(textGrob(l[[1]], gp = gpar(fontface = l[[2]], fontsize = l[[3]]))),
                    "mm", valueOnly = TRUE)
  data.frame(text = l[[1]], width_mm = round(w, 1), box_mm = b[[4]], fits = w <= b[[4]] - 4)
}))))
invisible(dev.off())
cat("Lines checked:", nrow(chk), " not fitting:", sum(!chk$fits),
    " widest fill:", round(max(chk$width_mm / chk$box_mm), 2), "\n")
if (any(!chk$fits)) print(chk[!chk$fits, ], row.names = FALSE)
