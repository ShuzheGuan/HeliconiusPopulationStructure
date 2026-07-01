# Purpose: plot best qpGraph models saved by find_graphs().
# Input: find_graphs .rds files and a matching precomputed F2 block directory.
# Output: one qpGraph figure PDF per .rds file.
# Software: R, admixtools, ggplot2, ggtext

suppressPackageStartupMessages({
  library(admixtools)
  library(ggplot2)
  library(ggtext)
})

QPGRAPH_RDS_DIR <- "<path/to/find_graphs_output_dir>"
F2_DIR <- "<path/to/f2_blocks_dir>"
OUTPUT_DIR <- "<path/to/qpgraph_figure_output_dir>"

WIDTH_IN <- 7
HEIGHT_IN <- 5
FONT_FAMILY <- "Helvetica"
TEXT_SIZE <- 2.5
TIP_LABEL_OFFSET <- 0.25

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

TIP_LABELS <- c(
  sample_or_population_id = "display label"
)

tip_label <- function(x) {
  x <- as.character(x)
  if (x %in% names(TIP_LABELS)) return(TIP_LABELS[[x]])
  x
}

plot_qpgraph <- function(edges) {
  edges <- as.data.frame(edges)
  plot_data <- admixtools:::graph_to_plotdat(edges, fix = TRUE, fix_down = TRUE)

  plot_data$eg$linetype_plot <- plot_data$eg$type
  admix_edges <- edges[edges$type == "admix", , drop = FALSE]

  for (node in unique(admix_edges$to)) {
    incoming <- admix_edges[admix_edges$to == node, , drop = FALSE]
    major_from <- incoming$from[which.max(incoming$weight)]
    is_major <- plot_data$eg$to == node & plot_data$eg$name == major_from
    plot_data$eg$linetype_plot[is_major] <- "normal"
  }

  plot_data$eg$label <- ifelse(plot_data$eg$type == "admix", plot_data$eg$label, "")

  tip_direction <- plot_data$eg[plot_data$eg$to %in% plot_data$nodes$name, , drop = FALSE]
  tip_direction$dx <- tip_direction$xend - tip_direction$x
  tip_direction$dy <- tip_direction$yend - tip_direction$y
  tip_direction <- tip_direction[!duplicated(tip_direction$to), c("to", "dx", "dy")]
  names(tip_direction)[1] <- "name"

  plot_data$nodes <- merge(plot_data$nodes, tip_direction, by = "name", all.x = TRUE)
  plot_data$nodes$label <- vapply(plot_data$nodes$name, tip_label, character(1))
  plot_data$nodes$dx[is.na(plot_data$nodes$dx)] <- 1
  plot_data$nodes$dy[is.na(plot_data$nodes$dy)] <- 0
  branch_len <- pmax(sqrt(plot_data$nodes$dx^2 + plot_data$nodes$dy^2), 1e-9)
  plot_data$nodes$label_x <- plot_data$nodes$x + TIP_LABEL_OFFSET * plot_data$nodes$dx / branch_len
  plot_data$nodes$label_y <- plot_data$nodes$y + TIP_LABEL_OFFSET * plot_data$nodes$dy / branch_len

  tree_edges <- plot_data$eg[plot_data$eg$linetype_plot == "normal", , drop = FALSE]
  dashed_edges <- plot_data$eg[plot_data$eg$linetype_plot == "admix", , drop = FALSE]
  dashed_edges$x_arrow <- dashed_edges$x + 0.98 * (dashed_edges$xend - dashed_edges$x)
  dashed_edges$y_arrow <- dashed_edges$y + 0.98 * (dashed_edges$yend - dashed_edges$y)

  ggplot(plot_data$eg, aes(x = x, xend = xend, y = y, yend = yend)) +
    geom_segment(data = tree_edges) +
    geom_segment(data = dashed_edges, aes(xend = x_arrow, yend = y_arrow), linetype = "dashed") +
    geom_segment(
      data = dashed_edges,
      aes(x = x_arrow, y = y_arrow, xend = xend, yend = yend),
      inherit.aes = FALSE,
      arrow = arrow(type = "open", angle = 22, length = grid::unit(0.05, "inches"))
    ) +
    geom_richtext(
      data = plot_data$nodes,
      aes(x = label_x, y = label_y, label = label),
      inherit.aes = FALSE,
      size = TEXT_SIZE,
      family = FONT_FAMILY,
      fill = NA,
      label.color = NA
    ) +
    geom_text(
      aes(x = (x + xend) / 2, y = (y + yend) / 2, label = label),
      size = TEXT_SIZE,
      nudge_y = -0.15
    ) +
    theme_void(base_family = FONT_FAMILY) +
    theme(legend.position = "none")
}

plot_rds <- function(rds_path) {
  prefix <- sub("[.]rds$", "", basename(rds_path))

  results <- readRDS(rds_path)
  best <- results[which.min(results$score), , drop = FALSE]
  graph <- best$graph[[1]]
  fit <- qpgraph(
    f2_from_precomp(F2_DIR),
    graph,
    return_fstats = FALSE,
    return_pvalue = FALSE,
    verbose = FALSE
  )

  figure <- plot_qpgraph(fit$edges)
  ggsave(file.path(OUTPUT_DIR, paste0(prefix, ".pdf")), figure,
         width = WIDTH_IN, height = HEIGHT_IN, units = "in")
}

rds_files <- sort(list.files(QPGRAPH_RDS_DIR, pattern = "[.]rds$", full.names = TRUE))
invisible(lapply(rds_files, plot_rds))
