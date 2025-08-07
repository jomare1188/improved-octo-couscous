# PCA plotting functions for RNA-seq analysis
# Author: Jorge Mario Muñoz Pérez


#' Plot PCA from raw data
#' 
#' @param data Input data for PCA
#' @param color Color variable for points
#' @param shape Shape variable for points
#' @param title Plot title
#' @param var Variance explained by each PC
#' @param file Output filename
#' @param lab_color Label for color legend
#' @return ggplot object
plot_pca <- function(data, color, shape, title, var, file, lab_color) {
  ggplot(df, aes(x = PC1, y = PC2, colour = color, shape = shape)) +
    geom_point(size = 4.5) +
    theme_bw(base_size = 20) +
    labs(title = title, shape = "Group", color = lab_color) +
    labs(x = paste0("PC1: ", round(var[1] * 100, 1), "%"),
         y = paste0("PC2: ", round(var[2] * 100, 1), "%")) +
    theme(text = element_text(size = 20),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(colour = "black", angle = 8, vjust = 0.7, hjust = 0.5),
          axis.text.y = element_text(colour = "black"),
          axis.line = element_line(colour = "black"))
  
  ggsave(filename = file, units = "cm", width = 25, height = 25, dpi = 320)
}

#' Plot PCA using DESeq2 plotPCA function with custom aesthetics
#' 
#' @param data DESeq2 object (vst or rlog transformed)
#' @param title Plot title
#' @param grouping_vars Vector of column names from colData to use for grouping
#' @param n_top Number of top genes to use for PCA
#' @param path Output file path
#' @param color_var Column name to map to color aesthetic (optional)
#' @param shape_var Column name to map to shape aesthetic (optional)
#' @param size_var Column name to map to size aesthetic (optional)
#' @return ggplot object
plot_pca_deseq_custom <- function(data, title, grouping_vars, n_top, path,
                                  color_var = NULL, shape_var = NULL,
                                  size_var = NULL) {
  
  # Generate PCA data using DESeq2
  pca_data <- plotPCA(data, intgroup = grouping_vars, returnData = TRUE, ntop = n_top)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  # Create aesthetic mapping
  aes_mapping <- aes(PC1, PC2)
  
  if (!is.null(color_var) && color_var %in% grouping_vars) {
    aes_mapping$colour <- sym(color_var)
    aes_mapping$fill <- sym(color_var)
  }
  
  if (!is.null(shape_var) && shape_var %in% grouping_vars) {
    aes_mapping$shape <- sym(shape_var)
  }
  
  if (!is.null(size_var) && size_var %in% grouping_vars) {
    aes_mapping$size <- sym(size_var)
  }
  
  # Base plot with dynamic aesthetics
  p <- ggplot(pca_data, aes_mapping) +
    theme_bw() +
    labs(title = title,
         x = paste0("PC1: ", percentVar[1], "% variance"),
         y = paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 22),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black"))
  
  # Add geom_point with conditional size
  if (!is.null(size_var) && size_var %in% grouping_vars) {
    p <- p + geom_point(stroke = 1)  # No fixed size when using size aesthetic
  } else {
    p <- p + geom_point(size = 4.5, stroke = 1)  # Fixed size when not using size aesthetic
  }
  
  # Add labels and scales
  if (!is.null(color_var) && color_var %in% grouping_vars) {
    p <- p + labs(color = color_var) +
      guides(fill = "none")
  }
  
  if (!is.null(shape_var) && shape_var %in% grouping_vars) {
    p <- p + labs(shape = shape_var) +
      scale_shape_manual(values = c(21, 22, 23, 24, 25, 8))
  }
  
  if (!is.null(size_var) && size_var %in% grouping_vars) {
    p <- p + labs(size = size_var) +
      scale_size_manual(values = c(4, 6))  # Customize sizes for better differentiation
  }
  
  # Save plot
  ggsave(plot = p, filename = path, units = "cm",
         width = 15 * 2, height = 15 * 2, dpi = 320)
  
  return(p)
}
