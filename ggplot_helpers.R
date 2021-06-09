# GGPLOT HELPERS
# ===============================================================================================

theme_custom <- function(legend = "none", x.text.angle = 0, margin = T, base_size = 12, border = T){

  theme_pubr(legend = legend, border = border, x.text.angle = x.text.angle, margin = margin, base_size = base_size) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic"),
          axis.title = element_text(face = "bold"),
          legend.title = element_blank())
}
