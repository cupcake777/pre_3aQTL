library(tidyverse)
library(ggplot2)
library(ggsci)
library(patchwork)
gw_labeller <- function(x) {
  x[is.na(x)] <- 0
  sapply(x, function(val) {
    if (val < -0.01) { 
      gw_val <- round((val * 365 / 7) + 40, 0)
      return(paste0(gw_val, "GW"))
    } else {
      return(as.character(val))
    }
  })
}

meta_file <- "../../raw_data/CHB.all.meta"
if(!file.exists(meta_file)) stop("Meta file not found!")

df <- read.table(meta_file, header = TRUE, stringsAsFactors = FALSE)

# Filter out samples with "_rep" suffix to keep only genotype data
plot_data <- df %>%
  filter(!grepl("_rep$", sample)) %>%  # Remove samples ending with "_rep"
  mutate(
    Sex_Label = factor(sex, levels = c(0, 1), labels = c("Male", "Female")),
    LifeStage = factor(LifeStage, levels = c("Stage1", "Stage2", "Stage3")),
    Batch = as.factor(batch)
  )

message(paste("Total Samples Loaded:", nrow(plot_data)))

stage_color <- c("Stage1" = "#56B4E9", 
                 "Stage2" = "#FFB482", 
                 "Stage3" = "#8DE5A1")


p1 <- ggplot(plot_data, aes(x = Age, y = RIN, color = LifeStage)) +
  geom_point(aes(shape = Sex_Label), alpha = 0.7, size = 2.5) + 
  geom_smooth(method = "loess", se = FALSE, color = "grey40", linetype = "dashed", linewidth = 0.8, alpha = 0.5) +
  facet_grid(~LifeStage, scales = "free_x") +
  scale_color_manual(values = stage_color) +
  scale_x_continuous(labels = gw_labeller) +
  scale_y_continuous(limits = c(min(plot_data$RIN) - 0.5, 10), breaks = seq(4, 10, 1)) +
  
  labs(
    title = "A. Sample Distribution & Quality (N=761)",
    subtitle = "Age span: 14GW to 106 Years",
    x = "Age (GW / Years)",
    y = "RNA Integrity Number",
    color = "Stage",
    shape = "Sex"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "#E8F5E9"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.box = "horizontal", 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )


stats_df <- plot_data %>%
  group_by(LifeStage, Sex_Label) %>%
  summarise(Count = n(), .groups = "drop")

p2 <- ggplot(stats_df, aes(x = LifeStage, y = Count, fill = Sex_Label)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", linewidth = 0.3) +
  
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), 
            size = 4, color = "white", fontface = "bold") +
  
  scale_fill_manual(values = c("Male" = "#FFD700", "Female" = "#7EC0EE")) +
  
  labs(
    title = "B. Demographics",
    subtitle = "Sample Count by Sex",
    x = NULL,
    y = "Number of Samples",
    fill = "Sex"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1), 
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )


final_plot <- p1 + p2 + 
  plot_layout(widths = c(2.5, 1))

outfile <- "Dataset_Overview_RealData.pdf"
ggsave(outfile, final_plot, width = 16, height = 7) 
message("Finish! check your plot~")
