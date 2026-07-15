# Clear environment
rm(list=ls())

## The R script for Žužić et al. manuscript

# Set path to the working directory containing prepared data set

library(tidyverse)
getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = 'right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file) == 0) {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

setwd(getCurrentFileLocation())

library(dplyr)
library(ggplot2)
library(readr)
library(extrafont) # Added extrafont

# (Optional: run loadfonts() if your PDF output isn't rendering the font correctly)
loadfonts(device = "pdf", quiet = TRUE)

# 1. Define the dictionary
strnn <- c(
  "all" = "All events",
  "1" = "Storm 1 (17.-19.03.2015)",
  "2" = "Storm 2 (27.-29.05.2017)",
  "3" = "Storm 3 (07.-09.09.2017)",
  "4" = "Storm 4 (26.-28.09.2017)",
  "1_0760" = "Storm 1 (17.03.2015)",
  "1_0770" = "Storm 1 (18.03.2015)",
  "1_0780" = "Storm 1 (19.03.2015)",
  "2_1470" = "Storm 2 (27.05.2017)",
  "2_1480" = "Storm 2 (28.05.2017)",
  "2_1490" = "Storm 2 (29.05.2017)",
  "3_2500" = "Storm 3 (07.09.2017)",
  "3_2510" = "Storm 3 (08.09.2017)",
  "3_2520" = "Storm 3 (09.09.2017)",
  "4_2690" = "Storm 4 (26.09.2017)",
  "4_2700" = "Storm 4 (27.09.2017)",
  "4_2710" = "Storm 4 (28.09.2017)"
)

var_start <- "total seconds"
var_end <- "horizontal(deg)"
num_groups_first <- 24
num_groups_second <- 5

# Helper function for plotmath expressions
to_pm <- function(v) {
  if (v == 0) return("0")
  e <- floor(log10(abs(v)))
  m <- v / 10^e
  if (e == 0) return(sprintf("%.2f", m))
  sprintf("%.2f %%*%% 10^%d", m, e) 
}

for (nn in names(strnn)) {
  
  file_prefix <- ifelse(nn == "all", "", paste0("_", nn))
  file_name <- paste0("afullframe", file_prefix, "_new.csv")
  
  if (!file.exists(file_name)) next
  df <- read_csv(file_name, show_col_types = FALSE)
  if (nrow(df) == 0) next
  
  print(length(df[["latitude(deg)"]]))
  
  # Format Title
  clean_title <- var_end
  clean_title <- paste0(toupper(substr(clean_title, 1, 1)), substr(clean_title, 2, nchar(clean_title)))
  clean_title <- gsub("(", " error (", clean_title, fixed = TRUE)
  clean_title <- gsub("deg", "m", clean_title, fixed = TRUE)
  clean_title <- gsub("(m)", "(m)", clean_title, fixed = TRUE)
  clean_title <- gsub("(", "[", clean_title, fixed = TRUE)
  clean_title <- gsub(")", "]", clean_title, fixed = TRUE)
  clean_title <- gsub("error", "positioning errors", clean_title, fixed = TRUE)
  plot_title <- paste0(clean_title, " by hour of day\nfor ", tolower(strnn[[nn]]))
  
  # Create limits and bins
  val_bins_start <- seq(0, 60 * 60 * 24, length.out = num_groups_first + 1)
  val_bins_end <- seq(min(df[[var_end]], na.rm = TRUE), max(df[[var_end]], na.rm = TRUE), length.out = num_groups_second + 1)
  
  # Format plotmath legend labels [a, b>
  pm_labels <- sapply(1:(length(val_bins_end) - 1), function(i) {
    v1 <- to_pm(val_bins_end[i])
    v2 <- to_pm(val_bins_end[i + 1])
    if (i == length(val_bins_end) - 1) {
      sprintf("paste('[', %s, ', ', %s, ']')", v1, v2)
    } else {
      sprintf("paste('[', %s, ', ', %s, '>')", v1, v2)
    }
  })
  
  # Process Data
  plot_data <- df %>%
    filter(!is.na(.data[[var_start]]), !is.na(.data[[var_end]])) %>%
    mutate(
      start_bin = cut(.data[[var_start]], breaks = val_bins_start, include.lowest = TRUE, right = FALSE, labels = 0:(num_groups_first - 1)),
      end_bin = cut(.data[[var_end]], breaks = val_bins_end, include.lowest = TRUE, right = FALSE, labels = 1:num_groups_second)
    ) %>%
    count(start_bin, end_bin, .drop = FALSE) %>%
    group_by(start_bin) %>%
    mutate(
      total = sum(n),
      prop = ifelse(total > 0, n / total, 0),
      label_text = ifelse(prop > 0.5, as.character(round(prop * 100)), ""),
      text_color = ifelse(as.numeric(end_bin) >= 4, "white", "black")
    ) %>%
    ungroup()
  
  # Insert the empty "gap" to fit the scale
  plot_data$start_bin <- factor(plot_data$start_bin, levels = c("gap", as.character(0:23)))
  
  # Calculate tangential rotation angles
  plot_data <- plot_data %>%
    mutate(
      theta_deg = (as.numeric(start_bin) - 1) * 360 / 25,
      text_angle = (-theta_deg) %% 360
    )
  
  # Create a separate dataset for the 50% to 100% scale printed inside the top gap
  axis_labels <- data.frame(
    start_bin = factor("gap", levels = c("gap", as.character(0:23))),
    y = seq(0.55, 1.05, 0.1), # Centering text
    label = paste0(seq(50, 100, 10), "%") # Added 0% here
  )
  
  # Create a separate dataset for the manually rotated outer hour labels
  hour_labels <- data.frame(
    start_bin = factor(as.character(0:23), levels = c("gap", as.character(0:23))),
    y = 1.15, # Position just outside the bars
    label = as.character(0:23)
  ) %>%
    mutate(
      theta_deg = (as.numeric(start_bin) - 1) * 360 / 25,
      text_angle = (-theta_deg) %% 360
    )
  
  # Plot
  p <- ggplot(plot_data, aes(x = start_bin, y = prop, fill = end_bin)) +
    
    # Draw radial grid rings (starts at 0)
    geom_hline(yintercept = seq(0, 1, 0.1), color = "gray40", linewidth = 0.4) +
    
    # Draw the wedges
    geom_col(width = 0.85, color = "black", linewidth = 0.5, position = position_stack(reverse = TRUE), family = "DejaVu Sans") +
    
    # Draw the rotating wedge percentage labels 
    geom_text(aes(label = label_text, color = text_color, angle = text_angle), 
              position = position_stack(vjust = 0.9, reverse = TRUE), size = 5, show.legend = FALSE, family = "DejaVu Sans") +
    scale_color_identity() +
    
    # Draw the 0% to 100% labels exactly centered in the gap
    geom_text(data = axis_labels, aes(x = start_bin, y = y, label = label), 
              inherit.aes = FALSE, size = 5, color = "black", hjust = 0.5, vjust = 0.5, family = "DejaVu Sans") +
              
    # Draw the rotating outer hour labels
    geom_text(data = hour_labels, aes(x = start_bin, y = y, label = label, angle = text_angle),
              inherit.aes = FALSE, size = 10, color = "black", hjust = 0.5, vjust = 0.5) +
    
    # Rotate backwards by exactly half a slice (pi/25) to center the gap at 12 o'clock
    coord_polar(theta = "x", start = -pi/25) +
    
    scale_fill_viridis_d(labels = parse(text = pm_labels), option = "plasma", direction = -1) +
    
    # Drop=FALSE forces the gap to be drawn. Hide the default axis text.
    scale_x_discrete(drop = FALSE) +
    
    # Expand limits into the negative to create the hollow center, and 1.2 to fit the outer labels
    scale_y_continuous(limits = c(-0.3, 1.2), breaks = seq(0, 1, 0.1)) +
    
    # Standard legend order
    guides(fill = guide_legend(ncol = 2, byrow = FALSE)) +
    
    theme_void(base_size = 12, base_family = "DejaVu Sans") +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 24, family = "DejaVu Sans"),
      # Pull the legend UP into the dead space by 30 pixels
      legend.margin = margin(t = -30, b = 0),
      
      # Keep the title centered, but remove the bottom padding
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(t = -20, b = 0), family = "DejaVu Sans"),
      
      legend.key.size = unit(1.2, "lines"),
      
      # Strip all outer plot margins entirely
      plot.margin = margin(0, 0, 0, 0),
      
      axis.text.x = element_blank(), 
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank()
    ) +
    labs(title = plot_title)
  
  save_name <- paste0("afullframe", file_prefix, "_", var_start, " ", var_end, "_R.pdf")
  ggsave(save_name, p, width = 6, height = 8, device = cairo_pdf)
  cat("Saved:", save_name, "\n")
}