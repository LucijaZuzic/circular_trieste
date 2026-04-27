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
library(tidyr)
library(lubridate)
library(geosphere)
library(ggplot2)
library(readr)
library(stringr)
library(extrafont) # Added extrafont

# (Optional: run loadfonts() if your PDF output isn't rendering the font correctly)
loadfonts(device = "pdf", quiet = TRUE)

# 1. Initialization and reference coordinates
strr <- c(
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

ref_coords <- list(
  lon = 13.763519885584342,
  lat = 45.70975671789818,
  height = 323.4
)

# Shared theme
bigger_theme <- theme_classic(base_size = 20, base_family = "DejaVu Sans") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),       # Pure black text
    axis.ticks = element_line(color = "black", linewidth = 1),  # Pure black tick marks
    axis.line = element_line(color = "black", linewidth = 1),   # Pure black axis lines
    legend.text = element_text(size = 20, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

full_data <- list()
case_data <- list()
date_data <- list()

# 2. File parsing and distance computation (vectorized)
for (stormnum in 1:4) {
  storm_dir <- paste0("Storm ", stormnum)
  if (!dir.exists(storm_dir)) next
  
  case_data[[as.character(stormnum)]] <- list()
  date_data[[as.character(stormnum)]] <- list()
  
  files <- list.files(storm_dir, pattern = "\\.pos$", full.names = TRUE)
  
  for (file_path in files) {
    filename <- basename(file_path)
    if (grepl("events|stat", filename)) next
    
    shortfile <- gsub(".pos|trie", "", filename)
    
    # Read the file header to extract column names
    raw_lines <- readLines(file_path, n = 20)
    header_line <- raw_lines[15] # skip = 14
    
    # 1. Clean the header and remove the "%" and "GPST" keywords
    header_line <- gsub("%|GPST", "", header_line) 
    # 2. Trim excess spaces
    header_line <- trimws(gsub("\\s+", " ", header_line)) 
    
    # 3. Split into remaining column names and drop any accidental empty strings
    data_cols <- unlist(strsplit(header_line, " "))
    data_cols <- data_cols[data_cols != ""] 
    
    # 4. Prepend date and time to perfectly match the 15 data columns
    final_col_names <- c("date", "time", data_cols)
    
    # Read the data block using the perfectly aligned names
    df_raw <- read_table(file_path, skip = 15, col_names = final_col_names, show_col_types = FALSE)
    
    if (nrow(df_raw) == 0) next
    
    # Vectorized timestamp parsing
    df <- df_raw %>%
      mutate(
        datetime = ymd_hms(paste(gsub("/", "-", date), time)),
        d = day(datetime),
        m = month(datetime),
        y = year(datetime),
        hr = hour(datetime),
        min = minute(datetime),
        sec = second(datetime),
        `total seconds` = hr * 3600 + min * 60 + sec
      )
    
    # Vectorized geodesic distance calculations
    df <- df %>%
      mutate(
        # Point 1 is a simple vector, Point 2 uses cbind() to create an N x 2 matrix
        dy = distGeo(c(ref_coords$lon, ref_coords$lat), 
                     cbind(ref_coords$lon, `latitude(deg)`)),
        
        dx = distGeo(c(ref_coords$lon, ref_coords$lat), 
                     cbind(`longitude(deg)`, ref_coords$lat)),
        
        `absolute height(m)` = abs(ref_coords$height - `height(m)`),
        `absolute latitude(deg)` = abs(dy),
        `absolute longitude(deg)` = abs(dx),
        
        `horizontal(deg)` = distGeo(c(ref_coords$lon, ref_coords$lat), 
                                    cbind(`longitude(deg)`, `latitude(deg)`)),
        
        ang = atan2(dy, dx) / pi * 180,
        `horizontal angle(deg)` = ang,
        `absolute horizontal angle(deg)` = abs(ang),
        `minimum horizontal angle(deg)` = ifelse(abs(ang) > 90, 180 - abs(ang), abs(ang))
      ) %>%
      select(-date, -time, -datetime, -dy, -dx, -ang)
    
    # Store in memory hierarchies
    full_data[[length(full_data) + 1]] <- df
    case_data[[as.character(stormnum)]][[length(case_data[[as.character(stormnum)]]) + 1]] <- df
    date_data[[as.character(stormnum)]][[shortfile]] <- df
    
    print(paste(stormnum, filename, ncol(df), nrow(df)))
  }
}

# Combine lists into dataframes
full_df <- bind_rows(full_data)
write_csv(full_df, "afullframe_new.csv")

if ("GPST" %in% names(full_df)) {
  print(length(full_df$GPST))
}

makeplot <- TRUE
printpart <- TRUE
printday <- FALSE
hrlab <- 0
daylab <- 10
monthlab <- 10

rdim <- list(hr = 0, len = 0, min = 3, Q1 = 2, Q2 = 2, avg = 2, Q3 = 2, max = 3, SD = 2)
nks <- c("hr", "len", "min", "Q1", "Q2", "avg", "Q3", "max", "SD")

# 3. Storm day generation (statistics and plots)
for (stormnum in names(case_data)) {
  storm_df <- bind_rows(case_data[[stormnum]])
  write_csv(storm_df, paste0("afullframe_", stormnum, "_new.csv"))
  
  for (shortfile in names(date_data[[stormnum]])) {
    day_df <- date_data[[stormnum]][[shortfile]]
    write_csv(day_df, paste0("afullframe_", stormnum, "_", shortfile, "_new.csv"))
    
    # Plotting variables
    maxofs <- 1
    stepv <- 2
    mind <- min(day_df$d)
    minm <- min(day_df$m)
    miny <- min(day_df$y)
    
    # Stats calculation
    stats_df <- day_df %>%
      group_by(hr) %>%
      summarise(
        len = n(),
        min = min(`horizontal(deg)`, na.rm = TRUE),
        Q1 = quantile(`horizontal(deg)`, 0.25, na.rm = TRUE),
        Q2 = median(`horizontal(deg)`, na.rm = TRUE),
        avg = mean(`horizontal(deg)`, na.rm = TRUE),
        Q3 = quantile(`horizontal(deg)`, 0.75, na.rm = TRUE),
        max = max(`horizontal(deg)`, na.rm = TRUE),
        SD = sd(`horizontal(deg)`, na.rm = TRUE)
      ) %>%
      complete(hr = 0:23, fill = list(len=0, min=0, Q1=0, Q2=0, avg=0, Q3=0, max=0, SD=0)) %>%
      mutate(SD = replace_na(SD, 0))
    
    # LaTeX table print loop (day)
    if (printday) {
      for (vix in 3:length(nks)) {
        col_name <- nks[vix]
        max_val <- max(stats_df[[col_name]])
        min_val <- min(stats_df[[col_name]])
        max_idx <- which.max(stats_df[[col_name]]) - 1
        min_idx <- which.min(stats_df[[col_name]]) - 1
        cat(sprintf("$%s$ $%s$ $%s$ $%s$ $%s$\n", col_name, max_idx, round(max_val, rdim[[col_name]]), min_idx, round(min_val, rdim[[col_name]])))
      }
      
      seenvals <- setNames(replicate(length(nks), character(0), simplify = FALSE), nks)
      
      for (h in 0:23) {
        h1 <- sprintf("%02d", h)
        pvals <- paste0("$", h1, "$")
        smz <- sum(stats_df[stats_df$hr == h, 3:9])
        
        for (vix in 3:length(nks)) {
          col_name <- nks[vix]
          val <- stats_df[[col_name]][h + 1]
          strv <- as.character(round(val, rdim[[col_name]]))
          
          nonzero_vals <- stats_df[[col_name]][stats_df[[col_name]] != 0]
          s1 <- if(length(nonzero_vals) > 0) as.character(round(max(nonzero_vals), rdim[[col_name]])) else "0"
          s2 <- if(length(nonzero_vals) > 0) as.character(round(min(nonzero_vals), rdim[[col_name]])) else "0"
          
          if (strv == s1) pvals <- c(pvals, paste0("$\\mathbf{", strv, "}$"))
          else if (strv == s2) pvals <- c(pvals, paste0("$\\underline{\\mathbf{", strv, "}}$"))
          else pvals <- c(pvals, paste0("$", strv, "$"))
          
          if (!(strv %in% seenvals[[col_name]])) {
            seenvals[[col_name]] <- c(seenvals[[col_name]], strv)
          } else if (strv == s1 || strv == s2) {
            cat("seen me", col_name, h, strv, "\n")
          }
        }
        if (smz > 0) cat(paste(pvals, collapse = " & "), " \\\\ \\hline\n")
      }
    }
    
    # Make plot (day)
    if (makeplot) {
      secticks <- seq(0, (maxofs * 25 - maxofs) * 3600, by = stepv * 3600)
      seclabs <- sprintf(ifelse(stepv < hrlab, "%02d", "%d"), (seq(0, maxofs * 25 - maxofs, stepv) %% 24))
      
      file_key <- paste0(stormnum, "_", shortfile)
      ttln <- ifelse(file_key %in% names(strr), strr[[file_key]], "")
      
      plot_title <- paste0("Time series of horizontal positioning errors [m]\nfor ", tolower(ttln))
      
      p <- ggplot(day_df, aes(x = `total seconds`, y = `horizontal(deg)`)) +
        geom_line(color = "blue", linewidth = 1) +
        scale_x_continuous(breaks = secticks, labels = seclabs, limits = c(min(secticks), max(secticks)), expand = c(0, 0)) +
        labs(title = plot_title, x = "Hour of day", y = "Horizontal\npositioning\nerrors [m]") +
        bigger_theme
      
      v1 <- min(day_df$`horizontal(deg)`); v2 <- max(day_df$`horizontal(deg)`)
      
      for (ofs in 0:maxofs) {
        txtt1 <- sprintf(ifelse(mind + ofs < daylab, " %02d.%02d.%d", " %d.%02d.%d"), mind + ofs, minm, miny)
        if (ofs > 0 && ofs < maxofs) {
          p <- p + geom_vline(xintercept = ofs * 24 * 3600, linetype = "dashed", color = "gray")
        }
        if (ofs < maxofs) {
          p <- p + annotate("text", x = ofs * 24 * 3600, y = v2, label = txtt1, vjust = 1, hjust = 0, size = 6, family = "DejaVu Sans")
        }
      }
      ggsave(paste0("afullframe_", stormnum, "_", shortfile, "_time_R.pdf"), plot = p, width = 10, height = 4.5, device = cairo_pdf)
    }
  }
  
  # 4. Storm aggregation generation (statistics and plots)
  maxofs <- 3
  stepv <- 4
  mind <- min(storm_df$d)
  minm <- min(storm_df$m)
  miny <- min(storm_df$y)
  
  storm_df <- storm_df %>%
    mutate(sec_adj = (d - mind) * 3600 * 24 + `total seconds`)
  
  stats_storm <- storm_df %>%
    group_by(hr) %>%
    summarise(
      len = n(),
      min = min(`horizontal(deg)`, na.rm = TRUE),
      Q1 = quantile(`horizontal(deg)`, 0.25, na.rm = TRUE),
      Q2 = median(`horizontal(deg)`, na.rm = TRUE),
      avg = mean(`horizontal(deg)`, na.rm = TRUE),
      Q3 = quantile(`horizontal(deg)`, 0.75, na.rm = TRUE),
      max = max(`horizontal(deg)`, na.rm = TRUE),
      SD = sd(`horizontal(deg)`, na.rm = TRUE)
    ) %>%
    complete(hr = 0:23, fill = list(len=0, min=0, Q1=0, Q2=0, avg=0, Q3=0, max=0, SD=0)) %>%
    mutate(SD = replace_na(SD, 0))
  
  # LaTeX table print loop (storm)
  if (printpart) {
    for (vix in 3:length(nks)) {
      col_name <- nks[vix]
      cat(sprintf("$%s$ $%s$ $%s$ $%s$ $%s$\n", col_name, which.max(stats_storm[[col_name]]) - 1, 
                  round(max(stats_storm[[col_name]]), rdim[[col_name]]), 
                  which.min(stats_storm[[col_name]]) - 1, round(min(stats_storm[[col_name]]), rdim[[col_name]])))
    }
    for (h in 0:23) {
      h1 <- sprintf("%02d", h)
      pvals <- paste0("$", h1, "$")
      for (vix in 3:length(nks)) {
        col_name <- nks[vix]
        val <- stats_storm[[col_name]][h + 1]
        strv <- as.character(round(val, rdim[[col_name]]))
        s1 <- as.character(round(max(stats_storm[[col_name]]), rdim[[col_name]]))
        s2 <- as.character(round(min(stats_storm[[col_name]]), rdim[[col_name]]))
        
        if (strv == s1) pvals <- c(pvals, paste0("$\\mathbf{", strv, "}$"))
        else if (strv == s2) pvals <- c(pvals, paste0("$\\underline{\\mathbf{", strv, "}}$"))
        else pvals <- c(pvals, paste0("$", strv, "$"))
      }
      cat(paste(pvals, collapse = " & "), " \\\\ \\hline\n")
    }
  }
  
  # Make plot (storm)
  if (makeplot) {
    secticks <- seq(0, (maxofs * 25 - maxofs) * 3600, by = stepv * 3600)
    seclabs <- sprintf(ifelse(stepv < hrlab, "%02d", "%d"), (seq(0, maxofs * 25 - maxofs, stepv) %% 24))
    
    ttln <- ifelse(as.character(stormnum) %in% names(strr), strr[[as.character(stormnum)]], "")
    plot_title <- paste0("Time series of horizontal positioning errors [m]\nfor ", tolower(ttln))
    
    p <- ggplot(storm_df, aes(x = sec_adj, y = `horizontal(deg)`)) +
      geom_line(color = "blue", linewidth = 1) +
      scale_x_continuous(breaks = secticks, labels = seclabs, limits = c(min(secticks), max(secticks)), expand = c(0, 0)) +
      labs(title = plot_title, x = "Hour of day", y = "Horizontal\npositioning\nerrors [m]") +
      bigger_theme
    
    v1 <- min(storm_df$`horizontal(deg)`); v2 <- max(storm_df$`horizontal(deg)`)
    
    for (ofs in 0:maxofs) {
      txtt1 <- sprintf(ifelse(mind + ofs < daylab, " %02d.%02d.%d", " %d.%02d.%d"), mind + ofs, minm, miny)
      if (ofs > 0 && ofs < maxofs) {
        p <- p + geom_vline(xintercept = ofs * 24 * 3600, linetype = "dashed", color = "gray")
      }
      if (ofs < maxofs) {
        p <- p + annotate("text", x = ofs * 24 * 3600, y = v2, label = txtt1, vjust = 1, hjust = 0, size = 6, family = "DejaVu Sans")
      }
    }
    ggsave(paste0("afullframe_", stormnum, "_time_R.pdf"), plot = p, width = 10, height = 4.5, device = cairo_pdf)
  }
}