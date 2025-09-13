library(stringr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(lubridate)

# Define a list of European countries
european_countries <- c(
  "Albania", "Andorra", "Austria", "Belarus", "Belgium", "Bosnia and Herzegovina",
  "Bulgaria", "Croatia", "Czech Republic", "Denmark", "Estonia", "Finland",
  "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy",
  "Latvia", "Liechtenstein", "Lithuania", "Luxembourg", "Malta", "Moldova",
  "Monaco", "Montenegro", "Netherlands", "North Macedonia", "Norway", "Poland",
  "Portugal", "Romania", "Russia", "San Marino", "Serbia", "Slovakia", "Slovenia",
  "Spain", "Sweden", "Switzerland", "Ukraine", "United Kingdom", "Vatican City"
)
# ------------------- Read clinical and wastewater data ----------------------

clinical_data <- read.csv(
  'metadata/rsv-a_metadata_2025-08-17T2200.tsv',
  sep = '\t'
)[, c('lineage', 'sampleCollectionDate', 'geoLocCountry')]
clinical_data$counts <- 1
clinical_data$sampleCollectionDate <- as.Date(clinical_data$sampleCollectionDate)

wastewater_data <- read.csv(
  '../../../preprint/data/lollipop/rsva/deconvolved.csv',
  sep = '\t'
)
wastewater_data$date <- as.Date(wastewater_data$date)




# Filter clinical data (filter by calendorical weeks (starting on Monday))
clinical_data <- clinical_data %>%
  filter(sampleCollectionDate > "2023-10-29",
         sampleCollectionDate < "2024-03-04",
         geoLocCountry %in% european_countries) ## available European countries


# Select the samples of only those lineages that were detected in Switzerland
swiss_sub_lineages <- c("A.D.1", "A.D.1.5", "A.D.1.6", "A.D.2.1",
                        "A.D.3", "A.D.3.1", "A.D.5.1", "A.D.5.2")
swiss_lineages <- c("A.D.1","A.D.2", "A.D.3", "A.D.5")


clinical_aggregated <- clinical_data

# Aggregate sublineages into lineages
clinical_aggregated$lineage <- gsub("^A\\.D\\.1.*", "A.D.1", clinical_aggregated$lineage)
clinical_aggregated$lineage <- gsub("^A\\.D\\.2.*", "A.D.2", clinical_aggregated$lineage)
clinical_aggregated$lineage <- gsub("^A\\.D\\.3.*", "A.D.3", clinical_aggregated$lineage)
clinical_aggregated$lineage <- gsub("^A\\.D\\.5.*", "A.D.5", clinical_aggregated$lineage)

clinical_aggregated <- clinical_aggregated %>%
  filter(lineage %in% swiss_lineages)

clinical_aggregated <- clinical_aggregated %>%
  group_by(lineage, sampleCollectionDate, geoLocCountry) %>%
  summarise(counts = sum(counts, na.rm = TRUE), .groups = "drop") %>%
  mutate(sampleCollectionDate = as.Date(sampleCollectionDate)) %>%
  mutate(week = floor_date(sampleCollectionDate, unit = "week")) %>%
  group_by(variant = lineage, week) %>%
  summarise(lineage_counts = sum(counts, na.rm = TRUE), .groups = "drop")


# When there is no clinical samples for a specific lineage at specific week -> impute zeros
weeks <- unique(clinical_aggregated$week)
clinical_aggregated <- clinical_aggregated %>%
  complete(week = unique(week),
           variant = unique(variant),
           fill = list(frequency = 0, lineage_counts = 0)) %>%
  arrange(variant, week) %>%
  group_by(week) %>%
  mutate(total_counts = sum(lineage_counts),
         frequency_clin = lineage_counts / total_counts) %>%
  ungroup()


# Step 1: Combine WW abundances from Zurich and Geneva with inverse variance weighting
# Assume ww_zh and ww_ge are data frames with week, variant, proportion, proportionLower, proportionUpper
ww_combined <- wastewater_data %>%  # Combine Zurich and Geneva
  mutate(week = floor_date(date, unit = "week")) %>%
  dplyr::group_by(variant, week, location) %>%
  #dplyr::select(variant, location, week, proportion) %>%
  group_by(variant, week, location) %>%
  mutate(
    variance = ((proportionUpper - proportionLower) / (2 * 1.96))^2,  # Variance from 95% CI
    weight = 1 / variance  # Inverse variance
  ) %>%
  summarise(
    combined_proportion = weighted.mean(proportion, weight, na.rm = TRUE),
    .groups = "drop"
  )  %>%
  rename(proportion = combined_proportion)
  # first average
  #summarise(proportion = mean(proportion, na.rm = TRUE), .groups = "drop")

  

# then aggregate sublineages into lineages
ww_combined$variant <- gsub("^A\\.D\\.1.*", "A.D.1", ww_combined$variant)
ww_combined$variant <- gsub("^A\\.D\\.2.*", "A.D.2", ww_combined$variant)
ww_combined$variant <- gsub("^A\\.D\\.3.*", "A.D.3", ww_combined$variant)
ww_combined$variant <- gsub("^A\\.D\\.5.*", "A.D.5", ww_combined$variant)
ww_combined <- ww_combined %>%
  filter(variant %in% c("A.D.1", "A.D.2", "A.D.3", "A.D.5")) %>%
  group_by(variant, week, location) %>%
# then sum
  summarise(proportion = sum(proportion, na.rm = TRUE), .groups = "drop") %>%
  ungroup() %>%
  group_by(variant, week, location)

# 

ww_combined <- ww_combined %>% 
  group_by(week, variant) %>%
  #summarise(freq_all_var = sum(proportion), .groups = "drop") %>%
  pivot_wider(names_from = variant, values_from = proportion) %>%
  arrange(week)

ww_combined_geneva <- ww_combined %>%
  filter(location == "Genève (GE)")

ww_combined_zurich <- ww_combined %>%
  filter(location == "Zürich (ZH)")
# Add week back if needed

# Step 2: Smooth clinical data (assuming clin_with_week has week and variants)
all_weeks <- sort(unique(clinical_aggregated$week))
clinical_aggregated <- clinical_aggregated %>%
  dplyr::select(week,variant, frequency_clin, total_counts) %>%
  pivot_wider(names_from = variant, values_from = frequency_clin) %>%
  arrange(week)

# Estimate the distance between clinical and wastewater estimates for different lineages
clin_ww_combined_zh <- merge(ww_combined_zurich, clinical_aggregated, by = "week", all = TRUE)
clin_ww_combined_ge <- merge(ww_combined_geneva, clinical_aggregated, by = "week", all = TRUE)
clin_ww_combined_zh <- clin_ww_combined_zh %>%
  rowwise() %>%
  mutate(
    delta_abs_AD1 = abs(`A.D.1.x` - `A.D.1.y`),
    delta_abs_AD2 = abs(`A.D.2.x` - `A.D.2.y`),
    delta_abs_AD3 = abs(`A.D.3.x` - `A.D.3.y`),
    delta_abs_AD5 = abs(`A.D.5.x` - `A.D.5.y`)
  )
clin_ww_combined_ge <- clin_ww_combined_ge %>%
  rowwise() %>%
  mutate(
    delta_abs_AD1 = abs(`A.D.1.x` - `A.D.1.y`),
    delta_abs_AD2 = abs(`A.D.2.x` - `A.D.2.y`),
    delta_abs_AD3 = abs(`A.D.3.x` - `A.D.3.y`),
    delta_abs_AD5 = abs(`A.D.5.x` - `A.D.5.y`)
  )


max_distance_zh <- max(clin_ww_combined_zh[,c("delta_abs_AD1","delta_abs_AD2","delta_abs_AD3","delta_abs_AD5")], na.rm = TRUE)
max_distance_ge <- max(clin_ww_combined_ge[,c("delta_abs_AD1","delta_abs_AD2","delta_abs_AD3","delta_abs_AD5")], na.rm = TRUE)

# Pivot longer

delta_long_zh <- clin_ww_combined_zh %>%
  dplyr::select(week, starts_with("delta_abs_AD"), total_counts) %>%
  pivot_longer(
    cols = starts_with("delta_abs_AD"),
    names_to = "lineage",
    values_to = "delta_abs"
  )
delta_long_ge <- clin_ww_combined_ge %>%
  dplyr::select(week, starts_with("delta_abs_AD"),total_counts) %>%
  pivot_longer(
    cols = starts_with("delta_abs_AD"),
    names_to = "lineage",
    values_to = "delta_abs"
  )
# Visualize

max_delta <- max_distance_zh # proportion/frequency goes 0-1
max_count <- max(delta_long_zh$total_counts, na.rm = TRUE)
custom_labels <- c(
  "delta_abs_AD1" = "A.D.1",
  "delta_abs_AD2" = "A.D.2",
  "delta_abs_AD3" = "A.D.3",
  "delta_abs_AD5" = "A.D.5"
  
)
scale_factor <-max_delta/ max_count

bars_df_zh <- delta_long_zh %>%
  dplyr::select(week, total_counts) %>%
  distinct()
p_zurich <- ggplot(delta_long_zh, aes(x = week, y = delta_abs)) + 
  facet_wrap(~ lineage, scales = "fixed", labeller = as_labeller(custom_labels))+
  geom_col(
    data = bars_df_zh,
    aes(x = week, y = total_counts * scale_factor),
    fill = "grey80", color = "black", alpha = 0.5
  )+
  geom_line(color = "#1E90FF", size = 1.2) +
  geom_point(color = "#1E90FF", size = 2.5, alpha = 0.9) +
  
  labs(
    title = "RSV-A Zurich",
    x = "Week"
  ) +
  
  scale_x_date(
    date_labels = "%b %Y",
    date_breaks = "1 month",
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    name =  "Absolute difference",
    sec.axis = sec_axis(~ . / scale_factor, name = "Number of clinical sequences")#,
    #limits = c(0, max_distance_zh)
  ) +
  
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6)
  )

print(p_zurich)




max_delta <- max_distance_ge # proportion/frequency goes 0-1
max_count <- max(delta_long_ge$total_counts, na.rm = TRUE)

scale_factor <-max_delta/ max_count
custom_labels <- c(
  "delta_abs_AD1" = "A.D.1",
  "delta_abs_AD2" = "A.D.2",
  "delta_abs_AD3" = "A.D.3",
  "delta_abs_AD5" = "A.D.5"
  
)
bars_df_ge <- delta_long_ge %>%
  dplyr::select(week, total_counts) %>%
  distinct()
p_geneva <- ggplot(delta_long_ge, aes(x = week, y = delta_abs)) + 
  facet_wrap(~ lineage, scales = "fixed", labeller = as_labeller(custom_labels))+
  geom_col(
    data = bars_df_ge,
    aes(x = week, y = total_counts * scale_factor),
    fill = "grey80", color = "black", alpha = 0.5
  )+
  geom_line(color = "#1E90FF", size = 1.2) +
  geom_point(color = "#1E90FF", size = 2.5, alpha = 0.9) +
  
  labs(
    title = "RSV-A Geneva",
    x = "Week") +
  
  scale_x_date(
    date_labels = "%b %Y",
    date_breaks = "1 month",
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    name =  "Absolute difference", 
    sec.axis = sec_axis(~ . / scale_factor, name = "Number of clinical sequences")#,
    #limits = c(0, max_distance_ge)
  ) +
  
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6)
  )

print(p_geneva)

summary_zh <- delta_long_zh %>%
  group_by(lineage) %>%
  summarise(mean_delta = mean(delta_abs, na.rm = TRUE), .groups = "drop")
print(summary_zh)

summary_ge <- delta_long_ge %>%
  group_by(lineage) %>%
  summarise(mean_delta = mean(delta_abs, na.rm = TRUE), .groups = "drop")
print(summary_ge)
