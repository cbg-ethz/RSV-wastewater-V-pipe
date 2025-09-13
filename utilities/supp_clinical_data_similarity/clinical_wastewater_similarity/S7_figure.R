library(stringr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(lubridate)
library(RColorBrewer)
library(purrr)
library(binom)
library(DescTools)


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



# Summarising wastewater frequencies weekly
wastewater_data <- wastewater_data %>%
  mutate(week = floor_date(date, unit = "week")) %>%
  dplyr::group_by(location, variant, week) %>%
  #  summarise weekly (separately per location and per variant)
  summarise(mean_freq = weighted.mean(
    proportion,
    weight = 1/(((proportionLower - proportionUpper)/(2*1.96))^2), # inverse variance weighting
    .groups = "drop"
  )) %>%
  dplyr::group_by(location, week) %>%
  mutate(variant_freq = mean_freq/sum(mean_freq)) %>%
  ungroup()


# aggregate sublineages into lineages
wastewater_data$variant <- gsub("^A\\.D\\.1.*", "A.D.1", wastewater_data$variant)
wastewater_data$variant <- gsub("^A\\.D\\.2.*", "A.D.2", wastewater_data$variant)
wastewater_data$variant <- gsub("^A\\.D\\.3.*", "A.D.3", wastewater_data$variant)
wastewater_data$variant <- gsub("^A\\.D\\.5.*", "A.D.5", wastewater_data$variant)
wastewater_data <- wastewater_data %>%
  group_by(variant, week, location) %>%
  summarise(proportion = sum(mean_freq, na.rm = TRUE), .groups = "drop")


# Filter clinical data (filter by calendorical weeks (starting on Monday))
clinical_data <- clinical_data %>%
  filter(sampleCollectionDate > "2023-10-29",
         sampleCollectionDate < "2024-03-04",
         geoLocCountry %in% european_countries) ## available European countries


# Select the samples of only those lineages that were detected in Switzerland
swiss_sub_lineages <- c("A.D.1", "A.D.1.5", "A.D.1.6", "A.D.2.1",
                        "A.D.3", "A.D.3.1", "A.D.5.1", "A.D.5.2")
swiss_lineages <- c("A.D.1","A.D.2", "A.D.3", "A.D.5")


clinical_aggregated <- clinical_data %>%
  filter(lineage %in% swiss_sub_lineages)

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


# Compute multinomial confidence intervals for each lineage
clinical_aggregated_ci <- clinical_aggregated %>%
  dplyr::select(week, variant, lineage_counts, total_counts) %>%
  pivot_wider(
    names_from = variant,
    values_from = lineage_counts
  )
clinical_aggregated_ci <- clinical_aggregated_ci[c(1, 3, 4, 5, 6, 2)]
# Function to compute multinomial CIs for each row
compute_multinom_ci <- function(v1, v2, v3, v4, v_total) {
  counts <- c(v1, v2, v3, v4)
  ci <- MultinomCI(counts, conf.level = 0.95, method = "sisonglaz")  # sisonglaz method
  return(data.frame(
    A.D.1_lower = ci[1, "lwr.ci"],
    A.D.1_upper = ci[1, "upr.ci"],
    A.D.2_lower = ci[2, "lwr.ci"],
    A.D.2_upper = ci[2, "upr.ci"],
    A.D.3_lower = ci[3, "lwr.ci"],
    A.D.3_upper = ci[3, "upr.ci"],
    A.D.5_lower = ci[4, "lwr.ci"],
    A.D.5_upper = ci[4, "upr.ci"]
  ))
}

# Apply row-wise using apply or dplyr
clinical_aggregated_ci <- clinical_aggregated_ci %>%
  mutate(ci_results = pmap(list(A.D.1, A.D.2, A.D.3, A.D.5, total_counts),
                           ~compute_multinom_ci(..1, ..2, ..3, ..4))) %>%
  unnest(ci_results)


clinical_aggregated_ci_pivot <- clinical_aggregated_ci %>%
  pivot_longer(
    cols = matches("^A\\.D\\.[0-9]+(_lower|_upper)?$"),
    names_to = c("variant", "ci_type"),
    names_pattern = "(A\\.D\\.[0-9]+)(?:_(lower|upper))?",
    values_to = "value"
  )
clinical_aggregated_ci_pivot <- clinical_aggregated_ci_pivot %>%
  pivot_wider(
    names_from = ci_type,
    values_from = value,
    names_prefix = "ci_"
  ) %>%
  rename(counts = ci_)  # because NA becomes "ci_NA" in names

clinical_aggregated_merged <- merge(clinical_aggregated, clinical_aggregated_ci_pivot, by = c("week", "variant", "total_counts"))
clinical_aggregated <- clinical_aggregated_merged

# -------------------- Plotting relative abundances ---------------------------

clinical_aggregated <- clinical_aggregated %>%
  mutate(type = "Clinical")  # for points
# Define colors
wastewater_colors <- c("#66C2A5", "#FC8D62")# brewer.pal(n = length(unique(wastewater_data$location)), "Set2")
clinical_color <- "blue"



# Max values for scaling
max_freq <- 1.01  # proportion/frequency goes 0-1
max_count <- max(clinical_aggregated$total_counts)

# Scale factor
scale_factor <- max_freq / max_count
# Change Swiss -> to English names

wastewater_data$location <- gsub("Zürich", "Zurich", wastewater_data$location)
wastewater_data$location <- gsub("Genève", "Geneva", wastewater_data$location)


p <- ggplot() +
  # Clinical sample counts as bars (scaled)
  geom_col(data = clinical_aggregated,
           aes(x = week, y =  total_counts* scale_factor),
           fill = "grey80", alpha = 0.6) +
  
  # Wastewater lines
  geom_line(data = wastewater_data,
            aes(x = week, y = proportion, color = location, group = location),
            size = 1) +
  
  # Clinical ribbons
  geom_ribbon(data = clinical_aggregated,
              aes(x = week, ymin = ci_lower, ymax = ci_upper),
              fill = "blue",alpha = 0.2, color = NA) +
  # Clinical points
  geom_point(data = clinical_aggregated,
             aes(x = week, y = frequency_clin),
             color = clinical_color,
             size = 0.8) +

  # Clinical line
  geom_line(data = clinical_aggregated,
            aes(x = week, y = frequency_clin, color=type),
            color = clinical_color,
            size = 0.8) +
  
  # Facet by variant
  facet_wrap(~ variant) +
  
  # Scales
  scale_color_manual(values = wastewater_colors) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  
  # Y-axis with secondary axis
  scale_y_continuous(
    name = "Relative abundances",
    limits = c(0, max_freq),
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "Number of clinical sequences")
  ) +
  scale_color_manual(
    name = "",
    values = c(setNames(wastewater_colors, c("Zurich (ZH)", "Geneva (GE)")),
               "Clinical" = clinical_color)
  ) +
  
  labs(x = "Sample Collection Date",
       color = "Wastewater Location",
       fill = "Wastewater Location") +

theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 17))

print(p)


rm(list=ls())
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# RSV-B data analysis
european_countries <- c(
  "Albania", "Andorra", "Austria", "Belarus", "Belgium", "Bosnia and Herzegovina",
  "Bulgaria", "Croatia", "Czech Republic", "Denmark", "Estonia", "Finland",
  "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy",
  "Latvia", "Liechtenstein", "Lithuania", "Luxembourg", "Malta", "Moldova",
  "Monaco", "Montenegro", "Netherlands", "North Macedonia", "Norway", "Poland",
  "Portugal", "Romania", "Russia", "San Marino", "Serbia", "Slovakia", "Slovenia",
  "Spain", "Sweden", "Switzerland", "Ukraine", "United Kingdom", "Vatican City"
)

clinical_data <- read.csv(
  'metadata/rsv-b_metadata_2025-08-17T2158.tsv',
  sep = '\t'
)[, c('lineage', 'sampleCollectionDate', 'geoLocCountry')]
clinical_data$counts <- 1
clinical_data$sampleCollectionDate <- as.Date(clinical_data$sampleCollectionDate)

wastewater_data <- read.csv(
  '../../../preprint/data/lollipop/rsvb/deconvolved.csv',
  sep = '\t'
)
wastewater_data$date <- as.Date(wastewater_data$date)



# Summarising wastewater frequencies weekly
wastewater_data <- wastewater_data %>%
  mutate(week = floor_date(date, unit = "week")) %>%
  dplyr::group_by(location, variant, week) %>%
  #  summarise weekly (separately per location and per variant)
  summarise(mean_freq = weighted.mean(
    proportion,
    weight = 1/(((proportionLower - proportionUpper)/(2*1.96))^2), # inverse variance weighting
    .groups = "drop"
  )) %>%
  dplyr::group_by(location, week) %>%
  mutate(variant_freq = mean_freq/sum(mean_freq)) %>%
  ungroup()

# Select the samples of only those lineages that were detected in nonzero frequencies in Swiss wastewater
swiss_sub_lineages <- c("B.D.4.1.1", "B.D.E.1", "B.D.E.4", "undetermined")
# 

wastewater_data <- wastewater_data %>%
  group_by(variant, week, location) %>%
  summarise(proportion = sum(mean_freq, na.rm = TRUE), .groups = "drop") %>%
  filter(variant %in% swiss_sub_lineages)


# Filter clinical data (filter by calendorical weeks (starting on Monday))
clinical_data <- clinical_data %>%
  filter(sampleCollectionDate > "2022-11-06",
         sampleCollectionDate < "2023-02-06",
         geoLocCountry %in% european_countries) ## available European countries




clinical_aggregated <- clinical_data %>%
  filter(lineage %in% swiss_sub_lineages)

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



# Compute multinomial confidence intervals for each lineage
clinical_aggregated_ci <- clinical_aggregated %>%
  dplyr::select(week, variant, lineage_counts, total_counts) %>%
  pivot_wider(
    names_from = variant,
    values_from = lineage_counts
  )
clinical_aggregated_ci <- clinical_aggregated_ci[c(1, 3, 4, 5, 2)]
# Function to compute multinomial CIs for each row
compute_multinom_ci <- function(v1, v2, v3, v_total) {
  counts <- c(v1, v2, v3)
  ci <- MultinomCI(counts, conf.level = 0.95, method = "sisonglaz")  # sisonglaz method
  return(data.frame(
    B.D.4.1.1_lower = ci[1, "lwr.ci"],
    B.D.4.1.1_upper = ci[1, "upr.ci"],
    B.D.E.1_lower = ci[2, "lwr.ci"],
    B.D.E.1_upper = ci[2, "upr.ci"],
    B.D.E.4_lower = ci[3, "lwr.ci"],
    B.D.E.4_upper = ci[3, "upr.ci"]
  ))
}

# Apply row-wise using apply or dplyr
clinical_aggregated_ci <- clinical_aggregated_ci %>%
  mutate(ci_results = pmap(list(B.D.4.1.1, B.D.E.1, B.D.E.4, total_counts),
                           ~compute_multinom_ci(..1, ..2, ..3, ..4))) %>%
  unnest(ci_results)


clinical_aggregated_ci_pivot <- clinical_aggregated_ci %>%
  pivot_longer(
    cols = matches("^B\\.D\\..+(_lower|_upper)?$"),
    names_to = c("variant", "ci_type"),
    names_pattern = "^(B\\.D\\..+?)(?:_(lower|upper))?$",
    values_to = "value"
  )
clinical_aggregated_ci_pivot <- clinical_aggregated_ci_pivot %>%
  pivot_wider(
    names_from = ci_type,
    values_from = value,
    names_prefix = "ci_"
  ) %>%
  rename(counts = ci_)  # because NA becomes "ci_NA" in names

clinical_aggregated_merged <- merge(clinical_aggregated, clinical_aggregated_ci_pivot, by = c("week", "variant", "total_counts"))
clinical_aggregated <- clinical_aggregated_merged


# # -----------------------------------------------------------------------------

# -------------------- Plotting relative abundances ---------------------------


clinical_aggregated <- clinical_aggregated %>%
  mutate(type = "Clinical")  # for points
# Define colors
wastewater_colors <- wastewater_colors <- c("#66C2A5", "#FC8D62")
clinical_color <- "blue"



# Max values for scaling
max_freq <- 1.01 #  # proportion/frequency goes 0-1
max_count <- max(clinical_aggregated$total_counts)

# Scale factor
scale_factor <- max_freq / max_count
# Change Swiss -> to English names

wastewater_data$location <- gsub("Zürich", "Zurich", wastewater_data$location)
wastewater_data$location <- gsub("Genève", "Geneva", wastewater_data$location)


p <- ggplot() +
  # Clinical sample counts as bars (scaled)
  geom_col(data = clinical_aggregated,
           aes(x = week, y =  total_counts* scale_factor),
           fill = "grey80", alpha = 0.6) +
  
  # Wastewater lines
  geom_line(data = wastewater_data,
            aes(x = week, y = proportion, color = location, group = location),
            size = 1) +
  
  # Clinical ribbons
  geom_ribbon(data = clinical_aggregated,
              aes(x = week, ymin = ci_lower, ymax = ci_upper),
              fill = "blue", alpha = 0.2, color = NA) +
  # Clinical points
  geom_point(data = clinical_aggregated,
             aes(x = week, y = frequency_clin),
             color = clinical_color,
             size = 0.8) +
  
  # Clinical line
  geom_line(data = clinical_aggregated,
            aes(x = week, y = frequency_clin, color=type),
            color = clinical_color,
            size = 0.8) +
  
  # Facet by variant
  facet_wrap(~ variant) +
  
  # Scales
  scale_color_manual(values = wastewater_colors) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  # Y-axis with secondary axis
  scale_y_continuous(
    name = "Relative abundances",
    limits = c(0, max_freq),
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "Number of clinical sequences")
  ) +
  scale_color_manual(
    name = "",
    values = c(setNames(wastewater_colors, c("Zurich (ZH)", "Geneva (GE)")),
               "Clinical" = clinical_color)
  ) +
  
  labs(x = "Sample Collection Date",
       color = "Wastewater Location",
       fill = "Wastewater Location") +

  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 17))

print(p)

