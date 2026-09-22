
#--------------------------------------------------------
#Title: Does Habitat Protection Status Matter? Insights From a Tropical Primate Community in Northeast Bangladesh 
# Authors: Atikul Islam Mithu, Sabit Hasan, Md. Sakhawat Hossain, Christian Roos and Tanvir Ahmed
# Corresponding Email: tahmed@dpz.eu (Tanvir Ahmed)
#--------------------------------------------------------


# 1. LOAD PACKAGES
library(dplyr)
library(ggplot2)
library(MASS)
library(glmmTMB)
library(DHARMa)
library(broom)
library(knitr)
library(tidyr)
library(forcats)
library(patchwork)
library(viridis)

#--------------------------------------------------------
# 2. DATA
raw <- read.csv(
  "Data.csv",
  stringsAsFactors = FALSE
)

raw <- raw %>%
  mutate(
    Survey_ID = paste(Transect.ID, Date, sep = "_")
  )

#--------------------------------------------------------
# 3. DATASET

survey_area <- raw %>%
  distinct(
    Survey_ID,
    Date,
    Transect.ID,
    PA.NPA,
    Transect.Length..km.
  )

#--------------------------------------------------------
# 4. COUNT NON-LORIS PRIMATE CLUSTERS
cluster_counts <- raw %>%
  filter(Species != "Bengal slow loris") %>%
  count(
    Survey_ID,
    Date,
    Transect.ID,
    PA.NPA,
    Transect.Length..km.,
    name = "Total_Clusters"
  )

#--------------------------------------------------------
# 5. COMBINE SURVEY FRAME WITH CLUSTER COUNTS
survey_cluster <- survey_area %>%
  left_join(
    cluster_counts,
    by = c(
      "Survey_ID",
      "Date",
      "Transect.ID",
      "PA.NPA",
      "Transect.Length..km."
    )
  ) %>%
  mutate(
    Total_Clusters = coalesce(Total_Clusters, 0L)
  ) %>%
  rename(
    Transect_ID = Transect.ID,
    Area_Type = PA.NPA,
    Length_km = Transect.Length..km.
  ) %>%
  mutate(
    Area_Type = factor(
      Area_Type,
      levels = c("NPA", "PA")
    ),
    Transect_ID = factor(Transect_ID),
    Survey_ID = factor(Survey_ID)
  )

#--------------------------------------------------------
# 6. CHECK
print(table(survey_cluster$Area_Type))
print(table(survey_cluster$Total_Clusters))

#--------------------------------------------------------
# 7. SURVEY EFFORT AND ENCOUNTER SUMMARY
summary_effort <- survey_cluster %>%
  group_by(Area_Type) %>%
  summarise(
    Surveys = n(),
    Unique_Transects = n_distinct(Transect_ID),
    Total_Transect_Length_km = sum(Length_km),
    Total_Clusters = sum(Total_Clusters),
    Mean_Clusters_per_Survey = mean(Total_Clusters),
    SD_Clusters_per_Survey = sd(Total_Clusters),
    .groups = "drop"
  )

knitr::kable(
  summary_effort,
  digits = 2,
  caption =
    ""
)

#--------------------------------------------------------
# 8. ENCOUNTER RATE
survey_cluster <- survey_cluster %>%
  mutate(
    Cluster_Rate = Total_Clusters / Length_km
  )

encounter_rate_summary <- survey_cluster %>%
  group_by(Area_Type) %>%
  summarise(
    N = n(),
    Mean_Rate = mean(Cluster_Rate),
    SD_Rate = sd(Cluster_Rate),
    Median_Rate = median(Cluster_Rate),
    Min_Rate = min(Cluster_Rate),
    Max_Rate = max(Cluster_Rate),
    .groups = "drop"
  )

knitr::kable(
  encounter_rate_summary,
  digits = 3,
  caption =
    ""
)


#--------------------------------------------------------
# 9. MIXED-EFFECTS MODEL
# Transect length is included as an offset.
# Transect_ID accounts for repeated visits to fixed transects.

primary_model <- glmmTMB(
  Total_Clusters ~
    Area_Type +
    offset(log(Length_km)) +
    (1 | Transect_ID),
  family = poisson,
  data = survey_cluster
)

summary(primary_model)

#--------------------------------------------------------
# 10. PA VS NPA INCIDENCE RATE RATIO

model_coef <- summary(primary_model)$coefficients$cond

pa_row <- which(
  rownames(model_coef) == "Area_TypePA"
)

primary_model_results <- data.frame(
  Term = "PA vs NPA",
  Estimate = model_coef[pa_row, "Estimate"],
  SE = model_coef[pa_row, "Std. Error"],
  Z = model_coef[pa_row, "z value"],
  P_value = model_coef[pa_row, "Pr(>|z|)"]
) %>%
  mutate(
    IRR = exp(Estimate),
    Lower_CI = exp(Estimate - 1.96 * SE),
    Upper_CI = exp(Estimate + 1.96 * SE)
  )

knitr::kable(
  primary_model_results,
  digits = 3,
  caption =
    "Mixed-effects model of primate cluster encounter rates."
)

#--------------------------------------------------------
# 11. DIAGNOSTICS

simulation <- DHARMa::simulateResiduals(
  primary_model,
  n = 1000
)

plot(simulation)

DHARMa::testDispersion(simulation)
DHARMa::testZeroInflation(simulation)
DHARMa::testUniformity(simulation)


#--------------------------------------------------------
# 12. DETECTION SENSITIVITY ANALYSIS
# Probability of detecting at least one non-loris cluster (secondary analysis).


survey_cluster <- survey_cluster %>%
  mutate(
    Detection = as.integer(Total_Clusters > 0)
  )

detection_model <- glmmTMB(
  Detection ~
    Area_Type +
    (1 | Transect_ID),
  family = binomial,
  data = survey_cluster
)

summary(detection_model)


#--------------------------------------------------------
# 13. DETECTION ODDS RATIO

detection_coef <- summary(
  detection_model
)$coefficients$cond

pa_detection_row <- which(
  rownames(detection_coef) == "Area_TypePA"
)

detection_results <- data.frame(
  Term = "PA vs NPA",
  Estimate = detection_coef[pa_detection_row, "Estimate"],
  SE = detection_coef[pa_detection_row, "Std. Error"],
  Z = detection_coef[pa_detection_row, "z value"],
  P_value = detection_coef[pa_detection_row, "Pr(>|z|)"]
) %>%
  mutate(
    OR = exp(Estimate),
    Lower_CI = exp(Estimate - 1.96 * SE),
    Upper_CI = exp(Estimate + 1.96 * SE)
  )

knitr::kable(
  detection_results,
  digits = 3,
  caption =
    "Mixed-effects logistic model of primate cluster detection probability."
)

#--------------------------------------------------------
# 14. TRANSECT-LEVEL DETECTION SUMMARY
detection_by_transect <- survey_cluster %>%
  group_by(Transect_ID) %>%
  summarise(
    N = n(),
    Detected = sum(Detection),
    No_Detection = sum(Detection == 0),
    Detection_Rate = mean(Detection),
    .groups = "drop"
  )

knitr::kable(
  detection_by_transect,
  digits = 3,
  caption =
    "Transect-level variation in primate cluster detection."
)


#--------------------------------------------------------
# 15. FIGURE S1

ggplot(
  survey_cluster,
  aes(
    x = Area_Type,
    y = Cluster_Rate
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    alpha = 0.6
  ) +
  labs(
    x = "Land-use type",
    y = "Primate cluster encounter rate (clusters/km)"
  ) +
  theme_classic(base_size = 13)

#--------------------------------------------------------
# 16. BENGAL SLOW LORIS ONLY

loris_data <- raw %>%
  filter(Species == "Bengal slow loris")

survey_loris <- loris_data %>%
  group_by(
    Survey_ID,
    Date,
    Transect.ID,
    PA.NPA,
    Transect.Length..km.
  ) %>%
  summarise(
    Presence = 1,
    .groups = "drop"
  ) %>%
  rename(
    Transect_ID = Transect.ID,
    Area_Type = PA.NPA,
    Length_km = Transect.Length..km.
  )

loris_summary <- survey_loris %>%
  group_by(Area_Type) %>%
  summarise(
    Surveys_with_Detections = n(),
    .groups = "drop"
  )

knitr::kable(
  loris_summary,
  caption =
)

#--------------------------------------------------------
# 17. SPECIES-WISE CLUSTER SIZE

cluster_data <- raw %>%
  filter(
    Species != "Bengal slow loris",
    Cluster.Size > 1
  ) %>%
  rename(PA_NPA = PA.NPA) %>%
  mutate(
    Species = factor(Species),
    PA_NPA = factor(PA_NPA),
    PA_NPA = relevel(PA_NPA, ref = "NPA")
  )

results <- list()

for (sp in levels(cluster_data$Species)) {
  
  sub <- cluster_data %>%
    filter(Species == sp)
  
  if (n_distinct(sub$PA_NPA) < 2) {
    
    results[[sp]] <- data.frame(
      Species = sp,
      Beta = NA,
      SE = NA,
      Z_value = NA,
      P_value = NA,
      Note = "Only one land-use type"
    )
    
    next
  }
  
  m <- tryCatch(
    glm.nb(
      Cluster.Size ~ PA_NPA,
      data = sub,
      control = glm.control(maxit = 50)
    ),
    error = function(e) NULL
  )
  
  if (is.null(m)) {
    
    results[[sp]] <- data.frame(
      Species = sp,
      Beta = NA,
      SE = NA,
      Z_value = NA,
      P_value = NA,
      Note = "Model failed to converge"
    )
    
    next
  }
  
  coef_m <- summary(m)$coefficients
  pa_row <- grep("^PA_NPA", rownames(coef_m))
  
  if (length(pa_row) == 0) {
    
    results[[sp]] <- data.frame(
      Species = sp,
      Beta = NA,
      SE = NA,
      Z_value = NA,
      P_value = NA,
      Note = "PA effect not estimable"
    )
    
  } else {
    
    results[[sp]] <- data.frame(
      Species = sp,
      Beta = round(coef_m[pa_row, "Estimate"], 3),
      SE = round(coef_m[pa_row, "Std. Error"], 3),
      Z_value = round(coef_m[pa_row, "z value"], 3),
      P_value = round(coef_m[pa_row, "Pr(>|z|)"], 3),
      Note = "Tested"
    )
  }
}

cluster_size_results <- bind_rows(results)

knitr::kable(
  cluster_size_results,
  caption =""
    
)

#--------------------------------------------------------
# 18. HABITAT USE

habitat_data <- data.frame(
  Species = c(
    "Western hoolock gibbon",
    "Capped langur",
    "Phayre's langur",
    "Rhesus macaque",
    "Northern pig-tailed macaque",
    "Bengal slow loris"
  ),
  Less.disturbed.forest = c(100, 33.33, 25, 6.25, 0, 10),
  Disturbed.forest = c(0, 9.52, 12.5, 25, 50, 10),
  Forest.edge = c(0, 42.86, 62.5, 43.75, 50, 45),
  Plantations = c(0, 14.29, 0, 25, 0, 15),
  Human.settlement = c(0, 0, 0, 0, 0, 20)
)

habitat_long <- habitat_data %>%
  pivot_longer(
    -Species,
    names_to = "Habitat",
    values_to = "Percentage"
  ) %>%
  mutate(
    Habitat = gsub("\\.", " ", Habitat),
    Species = factor(
      Species,
      levels = c(
        "Western hoolock gibbon",
        "Capped langur",
        "Phayre's langur",
        "Rhesus macaque",
        "Northern pig-tailed macaque",
        "Bengal slow loris"
      )
    ),
    Habitat = factor(
      Habitat,
      levels = c(
        "Less disturbed forest",
        "Disturbed forest",
        "Forest edge",
        "Plantations",
        "Human settlement"
      )
    )
  )

# Habitat-use plot with species on the Y-axis
habitat_plot <- ggplot(
  habitat_long,
  aes(x = Percentage, y = Species, fill = Habitat)
) +
  geom_col() +
  scale_fill_viridis_d(name = "Habitat type") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  ) +
  labs(
    x = "Habitat use (%)",
    y = NULL
  )

habitat_plot

#Figure 3

ggsave(
  "Habitat use.pdf",
  habitat_plot,
  width = 9,
  height = 5,
  dpi = 600
)

#--------------------------------------------------------
# 19. DENSITY ESTIMATION

species_order <- c(
  "Western hoolock gibbon",
  "Capped langur",
  "Phayre’s langur",
  "Rhesus macaque",
  "Northern pig-tailed macaque",
  "Bengal slow loris"
)

density_data <- data.frame(
  Species = species_order,
  PA_clusters = c(
    0.23, 1.14, 0.65, 1.53, NA, NA
  ),
  NonPA_clusters = c(
    NA, 2.69, 1.16, 1.32, 0.46, NA
  ),
  PA_individuals = c(
    0.47, 8.89, 9.93, 42.13, NA, 1.61
  ),
  NonPA_individuals = c(
    NA, 21.3, 13.01, 36.04, 15.7, 1.01
  )
)


# CLUSTER
clusters <- density_data %>%
  dplyr::select(
    Species,
    PA_clusters,
    NonPA_clusters
  ) %>%
  pivot_longer(
    cols = -Species,
    names_to = "Area",
    values_to = "Density"
  ) %>%
  filter(!is.na(Density)) %>%
  mutate(
    Area = recode(
      Area,
      PA_clusters = "Protected Area",
      NonPA_clusters = "Non-Protected Area"
    ),
    Density = ifelse(
      Area == "Non-Protected Area",
      -Density,
      Density
    ),
    Species = factor(
      Species,
      levels = rev(species_order)
    )
  )

p_clusters <- ggplot(
  clusters,
  aes(Species, Density, fill = Area)
) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous(
    limits = c(-3, 3),
    breaks = seq(-3, 3, 1),
    labels = abs
  ) +
  scale_fill_viridis_d(
    option = "viridis",
    begin = 0.2,
    end = 0.8
  ) +
  labs(
    y = "Clusters/km²",
    x = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top"
  )

p_clusters

# INDIVIDUAL DENSITY

individuals <- density_data %>%
  dplyr::select(
    Species,
    PA_individuals,
    NonPA_individuals
  ) %>%
  pivot_longer(
    cols = -Species,
    names_to = "Area",
    values_to = "Density"
  ) %>%
  filter(!is.na(Density)) %>%
  mutate(
    Area = recode(
      Area,
      PA_individuals = "Protected Area",
      NonPA_individuals = "Non-Protected Area"
    ),
    Density = ifelse(
      Area == "Non-Protected Area",
      -Density,
      Density
    ),
    Species = factor(
      Species,
      levels = rev(species_order)
    )
  )

p_individuals <- ggplot(
  individuals,
  aes(Species, Density, fill = Area)
) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous(
    limits = c(-45, 45),
    breaks = seq(-45, 45, 15),
    labels = abs
  ) +
  scale_fill_viridis_d(
    option = "viridis",
    begin = 0.2,
    end = 0.8
  ) +
  labs(
    y = "Individuals/km²",
    x = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none"
  )

individuals
#--------------------------------------------------------
# 20. COMBINE

Density <- p_clusters / p_individuals
Density

#Figure 4
ggsave("Density.pdf", Density, width=7.5, height=6, dpi=600)
#--------------------------------------------------------
#END