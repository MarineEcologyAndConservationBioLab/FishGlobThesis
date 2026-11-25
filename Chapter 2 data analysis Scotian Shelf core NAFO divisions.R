#####################################
######## Initialize FISHGLOB ########
#####################################
library(readr)
FishGlob_public_clean <- read_csv("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/FishGlob_data-v2.0.1/AquaAuma-FishGlob_data-508ea5f/outputs/Compiled_data/FishGlob_public_clean.csv/FishGlob_public_clean.csv")

#########################################################
######## Initialize SCS #################################
#########################################################
scotian_shelf_year_round <- FishGlob_public_clean %>%
  filter(survey =="SCS")
scotian_shelf_summer <- scotian_shelf_year_round %>%
  filter(quarter == 3)

scotian_shelf_summer$body_weight <- scotian_shelf_summer$wgt / scotian_shelf_summer$num
scotian_shelf_summer <- bin_coordinates_fn(scotian_shelf_summer)

# add NAFO divisions
nafo_sf <- sf::st_read("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/DFO advice on SCS subregioning/NAFO divisions/Divisions.shp")
scotian_shelf_summer <- assign_NAFO_fn(scotian_shelf_summer, nafo_sf)
scotian_shelf_summer <- scotian_shelf_summer %>%
  select(-geometry)

# add Federal Marine Bioregions
bioregions_sf <- sf::st_read("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/DFO advice on SCS subregioning/Canadian FederalMarineBioregions_SHP/FederalMarineBioregions_SHP/FederalMarineBioregions.shp")

############################################################
##### Main df for Chapter 1 scotian_shelf_summer_core ######
############################################################
scotian_shelf_summer_core <- scotian_shelf_summer %>%
  filter(Nafo_division %in% c("4Vn", "4Vs", "4W","4X"))

############################################################
######## check strata coverage #############################
############################################################
Strata_areas <- read_csv("Strata.csv")
Strata_areas <- Strata_areas %>%
  select(-`...1`)%>%
  rename(Stratum = STRAT, `Area (km2)` = AREA)

SCS_strata_with_areas <- data.frame(
  Stratum = unique(scotian_shelf_summer_core$stratum))

SCS_strata_with_areas <- left_join(SCS_strata_with_areas, Strata_areas, by = "Stratum")


num_total_catches <- length(scotian_shelf_summer_core$survey)
num_total_taxa <- length(unique(scotian_shelf_summer_core$accepted_name))
num_total_hauls <- length(unique(scotian_shelf_summer_core$haul_id))

#####################################################
############ Section 1: Logistics ###################
#####################################################

########### Mapping of SCS data #####################

coordinates <- data.frame(longitude = scotian_shelf_summer$longitude, latitude = scotian_shelf_summer$latitude)

# Get the map data for Canada and the United States
na_map <- map_data("world") %>% 
  filter(region %in% c("Canada", "USA", "Saint Pierre and Miquelon"))

na_map <- na_map %>%
  mutate(region = recode(region, "Saint Pierre and Miquelon" = "France"))

# Define the zoom bounding box around SCS
zoom_bbox <- c(left = -71.0, bottom = 39.0, right = -55.0, top = 50.0)

# 1. CONVERT COORDINATE SYSTEMS (Crucial Step)
# This moves your bioregions from "meters" to "decimal degrees"
bioregions_sf <- st_transform(bioregions_sf, 4326)
nafo_sf <- st_transform(nafo_sf, 4326)

# 2. THE PLOT
SCS_map <- ggplot() +
  # Landmass (WGS84)
  geom_polygon(data = na_map, aes(x = long, y = lat, group = group, fill = region), color = "black") +
  scale_fill_manual(values = c("Canada" = "lightcoral", "USA" = "lightblue", "France" = "darkblue"), name = "Country") +
  
  new_scale_fill() +
  
  # Bioregions (Now converted to WGS84)
  geom_sf(data = bioregions_sf, aes(fill = NAME_E), color = "black", alpha = 0.6) +
  scale_fill_manual(
    values = c(
      "Gulf of Saint Lawrence" = "yellow",
      "Scotian Shelf" = "lightgreen",
      "Saint-Pierre et Miquelon" = "chocolate4",
      "Newfoundland-Labrador Shelves" = "orange"
    ),
    name = "Federal Marine Bioregion"
  ) +
  
  # Survey Points
  geom_point(data = coordinates, aes(x = longitude, y = latitude, color = "DFO R/V Survey Hauls")) +
  scale_color_manual(values = c("DFO R/V Survey Hauls" = "dimgrey"), name = "Data") +
  
  # NAFO Zones and Labels
  geom_sf(data = nafo_sf, aes(linetype = "NAFO Divisions"), fill = NA, color = "black", linewidth = 1.5) + 
  scale_linetype_manual(values = c("NAFO Divisions" = "dashed"), name = "Boundaries") +
  geom_sf_text(data = nafo_sf, aes(label = ZONE), color = "black", size = 7) +
  
  # Crop to your zoom window
  coord_sf(xlim = c(zoom_bbox["left"], zoom_bbox["right"]),
           ylim = c(zoom_bbox["bottom"], zoom_bbox["top"]),
           expand = FALSE) + # expand=FALSE prevents extra padding
  theme_minimal() +
  theme(
    axis.text = element_text(size = 24, color = "black"),    # Sizes for the numbers (Long/Lat)
    axis.title = element_text(size = 36, face = "bold"),    # Sizes for the "Longitude/Latitude" labels
    legend.title = element_text(size = 36),                 # Optional: makes legend headers larger
    legend.text = element_text(size = 24)                   # Optional: makes legend items larger
  ) +
  labs(x = "Longitude", y = "Latitude")
print(SCS_map)

###### Summary statistics for SCS Survey Region #########

SCS_bio_summary <- summarize_biological(scotian_shelf_summer_core) # use barplot from here for survey length biodiversity by rank
# returns df with list and count of num taxa present in survey region at each taxonomic rank and plots

SCS_operational_summary <- summarize_operational(scotian_shelf_summer_core)
SCS_haul_prevalence <- haul_prevalence(scotian_shelf_summer_core, SCS_bio_summary)
SCS_haul_prevalence$haul_prevalence_list <- as.numeric(SCS_haul_prevalence$haul_prevalence_list)
SCS_change_in_total_catch <- total_catch_change(scotian_shelf_summer_core)
SCS_common_change <- most_common_change(scotian_shelf_summer_core)
#lat_lon_frequency(scotian_shelf_summer) # generate histograms of haul occurrence by latitude and longitude from SCS data


############ Gears used over time #######################
scotian_shelf_western <- scotian_shelf_summer_core %>%
  filter(gear == "Western IIA trawl")

'scotian_shelf_bridle <- scotian_shelf_summer %>%
  filter(gear == "US 4 seam 3 bridle survey trawl")
not used in summertime!
'

scotian_shelf_yankee <- scotian_shelf_summer_core %>%
  filter(gear == "Yankee #36 otter trawl")

phrase1 <- paste("The Western IIA trawl gear was used between", min(scotian_shelf_western$year), "and", max(scotian_shelf_western$year))
#phrase2 <- paste("US 4 seam 3 bridle survey trawl", min(scotian_shelf_bridle$year), "and", max(scotian_shelf_bridle$year))
phrase3 <- paste("Yankee #36 otter trawl", min(scotian_shelf_yankee$year), "and", max(scotian_shelf_yankee$year))
print(paste(phrase1, phrase3, sep = ". "))

########### plot area surveyed and gear used ########### 
# Aggregate total area_swept by year
total_area_swept <- aggregate(area_swept ~ year, data = scotian_shelf_summer_core, sum)

# Cutoff year
cutoff_x <- 1981

# split swept area into two dfs to plot shaded regions separately
df_before <- total_area_swept %>% filter(year <= cutoff_x)
df_after <- total_area_swept %>% filter(year >= cutoff_x)

ggplot(total_area_swept, aes(x = year, y = area_swept)) +
  geom_area(data = df_before, fill = "darkgoldenrod", alpha = 0.6) +
  geom_area(data = df_after, fill = "blue", alpha = 0.6) +
  geom_line() +
  geom_point() +
  annotate("text", x = 1975.5, y = 50, label = "Western IIA Trawl Gear", size = 5, fontface = "bold") +
  annotate("text", x = 2000, y = 85, label = "Yankee #36 Otter Trawl Gear", size = 5, fontface = "bold") +
  labs(x = "Year", y = expression(paste("Total area swept (", km^2, ")"))) + 
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )

#### Generate data subset at species rank# ####
scotian_shelf_summer_core_by_species <- scotian_shelf_summer_core %>%
  filter(rank == "Species")

####################################################################################################
##### Make rarefaction curves for haul sampling effort and number of unique taxa observed ##########
####################################################################################################
# try rarefaction by area
SCS_summer_rarefied_by_area <- scotian_shelf_summer_core %>%
  group_by(year) %>%
  summarise(
    total_area_swept = sum(area_swept),
    unique_accepted_names = n_distinct(accepted_name),
    .groups = 'drop' # Remove grouping after summarization
  )

rarefaction_model <- lm(unique_accepted_names ~ total_area_swept, data = SCS_summer_rarefied_by_area)
r_squared <- summary(rarefaction_model)$r.squared

# Create the plot using ggplot2
ggplot(SCS_summer_rarefied_by_area, aes(x = total_area_swept, y = unique_accepted_names)) +
  # Add the points (scatter plot)
  geom_point() +
  
  # Add the linear trendline. The `method = "lm"` specifies a linear model.
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  
  # Set the plot limits and labels
  scale_x_continuous(name = "Total area swept (km^2)") +
  scale_y_continuous(name = "Number of uniquely identified taxa") +
  ggtitle("Species accumulation curve of Scotian Shelf") +
  
  # Add the R-squared value as text on the plot. Adjust the coordinates (x, y)
  # to position the text where you want it.
  annotate("text", x = 250, y = 50, label = paste("R^2 =", round(r_squared, 3))) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
  )

# note lack of sampling plateau, suggesting that number of hauls per year increases number of possible taxa identified 
# bad, rarefaction by areal sampling effort even more consistent with more sampling means more species

#########################################
###### Section 1: Stability #############
#########################################

######## attempt community matrix investigation following Britten et al. 2014 ------------------------------------
######## begin by finding the trophic level of all sampled species -----------------------------------------------

species_list <- unique(scotian_shelf_summer_core_by_species$accepted_name) 
species_list <- validate_names(species_list)

species_info <- species(species_list)

diet_table <- diet(species_list) # diet breakdown of SCS summer species including by lifestage

trophic_key <- data.frame(Species = diet_table$Species,Troph = diet_table$Troph, seTroph = diet_table$seTroph,  Remark = diet_table$Remark, Source = "Diet table") # trophic breakdown of SCS summer species with only relevant columns kept from diet info

hist(trophic_key$Troph) # note presence fo sub 2 values for troph which doesnt make sense, so remove these records
# Define the trophic threshold for demersal fish
trophic_threshold <- 2.0
# Remove rows where any value in a row is below the threshold
trophic_key <- trophic_key[trophic_key[, "Troph"] >= trophic_threshold, ] 

hist(trophic_key$Troph) # note loss of sub 2 trophic levels

# remove any rows with missing values
trophic_key<- trophic_key %>% filter(!is.na(Troph))

# add trophic estimates from ecology table
eco_table <- ecology(species_list)
eco_trophic_key <- data.frame(Species = eco_table$Species,Troph = eco_table$FoodTroph, seTroph = eco_table$FoodSeTroph,  Remark = eco_table$FoodRemark, Source = "Ecology table, food estimate") # trophic breakdown of SCS summer species with only relevant columns kept from diet info
# remove any rows with missing values
eco_trophic_key <- eco_trophic_key %>% filter(!is.na(Troph))
trophic_key <- rbind(trophic_key, eco_trophic_key)

eco_trophic_key_diet <- data.frame(Species = eco_table$Species,Troph = eco_table$DietTroph, seTroph = eco_table$DietSeTroph,  Remark = eco_table$DietRemark, Source = "Ecology table, diet estimate") # trophic breakdown of SCS summer species with only relevant columns kept from diet info
eco_trophic_key_diet <- eco_trophic_key_diet %>% filter(!is.na(Troph))
eco_trophic_key_diet <- eco_trophic_key_diet %>% filter(!is.na(seTroph))
trophic_key <- rbind(trophic_key, eco_trophic_key_diet) # now trophic key incorporates trophic estimates from 3 distinct sources
# remove 0 value SeTroph to allow for weighted average of trophic level
trophic_key <- trophic_key %>% filter(seTroph != 0)

test <- weighted_average_troph("Alepisaurus ferox", trophic_key)

# now generate weighted average trophic estimates using troph and seTroph from trophic_key for all SCS species with good trophic estimates

trophic_species_list <- as.list(unique(trophic_key$Species))
avg_trophic_key <- data.frame(Species = list(), WE_trophic_estimate = list(), SE_WE_trophic_estimate = list())
for (species in trophic_species_list){
  new_we <- weighted_average_troph(species, trophic_key)
  avg_trophic_key <- rbind(avg_trophic_key, new_we)
}
dev.new()
T_level_hist <- ggplot(avg_trophic_key, aes(x = WE_trophic_estimate)) +
  geom_histogram(fill = "orange", color = "black", binwidth = 0.1) +  # adjust binwidth as needed
  labs(
    title = "Distribution of weighted average trophic estimates",
    x = "Weighted average trophic level estimate",
    y = "Frequency (number of species)"
  ) +
  scale_x_continuous(limits = c(2.5, 4.75)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
  )
print(T_level_hist)

# attempt to assign trophic levels to species with good trophic information
keep_values <- avg_trophic_key$Species

# generate dictionary to find trophic levels given species key
values <- avg_trophic_key$WE_trophic_estimate
keys <- keep_values
troph_dict <- setNames(values, keys) 

# generate another dictionary to find SE trophic level given species key
values <- avg_trophic_key$SE_WE_trophic_estimate
keys <- keep_values
SE_troph_dict <- setNames(values, keys)

# Filter the scotian shelf year round dataframe to analyze only species with trophic information
troph_filtered_SCS <- scotian_shelf_summer_core_by_species %>% # generates new df of SCS year round only including species whose trophic information could be found
  filter(accepted_name %in% keep_values)

scotian_with_troph <- assign_SEtroph_fn(troph_filtered_SCS, SE_troph_dict) # all species get trophic level estimate
scotian_with_troph <- assign_troph_fn(scotian_with_troph, troph_dict)

#### check coverage of trophic assignment and difference between possible and caught trophic distribution #########

NA_SCS_with_troph <-  assign_troph_fn(scotian_shelf_summer_core_by_species, troph_dict) # has empty trophic levels
num_missing_troph_records <- length(scotian_shelf_summer_core_by_species$survey) - length(scotian_with_troph$survey) # reflects number of species rank observations we could not assign trophic level

missing_troph_species <- species_list[!species_list %in% keep_values]
num_missing_troph_species <- length(missing_troph_species)

print(c(num_missing_troph_records, num_missing_troph_species))

hist(troph_dict, main = "Distribution of possible trophic level estimates in SCS", xlab = "Trophic level estimate", ylab = "Frequency (number of unique species)")
boxplot(troph_dict, main = "Boxplot of possible trophic levels in SCS Community", xlab = "Note lower average than catch average", ylab = "Weighted average trophic level estimate", fill = "orange")

hist(scotian_with_troph$Troph, main = "Distribution of trophic level estimates in SCS", xlab = "Trophic level estimate", ylab = "Frequency (number of catches)",  col = "orange")
boxplot(scotian_with_troph$Troph, main = "Boxplot of trophic levels by catch in SCS Community", xlab = "Note higher average than community average", ylab = "Weighted average trophic level estimate")

troph_level_boxplot_data <- bind_rows(
  tibble(Category = "Possible Community", Measurement = troph_dict),
  tibble(Category = "Caught Community", Measurement = scotian_with_troph$Troph)
)


# The 'data_df' now contains all values in 'Measurement' and labels in 'Category', which is the correct structure for ggplot's boxplot.
boxplot_fill <- ggplot(troph_level_boxplot_data, aes(x = Category, y = Measurement)) +
  # Set the fill color for the box bodies to light green
  geom_boxplot(fill = "orange") +
  labs(title = "High trophic level bias in SCS catch") +
  theme_minimal()
print(boxplot_fill)

# visualize number of identifications with trophic assignment
total_taxa <- unique(scotian_shelf_year_round$accepted_name)
num_total_taxa <- length(total_taxa)
#num_total_species <- length(unique(SCS_by_species$accepted_name))
#num_trophic_identified_species <- length(unique(scotian_with_troph$accepted_name))
rank_names <- c("Total taxonomic identifications", "Total species rank identifications", "Total trophic assigned species rank identifications")
#taxa_numbers <- c(num_total_taxa, num_total_species, num_trophic_identified_species)
# Create a data frame from the lists
#taxa_barplot <- data.frame(Name = rank_names, Value = taxa_numbers)

# Create the barplot
#barplot(taxa_barplot$Value, 
#names.arg = taxa_barplot$Name, 
# col = "blue", 
#main = "Breakdown of taxonomic identifications and capture of trophic assignment", 
#ylab = "Number of unique taxonomic identifications")

# generate CPUA for trophic groups as time series after defining trophic groups considered'


# suggest trophic groupings of <3.5, 3.5-4.0, and >4.0
scotian_with_troph <- scotian_with_troph %>%
  mutate(
    trophic_group = case_when(
      Troph < 3.5 ~ "TL3",
      Troph >= 3.5 & Troph < 4.0 ~ "TL3.5",
      TRUE ~ "TL4" 
    )
  )
# now trophic groups are assigned to all observations in scotian_with_troph

TL3_all <- subset(scotian_with_troph, trophic_group == 'TL3')
TL3t <- aggregate(wgt_cpua ~ year, data = TL3_all, mean)
plot(TL3t$year, TL3t$wgt_cpua, type = "l", main = "Mean annual weight CPUA for trophic group TL3: T < 3.5", xlab = "Year", ylab = "Mean annual weight CPUA (kg/km^2)", ylim = c(0, 500))
TL3.5_all <- subset(scotian_with_troph, trophic_group == 'TL3.5')
TL3.5t <- aggregate(wgt_cpua ~ year, data = TL3.5_all, mean)
plot(TL3.5t$year, TL3.5t$wgt_cpua, type = "l", main = "Mean annual weight CPUA for trophic group TL3.5: 3.5 =< T < 4.0", xlab = "Year", ylab = "Mean annual weight CPUA (kg/km^2)", ylim = c(0, 500))
TL4_all <- subset(scotian_with_troph, trophic_group == 'TL4')
TL4t <- aggregate(wgt_cpua ~ year, data = TL4_all, mean)
plot(TL4t$year, TL4t$wgt_cpua, type = "l", main = "Mean annual weight CPUA for  group TL4: T =< 4.0", xlab = "Year", ylab = "Mean annual weight CPUA (kg/km^2)", ylim = c(0, 500))

####################################################################
############ Build haul level average trophic level ################
####################################################################
# Calculate weighted average trophic level for each haul, accounting for SE_Troph, but not biomass
SCS_hauls_with_avg_troph <- scotian_with_troph %>%
  group_by(haul_id, year) %>%
  summarise(
    # The weight is the inverse of the variance (SE_Troph^2)
    weight = sum(1 / (SE_Troph^2), na.rm = TRUE),
    
    # Weighted average formula: sum(Troph * weight) / sum(weight)
    avg_trophic_level = sum(Troph / (SE_Troph^2), na.rm = TRUE) / weight,
    
    # Calculate the standard error of the weighted average for each haul
    # Formula: 1 / sqrt(sum(1 / SE_Troph^2))
    SE_avg_trophic_level = 1 / sqrt(weight),
    .groups = 'drop'
  )
# Calculate weighted average trophic level for each haul, accounting for SE_Troph AND biomass V3!!
SCS_hauls_with_avg_troph <- scotian_with_troph %>%
  group_by(haul_id, year) %>%
  # Filter out rows missing critical weighting data to avoid NaN results
  filter(!is.na(wgt_cpua), !is.na(Troph), !is.na(SE_Troph), SE_Troph > 0) %>%
  summarise(
    # 1. Calculate the composite weight: (Biomass) * (1 / Variance)
    # This gives more power to high-biomass species with certain TLs
    comp_weight = sum(wgt_cpua * (1 / SE_Troph^2), na.rm = TRUE),
    
    # 2. Weighted average: sum(Value * CompositeWeight) / sum(CompositeWeight)
    avg_trophic_level = sum(Troph * (wgt_cpua * (1 / SE_Troph^2)), na.rm = TRUE) / comp_weight,
    
    # 3. Calculate the standard error of the weighted mean
    # Note: This reflects the uncertainty of the estimate within that haul
    SE_avg_trophic_level = 1 / sqrt(sum(1 / SE_Troph^2, na.rm = TRUE)),
    
    # 4. Keep track of total biomass for the haul
    total_biomass_cpua = sum(wgt_cpua, na.rm = TRUE),
    
    .groups = 'drop'
  )

# Calculate yearly average, standard deviation, and standard error
yearly_average_trophic_level <- SCS_hauls_with_avg_troph %>%
  group_by(year) %>%
  summarise(
    overall_yearly_avg_trophic_level = mean(avg_trophic_level, na.rm = TRUE),
    yearly_sd = sd(avg_trophic_level, na.rm = TRUE),
    n_hauls = n(),
    .groups = 'drop'
  ) %>%
  # Calculate the standard error and the 95% confidence interval
  mutate(
    yearly_se = yearly_sd / sqrt(n_hauls),
    ci_lower = overall_yearly_avg_trophic_level - 1.96 * yearly_se,
    ci_upper = overall_yearly_avg_trophic_level + 1.96 * yearly_se,
    # Add a 5-year rolling average. `fill = NA` ensures the first and last
    # values (where a full 5-year window isn't available) are NA.
    rolling_avg = rollmean(overall_yearly_avg_trophic_level, k = 5, fill = NA, align = "center")
  )

yearly_trophic_level_plot <- ggplot(yearly_average_trophic_level, 
                                    aes(x = as.numeric(year), y = overall_yearly_avg_trophic_level)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "orange", alpha = 0.4) +
  geom_line(linetype = "dashed") +
  # Add a line for the 5-year rolling average
  geom_line(aes(y = rolling_avg), color = "red", linewidth = 1.2) +
  # Add points for each year
  geom_point(size = 2) +
  labs(
    title = "Shelf average trophic level",
    x = "Year",
    y = "Yearly weighted average trophic level"
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
  )

print(yearly_trophic_level_plot)

############################################################################
######## analyze body weight trends in trophic groups over time ############
############################################################################

clean_body_weight_scotian_with_troph <- scotian_with_troph %>%
  filter(body_weight > 0 & body_weight < 1000)

# 1. Calculate the yearly average body weight for each trophic group
yearly_avg_weight <- clean_body_weight_scotian_with_troph %>%
  group_by(trophic_group, year) %>%
  summarise(average_body_weight = mean(body_weight, na.rm = TRUE),
            .groups = 'drop') # Drop grouping after summarization

# 2. Create the plot using ggplot2
ggplot(yearly_avg_weight, aes(x = year, y = average_body_weight, color = trophic_group)) +
  geom_line() +
  geom_point() + # Optional: Add points for each year
  labs(
    title = "Yearly Average Body Weight by Trophic Group",
    x = "Year",
    y = "Summertime Yearly Average Body Weight (kg/km^2/individual) of Trophic Group",
    color = "Trophic Group"
  ) +
  theme_minimal() + # You can choose other themes as well
  scale_x_continuous(breaks = unique(yearly_avg_weight$year)) # Ensure all years are shown on x-axis

# 3. Now plot proportional change in body weight as percent difference from initial value
TL3_body_weight_1970 <- yearly_avg_weight %>%
  filter(trophic_group == 'TL3' & year == 1970) %>%
  pull(average_body_weight) # Use pull() to extract the single value
TL3.5_body_weight_1970 <- yearly_avg_weight %>%
  filter(trophic_group == 'TL3.5' & year == 1970) %>%
  pull(average_body_weight) # Use pull() to extract the single value
TL4_body_weight_1970 <- yearly_avg_weight %>%
  filter(trophic_group == 'TL4' & year == 1970) %>%
  pull(average_body_weight) # Use pull() to extract the single value

TL3_body_weight <- yearly_avg_weight %>%
  filter(trophic_group == 'TL3')
TL3_body_weight$average_body_weight <- TL3_body_weight$average_body_weight / TL3_body_weight_1970

TL3.5_body_weight <- yearly_avg_weight %>%
  filter(trophic_group == 'TL3.5')
TL3.5_body_weight$average_body_weight <- TL3.5_body_weight$average_body_weight / TL3.5_body_weight_1970

TL4_body_weight <- yearly_avg_weight %>%
  filter(trophic_group == 'TL4')
TL4_body_weight$average_body_weight <- TL4_body_weight$average_body_weight / TL4_body_weight_1970

year_list <- sort(unique(scotian_with_troph$year))
proportional_change_in_body_weight <- data.frame(Year = year_list, TL3 = TL3_body_weight$average_body_weight, TL3.5 = TL3.5_body_weight$average_body_weight, TL4 = TL4_body_weight$average_body_weight)

# 1. Reshape the data into a long format
proportional_change_long <- proportional_change_in_body_weight %>%
  pivot_longer(cols = c(TL3, TL3.5, TL4),
               names_to = "Trophic_Group",
               values_to = "Proportional_Change")

# 2. Create the ggplot
ggplot(proportional_change_long, aes(x = Year, y = Proportional_Change, color = Trophic_Group)) +
  geom_line() +
  geom_point() + # Optional: Add points
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") + # line at 100% for reference
  labs(
    title = "Proportional Change in Body Weight by Trophic Group",
    x = "Year",
    y = "Proportion of initial body weight (kg/km^2/individual)",
    color = "Trophic Group"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = unique(proportional_change_long$Year)) # Ensure all years are shown

trophic_group_colors <- c(
  "TL3" = "red", 
  "TL3.5" = "green",  
  "TL4" = "blue"
)
# -------------------------


# 2. Create the ggplot
proportional_change_plot <- ggplot(proportional_change_long, 
                                   aes(x = Year, y = Proportional_Change, color = Trophic_Group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) + # Added size for better visibility
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 1) + # line at 100% for reference
  
  # MODIFICATION: Use scale_color_manual to map specific colors to specific groups
  scale_color_manual(values = trophic_group_colors) +
  
  labs(
    title = "Proportional Change in Body Weight by Trophic Group",
    x = "Year",
    y = "Proportion of initial body weight",
    color = "Trophic Group"
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 16, face = "bold")
  )

print(proportional_change_plot)


##############################################################################


X <- cbind(TL3t$wgt_cpua, TL3.5t$wgt_cpua, TL4t$wgt_cpua)

# alternative plot -> Stacked area
# 1. Prepare the data for ggplot2 (long format)
# This transforms the matrix into a data frame suitable for ggplot2
X_df <- data.frame(Time = 1:nrow(X), X) # Add a Time column
X_df_long <- reshape2::melt(X_df, id.vars = "Time", variable.name = "Series", value.name = "Value")
X_df_long$year <- c(1970:2020)

########################################################################
############### Build community datatable for SCS summer ###############
########################################################################

trophic_source_table <- avg_trophic_key
trophic_source_table <- trophic_source_table %>%
  mutate(
    trophic_group = case_when(
      WE_trophic_estimate < 3.5 ~ "TL3",
      WE_trophic_estimate >= 3.5 & WE_trophic_estimate < 4.0 ~ "TL3.5",
      TRUE ~ "TL4" 
    )
  )

source <- data.frame(Diet_table = logical(), Food_troph_ecology_table = logical(), Diet_troph_ecology_table = logical())

for (species in trophic_source_table$Species){
  print(species)
  
  new_source <- data.frame(Diet_table = FALSE, Food_troph_ecology_table = FALSE, Diet_troph_ecology_table = FALSE)
  if (species %in% as.list(diet_table$Species)){
    new_source$Diet_table <- TRUE
  }
  if(species %in% as.list(eco_trophic_key$Species)){
    new_source$Food_troph_ecology_table <- TRUE
  }
  if(species %in% as.list(eco_trophic_key_diet$Species)){
    new_source$Diet_troph_ecology_table <- TRUE
  }
  source<- bind_rows(source, new_source)
}
trophic_source_table <- cbind(trophic_source_table, source)


################## Examine rarefaction by trophic group #########################################

scotian_summer_with_troph_rarefied_by_area <- scotian_with_troph %>%
  group_by(year, trophic_group) %>%
  summarise(
    total_area_swept = sum(area_swept),
    unique_accepted_names = n_distinct(accepted_name),
    .groups = 'drop' # Remove grouping after summarization
  )

r_squared_data <- scotian_summer_with_troph_rarefied_by_area %>%
  group_by(trophic_group) %>%
  summarise(
    # Fit a linear model for each group
    mod = list(lm(unique_accepted_names ~ total_area_swept, data = .)),
    .groups = 'drop'
  ) %>%
  # Extract the R-squared from each model's summary
  mutate(
    r_squared = sapply(mod, function(x) summary(x)$r.squared),
    # Create the label for the plot
    r2_label = paste("R^2 =", round(r_squared, 3))
  )

ggplot(scotian_summer_with_troph_rarefied_by_area, aes(x = total_area_swept, y = unique_accepted_names, color = trophic_group)) +
  geom_point() +
  
  # Add a linear trendline for each facet
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  
  # Add the R-squared value as text to each facet
  # We use a separate data frame for the text and specify its position
  geom_text(data = r_squared_data, aes(label = r2_label), 
            x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5, 
            inherit.aes = FALSE, size = 3) +
  scale_color_manual(values = trophic_group_colors) +
  facet_wrap(~trophic_group) +
  labs(
    x = "Yearly Area Swept (km^2)",
    y = "Yearly Number of Unique Accepted Names",
    color = "Trophic Group",
    title = "Rarefaction curves for yearly species richness by haul area sampling effort"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )


# looks like no diminishing return here either, maybe evidence of critical value, beyond which less benefit, especially in TL4


#####################################################################################
###### Synchrony amongst trophic groups following Loreau and Mazancourt 2008 ########
#####################################################################################

SCS_trophic_compensation_index <- generate_compensation_index(X_df)
SCS_trophic_synchrony <- generate_synchrony_index(X_df)

# Define window size
window_size <- 5
years <- sort(X_df$Time)
min_year <- min(years)
max_year <- max(years)
window_starts <- seq(min_year, max_year - window_size + 1, by = 1)
window_center_years <- window_starts + floor(window_size / 2)

# Initialize vectors to store results
trophic_comp_indices <- numeric(length(window_starts))
trophic_synch_indices <- numeric(length(window_starts))

# Loop through each sliding window
for (i in seq_along(window_starts)) {
  # Filter data for the current window
  window_data <- X_df %>%
    filter(Time >= window_starts[i], Time < window_starts[i] + window_size)
  
  # Check if there is enough data in the window
  if (nrow(window_data) >= window_size) {
    # Calculate compensation for the current window
    comp_result <- generate_compensation_index(window_data)
    trophic_comp_indices[i] <- comp_result$compensation_index
    
    # Calculate synchrony for the current window
    synch_result <- generate_synchrony_index(window_data)
    trophic_synch_indices[i] <- synch_result
  } else {
    trophic_comp_indices[i] <- NA
    trophic_synch_indices[i] <- NA
  }
}

# Create a data frame for plotting the results
trophic_comp_sync_indices <- tibble(
  Window = window_center_years,
  TrophicCompensation = trophic_comp_indices,
  TrophicSynchrony = trophic_synch_indices
)

# Plot 1: 5-Year Rolling Trophic Group Compensation

ggplot(trophic_comp_sync_indices, aes(x = Window, y = TrophicCompensation, group = 1)) +
  geom_line() +
  geom_point() +
  labs(title = "5-Year Rolling Compensation Index between Trophic Groups",
       x = "Center Year of Window",
       y = "Compensation Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: 5-Year Rolling Trophic Group Synchrony
ggplot(trophic_comp_sync_indices, aes(x = Window, y = TrophicSynchrony, group = 1)) +
  geom_line() +
  geom_point() +
  labs(title = "5-Year Rolling Synchrony Index between Trophic Groups",
       x = "Center Year of Window",
       y = "Synchrony Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##################################################
##### make combined trophic facet figure #########
##################################################

combined_trophic_plot <- (T_level_hist + yearly_trophic_level_plot) +
  plot_annotation(tag_levels = 'A') 

print(combined_trophic_plot)

##########################################################################
###### Section 2: Biodiversity haul and survey level biodiversity ########
##########################################################################
SCS_bio_summary <- summarize_biological(scotian_shelf_summer_core) #use barplot from here for survey length biodiversity by rank
barplot_richness <- ggplot2::last_plot()
# make species and supraspecies richness plot
SCS_haul_supraspecies_richness <- haul_level_supra_species_richness(scotian_shelf_summer_core)
SCS_haul_species_richness <- haul_level_species_richness(scotian_shelf_summer_core_by_species)
species_richness_plot <- ggplot2::last_plot()  
SCS_haul_biodiversity <- right_join(SCS_haul_species_richness, SCS_haul_supraspecies_richness, by = "haul_id")
SCS_haul_evenness <- haul_level_evenness(scotian_shelf_summer_core)
evenness_plot <- ggplot2::last_plot()
SCS_haul_biodiversity <- left_join(SCS_haul_biodiversity, SCS_haul_evenness, by = "haul_id")   #adds evenness as another column to the haul_level dataframe of biodiversity composed of richness and evenness estimates with haul data

#################################################################
########### Beta diversity between NAFO divisions ###############
#################################################################

SCS_dissimilarity_by_NAFO <- dissimilarity_matrices_by_NAFO_division(scotian_shelf_summer_core)
beta_div_plot <- ggplot2::last_plot()
###################################################################################
########## Change in body weight and trophic level by NAFO division ###############
###################################################################################
"
# 1. Calculate the yearly average body weight for each NAFO division
NAFO_yearly_avg_weight <- scotian_with_troph %>%
  group_by(Nafo_division, year) %>%
  summarise(average_body_weight = mean(body_weight, na.rm = TRUE),
            .groups = 'drop') # Drop grouping after summarization

# 2. Create the plot using ggplot2
ggplot(NAFO_yearly_avg_weight, aes(x = year, y = average_body_weight, color = Nafo_division)) +
  geom_line() +
  geom_point() + # Optional: Add points for each year
  labs(
    title = Yearly Average Body Weight by NAFO Division,
    x = Year,
    y = Yearly Average Body Weight (kg/km^2/individual) of NAFO division,
    color = NAFO Division
  ) +
  theme_minimal() + # You can choose other themes as well
  scale_x_continuous(breaks = unique(yearly_avg_weight$year)) # Ensure all years are shown on x-axis

# 1. Calculate the yearly summertime average trophic level for each NAFO division
NAFO_yearly_avg_trophic_level <- scotian_with_troph %>%
  group_by(Nafo_division, year) %>%
  summarise(average_trophic_level = mean(Troph, na.rm = TRUE),
            .groups = 'drop') # Drop grouping after summarization

# 2. Create the plot using ggplot2
ggplot(NAFO_yearly_avg_trophic_level, aes(x = year, y = average_trophic_level, color = Nafo_division)) +
  geom_line() +
  geom_point() + # Optional: Add points for each year
  labs(
    title = Yearly Summertime Average Trophic Level by NAFO Division,
    x = Year,
    y = Yearly Average Trophic Level of NAFO division,
    color = NAFO Division
  ) +
  theme_minimal() + # You can choose other themes as well
  scale_x_continuous(breaks = unique(yearly_avg_weight$year)) # Ensure all years are shown on x-axis
"

######## construct faceted biodiversity plot ###########

# You can arrange them in various ways, for example 2TL3.5:
combined_biodiversity_plot <- (barplot_richness + species_richness_plot) / (evenness_plot + beta_div_plot) +
  plot_annotation(tag_levels = 'A') +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.tag = element_text(face = "bold", size = 16)
  )

print(combined_biodiversity_plot)


##################################################################
############### Investigation of rarity on SCS ###################
##################################################################
SCS_haul_prevalence <- SCS_haul_prevalence %>%
  mutate(
    rarity_group = case_when(
      haul_prevalence_list >= 1000 ~ "R1",
      haul_prevalence_list < 1000 & haul_prevalence_list >= 100 ~ "R2",
      haul_prevalence_list < 100 & haul_prevalence_list >= 10 ~ "R3",
      TRUE ~ "R4" 
    )
  )
# now add rarity group to SCS
names(SCS_haul_prevalence)[1] <- "accepted_name" # so species names match up with scotian_shelf_year_round
scotian_shelf_summer_core <- merge(scotian_shelf_summer_core, SCS_haul_prevalence, by = "accepted_name", all.x = TRUE)

# now subset SCS by rarity group

scotian_shelf_year_round_R1 <- scotian_shelf_summer_core %>%
  filter(rarity_group == "R1")
scotian_shelf_year_round_R2 <- scotian_shelf_summer_core %>%
  filter(rarity_group == "R2")
scotian_shelf_year_round_R3 <- scotian_shelf_summer_core %>%
  filter(rarity_group == "R3")
scotian_shelf_year_round_R4 <- scotian_shelf_summer_core %>%
  filter(rarity_group == "R4")

dev.new()
# and has columns 'value' (numerical) and 'group' (categorical with 4 levels)

ggplot(scotian_shelf_summer_core, aes(x = year, fill = rarity_group)) +
  geom_histogram(binwidth = 1, color = "black") +
  facet_wrap(~ rarity_group, nrow = 2, scales = "free_y") +
  labs(
    title = "Yearly haul prevalence by rarity group",
    x = "Year",
    y = "Yearly haul frequency",
    fill = "Rarity Group"  # Add a legend title
  ) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "R1" = "darkblue",
    "R2" = "lightblue",
    "R3" = "purple",
    "R4" = "red"
  ))


# check body weight distribution of rarest taxa
scotian_shelf_year_round_R4 <- scotian_shelf_summer_core %>%
  filter(rarity_group == "R4" & body_weight > 0)

hist(log(scotian_shelf_year_round_R4$body_weight), main = "Distribution of body weights among rarest R4 taxa", xlab = "Log(body weight)", ylab = "Frequency (number of catches)") # note many small taxa and some rare large ones

# Define the custom color vector outside the plot, ensuring factor levels match
custom_colors <- c(
  "All Taxa" = "black",
  "R1" = "red",
  "R2" = "lightblue",
  "R3" = "lightgreen",
  "R4" = "purple"
)

# ----------------------------------------------------------------------
# 1. Calculate Average Body Weight for each Rarity Group 
# ----------------------------------------------------------------------
average_weights_R4 <- scotian_shelf_year_round_R4 %>%
  group_by(accepted_name) %>%
  filter(body_weight <= 1000) %>%
  summarise(average_weight = log(mean(body_weight, na.rm = TRUE)))

average_weights_R3 <- scotian_shelf_year_round_R3 %>%
  group_by(accepted_name) %>%
  filter(body_weight <= 1000) %>%
  summarise(average_weight = log(mean(body_weight, na.rm = TRUE)))

average_weights_R2 <- scotian_shelf_year_round_R2 %>%
  group_by(accepted_name) %>%
  filter(body_weight <= 1000) %>%
  summarise(average_weight = log(mean(body_weight, na.rm = TRUE)))

average_weights_R1 <- scotian_shelf_year_round_R1 %>%
  group_by(accepted_name) %>%
  filter(body_weight <= 1000) %>%
  summarise(average_weight = log(mean(body_weight, na.rm = TRUE)))

# ----------------------------------------------------------------------
# 2. Calculate Average Body Weight for ALL Taxa (NEW)
# ----------------------------------------------------------------------

# We need a data frame containing ALL unique taxa from R1, R2, R3, R4 combined.
# The original data (scotian_shelf_year_round) should contain all records.
# Assuming 'scotian_shelf_year_round' is the master dataframe.
all_taxa_data <- bind_rows(scotian_shelf_year_round_R1, 
                           scotian_shelf_year_round_R2, 
                           scotian_shelf_year_round_R3, 
                           scotian_shelf_year_round_R4) %>%
  # Filter for unique names to get the full list of species considered
  distinct(accepted_name, .keep_all = TRUE) 


average_weights_All <- all_taxa_data %>%
  group_by(accepted_name) %>%
  filter(body_weight <= 1000) %>%
  summarise(average_weight = log(mean(body_weight, na.rm = TRUE))) %>%
  mutate(Series = "All Taxa") # Label the new series

# ----------------------------------------------------------------------
# 3. Combine ALL results and set factor order
# ----------------------------------------------------------------------

combined_R_wgt_df <- bind_rows(
  average_weights_All, # NEW: Include the 'All Taxa' data first
  average_weights_R1 %>% mutate(Series = "R1"),
  average_weights_R2 %>% mutate(Series = "R2"),
  average_weights_R3 %>% mutate(Series = "R3"),
  average_weights_R4 %>% mutate(Series = "R4")
)

combined_R_wgt_df <- combined_R_wgt_df %>%
  filter(is.finite(average_weight)) # Remove any rows with NA average weights
# Set the factor levels to ensure correct plotting order (All Taxa first)
combined_R_wgt_df$Series <- factor(combined_R_wgt_df$Series, 
                                   levels = c("All Taxa", "R1", "R2", "R3", "R4"))

# ----------------------------------------------------------------------
# 4. Generate the Boxplot with Custom Colors and Attempt more robust test for stat sig diff between rarity group body weight
# ----------------------------------------------------------------------

# 2. Run the Linear Model (Equivalent to ANOVA)
# Formula: MTC_value ~ Series (Rarity Group)
w_anova_model <- lm(average_weight ~ Series, data = combined_R_wgt_df)

# 3. View the ANOVA Table
summary(w_anova_model)
anova(w_anova_model)

# 4. Post-hoc test (Tukey) to see which specific groups differ
# Note: We use the 'emmeans' package as it handles weights correctly
library(emmeans)
post_hoc <- emmeans(w_anova_model, list(pairwise ~ Series), adjust = "tukey")
print(post_hoc)

# 1. Create a data frame for the labels based on your Tukey output
# Group 'a' is significantly different from 'b'
label_data <- data.frame(
  Series = c("R1", "All Taxa", "R2", "R3", "R4"),
  Letters = c("b", "a", "a", "a", "a"),
  y_pos = 6 # This puts the letters just above your max residual/value
)

# 2. Build the plot
boxplot_bwgt_R <- ggplot(combined_R_wgt_df, aes(x = Series, y = average_weight, fill = Series)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  geom_text(data = label_data, aes(x = Series, y = y_pos, label = Letters), 
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Average body weight by rarity group",
    subtitle = "Groups sharing a letter are not significantly different (Tukey HSD, p < 0.05)",
    x = "Taxa Group",
    y = "Ln Average Body Weight (kg)",
    fill = "Taxa Group"
  ) +
  # Apply a clean theme and optionally adjust plot elements
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
  )

print(boxplot_bwgt_R)

##############################################################
########## Time of emergence / disappearance in R4 ##########
##############################################################

# note difference in emergence_disappearance between R3 and R4 with less agreement in emergence_disappearance in R3

# --- 1. Calculate Emergence/Disappearance Summaries for R4 and R3 ---
R4_summary <- scotian_shelf_year_round_R4 %>%
  group_by(accepted_name) %>%
  summarize(
    first_year = min(year),
    last_year = max(year)
  )

R3_summary <- scotian_shelf_year_round_R3 %>%
  group_by(accepted_name) %>%
  summarize(
    first_year = min(year),
    last_year = max(year)
  )

# Define shared parameters
custom_bins <- 51
x_lims <- c(1970, 2020)

# 3. R4 - Emergence
p3 <- ggplot(R4_summary, aes(x = first_year)) +
  geom_histogram(bins = custom_bins, fill = "dodgerblue", color = "black", alpha = 0.8) +
  coord_cartesian(xlim = x_lims) +
  labs(title = "R4 - Emergence (First Year)", x = "Year", y = "Unique taxa") +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
  )

# 4. R4 - Disappearance
p4 <- ggplot(R4_summary, aes(x = last_year)) +
  geom_histogram(bins = custom_bins, fill = "firebrick", color = "black", alpha = 0.8) +
  coord_cartesian(xlim = x_lims) +
  labs(title = "R4 - Disappearance (Last Year)", x = "Year", y = "Unique Taxa") +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
  )



# To view them, simply call the object name:
p3
p4

# --- 3. Create the Faceted Histogram Plot ---
combined_ToE_ToD_plot <- (p3 + p4) +
  plot_annotation(tag_levels = 'A') + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.tag = element_text(face = "bold", size = 16)
  )

print(combined_ToE_ToD_plot)


#####################################################
### Change in catch of rare species ############### REVIST AND DECIDE TO DO YEARLY AVERAGE OF TOTAL CPUA PER HAUL OR TAXON CPUA PER HAUL
######################################################

R1t <- aggregate(wgt_cpua ~ year, data = scotian_shelf_year_round_R1, mean)
R2t <- aggregate(wgt_cpua ~ year, data = scotian_shelf_year_round_R2, mean)
R3t <- aggregate(wgt_cpua ~ year, data = scotian_shelf_year_round_R3, mean)
R4t <- aggregate(wgt_cpua ~ year, data = scotian_shelf_year_round_R4, mean) # note no rare species in 1970 so gotta add place holder row
# Placeholder for R4 1970
if(!any(R4t$year == 1970)) {
  R4t <- rbind(data.frame(year = 1970, wgt_cpua = 0), R4t)
  R4t <- R4t[order(R4t$year), ] # Ensure it's sorted by year
}

# PLOT
R <- cbind(R1t$wgt_cpua, R2t$wgt_cpua, R3t$wgt_cpua, R4t$wgt_cpua)

# Vector of actual years
actual_years <- R1t$year

# 1. Prepare the data for ggplot2 (long format)
# Use actual_years instead of 1:nrow
R_df <- data.frame(Year = R1t$year, R) 
names(R_df) <- c("Year", "R1", "R2", "R3", "R4")

R_df_long <- reshape2::melt(R_df, id.vars = "Year", variable.name = "Series", value.name = "Value")

#########################################################
##### Compensation index between rarity groups ##########
#########################################################
R_df$Time <- sort(unique(scotian_shelf_summer_core$year))
SCS_rarity_compensation_index <- generate_compensation_index(R_df)
'
Coefficient of Variation for each group:
     R1_CV    R2_CV    R3_CV    R4_CV
1 22.78739 47.37386 123.8912 152.9032

Coefficient of Variation for total wgt_cpua: 23.55162 

Compensation Index (based on CV of total):
[1] 0.7820427

suggests compensation is stabilizing community total wgt_cpua
' 
# now try 5 year window compensation index
window_size <- 5
years <- sort(R_df$Time)
min_year <- min(years)
max_year <- max(years)
window_starts <- seq(min_year, max_year - window_size + 1, by = 1)
comp_indices <- numeric(length(window_starts))
window_labels <- character(length(window_starts))

for (i in seq_along(window_starts)) {
  window_data <- R_df %>%
    filter(Time >= window_starts[i], Time < window_starts[i] + window_size)
  if (nrow(window_data) >= window_size) {
    result <- generate_compensation_index(window_data)
    comp_indices[i] <- result$compensation_index
  } else {
    comp_indices[i] <- NA
  }
  window_labels[i] <- paste0(window_starts[i], "-", window_starts[i] + window_size - 1)
}

# Plot compensation index over time
rare_comp_indices <- tibble(Window = window_labels, CompensationIndex = comp_indices)
ggplot(rare_comp_indices, aes(x = Window, y = CompensationIndex, group = 1)) +
  geom_line() +
  geom_point() +
  labs(title = "5-Year Rolling Compensation Index among Rarity Groups", x = "Year Window", y = "Compensation Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#now try synchrony
SCS_rarity_synchrony_index <- generate_synchrony_index(R_df)

# now try 5 year window synchrony index
window_size <- 5
years <- sort(R_df$Time)
min_year <- min(years)
max_year <- max(years)
window_starts <- seq(min_year, max_year - window_size + 1, by = 1)
synch_indices <- numeric(length(window_starts))
window_labels <- character(length(window_starts))

for (i in seq_along(window_starts)) {
  window_data <- R_df %>%
    filter(Time >= window_starts[i], Time < window_starts[i] + window_size)
  if (nrow(window_data) >= window_size) {
    result <- generate_synchrony_index(window_data)
    synch_indices[i] <- result
  } else {
    synch_indices[i] <- NA
  }
  window_labels[i] <- paste0(window_starts[i], "-", window_starts[i] + window_size - 1)
}

# Plot synchrony index over time
rare_sync_indices <- tibble(Window = window_labels, SynchIndex = synch_indices)
ggplot(rare_sync_indices, aes(x = Window, y = SynchIndex, group = 1)) +
  geom_line() +
  geom_point() +
  labs(title = "5-Year Rolling Synchrony Index among Rarity Groups", x = "Year Window", y = "Synchrony Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



######################################################
########## MTC using new thermal data ################
######################################################
'
# V2 MTC

load("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/FISHGLOB Data Project/thermal affinity data/V2 thermal affinity data/Archive/META_02122025.Rdata")
length(unique(scotian_shelf_summer_core_by_species$aphia_id))

id_list_for_MTC2 <- unique(scotian_shelf_summer_core_by_species$aphia_id)
META_SCS <- META %>%
  filter(AphiaID %in% id_list_for_MTC2)

length(unique(META_SCS$AphiaID))
'

#####################################################

# V3 MTC
load("C:/Users/besch/Desktop/New Niche Files/META_02122025.Rdata")

niche_fields <- colnames(META)
print(niche_fields)

# Strategy: use seafloor occ niche files when available, then use pelagic occ niche files if unavailable

id_list_for_MTC3 <- unique(scotian_shelf_summer_core$aphia_id) # 353 things
META_SCS <- META %>%
  filter(AphiaID %in% id_list_for_MTC3)

length(unique(META_SCS$AphiaID)) # 279 things with niche data

############# add thermal data from individual r data files ########

# 1. Setup paths
data_dir <- "C:/Users/besch/Desktop/New Niche Files/NICHE_SEAFLOOR_OCC/"
AphiaID_list <- sort(unique(META_SCS$AphiaID))

# 2. Run the loop
seafloor_occ_thermals <- map_df(AphiaID_list, function(id) {
  
  # Construct path (trying .Rdata based on your successful manual load)
  file_path <- paste0(data_dir, "NICHE_MOD_SP_", id, ".Rdata")
  
  # If that doesn't exist, try the capital D version
  if (!file.exists(file_path)) {
    file_path <- paste0(data_dir, "NICHE_MOD_SP_", id, "_SEAFLOOR_OCC_FLAG1_NOOUTLIERS_ALLTEMP.RData")
  }
  
  if (file.exists(file_path)) {
    temp_env <- new.env()
    load(file_path, envir = temp_env)
    
    # Grab the object inside (the one you saw as 'my_data')
    data <- as.list(temp_env)[[1]]
    
    # Ensure it's a dataframe and add the AphiaID
    data <- as.data.frame(data) %>%
      mutate(AphiaID = id)
    
    return(data)
  } else {
    # If file is missing, return NULL (map_df will just skip it)
    return(NULL)
  }
})

# 3. Check the result
print(paste("Successfully loaded data for", length(unique(seafloor_occ_thermals$AphiaID)), "species."))
head(seafloor_occ_thermals)


# transform to wide data

seafloor_occ_thermals_wide <- seafloor_occ_thermals %>%
  pivot_wider(
    id_cols = AphiaID,
    names_from = variable, 
    values_from = c(mean, min, max, q01, q05, q25, q50, q75, q95, q99),
    names_glue = "{variable}_{.value}"
  )

# Preview the wide data
seafloor_occ_thermals_wide_MEANS_ONLY <- seafloor_occ_thermals_wide %>%
  select(
    AphiaID, 
    starts_with("thetao_mean")
  )

# Join by AphiaID
SCS_META_and_THERMALS <- META_SCS %>%
  mutate(AphiaID = as.double(AphiaID)) %>% # Ensure types match
  left_join(seafloor_occ_thermals_wide_MEANS_ONLY, by = "AphiaID")

# Check how many species now have thermal data
sum(!is.na(SCS_META_and_THERMALS$thetao_mean_mean)) # 186 aphiaIDs are benthic

reduced_SCS_META_and_THERMALS <- SCS_META_and_THERMALS %>%
  select(c(AphiaID, thetao_mean_mean, thetao_mean_min,thetao_mean_max,thetao_mean_q01,thetao_mean_q05,thetao_mean_q25,thetao_mean_q50,thetao_mean_q75,thetao_mean_q95, thetao_mean_q99))

# Now need thermal niches for species without benthic models
# 1. Setup paths
data_dir <- "C:/Users/besch/Desktop/New Niche Files/NICHE_PELAGIC_OCC/"
AphiaID_list <- sort(unique(META_SCS$AphiaID))

pelagic_occ_thermals <- map_df(AphiaID_list, function(id) {
  
  # Construct path (trying .Rdata based on your successful manual load)
  file_path <- paste0(data_dir, "NICHE_MOD_SP_", id, ".Rdata")
  
  # If that doesn't exist, try the capital D version
  if (!file.exists(file_path)) {
    file_path <- paste0(data_dir, "NICHE_MOD_SP_", id, "_PELAGIC_OCC_FLAG1_NOOUTLIERS_ALLTEMP.RData")
  }
  
  if (file.exists(file_path)) {
    temp_env <- new.env()
    load(file_path, envir = temp_env)
    
    # Grab the object inside (the one you saw as 'my_data')
    data <- as.list(temp_env)[[1]]
    
    # Ensure it's a dataframe and add the AphiaID
    data <- as.data.frame(data) %>%
      mutate(AphiaID = id)
    
    return(data)
  } else {
    # If file is missing, return NULL (map_df will just skip it)
    return(NULL)
  }
})

# 3. Check the result
print(paste("Successfully loaded data for", length(unique(pelagic_occ_thermals$AphiaID)), "species."))

# transform to wide data

pelagic_occ_thermals_wide <- pelagic_occ_thermals %>%
  pivot_wider(
    id_cols = AphiaID,
    names_from = variable, 
    values_from = c(mean, min, max, q01, q05, q25, q50, q75, q95, q99),
    names_glue = "{variable}_{.value}"
  )

# Preview the wide data
pelagic_occ_thermals_wide_MEANS_ONLY <- pelagic_occ_thermals_wide %>% # 278 species
  select(
    AphiaID, 
    starts_with("thetao_mean")
  )

pelagic_ids <- unique(pelagic_occ_thermals_wide_MEANS_ONLY$AphiaID)
seafloor_ids <- unique(seafloor_occ_thermals_wide_MEANS_ONLY$AphiaID)

head(pelagic_occ_thermals_wide_MEANS_ONLY)
head(seafloor_occ_thermals_wide_MEANS_ONLY)

# 1. Identify species in pelagic that are NOT in seafloor
pelagic_unique <- pelagic_occ_thermals_wide_MEANS_ONLY %>%
  anti_join(seafloor_occ_thermals_wide_MEANS_ONLY, by = "AphiaID")

# 2. Combine the preferred seafloor data with the unique pelagic rows
combined_pelagic_seafloor_occ_data <- bind_rows(seafloor_occ_thermals_wide_MEANS_ONLY, pelagic_unique)

# Join by AphiaID
SCS_META_and_THERMALS <- META_SCS %>%
  mutate(AphiaID = as.double(AphiaID)) %>% # Ensure types match
  left_join(combined_pelagic_seafloor_occ_data, by = "AphiaID")

# Check how many species now have thermal data
sum(!is.na(SCS_META_and_THERMALS$thetao_mean_mean))

missing_AphiaID <- setdiff(as.character(unique(META_SCS$AphiaID)), 
                           unique(combined_pelagic_seafloor_occ_data$AphiaID))
print(missing_AphiaID) # NO SDM FOR THIS GUY Cusk

colnames(SCS_META_and_THERMALS)

reduced_SCS_META_and_THERMALS <- SCS_META_and_THERMALS %>%
  select(c(AphiaID, thetao_mean_mean, thetao_mean_min,thetao_mean_max,thetao_mean_q01,thetao_mean_q05,thetao_mean_q25,thetao_mean_q50,thetao_mean_q75,thetao_mean_q95, thetao_mean_q99))

##############################################################
#### calculate coverage of MTC3 for SCS and Fishglob #########
##############################################################
SCS_temp3_hauls <- left_join(
  scotian_shelf_summer_core, 
  reduced_SCS_META_and_THERMALS, 
  by = c("aphia_id" = "AphiaID")
)


SCS_temp3_hauls <- SCS_temp3_hauls %>%
  filter(!is.na(thetao_mean_mean))

nhaul_SCS_MTC3 <- length(unique(SCS_temp3_hauls$haul_id)) # 8941
ntaxa_SCS_MTC3 <- length(unique(SCS_temp3_hauls$accepted_name)) # 278


#################################################################
####### V3 MTC using AquaX ######################################
'Citation: Reygondeau G, Egorova Y, Boerder K, Tittensor DP, Kaschner K, Kesner-Reyes K, et al. (2026) AquaX: An enhanced and revised AquaMaps framework to model marine species distributions and biodiversity. PLoS One 21(2): e0335823. https://doi.org/10.1371/journal.pone.0335823
'
#################################################################

############################# MTC Derivation (Haul Level) ############################

# 1. Enhanced Thermal Dictionary Preparation
# This uses multiple points of the distribution to reach a consensus SD
simplified_thermal_dictionary <- reduced_SCS_META_and_THERMALS %>%
  mutate(
    # Estimate SD using three different quantile spreads
    sd_98 = (thetao_mean_q99 - thetao_mean_q01) / 4.645,
    sd_90 = (thetao_mean_q95 - thetao_mean_q05) / 3.290,
    sd_iqr = (thetao_mean_q75 - thetao_mean_q25) / 1.349,
    
    # Take the mean of these estimates for a more robust SD_Ti
    # This uses the breadth of the entire distribution provided
    SD_Ti = (sd_98 + sd_90 + sd_iqr) / 3,
    
    # Calculate Variance
    Variance_Ti = SD_Ti^2,
    
    # Rename mean for the join
    estim_thermal_mean = thetao_mean_mean
  ) %>%
  select(AphiaID, estim_thermal_mean, SD_Ti, Variance_Ti)

# 2. Join and MTC Calculation (Logic remains consistent)
haul_species_data <- SCS_temp3_hauls %>%
  mutate(wgt_cpua = replace_na(wgt_cpua, 0)) %>%
  left_join(simplified_thermal_dictionary, by = c("aphia_id" = "AphiaID")) %>%
  filter(!is.na(estim_thermal_mean))

# 3. Propagate uncertainty to the Haul Level
haul_mtc_results <- haul_species_data %>%
  group_by(haul_id, year) %>%
  summarise(
    total_wgt = sum(wgt_cpua),
    # Weighted average for MTC
    MTC_haul = sum(estim_thermal_mean * wgt_cpua) / total_wgt,
    # Variance propagation: (Sum of W_i^2 * Var_i) / (Sum of W_i)^2
    MTC_variance_Ti = sum((wgt_cpua^2 * Variance_Ti)) / (total_wgt^2),
    .groups = "drop"
  ) %>%
  filter(total_wgt > 0)

############################# MTC Derivation (Year Level) ############################

# 4. Yearly Summary using the Law of Total Variance
# Total Variance = Var(Haul Means) + Mean(Haul Variances)
MTC_yearly_summary <- haul_mtc_results %>%
  group_by(year) %>%
  summarise(
    MTC_year = mean(MTC_haul, na.rm = TRUE),
    n_hauls = n(),
    # Spatial/Sampling variance between hauls
    spatial_var = var(MTC_haul, na.rm = TRUE),
    # Average "internal" uncertainty from species niche widths
    niche_uncertainty = mean(MTC_variance_Ti, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # The SE of the year mean accounts for both types of noise
    MTC_SE = sqrt((spatial_var / n_hauls) + (niche_uncertainty / n_hauls)),
    CI_lower = MTC_year - 1.96 * MTC_SE,
    CI_upper = MTC_year + 1.96 * MTC_SE
  )
# 5. Final Formatting and CI Calculation
MTC <- MTC_yearly_summary %>%
  mutate(
    CI_lower = MTC_year - 1.96 * MTC_SE,
    CI_upper = MTC_year + 1.96 * MTC_SE
  ) %>%
  arrange(year)

# 5-year running average
MTC$MTC_5yr_avg <- rollmean(MTC$MTC_year, k = 5, fill = NA, align = "center")

############################# Visualization ############################

MTC_plot <- ggplot(MTC, aes(x = year)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), fill = "lightgreen", alpha = 0.5) + #Add 95% CI shaded region (light green)
  geom_line(aes(y = MTC_year), color = "black", linetype = "dashed") +
  geom_point(aes(y = MTC_year), color = "black", size = 2) + #Add yearly MTC points (black)
  geom_line(aes(y = MTC_5yr_avg), color = "red", size = 1.2) + #Add 5-year running average line (red)
  labs(
    y = "MTC (°C)", 
    x = "Year", 
    title = "Shelf average mean temperature of the catch",
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
  )

print(MTC_plot)

########################### Combine thermal data into main dataset #################################

scotian_shelf_summer_core <- left_join( # note NA's present thrughout!
  scotian_shelf_summer_core, 
  simplified_thermal_dictionary, 
  by = c("aphia_id" = "AphiaID")
)

##############################################################
######### Compare thermal affinity of Rarity groups ##########
##############################################################

# 1. Function to extract MTC for a given rarity group's taxa
# Assumes SCS_thermal_dictionary contains 'accepted_name' and 'estim_thermal_mean'
extract_MTC_data <- function(rarity_df, group_name) {
  # Get the unique accepted names for the current rarity group
  taxa <- unique(rarity_df$aphia_id)
  
  # Filter the MTC dictionary (assuming it's loaded as SCS_thermal_dictionary)
  # and get the estimated thermal mean for those species
  MTC_data <- simplified_thermal_dictionary %>%
    filter(AphiaID %in% taxa) %>%
    # Rename MTC column for clearer use in the plot
    rename(MTC_value = estim_thermal_mean) %>% 
    # Add the rarity group identifier
    mutate(Series = group_name)
  
  return(MTC_data)
}

# 2. Extract MTC data for all four rarity groups
MTC_R1 <- extract_MTC_data(scotian_shelf_year_round_R1, "R1")
MTC_R2 <- extract_MTC_data(scotian_shelf_year_round_R2, "R2")
MTC_R3 <- extract_MTC_data(scotian_shelf_year_round_R3, "R3")
MTC_R4 <- extract_MTC_data(scotian_shelf_year_round_R4, "R4")

# 3. Create the "All Taxa" dataset
MTC_All_Taxa <- simplified_thermal_dictionary %>%
  rename(MTC_value = estim_thermal_mean) %>% 
  mutate(Series = "All Taxa")

# 4. Combine ALL results into a single data frame
combined_R_MTC_df <- bind_rows(MTC_All_Taxa, MTC_R1, MTC_R2, MTC_R3, MTC_R4)

# 5. Optional: Ensure the Series factor levels are ordered correctly for the plot
combined_R_MTC_df$Series <- factor(combined_R_MTC_df$Series, 
                                   # Order: All Taxa first, then R1 through R4
                                   levels = c("All Taxa", "R1", "R2", "R3", "R4"))
# Define the custom color vector, ensuring the names match the 'Series' factor levels
custom_colors <- c(
  "All Taxa" = "black",
  "R1" = "red",
  "R2" = "lightblue",
  "R3" = "lightgreen",
  "R4" = "purple"
)

boxplot_MTC_R <- ggplot(combined_R_MTC_df, aes(x = Series, y = MTC_value, fill = Series)) +
  # Use geom_boxplot to create the boxplots
  geom_boxplot(outlier.shape = 1) +
  
  # 🔑 ADD THE CUSTOM COLOR SCALE HERE
  scale_fill_manual(values = custom_colors) +
  
  # Customize labels and title for MTC
  labs(
    title = "Preferred temperature by rarity groups",
    x = "Taxa Group",
    y = "Estimated Thermal Mean (°C)",
    fill = "Taxa Group"
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
  )
print(boxplot_MTC_R)

# note obviously higher T optimum for R4 than R1

########### Attempt more robust test for stat sig diff between rarity group MTC ############

# 1. Define Weights
# Inverse-variance weighting: species with smaller SD_Ti get higher weights
combined_R_MTC_df <- combined_R_MTC_df %>%
  mutate(weight = 1 / (SD_Ti^2))

# 2. Run the Weighted Linear Model (Equivalent to Weighted ANOVA)
# Formula: MTC_value ~ Series (Rarity Group)
w_anova_model <- lm(MTC_value ~ Series, data = combined_R_MTC_df, weights = weight)

# 3. View the ANOVA Table
summary(w_anova_model)
anova(w_anova_model)

# 4. Post-hoc test (Tukey) to see which specific groups differ
# Note: We use the 'emmeans' package as it handles weights correctly
library(emmeans)
post_hoc <- emmeans(w_anova_model, list(pairwise ~ Series), adjust = "tukey")
print(post_hoc)

boxplot_MTC_R <- ggplot(combined_R_MTC_df, aes(x = Series, y = MTC_value, fill = Series)) +
  # Use geom_boxplot to create the boxplots
  geom_boxplot(outlier.shape = 1) +
  
  # 🔑 ADD THE CUSTOM COLOR SCALE HERE
  scale_fill_manual(values = custom_colors) +
  
  # Customize labels and title for MTC
  labs(
    title = "Preferred temperature by rarity groups",
    x = "Taxa Group",
    y = "Estimated Thermal Mean (°C)",
    fill = "Taxa Group"
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
  )
print(boxplot_MTC_R)

# 1. Create the label dataframe for MTC
# Based on Tukey: Only R3 and R4 are significantly different.
mtc_label_data <- data.frame(
  Series = c("R4", "All Taxa", "R2", "R1", "R3"),
  Letters = c("b", "ab", "ab", "ab", "a"),
  y_pos = 25 # Adjusted to sit above the MTC scale
)

# 2. Updated Boxplot
boxplot_MTC_R <- ggplot(combined_R_MTC_df, aes(x = Series, y = MTC_value, fill = Series)) +
  geom_boxplot(outlier.shape = 1) +
  
  # Add the significance letters
  geom_text(data = mtc_label_data, aes(x = Series, y = y_pos, label = Letters), 
            size = 5, fontface = "bold", vjust = 0) +
  
  # Styling
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Preferred temperature by rarity group",
    subtitle = "Groups sharing a letter are not significantly different (Tukey HSD, p < 0.05)",
    x = "Taxa Group",
    y = "Estimated Thermal Mean (°C)",
    fill = "Taxa Group"
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
  )
print(boxplot_MTC_R)


############## Build faceted plot of rarity ###################

combined_rarity_plot <-   (boxplot_bwgt_R + boxplot_MTC_R) +
  plot_annotation(tag_levels = 'A') +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.tag = element_text(face = "bold", size = 16)
  )

print(combined_rarity_plot)


##################################################
###### Section 3: Structure/productivity #########
##################################################
#############################################################################################
###### Interspecific synchrony and compensation following Loreau and Mazancourt 2008 ########
#############################################################################################
yearly_average_community <- scotian_shelf_summer_core %>%
  group_by(year, accepted_name) %>%
  summarize(mean_wgt_cpua = mean(wgt_cpua, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(
    names_from = accepted_name,
    values_from = mean_wgt_cpua,
    values_fill = list(mean_wgt_cpua = 0)
  ) %>%
  rename(Time = year)

SCS_interspecific_compensation_total <- generate_compensation_index(yearly_average_community)
#now try synchrony
SCS_interspecific_synchrony_total <- generate_synchrony_index(yearly_average_community)

# now try 5 year window compensation index
window_size <- 5
years <- sort(yearly_average_community$Time)
min_year <- min(years)
max_year <- max(years)
window_starts <- seq(min_year, max_year - window_size + 1, by = 1)
comp_indices <- numeric(length(window_starts))
window_labels <- character(length(window_starts))

for (i in seq_along(window_starts)) {
  window_data <- yearly_average_community %>%
    filter(Time >= window_starts[i], Time < window_starts[i] + window_size)
  if (nrow(window_data) >= window_size) {
    result <- generate_compensation_index(window_data)
    comp_indices[i] <- result$compensation_index
  } else {
    comp_indices[i] <- NA
  }
  window_labels[i] <- paste0(window_starts[i], "-", window_starts[i] + window_size - 1)
}

# Plot compensation index over time
spec_comp_indices <- tibble(Window = window_labels, CompensationIndex = comp_indices)
ggplot(spec_comp_indices, aes(x = Window, y = CompensationIndex, group = 1)) +
  geom_line() +
  geom_point() +
  labs(title = "5-Year Rolling Compensation Index between all SCS Taxa", x = "Year Window", y = "Compensation Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# now try 5 year window synchrony index
window_size <- 5
years <- sort(yearly_average_community$Time)
min_year <- min(years)
max_year <- max(years)
window_starts <- seq(min_year, max_year - window_size + 1, by = 1)
synch_indices <- numeric(length(window_starts))
window_labels <- character(length(window_starts))

for (i in seq_along(window_starts)) {
  window_data <- yearly_average_community %>%
    filter(Time >= window_starts[i], Time < window_starts[i] + window_size)
  if (nrow(window_data) >= window_size) {
    result <- generate_synchrony_index(window_data)
    synch_indices[i] <- result
  } else {
    synch_indices[i] <- NA
  }
  window_labels[i] <- paste0(window_starts[i], "-", window_starts[i] + window_size - 1)
}

# Plot synchrony index over time
spec_sync_indices <- tibble(Window = window_labels, SynchIndex = synch_indices)
ggplot(spec_sync_indices, aes(x = Window, y = SynchIndex, group = 1)) +
  geom_line() +
  geom_point() +
  labs(title = "5-Year Rolling Synchrony Index between all SCS Taxa", x = "Year Window", y = "Synchrony Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# test relationship between synchrony and compensation, note when all group means and variances are equal, Comp = 1 - sqrt(Synch)

test_relation <- SCS_interspecific_compensation_total$compensation_index - 1 + sqrt(SCS_interspecific_synchrony_total)
test_relation_windowed <- comp_indices - 1 + sqrt(synch_indices)

'
Note following definitions of catch statistics
1. num_cpua, number of individuals (abundance) in individuals/km2
2. num_cpue, number of individuals (abundance) in individuals/h 
3. num, number of individuals (abundance)
4. wgt_cpua, weight in kg/km2
5. wgt_cpue, weight in kg/h 
6. wgt, weight in kg
'

########################################################
######### Body weight investigation ####################
########################################################

clean_body_weight_SCS <- scotian_shelf_summer_core %>%
  filter(body_weight > 0 & num > 0 & wgt > 0 & body_weight < 1000)

CPUA_body_weight_SCS <- scotian_shelf_summer_core %>%
  filter(num_cpua > 0 & wgt_cpua > 0)
CPUA_body_weight_SCS$body_weight_cpua <- CPUA_body_weight_SCS$wgt_cpua / CPUA_body_weight_SCS$num_cpua # the equivalence of body weights reported by clean_body_weight_SCS and CPUA_body_weight_SCS demonstrates that calculating this quantity with or without cpua results in same estimate

# check distribution of body weights on SCS

hist(clean_body_weight_SCS$body_weight, main = "Distribution of haul recorded body weights across all taxa in SCS region", xlab = "Body weight (kg)", ylab = "Frequency", breaks = 50)
hist(log10(clean_body_weight_SCS$body_weight), main = "Distribution of haul recorded body weights across all taxa in SCS region on log10 scale", xlab = "Log10(Body weight (kg))", ylab = "Frequency", breaks = 50)

####################################################################
######## See how community average body weight has changed #########
####################################################################

# 1. Calculate the average body weight for each Haul
weight_haul_summary <- clean_body_weight_SCS %>%
  group_by(year, haul_id) %>%
  summarise(
    # The mean body weight of all fish within that specific haul
    haul_avg_weight = mean(body_weight, na.rm = TRUE),
    n_haul_observations = n(), # Number of fish in the haul
    .groups = 'drop'
  )

# 2. Calculate the yearly average (of the haul means) and statistics
yearly_summary_haul_avg <- weight_haul_summary %>%
  group_by(year) %>%
  summarise(
    # The new yearly average is the mean of all haul averages in that year
    avg_weight = mean(haul_avg_weight, na.rm = TRUE),
    n_hauls = n(), # Number of hauls in the year (N for statistical calculations)
    sd_haul_avg = sd(haul_avg_weight, na.rm = TRUE), # Standard Deviation of haul averages
    # Standard Error (SE) for the 95% CI is based on the number of hauls (n_hauls)
    se_weight = sd_haul_avg / sqrt(n_hauls), 
    .groups = 'drop'
  ) %>%
  # Calculate the sliding 5-year average (moving average)
  mutate(
    sliding_5yr_avg = zoo::rollmean(avg_weight, k = 5, fill = NA, align = "center")
  )

# 3. Create the plot using the new yearly_summary_haul_avg data
community_body_weight_plot_haul_avg <- ggplot(yearly_summary_haul_avg, aes(x = year, y = avg_weight)) +
  
  # A. Shade the 95% Confidence Interval region (CI = +/- 1.96 * SE)
  geom_ribbon(
    aes(ymin = avg_weight - 1.96 * se_weight, 
        ymax = avg_weight + 1.96 * se_weight),
    fill = "lightgreen", 
    alpha = 0.4
  ) +
  
  # B. Add the main yearly average line (based on haul means) using a dashed linetype
  geom_line(linetype = "dashed", color = "black") +
  
  # C. Keep the yearly points
  geom_point(color = "black") +
  
  # D. Add the Sliding 5-Year Average line
  geom_line(aes(y = sliding_5yr_avg), 
            color = "red", 
            linewidth = 1.2) +
  
  # E. Customize Labels and Theme
  labs(
    title = "Shelf average bodyweight",
    x = "Year",
    y = "Average haul body weight (kg)" # Clarify the y-axis represents the mean of haul means
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
  )
print(community_body_weight_plot_haul_avg)

weight_haul_summary <- weight_haul_summary[, -1] # remove year column from weight_haul_summary to add it to haul stats
weight_haul_summary <- weight_haul_summary[, -3] # remove num_obs column from weight_haul_summary to add it to haul stats

SCS_haul_biodiversity <- left_join(SCS_haul_biodiversity, weight_haul_summary, by = "haul_id")

'
Below average trophic level does not incorporate uncertainty of trophic estimates, above attempt does
'

####################################################################
######## See how community average trophic level has changed #######
####################################################################

# 1. Calculate the average Trophic Level for each Haul
haul_trophic_summary <- scotian_with_troph %>%
  group_by(year, haul_id) %>%
  summarise(
    # The mean Trophic Level (Troph) of all individuals within that specific haul
    haul_avg_troph = mean(Troph, na.rm = TRUE),
    .groups = 'drop'
  )

# 2. Calculate the yearly average (of the haul means) and statistics
yearly_trophic_summary <- haul_trophic_summary %>%
  group_by(year) %>%
  summarise(
    # The new yearly average is the mean of all haul averages in that year
    avg_troph = mean(haul_avg_troph, na.rm = TRUE),
    n_hauls = n(), # Number of hauls in the year (N for statistical calculations)
    sd_haul_avg = sd(haul_avg_troph, na.rm = TRUE), # Standard Deviation of haul averages
    # Standard Error (SE) for the 95% CI is based on the number of hauls (n_hauls)
    se_troph = sd_haul_avg / sqrt(n_hauls), 
    .groups = 'drop'
  ) %>%
  # Calculate the sliding 5-year average (moving average)
  mutate(
    sliding_5yr_avg = zoo::rollmean(avg_troph, k = 5, fill = NA, align = "center")
  )

# 3. Create the plot with Confidence Intervals and Sliding Average
community_average_trophic_plot <- ggplot(yearly_trophic_summary, aes(x = year, y = avg_troph)) +
  
  # A. Shade the 95% Confidence Interval region (CI = +/- 1.96 * SE)
  geom_ribbon(
    aes(ymin = avg_troph - 1.96 * se_troph, 
        ymax = avg_troph + 1.96 * se_troph),
    fill = "lightblue", # Use a light color for the shading
    alpha = 0.4 # Make it semi-transparent
  ) +
  
  # B. Add the main yearly average line using a DASHED linetype
  geom_line(linetype = "dashed", color = "black") +
  
  # C. Keep the yearly points
  geom_point(color = "black") +
  
  # D. Add the Sliding 5-Year Average line (solid and bold trendline)
  geom_line(aes(y = sliding_5yr_avg), 
            color = "red", 
            linewidth = 1.2) +
  
  # E. Customize Labels and Theme
  labs(
    title = "Temporal changes in average fish trophic level in SCS region",
    subtitle = "95% Confidence Interval Shaded and 5-Year Sliding Average Trendline",
    x = "Year",
    y = "Average Haul Level Trophic Level"
  ) +
  theme_minimal() 

print(community_average_trophic_plot)
haul_trophic_summary <- haul_trophic_summary[, -1] # remove year column from community average trophic summary to add it to haul
SCS_haul_biodiversity <- left_join(SCS_haul_biodiversity, haul_trophic_summary, by = "haul_id")


################################################################
####### add total biomass to haul indexed dataframe ############
################################################################

# 1. Calculate the TOTAL wgt_cpua for each Haul
# The primary sampling unit statistic is the *sum* of wgt_cpua in the haul.
haul_cpua_summary <- scotian_shelf_summer_core %>%
  group_by(year, haul_id) %>%
  summarise(
    # Sum the wgt_cpua for all entries within that specific haul
    haul_total_cpua = sum(wgt_cpua, na.rm = TRUE),
    .groups = 'drop'
  )

# 2. Calculate the yearly average (of the haul totals) and statistics
yearly_cpua_summary <- haul_cpua_summary %>%
  group_by(year) %>%
  summarise(
    # The yearly average is the mean of all haul total CPUA values in that year
    avg_cpua = mean(haul_total_cpua, na.rm = TRUE),
    n_hauls = n(), # Number of hauls in the year (N for statistical calculations)
    sd_haul_avg = sd(haul_total_cpua, na.rm = TRUE), # Standard Deviation of haul totals
    # Standard Error (SE) for the 95% CI is based on the number of hauls (n_hauls)
    se_cpua = sd_haul_avg / sqrt(n_hauls), 
    .groups = 'drop'
  ) %>%
  # Calculate the sliding 5-year average (moving average)
  mutate(
    sliding_5yr_avg = zoo::rollmean(avg_cpua, k = 5, fill = NA, align = "center")
  )

# 3. Create the plot with Confidence Intervals and Sliding Average
community_cpua_plot <- ggplot(yearly_cpua_summary, aes(x = year, y = avg_cpua)) +
  
  # A. Shade the 95% Confidence Interval region (CI = +/- 1.96 * SE)
  geom_ribbon(
    aes(ymin = avg_cpua - 1.96 * se_cpua, 
        ymax = avg_cpua + 1.96 * se_cpua),
    fill = "pink", # Using a different color for variety
    alpha = 0.4
  ) +
  
  # B. Add the main yearly average line using a DASHED linetype
  geom_line(linetype = "dashed", color = "black") +
  
  # C. Keep the yearly points
  geom_point(color = "black") +
  
  # D. Add the Sliding 5-Year Average line (solid and bold trendline)
  geom_line(aes(y = sliding_5yr_avg), 
            color = "red", 
            linewidth = 1.2) +
  # E. Customize Labels and Theme
  labs(
    x = "Year",
    y = expression(paste("Yearly average CPUA (kg ", km^-2, ")"))) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
  )

print(community_cpua_plot)

# attach haul cpua data to haul summary of SCS
haul_cpua_summary <- haul_cpua_summary[, -1]
SCS_haul_biodiversity <- left_join(SCS_haul_biodiversity, haul_cpua_summary, by = "haul_id")


###########----------------------------------------------------------------------------########
########## plot Variance in haul community biomass over time, both logged and unlogged ########
###########----------------------------------------------------------------------------########

# Calculate Coefficient of Variation (Relative Volatility)
yearly_biomass_CV <- SCS_haul_biodiversity %>%
  group_by(year) %>%
  summarize(CV_biomass = (sd(haul_total_cpua, na.rm = TRUE) / mean(haul_total_cpua, na.rm = TRUE)) * 100)

# Calculate 5-year windowed CV (centered)
window_size <- 5
windowed_CV <- yearly_biomass_CV %>%
  mutate(windowed_CV = zoo::rollapply(CV_biomass, width = window_size, FUN = mean, align = "center", fill = NA))

# Fit linear model for yearly CV
CV_lm <- lm(CV_biomass ~ year, data = yearly_biomass_CV)
summary_stats <- broom::glance(CV_lm)
annot_text <- paste0(
  "R² = ", signif(summary_stats$r.squared, 3),
  ", p = ", signif(summary_stats$p.value, 3)
)

# Plot with trendline, annotation, and 5-year windowed CV
ggplot(windowed_CV, aes(x = year, y = CV_biomass)) +
  geom_line(color = "red", size = 1) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_line(aes(y = windowed_CV), color = "darkgreen", size = 1.2, linetype = "dashed") +
  annotate("text", 
           x = min(windowed_CV$year), 
           y = max(windowed_CV$CV_biomass, na.rm = TRUE), 
           label = annot_text, 
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(
    title = "Yearly and 5-Year Windowed cv in Haul Level Total Biomass",
    x = "Year",
    y = "cv of Total Biomass (CPUA)",
    caption = "Red: yearly CV biomass, Green dashed: 5-year windowed CV"
  ) +
  theme_minimal()
'
# Calculate yearly CV now using logged values
yearly_biomass_CV <- SCS_haul_biodiversity %>%
  group_by(year) %>%
  summarize(CV_biomass = (sd(log(haul_total_cpua), na.rm = TRUE) / mean(log(haul_total_cpua), na.rm = TRUE)) * 100)

# Calculate 5-year windowed CV (centered)
window_size <- 5
windowed_CV <- yearly_biomass_CV %>%
  mutate(windowed_CV = zoo::rollapply(CV_biomass, width = window_size, FUN = mean, align = "center", fill = NA))

# Fit linear model for yearly CV
CV_lm <- lm(CV_biomass ~ year, data = yearly_biomass_CV)
summary_stats <- broom::glance(CV_lm)

# Format the R^2 and p-value text
annot_text <- paste0(
  "R² = ", signif(summary_stats$r.squared, 3),
  ", p = ", signif(summary_stats$p.value, 3)
)
'




# --- Calculate Annotation Position Buffer ---

# Find the maximum and minimum y-value and the range of the biomass variance
max_y <- max(windowed_CV$CV_biomass, na.rm = TRUE)
min_y <- min(windowed_CV$CV_biomass, na.rm = TRUE)
y_range <- max_y - min_y

# Calculate the final Y position by adding a 10% buffer (increased from 5%) 
# to ensure the text is visible above the data.
y_annot_pos <- max_y + (y_range * 0.10) 
# Note: This adjustment ensures the text is placed just above the plot data boundary.

# Find the minimum x-value and the range of the year
min_x <- min(windowed_CV$year, na.rm = TRUE)
max_x <- max(windowed_CV$year, na.rm = TRUE)
x_range <- max_x - min_x

# Calculate the final X position by adding a small 1% buffer from the left edge
x_annot_pos <- min_x + (x_range * 0.01)

# --- Plot with Annotations ---

CV_community_cpua_plot <- ggplot(windowed_CV, aes(x = year, y = CV_biomass)) +
  
  # Yearly Variance (Black Dashed)
  geom_line(color = "black", size = 0.9, linetype = "dashed") +
  geom_point(size = 2) +
  
  # Trendline (Linear Model) with SE ribbon
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "pink", alpha = 0.5) +
  
  # 5-Year Windowed Variance (Red Line)
  geom_line(aes(y = windowed_CV), color = "red", size = 1.2) +
  
  # Add the R^2 and p-value annotation
  annotate("text", 
           x = x_annot_pos, # Use the buffered X position (1% from left)
           y = y_annot_pos, # Use the buffered Y position (10% from top)
           label = annot_text, 
           hjust = 0, vjust = 0, # vjust = 0 aligns the bottom of the text to y_annot_pos
           size = 4, 
           color = "black") +
  
  labs(
    x = "Year",
    y = expression(paste("CV yearly average CPUA (kg ", km^-2, ")"))
  ) +
  
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
  )


print(CV_community_cpua_plot)

############## develop table of population level change in total biomass ##############

SCS_for_density <- scotian_shelf_summer_core %>%
  filter(wgt_cpua > 0)

# 1. Create annual abundance summary
SCS_density_summary <- SCS_for_density %>%
  group_by(`accepted_name`, year) %>%
  summarize(average_wgt_cpua = mean(wgt_cpua, na.rm = TRUE), .groups = "drop")

# 2. Pivot wide and FILL implicit missing years with 0
SCS_density_summary_wide <- SCS_density_summary %>%
  pivot_wider(
    names_from = `accepted_name`, 
    values_from = average_wgt_cpua,
    values_fill = 0  # <--- This fills years where a species was completely absent
  )

# 3. Apply normalization logic
SCS_relative_density_df <- SCS_density_summary_wide %>%
  mutate(across(-year, ~ {
    # Find the first non-zero value as the baseline
    # (Since we filled NAs with 0, we only need to check for > 0)
    first_val <- .[which(. > 0)[1]] 
    
    # Divide by baseline; if species never appeared, keep as 0 or NA
    if(is.na(first_val) || first_val == 0) NA else . / first_val
  }))

# 4. Final safety check: Fill any remaining NAs resulting from the math
SCS_relative_density_df <- SCS_relative_density_df %>%
  mutate(across(-year, ~replace_na(., 0)))

## now with relative annual density info for every species, plot that for top 20 most dominant taxa ##

# 1. Identify the Top 20 Taxa by average CPUA
top_20_taxa <- SCS_density_summary_wide %>%
  # Pivot to long so we can calculate the mean per species
  pivot_longer(cols = -year, names_to = "accepted_name", values_to = "wgt_cpua") %>%
  group_by(accepted_name) %>%
  # This mean now includes the 0s for missing years
  summarize(overall_avg_density = mean(wgt_cpua, na.rm = TRUE)) %>%
  slice_max(overall_avg_density, n = 20) %>%
  pull(accepted_name)

# 2. Prepare the data for plotting
# Pivot the relative density dataframe and filter for our top taxa
plot_data <- SCS_relative_density_df %>%
  pivot_longer(cols = -year, names_to = "accepted_name", values_to = "relative_wgt") %>%
  filter(accepted_name %in% top_20_taxa)

# 3. Plot the trends
ggplot(plot_data, aes(x = year, y = relative_wgt, color = accepted_name)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  # Using facet_wrap makes it easier to see individual trends since 
  # relative scales start at 1.0 for everyone
  facet_wrap(~accepted_name, scales = "free_y", ncol = 4) + 
  scale_color_viridis_d(option = "turbo") +
  labs(
    title = "Relative annual density of 20 most dominant taxa",
    subtitle = "Relative to first year of occurrence",
    x = "Year",
    y = "Relative density (kg/km2 / Baseline)"
  ) +  
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    legend.position = "none" # Hide legend since facets are labeled
  ) 

############ make visualization for population level change in abundance ##############
"
> unique(is.na(scotian_shelf_summer_core$num_cpua))
[1] FALSE
"
SCS_for_abundance <- scotian_shelf_summer_core %>%
  filter(num_cpua > 0)

# 1. Create annual abundance summary for observed years
SCS_abundance_summary <- SCS_for_abundance %>%
  group_by(`accepted_name`, year) %>%
  summarize(average_num_cpua = mean(num_cpua, na.rm = TRUE), .groups = "drop")

# 1. Pivot wide so each species is a column
SCS_abundance_summary_wide <- SCS_abundance_summary %>%
  pivot_wider(names_from = `accepted_name`, values_from = average_num_cpua)

# 2. Apply the normalization logic
# This sets the first recorded year for each species to 1.0
SCS_relative_abundance_df <- SCS_abundance_summary_wide %>%
  mutate(across(-year, ~ {
    # Find the first non-NA, non-zero value as the baseline
    first_val <- .[which(!is.na(.) & . > 0)[1]] 
    
    # Divide by that baseline (handle cases where a species is all zeros)
    if(is.na(first_val)) NA else . / first_val
  }))

# 3. Pivot back to 'long' format (best for ggplotting later)
SCS_relative_abundance_df_long <- SCS_relative_abundance_df %>%
  pivot_longer(-year, names_to = "accepted_name", values_to = "relative_abundance")


# 2. Fill in the gaps for individual taxa during missing years
SCS_relative_abundance_summary_complete <- SCS_relative_abundance_df_long %>%
  group_by(`accepted_name`) %>%
  complete(year = min(SCS_abundance_summary$year):max(SCS_abundance_summary$year), fill = list(average_num_cpua = 0)) %>%
  ungroup()

SCS_relative_abundance_trends <- SCS_relative_abundance_summary_complete %>%
  group_by(`accepted_name`) %>%
  summarize(
    # 1. Count years where average abundance is greater than 0
    num_obs_years = sum(relative_abundance > 0),
    
    # 2. Calculate the slope of the linear regression
    # We extract the second coefficient [2], which is the 'year' slope
    total_slope = if(n() > 1) {
      coef(lm(relative_abundance ~ year))[2]
    } else {
      NA
    },
    .groups = "drop"
  )

population_relative_abundance_change_plot <-ggplot(SCS_relative_abundance_trends, aes(x = total_slope)) +
  geom_histogram(bins = 50, fill = "pink", color = "black") + 
  labs(
    title = "Population level changes in relative abundance",
    x = "Slope of Relative Abundance Linear Regression", 
    y = "Number of taxa"
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
  ) +
  coord_cartesian(xlim = c(-2.25, 2), ylim = c(-2, 175))
print(population_relative_abundance_change_plot)

# make combined biomass plot

combined_biomass_plot <- (community_cpua_plot / CV_community_cpua_plot) + 
  plot_annotation(tag_levels = 'A')

print(combined_biomass_plot)

# Save the combined plot
ggsave("combined_biomass_plot.png", combined_biomass_plot, 
       width = 10, height = 12,  # adjust dimensions as needed
       dpi = 300)

########################################################################
########### Generate size spectra for community ########################
########################################################################

SCS_haul_biodiversity_and_size_spectra <- generate_linear_size_spectra_slope(clean_body_weight_SCS)

############# Now track linear size spectra slope changing over time #################

# Function to calculate 95% Confidence Interval (No change needed here)
# CI = mean +/- t * (SD / sqrt(n))
calculate_ci <- function(x, confidence_level = 0.95) {
  n <- length(x)
  if (n < 2) { # Need at least 2 data points to calculate CI
    return(data.frame(mean = mean(x, na.rm = TRUE),
                      lower_ci = NA,
                      upper_ci = NA))
  }
  mean_val <- mean(x, na.rm = TRUE)
  std_dev <- sd(x, na.rm = TRUE)
  std_err <- std_dev / sqrt(n)
  
  # For 95% CI, alpha = 0.05. We need alpha/2 in each tail.
  alpha <- 1 - confidence_level
  t_critical <- qt(1 - alpha / 2, df = n - 1)
  
  margin_error <- t_critical * std_err
  
  lower_ci <- mean_val - margin_error
  upper_ci <- mean_val + margin_error
  
  return(data.frame(mean = mean_val,
                    lower_ci = lower_ci,
                    upper_ci = upper_ci))
}


# --- Step 1, 2, and 3: Calculate yearly averages, CI, and 5-Year Sliding Average ---
yearly_summary_size_spectrum <- SCS_haul_biodiversity_and_size_spectra %>%
  group_by(year) %>%
  summarise(
    ci_data = list(calculate_ci(body_weight_coeff)) # Apply the CI function to each group
  ) %>%
  tidyr::unnest(ci_data) %>% # Unnest the list-column into separate columns
  ungroup() %>%
  # ADDED: Calculate the sliding 5-year average (moving average)
  mutate(
    sliding_5yr_avg = zoo::rollmean(mean, k = 5, fill = NA, align = "center")
  )

print(yearly_summary_size_spectrum)

# --- Step 4: Plotting with ggplot2 (Modified) ---

community_size_spectra_plot <- ggplot(yearly_summary_size_spectrum, aes(x = year, y = mean)) +
  geom_ribbon(
    aes(ymin = lower_ci, ymax = upper_ci),
    fill = "lightgreen",
    alpha = 0.3
  ) +
  geom_line(linetype = "dashed", color = "black") + 
  geom_point(color = "black", size = 2) + 
  geom_line(aes(y = sliding_5yr_avg), 
            color = "red", 
            linewidth = 1.2) +
  labs(
    title = "Shelf average linear size spectrum coefficient",
    x = "Year",
    y = "Average log body weight coefficient: kappa"
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black")
  ) + 
  coord_cartesian(
    ylim = c(min(yearly_summary_size_spectrum$mean) - 0.2,
             max(yearly_summary_size_spectrum$mean) + 0.2)
  )

print(community_size_spectra_plot)

# note strange reversal of body weight influence in 1995, maybe cod moratorium? and high variability in size spectrum coefficient in 1976

#########################################################################
###### Systematic review of body weight of all taxa over time ###########
#########################################################################
########### Something whacky going on with plot and relative body condition ##########

weight_listed_taxa <- unique(CPUA_body_weight_SCS$accepted_name)
year_list <- sort(unique(CPUA_body_weight_SCS$year))

body_condition_matrix <- c(year_list)
for (taxon in weight_listed_taxa){
  body_condition_over_time <- body_condition_track(CPUA_body_weight_SCS, taxon)
  body_condition_series <- body_condition_over_time$body_condition_over_time
  body_condition_matrix <- cbind(body_condition_matrix, body_condition_series)
}

body_condition_df <- data.frame(body_condition_matrix) # body weight CPUA of each taxon over time
names(body_condition_df) <- c("Year", weight_listed_taxa)

relative_body_condition_df <- body_condition_df %>%
  mutate(across(-Year, ~ . / .[which(!is.na(.))[1]]))

matplot(body_condition_df[, 1], body_condition_df[, 2:ncol(body_condition_df)], 
        type = "l",  # "l" for lines
        xlab = names(body_condition_df)[1],
        ylab = "Body weight (kg)",
        main = "Body weight over time of all taxa on SCS",
        #ylim = c(0,5),
        col = 1:(ncol(body_condition_df)-1), # Different colors for each line
        lty = 1:(ncol(body_condition_df)-1), # Different line types for each line
        lwd = 1.5 # Line width
)
legend("topright", legend = names(body_condition_df)[2:ncol(body_condition_df)], col = 1:(ncol(body_condition_df)-1), lty = 1:(ncol(body_condition_df)-1), lwd = 1.5) # Add a legend

matplot(relative_body_condition_df[, 1], relative_body_condition_df[, 2:ncol(relative_body_condition_df)], 
        type = "l",  # "l" for lines
        xlab = names(relative_body_condition_df)[1],
        ylab = "Body weight (kg) relative to 1970 bodyweight",
        main = "Relative body weight over time of all taxa on SCS",
        #ylim = c(0,5),
        col = 1:(ncol(relative_body_condition_df)-1), # Different colors for each line
        lty = 1:(ncol(relative_body_condition_df)-1), # Different line types for each line
        lwd = 1.5 # Line width
)
legend("topright", legend = names(relative_body_condition_df)[2:ncol(relative_body_condition_df)], col = 1:(ncol(relative_body_condition_df)-1), lty = 1:(ncol(relative_body_condition_df)-1), lwd = 1.5) # Add a legend

max_diffs <- sapply(body_condition_df, function(x) {
  if (all(is.na(x))) { # Check if all values are NA
    NA # return NA if all are NA
  } else {
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE) # Calculate range
  }
})
max_diffs <- max_diffs[-1]

max_change_in_body_condition <- as.data.frame(t(max_diffs))
max_col_indices <- apply(max_change_in_body_condition, 1, which.max)  # Find indices of max values in each row
max_cols <- colnames(max_change_in_body_condition)[max_col_indices]   # Get column names
print(max_cols)
# swordfish changed the most over time!!!

plot(body_condition_df$Year, body_condition_df$`Gadus morhua`, main = "Average body weight of Cod over time on SCS", xlab = "Year", ylab = "Body weight (kg)")

max_diffs <- sapply(relative_body_condition_df, function(x) {
  if (all(is.na(x))) { # Check if all values are NA
    NA # return NA if all are NA
  } else {
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE) # Calculate range
  }
})
max_diffs <- max_diffs[-1]

max_change_in_relative_body_condition <- as.data.frame(t(max_diffs))
max_col_indices <- apply(max_change_in_relative_body_condition, 1, which.max)  # Find indices of max values in each row
max_cols <- colnames(max_change_in_body_condition)[max_col_indices]   # Get column names
print(max_cols)
# chiasmodon changed the most over time!!!


# body weight models
MbodyweightFull <- glm(body_weight ~ year + gear + quarter + latitude + longitude + depth + sst*sbt, family = Gamma(link = "log"), data = clean_body_weight_SCS)
summary(MbodyweightFull)
MbodyweightStep <- stepAIC(MbodyweightFull, direction = "backward")
summary(MbodyweightStep)
MbodyweightYear <- glm(body_weight ~ year, family = Gamma(link = "log"), data = clean_body_weight_SCS)
summary(MbodyweightYear)


############ make visualization for population level change in relative body weight ##############
# Identify species columns (everything except 'Year')
species_names <- setdiff(names(relative_body_condition_df), "Year")
# Apply a function to each column
slopes <- sapply(species_names, function(sp) {
  model <- lm(relative_body_condition_df[[sp]] ~ relative_body_condition_df$Year)
  coef(model)["relative_body_condition_df$Year"] # Extract the slope
})
# Convert to a named data frame
rates_of_relative_bodyweight_change <- data.frame(Species = species_names, Slope = as.numeric(slopes))
population_bodyweight_change_plot <-ggplot(rates_of_relative_bodyweight_change, aes(x = Slope)) +
  geom_histogram(bins = 50, fill = "lightgreen", color = "black") + 
  labs(
    title = "Population level changes in relative body weight",
    x = "Slope of relative body condition linear regression", 
    y = "Number of taxa"
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
  )+
  coord_cartesian(xlim = c(-2.5, 11), ylim = c(-2, 200))
print(population_bodyweight_change_plot)


################################################
####### now investigate taxa of interest #######
################################################

#taxa_of_interest <- c("Gadus morhua", "Merluccius bilinearis", "Sebastes", "Melanogrammus aeglefinus", "Pollachius virens", "Squalus acanthius")
'
clean_scotian_shelf_year_round <- scotian_shelf_year_round[!(is.na(scotian_shelf_year_round$year) | is.na(scotian_shelf_year_round$latitude) | is.na(scotian_shelf_year_round$longitude) | is.na(scotian_shelf_year_round$season) | is.na(scotian_shelf_year_round$gear) | is.na(scotian_shelf_year_round$depth) | is.na(scotian_shelf_year_round$sst) | is.na(scotian_shelf_year_round$sst) | is.na(scotian_shelf_year_round$sbt)), ]

"Gadus morhua"
cod_SCS <- clean_scotian_shelf_year_round %>% # find all hauls which catch a non-negligible amount of cod
  filter(accepted_name == "Gadus morhua", num > 0, wgt > 0)
attach(cod_SCS)

Cod_Mlat <- glm(cod_SCS$body_weight ~ cod_SCS$latitude) # haul average body weight as a function of latitude
summary(Cod_Mlat)

Cod_Myeargear <- glm(cod_SCS$body_weight ~ cod_SCS$year + cod_SCS$gear)
summary(Cod_Myeargear)

Cod_Mfull <- glm(cod_SCS$body_weight~cod_SCS$year + cod_SCS$latitude + cod_SCS$longitude + cod_SCS$season + cod_SCS$gear + cod_SCS$depth + cod_SCS$sst*cod_SCS$sbt) # expect body weight to follow relatively normal gaussian distribution
summary(Cod_Mfull)

Cod_Mstep <- stepAIC(Cod_Mfull, direction = "backward")
summary(Cod_Mstep)

example_species <- weight_listed_taxa[1]
cod_over_time <- body_condition_track(CPUA_body_weight_SCS, example_species)

# Calculate the standard error
cod_over_time$se <- cod_over_time$body_condition_over_time / sqrt(cod_over_time$n_observations)

# Create the scatterplot with error bars
cod_body_weight_plot <- ggplot(cod_over_time, aes(x = year, y = body_condition_over_time)) +
  geom_point() +
  geom_errorbar(aes(ymin = body_condition_over_time - se, ymax = body_condition_over_time + se),
                width = 0.2) + # Adjust width of error bars as needed
  labs(title = "Yearly Average Body Weight of Cod Over Time",
       x = "Year",
       y = "Body Weight (kg)") +
  theme_minimal() # A clean and modern theme
print(cod_body_weight_plot)

"Merluccius bilinearis"

hake_SCS <- scotian_shelf_year_round %>% # find all hauls which catch a non-negligible amount of hake
  filter(accepted_name == "Merluccius bilinearis", num > 0, wgt > 0)

attach(hake_SCS)
Hake_Mlat <- glm(body_weight ~ latitude) # haul average body weight as a function of latitude
summary(Hake_Mlat)

Hake_Myeargear <- glm(hake_SCS$body_weight ~ hake_SCS$year + hake_SCS$gear)
summary(Hake_Myeargear)

Hake_Mfull <- glm(hake_SCS$body_weight~hake_SCS$year + hake_SCS$latitude + hake_SCS$longitude + hake_SCS$season + hake_SCS$gear + hake_SCS$depth + hake_SCS$sst*hake_SCS$sbt) # expect body weight to follow relatively normal gaussian distribution
summary(Hake_Mfull)

Hake_Mstep <- stepAIC(Hake_Mfull, direction = "backward")
summary(Hake_Mstep)

example_species <- "Merluccius bilinearis"
hake_over_time <- body_condition_track(CPUA_body_weight_SCS, example_species)

# Calculate the standard error
hake_over_time$se <- hake_over_time$body_condition_over_time / sqrt(hake_over_time$n_observations)

# Create the scatterplot with error bars
hake_body_weight_plot <- ggplot(hake_over_time, aes(x = year, y = body_condition_over_time)) +
  geom_point() +
  geom_errorbar(aes(ymin = body_condition_over_time - se, ymax = body_condition_over_time + se),
                width = 0.2) + # Adjust width of error bars as needed
  labs(title = "Yearly Average Body Weight of Hake Over Time",
       x = "Year",
       y = "Body Weight (kg)") +
  theme_minimal() 
print(hake_body_weight_plot)
# notice outlier in body weight making it difficult to analyze trends, reproduce graph with outlier removed
z_scores <- scale(hake_SCS$body_weight)
outliers <- which(abs(z_scores) > 3)
hake_SCS <- hake_SCS[-outliers,]
boxplot(body_weight ~ year, data=hake_SCS, main = "Body weight of Hake over time, outlier removed", xlab = "Year", ylab = "Average body weight in haul (kg)")
# evidence of outlier in high weight of hake in 1998

"Sebastes"

rockfish_SCS <- scotian_shelf_year_round %>% # find all hauls which catch a non-negligible amount of cod
  filter(accepted_name == "Sebastes", num > 0, wgt > 0)
rockfish_SCS$body_weight <- rockfish_SCS$wgt / rockfish_SCS$num
attach(rockfish_SCS)

Rockfish_Mlat <- glm(body_weight ~ latitude) # haul average body weight as a function of latitude
summary(Rockfish_Mlat)

Rockfish_Myeargear <- glm(body_weight ~ year + gear, data = rockfish_SCS)
summary(Rockfish_Myeargear)

Rockfish_Mfull <- glm(body_weight~year + latitude + longitude + season + gear + depth + sst*sbt, data = rockfish_SCS) # expect body weight to follow relatively normal gaussian distribution
summary(Rockfish_Mfull)

Rockfish_Mstep <- stepAIC(Rockfish_Mfull, direction = "forward")
summary(Rockfish_Mstep)

example_species <- "Sebastes"
rockfish_over_time <- body_condition_track(CPUA_body_weight_SCS, example_species)

# Calculate the standard error
rockfish_over_time$se <- rockfish_over_time$body_condition_over_time / sqrt(rockfish_over_time$n_observations)

# Create the scatterplot with error bars
rockfish_body_weight_plot <- ggplot(rockfish_over_time, aes(x = year, y = body_condition_over_time)) +
  geom_point() +
  geom_errorbar(aes(ymin = body_condition_over_time - se, ymax = body_condition_over_time + se),
                width = 0.2) + # Adjust width of error bars as needed
  labs(title = "Yearly Average Body Weight of Rockfish Over Time",
       x = "Year",
       y = "Body Weight (kg)") +
  theme_minimal() 
print(rockfish_body_weight_plot)

"Melanogrammus aeglefinus"

haddock_SCS <- scotian_shelf_year_round %>% # find all hauls which catch a non-negligible amount of cod
  filter(accepted_name == "Melanogrammus aeglefinus", num > 0, wgt > 0)
haddock_SCS$body_weight <- haddock_SCS$wgt / haddock_SCS$num
attach(haddock_SCS)

Haddock_Mlat <- glm(body_weight ~ latitude) # haul average body weight as a function of latitude
summary(Haddock_Mlat)

Haddock_Myeargear <- glm(body_weight ~ year + gear, data = haddock_SCS)
summary(Haddock_Myeargear)

Haddock_Mfull <- glm(body_weight~year + latitude + longitude + season + gear + depth + sst*sbt, data = haddock_SCS) # expect body weight to follow relatively normal gaussian distribution
summary(Haddock_Mfull)

Haddock_Mstep <- stepAIC(Haddock_Mfull, direction = "forward")
summary(Haddock_Mstep)

example_species <- "Melanogrammus aeglefinus"
haddock_over_time <- body_condition_track(CPUA_body_weight_SCS, example_species)

# Calculate the standard error
haddock_over_time$se <- haddock_over_time$body_condition_over_time / sqrt(haddock_over_time$n_observations)

# Create the scatterplot with error bars
haddock_body_weight_plot <- ggplot(haddock_over_time, aes(x = year, y = body_condition_over_time)) +
  geom_point() +
  geom_errorbar(aes(ymin = body_condition_over_time - se, ymax = body_condition_over_time + se),
                width = 0.2) + # Adjust width of error bars as needed
  labs(title = "Yearly Average Body Weight of Haddock Over Time",
       x = "Year",
       y = "Body Weight (kg)") +
  theme_minimal() 
print(haddock_body_weight_plot)

"Pollachius virens"

pollock_SCS <- scotian_shelf_year_round %>% # find all hauls which catch a non-negligible amount of cod
  filter(accepted_name == "Pollachius virens", num > 0, wgt > 0)
pollock_SCS$body_weight <- pollock_SCS$wgt / pollock_SCS$num
attach(pollock_SCS)

Pollock_Mlat <- glm(body_weight ~ latitude) # haul average body weight as a function of latitude
summary(Pollock_Mlat)

Pollock_Myeargear <- glm(body_weight ~ year + gear, data = pollock_SCS)
summary(Pollock_Myeargear)

Pollock_Mfull <- glm(body_weight~year + latitude + longitude + season + gear + depth + sst*sbt, data = pollock_SCS) # expect body weight to follow relatively normal gaussian distribution
summary(Pollock_Mfull)

Pollock_Mstep <- stepAIC(Pollock_Mfull, direction = "forward")
summary(Pollock_Mstep)

example_species <- "Pollachius virens"
pollock_over_time <- body_condition_track(CPUA_body_weight_SCS, example_species)

# Calculate the standard error
pollock_over_time$se <- pollock_over_time$body_condition_over_time / sqrt(pollock_over_time$n_observations)

# Create the scatterplot with error bars
pollock_body_weight_plot <- ggplot(pollock_over_time, aes(x = year, y = body_condition_over_time)) +
  geom_point() +
  geom_errorbar(aes(ymin = body_condition_over_time - se, ymax = body_condition_over_time + se),
                width = 0.2) + # Adjust width of error bars as needed
  labs(title = "Yearly Average Body Weight of Pollock Over Time",
       x = "Year",
       y = "Body Weight (kg)") +
  theme_minimal() 
print(pollock_body_weight_plot)

"Squalus acanthius" ########### Issue come back to later!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

spurdog_SCS <- scotian_shelf_year_round %>% # find all hauls which catch a non-negligible amount of cod
  filter(accepted_name == "Squalus acanthias")
spurdog_SCS$body_weight <- spurdog_SCS$wgt / spurdog_SCS$num
attach(spurdog_SCS)

'

######### FIX SPURDOG ####################
'
Spurdog_Mlat <- glm(body_weight ~ latitude, data = spurdog_SCS) # haul average body weight as a function of latitude
summary(Spurdog_Mlat)

Spurdog_Myeargear <- glm(body_weight ~ year + gear, data = spurdog_SCS)
summary(Spurdog_Myeargear)

Spurdog_Mfull <- glm(body_weight~year + latitude + longitude + season + gear + depth + sst*sbt, data = spurdog_SCS) # expect body weight to follow relatively normal gaussian distribution
summary(Spurdog_Mfull)

Spurdog_Mstep <- stepAIC(Spurdog_Mfull, direction = "forward")
summary(Spurdog_Mstep)

example_species <- "Squalus acanthius"
spurdog_over_time <- body_condition_track(CPUA_body_weight_SCS, example_species)

# Calculate the standard error
spurdog_over_time$se <- spurdog_over_time$body_condition_over_time / sqrt(spurdog_over_time$n_observations)

# Create the scatterplot with error bars
ggplot(spurdog_over_time, aes(x = year, y = body_condition_over_time)) +
  geom_point() +
  geom_errorbar(aes(ymin = body_condition_over_time - se, ymax = body_condition_over_time + se),
                width = 0.2) + # Adjust width of error bars as needed
  labs(title = "Yearly Average Body Weight of Spurdog Over Time",
       x = "Year",
       y = "Body Weight (kg)") +
  theme_minimal() 

################# Facet wrap six body weight plots (5 species, 1 community) #####################

# Now, combine them using patchwork
combined_plot <- community_body_weight_plot + cod_body_weight_plot + hake_body_weight_plot + rockfish_body_weight_plot + haddock_body_weight_plot + pollock_body_weight_plot

# To arrange them in a 2TL4 grid (or any arrangement you prefer):
final_plot <- combined_plot + plot_layout(ncol = 3, nrow = 2)

# To view the combined plot:
print(final_plot)

'


##############################################
######## Section 4: Subregions ###############
##############################################


##################################################
##### add sampling info back into haul stats #####
##################################################
compact_SCS <- scotian_shelf_summer_core %>%
  # 1. Select the columns you want to keep: 
  #    The unique identifier (haul_id) and all constant environmental variables.
  #    Replace 'env_var1', 'env_var2', etc., with your actual column names.
  select(haul_id, month, day, area_swept, gear, depth, sbt, sst, latitude, longitude, Nafo_division) %>%
  
  # 2. Use distinct() to keep only the unique combinations of these columns.
  #    Since the environmental variables are constant for each haul_id, 
  #    this effectively leaves one row per haul.
  distinct()

# The resulting dataframe 'compact_dataframe' will now have one row per unique haul_id,
# containing only the constant environmental variables.

SCS_haul_stats <- left_join(compact_SCS, SCS_haul_biodiversity_and_size_spectra)

'
# now haul average biomass by nafo division

average_biomass_by_nafo <-  SCS_haul_biodiversity_and_size_spectra %>% 
  group_by(year, Nafo_division) %>%
  summarise(
    avg_biomass = mean(total_biomass_cpua, na.rm = TRUE),
  )

ggplot(average_biomass_by_nafo, aes(x = year, y = avg_biomass, color = Nafo_division)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Nafo_division, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Average Haul Biomass by Year for Each NAFO Division",
    x = "Year",
    y = "Yearly Average Total Haul Biomass (wgt_cpua)"
  ) +
  theme_minimal()
'

#################################################################
######### Standardized graphing of prelim results ###############
#################################################################

####### Yearly single stats across region #######################
# average haul level: trophic level, k, biomass, evenness

vars <- c("evenness", "haul_avg_troph", "haul_total_cpua", "body_weight_coeff")

# Custom axis titles for each variable
x_axis_titles <- c(
  evenness = "Yearly average (Haul Evenness)",
  avg_troph = "Yearly average (haul average (Trophic Level))",
  total_biomass_cpua = "Yearly average (Haul Total Biomass CPUA)",
  body_weight_coeff = "Year average (Haul Size Spectrum Coefficient)"
)

# Calculate yearly means and confidence intervals
summary_df <- SCS_haul_stats %>%
  pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value") %>%
  group_by(year, variable) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    n = sum(!is.na(value)),
    se = sd(value, na.rm = TRUE) / sqrt(n),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop"
  ) %>%
  arrange(variable, year) %>%
  group_by(variable) %>%
  mutate(
    mean_5yr = zoo::rollapply(mean, width = 5, FUN = mean, align = "center", fill = NA, na.rm = TRUE)
  ) %>%
  ungroup()

# Custom labeller function
facet_labeller <- as_labeller(x_axis_titles)

# Plot
ggplot(summary_df, aes(x = year, y = mean, color = variable, fill = variable)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = mean_5yr), color = "black", size = 1.2, linetype = "solid", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free_y", labeller = facet_labeller) +
  labs(
    title = "Yearly Average and 95% CI for SCS Haul Stats\nwith 5-Year Centered Moving Average",
    y = "Mean (±95% CI)"
  ) +
  theme_minimal() +
  theme(
    strip.text.x = element_text(size = 12)
  )

# Custom labeller function
facet_labeller <- as_labeller(x_axis_titles)

# Plot
ggplot(summary_df, aes(x = year, y = mean, color = variable, fill = variable)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = mean_5yr), color = "black", size = 1.2, linetype = "solid", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free_y", labeller = facet_labeller) +
  labs(
    title = "Summertime Yearly Average and 95% CI for SCS Haul Stats\nwith 5-Year Centered Moving Average",
    y = "Mean (±95% CI)"
  ) +
  theme_minimal() +
  theme(
    strip.text.x = element_text(size = 12)
  )

# distinct metrics: average fish bodyweight across year community, average main NAFO division dissimilarity, average error in stability approximation
#SCS_dissimilarity_by_NAFO <- dissimilarity_matrices_by_NAFO_division(scotian_shelf_summer_core)

######################################################
## Community dissimilarity over time in core NAFO #######
######################################################
'
core_NAFO_divisions <- c("4X", "4Vn", "4Vs", "4W")
SCS_dissimilarity_between_years <- dissimilarity_matrices_for_shelf_community(scotian_shelf_summer_core)
PCoA_plot <- ggplot2::last_plot()

axes_communities <- enframe(SCS_dissimilarity_between_years$axes_cor_species, name = "Taxa", value = "Correlation")

# 2. Join with your core data to get aphia_ids
# This replaces the entire for loop
axes_communities_with_ids <- axes_communities %>%
  left_join(
    scotian_shelf_summer_core %>% 
      select(accepted_name, aphia_id) %>% 
      distinct(), 
    by = c("Taxa" = "accepted_name")
  )
my_ids <- axes_communities_with_ids$aphia_id
common_names_list <- wm_common_id_(id = my_ids)
english_names <- common_names_list %>%
  bind_rows(.id = "aphia_id") %>%           # Flatten list to a dataframe
  filter(language == "English") %>%        # Filter for English
  group_by(aphia_id) 

unique_english_names <- english_names %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup()

unique_english_names <- unique_english_names %>%
  select(c(id, vernacular))
unique_english_names <- unique_english_names %>%
  mutate(id = as.numeric(id))

axes_communities_with_ids <- axes_communities_with_ids %>%
  left_join(unique_english_names, by = c("aphia_id" = "id"))

axes_communities_with_ids <- axes_communities_with_ids %>%
  mutate(vernacular = case_when(
    Taxa == "Sebastes" ~ "Redfish",
    TRUE ~ vernacular # This keeps the original value for everything else
  ))

# 1. Your list of search terms
search_list <- unique(axes_communities_with_ids$Taxa)

# 2. Collapse the list into a regex pattern separated by | (the operator)
pattern <- paste(search_list, collapse = "|")

# 3. Filter the dataframe
filtered_rates_of_relative_bodyweight_change <- rates_of_relative_bodyweight_change %>%
  filter(str_detect(Species, pattern)) %>%
  # Extract the specific word that matched and put it in a new column
  mutate(Taxa = str_extract(Species, pattern))

filtered_rates_of_relative_bodyweight_change <- filtered_rates_of_relative_bodyweight_change %>%
  select(c(Taxa, Slope))

axes_communities_with_ids <- left_join(axes_communities_with_ids, filtered_rates_of_relative_bodyweight_change, by = "Taxa")

####### add total bodyweight average to this axis community ##########
species_total_avg_CPUA <- yearly_average_community %>%
  # 1. Pivot from wide to long (gather species columns)
  pivot_longer(
    cols = -Time, 
    names_to = "Species", 
    values_to = "CPUA"
  ) %>%
  # 2. Group by Species to calculate the overall average
  group_by(Species) %>%
  summarize(Average_CPUA = mean(CPUA, na.rm = TRUE)) %>%
  ungroup()

axis_community_bodyweights <- species_total_avg_CPUA %>%
  filter(Species %in% search_list)

axes_communities_with_ids <- left_join(axes_communities_with_ids, axis_community_bodyweights, by = c("Taxa" = "Species"))
axes_communities_with_ids$Axis <- 0
axes_communities_with_ids[1:6,]$Axis <- 1
axes_communities_with_ids[6:10,]$Axis <- 2
axes_communities_with_ids <- axes_communities_with_ids[, c("Taxa", "vernacular", "aphia_id", "Axis", "Correlation", "Average_CPUA", "Slope")]
'


################# Attempt community association by clusters instead of axes ##############

# first note optimal clustering is 2 with 1972 as lone year or 3 without
SCS_dissimilarity_between_years_with_optimal_k <- dissimilarity_matrices_for_shelf_community_with_optimal_k(scotian_shelf_summer_core)
# cut 1972 and check if things change
scotian_shelf_summer_core_without_lone_year <- scotian_shelf_summer_core %>%
  filter(year != 1972)
dissimilarity_matrices_for_shelf_community_with_optimal_k <- dissimilarity_matrices_for_shelf_community_with_optimal_k(scotian_shelf_summer_core_without_lone_year)
# now apparent that 2 most important clusters are before and after collapse after excluding 1972
SCS_dissimilarity_between_years <- dissimilarity_analysis_of_yearly_communities(scotian_shelf_summer_core_without_lone_year)
PCoA_plot <- ggplot2::last_plot()


# Now characterize the PCoA cluster communities 


# 1. Extract the Mean Abundance table
dominant_cluster_comm <- as.data.frame(SCS_dissimilarity_between_years$mean_abundance)
dominant_cluster_comm$Representation_class <- "Dominant"
dominant_cluster_comm$p_val <- NA

# 2. Extract the Indicator Species table
indicator_cluster_comm <- as.data.frame(SCS_dissimilarity_between_years$indicators)
indicator_cluster_comm$Representation_class <- "Indicator"

# merge communities 
indicator_cluster_comm <- indicator_cluster_comm %>%
  select(c(Cluster, Species, p_val, Representation_class))
dominant_cluster_comm <- dominant_cluster_comm %>%
  select(c(Cluster, Species, p_val, Representation_class))
PCoA_cluster_communities_df <- bind_rows(dominant_cluster_comm, indicator_cluster_comm)
PCoA_cluster_communities_df <- PCoA_cluster_communities_df %>%
  rename(Taxa = Species)

# 2. Join with your core data to get aphia_ids
# This replaces the entire for loop
cluster_communities_with_ids <- PCoA_cluster_communities_df %>%
  left_join(
    scotian_shelf_summer_core %>% 
      select(accepted_name, aphia_id) %>% 
      distinct(), 
    by = c("Taxa" = "accepted_name")
  )

my_ids <- cluster_communities_with_ids$aphia_id
common_names_list <- wm_common_id_(id = my_ids)
english_names <- common_names_list %>%
  bind_rows(.id = "aphia_id") %>%           # Flatten list to a dataframe
  filter(language == "English") %>%        # Filter for English
  group_by(aphia_id) 

unique_english_names <- english_names %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup()

unique_english_names <- unique_english_names %>%
  select(c(id, vernacular))
unique_english_names <- unique_english_names %>%
  mutate(id = as.numeric(id))

cluster_communities_with_ids <- cluster_communities_with_ids %>%
  left_join(unique_english_names, by = c("aphia_id" = "id"))

cluster_communities_with_ids <- cluster_communities_with_ids %>%
  mutate(vernacular = case_when(
    Taxa == "Sebastes" ~ "redfish",
    Taxa == "Myzopsetta ferruginea" ~ "yellowtail flounder",
    Taxa == "Gaidropsarus" ~ "rockling",
    Taxa == "Lycodes" ~ "eelpouts",
    Taxa == "Gadiformes" ~ "anacanthini",
    TRUE ~ vernacular # This keeps the original value for everything else
  ))

# 1. Your list of search terms
search_list <- unique(cluster_communities_with_ids$Taxa)

# 3. Filter the dataframe
filtered_rates_of_relative_bodyweight_change <- rates_of_relative_bodyweight_change %>%
  filter(Species %in% search_list)

cluster_communities_with_ids <- left_join(
  cluster_communities_with_ids, 
  filtered_rates_of_relative_bodyweight_change, 
  by = c("Taxa" = "Species")
)

####### add total density average to this cluster community ##########
species_total_avg_CPUA <- yearly_average_community %>%
  # 1. Pivot from wide to long (gather species columns)
  pivot_longer(
    cols = -Time, 
    names_to = "Species", 
    values_to = "CPUA"
  ) %>%
  # 2. Group by Species to calculate the overall average
  group_by(Species) %>%
  summarize(Average_CPUA = mean(CPUA, na.rm = TRUE)) %>%
  ungroup()

cluster_community_densities <- species_total_avg_CPUA %>%
  filter(Species %in% search_list)

cluster_communities_with_ids <- left_join(cluster_communities_with_ids, cluster_community_densities, by = c("Taxa" = "Species"))
cluster_communities_with_ids <- cluster_communities_with_ids[, c("Cluster", "Taxa", "vernacular", "aphia_id", "Representation_class", "p_val", "Average_CPUA", "Slope")]
cluster_communities_with_ids <- cluster_communities_with_ids %>%
  rename(Vernacular = vernacular)

### format into MEPS compatible table ####
colnames(cluster_communities_with_ids)[colnames(cluster_communities_with_ids) == "aphia_id"] <- "aphia ID"
colnames(cluster_communities_with_ids)[colnames(cluster_communities_with_ids) == "Representation_class"] <- "Representation class"
colnames(cluster_communities_with_ids)[colnames(cluster_communities_with_ids) == "p_val"] <- "p value"
colnames(cluster_communities_with_ids)[colnames(cluster_communities_with_ids) == "Average_CPUA"] <- "Avg CPUA"
colnames(cluster_communities_with_ids)[colnames(cluster_communities_with_ids) == "Slope"] <- "Wgt slope"
colnames(cluster_communities_with_ids)[colnames(cluster_communities_with_ids) == "Taxa"] <- "Taxon"
cluster_communities_with_ids$`Avg CPUA` <- signif(cluster_communities_with_ids$`Avg CPUA`, digits = 4)
cluster_communities_with_ids$`Wgt slope` <- signif(cluster_communities_with_ids$`Wgt slope`, digits = 3)

######## add abundance trends ######
unique_names <- unique(cluster_communities_with_ids$Taxon)

cluster_comm_abundances <- SCS_relative_abundance_trends %>%
  filter(accepted_name %in% unique_names) %>%
  select(accepted_name, total_slope)

cluster_communities_with_ids <- left_join(cluster_communities_with_ids, cluster_comm_abundances, by = c("Taxon" = "accepted_name"))
colnames(cluster_communities_with_ids)[colnames(cluster_communities_with_ids) == "total_slope"] <- "Num slope"

##### format for meps ########

cluster_community_table <- cluster_communities_with_ids %>% 
  meps_theme()
print(cluster_community_table)

#################################################################################################
######## Community cluster table V3 #############################################################
#################################################################################################
# Add before and after average weight trends to table instead of body weight and abundance trends

cluster_communities_with_ids_V3 <- cluster_communities_with_ids %>%
  select(c(Cluster, Taxon, Vernacular, `aphia ID`, `Representation class`, `p value`)) 

search_list <- unique(cluster_communities_with_ids_V3$Taxon)

pre90_yearly_average_community <- yearly_average_community %>%
  filter(Time <= 1990)
pre90_avg_CPUA <- pre90_yearly_average_community %>%
  # 1. Pivot from wide to long (gather species columns)
  pivot_longer(
    cols = -Time, 
    names_to = "Species", 
    values_to = "CPUA"
  ) %>%
  # 2. Group by Species to calculate the overall average
  group_by(Species) %>%
  summarize(Average_CPUA = mean(CPUA, na.rm = TRUE)) %>%
  ungroup()

pre90_cluster_community_densities <- pre90_avg_CPUA %>%
  filter(Species %in% search_list)

# now do post 90

post90_yearly_average_community <- yearly_average_community %>%
  filter(Time >= 1990)
post90_avg_CPUA <- post90_yearly_average_community %>%
  # 1. Pivot from wide to long (gather species columns)
  pivot_longer(
    cols = -Time, 
    names_to = "Species", 
    values_to = "CPUA"
  ) %>%
  # 2. Group by Species to calculate the overall average
  group_by(Species) %>%
  summarize(Average_CPUA = mean(CPUA, na.rm = TRUE)) %>%
  ungroup()

post90_cluster_community_densities <- post90_avg_CPUA %>%
  filter(Species %in% search_list)

cluster_communities_with_ids_V3$`Avg CPUA in cluster` <- NA

library(dplyr)

cluster_communities_with_ids_V3 <- cluster_communities_with_ids_V3 %>%
  rowwise() %>%
  mutate(
    `Avg CPUA in cluster` = if (Cluster == 1) {
      post90_cluster_community_densities %>%
        filter(Species == Taxon) %>%
        pull(Average_CPUA)
    } else {
      pre90_cluster_community_densities %>%
        filter(Species == Taxon) %>%
        pull(Average_CPUA)
    }
  )
cluster_communities_with_ids_V3 <- cluster_communities_with_ids_V3 %>%
  filter(`p value` <= 0.05 | is.na(`p value`))
cluster_communities_with_ids_V3 <- cluster_communities_with_ids_V3 %>%
  select(-`p value`)

cluster_community_table <- cluster_communities_with_ids_V3 
print(cluster_community_table)


############ combined community structure figure ##############
combined_comm_struc_plot <- (community_body_weight_plot_haul_avg + community_size_spectra_plot + proportional_change_plot) / (MTC_plot + PCoA_plot) +
  plot_annotation(tag_levels = 'A') 
print(combined_comm_struc_plot)

############################################################################
##### Make panelled compensation / synchrony figure for Article 1 ##########
############################################################################
#view(rare_comp_indices)
#view(rare_sync_indices)
#view(spec_comp_indices)
#view(spec_sync_indices)
#view(trophic_comp_sync_indices)
Synchrony_df <- data.frame(Window = spec_sync_indices$Window, Interspecific = spec_sync_indices$SynchIndex, Trophic = trophic_comp_sync_indices$TrophicSynchrony, Rarity = rare_sync_indices$SynchIndex)
Compensation_df <- data.frame(Window = spec_comp_indices$Window, Interspecific = spec_comp_indices$CompensationIndex, Trophic = trophic_comp_sync_indices$TrophicCompensation, Rarity = rare_comp_indices$CompensationIndex)

# --- 1. Combine and Prepare the Data (Same as before) ---
Synchrony_df$Type <- "Synchrony"
Compensation_df$Type <- "Compensation"
combined_df <- bind_rows(Synchrony_df, Compensation_df)

long_df <- combined_df %>%
  pivot_longer(
    cols = c(Interspecific, Trophic, Rarity), 
    names_to = "Index_Source",               
    values_to = "Value"                      
  )

library(ggplot2)
library(patchwork)

# 1. Create the Compensation Plot
p_comp <- ggplot(subset(long_df, Type == "Compensation"), 
                 aes(x = Window, y = Value, color = Index_Source)) +
  geom_line(aes(group = Index_Source), linewidth = 1, linetype = "dashed") +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(-0.25, 1)) +
  scale_color_manual(values = c(
    "Interspecific" = "lightgreen", 
    "Trophic" = "orange", 
    "Rarity" = "cornflowerblue"
  )) +
  labs(x = NULL, y = "Compensation Index", color = "Ecological unit") + 
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_blank(), # Recommended for MEPS,
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )

# 2. Create the Synchrony Plot
p_sync <- ggplot(subset(long_df, Type == "Synchrony"), 
                 aes(x = Window, y = Value, color = Index_Source)) +
  geom_line(aes(group = Index_Source), linewidth = 1, linetype = "dashed") +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(-0.25, 1)) +
  scale_color_manual(values = c(
    "Interspecific" = "lightgreen", 
    "Trophic" = "orange", 
    "Rarity" = "cornflowerblue"
  )) +
  labs(x = "Year Window", y = "Synchrony Index", color = "Ecological unit") + 
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
  )

# 3. Combine them with Patchwork
# The "/" stacks them vertically
comp_sync_index_plot <- p_comp / p_sync  + 
  plot_annotation(tag_levels = 'A')

print(comp_sync_index_plot)

###############################################################
#### Dynamical relationships between community indicators #####
###############################################################

########## Covariance and correlation between haul level statistics  #########

simplified_mtc <- haul_mtc_results[, c("haul_id", "MTC_haul")]
SCS_haul_stats <- left_join(SCS_haul_stats, simplified_mtc, by = "haul_id")
haul_indicator_cor_matrix <- cor(SCS_haul_stats[, c("haul_avg_troph", "haul_total_cpua", "species_richness", "evenness", "haul_avg_weight", "body_weight_coeff", "MTC_haul")],
                                 use = "pairwise.complete.obs")
print(haul_indicator_cor_matrix)
plot.new()
corrplot(haul_indicator_cor_matrix, method = "color",
         type = "upper", # Only show the upper triangle
         addCoef.col = "black", # Add the coefficient values
         tl.col = "black", tl.srt = 45, # Customize text labels
         col = COL2('RdBu', 20)) # Use a diverging color palette (Red/Blue)


simplified_haul_indicator_cor_matrix <- cor(SCS_haul_stats[, c("haul_avg_troph", "haul_total_cpua", "species_richness", "evenness", "haul_avg_weight", "body_weight_coeff", "MTC_haul")],
                                            use = "pairwise.complete.obs")
print(simplified_haul_indicator_cor_matrix)

plot.new()
corrplot(simplified_haul_indicator_cor_matrix, method = "color",
         type = "upper", # Only show the upper triangle
         addCoef.col = "black", # Add the coefficient values
         tl.col = "black", tl.srt = 45, # Customize text labels
         col = COL2('RdBu', 20)) # Use a diverging color palette (Red/Blue)


SCS_haul_stats$month<- as.numeric(as.character(SCS_haul_stats$month))
complex_haul_indicator_cor_matrix <- cor(SCS_haul_stats[, c("species_richness", "evenness", "haul_total_cpua", "haul_avg_weight", "body_weight_coeff", "haul_avg_troph", "MTC_haul", "sbt", "sst", "latitude", "depth", "area_swept", "year", "month")],
                                         use = "pairwise.complete.obs")
print(complex_haul_indicator_cor_matrix)

plot.new()

# Convert corrplot to a ggplot
indicator_corr_plot <- ggcorrplot(
  complex_haul_indicator_cor_matrix,
  type = "upper",
  hc.order = TRUE,
  title = "Correlation between haul metrics", 
  lab = TRUE,  # Show correlation coefficients
  lab_size = 3,  # Adjust size of correlation labels
  colors = c("#2166AC", "white", "#B2182B"),  # Blue-White-Red palette
  tl.srt = 45,  # Rotate text labels
  show.legend = TRUE) + 
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
    panel.spacing.y = unit(1.5, "lines"),
    panel.spacing.x = unit(1.5, "lines")
  )
print(indicator_corr_plot)

'
dont use covariance because that is sensitive to unit differences and requires knowledge of standard deviation to compare
complex_haul_indicator_cov_matrix <- cov(SCS_haul_stats[, c("haul_avg_troph", "haul_total_cpua", "species_richness", "evenness", "haul_avg_weight", "body_weight_coeff", "MTC_haul", "year", "month", "area_swept", "depth", "sbt", "sst", "latitude", "longitude")],
                                         use = "pairwise.complete.obs")

print(complex_haul_indicator_cov_matrix)
'


########### Construct yearly SCS stats for Structural break analysis and Time-varying vector autoregression #########
MTC_yearly_data <- MTC[, c("year", "MTC_year")]
troph_yearly_data <- yearly_average_trophic_level[, c("year", "overall_yearly_avg_trophic_level")]
cpua_yearly_data <- yearly_cpua_summary [, c("year", "avg_cpua")]
bodywgt_yearly_data <- yearly_summary_haul_avg [, c("year", "avg_weight")]
size_spec_yearly_data <- yearly_summary_size_spectrum [, c("year", "mean")]
size_spec_yearly_data$avg_kappa <- size_spec_yearly_data$mean # fix naming scheme of years for ToE/ToD
size_spec_yearly_data <- size_spec_yearly_data[,-2]

beta_div_vector <- as.list(SCS_dissimilarity_by_NAFO$yearly_avg_pairwise$Mean) # take the yearly beta div pairwise dissimilarity data
beta_div_vector <- c(
  beta_div_vector[1:48],
  NA,                 # add a NA for 2018 year because not enough info      
  beta_div_vector[49:50] 
)
beta_div_yearly_data <- data.frame(year = 1970:2020, dissimilarity = unlist(beta_div_vector))

evenness_yearly_data <- SCS_haul_stats %>%
  group_by(year) %>%
  summarize(ave_evenness = mean(evenness, na.rm = TRUE), .groups = 'drop')
evenness_yearly_data <- evenness_yearly_data[1:51,]

species_richness_yearly_data <- SCS_haul_stats %>%
  group_by(year) %>%
  summarize(ave_species_richness = mean(species_richness, na.rm = TRUE), .groups = 'drop')
species_richness_yearly_data <- species_richness_yearly_data[1:51,]

'
genus_richness_yearly_data <- SCS_haul_stats %>%
  group_by(year) %>%
  summarize(ave_genus_richness = mean(genus_richness, na.rm = TRUE), .groups = drop)
genus_richness_yearly_data <- genus_richness_yearly_data[1:51,]

family_richness_yearly_data <- SCS_haul_stats %>%
  group_by(year) %>%
  summarize(ave_family_richness = mean(family_richness, na.rm = TRUE), .groups = drop)
family_richness_yearly_data <- family_richness_yearly_data[1:51,]
'

CV_biomass_yearly_data <- yearly_biomass_CV

R4_ToE_count <- R4_summary %>%
  # 1. Group the data by the year the species was first observed
  group_by(first_year) %>%
  # 2. Count the number of unique accepted names in each group
  summarise(
    R4_appearances = n_distinct(accepted_name)
  ) %>%
  # 3. Ungroup the data (good practice after summarise)
  ungroup()
R4_ToE_count$year <- R4_ToE_count$first_year # fix naming scheme of years for ToE/ToD
R4_ToE_count <- R4_ToE_count[,-1]

R4_ToD_count <- R4_summary %>%
  # 1. Group the data by the year the species was first observed
  group_by(last_year) %>%
  # 2. Count the number of unique accepted names in each group
  summarise(
    R4_disappearances = n_distinct(accepted_name)
  ) %>%
  ungroup()
R4_ToD_count$year <- R4_ToD_count$last_year # fix naming scheme of years for ToE/ToD
R4_ToD_count <- R4_ToD_count[,-1]


R3_ToE_count <- R3_summary %>%
  # 1. Group the data by the year the species was first observed
  group_by(first_year) %>%
  # 2. Count the number of unique accepted names in each group
  summarise(
    R3_appearances = n_distinct(accepted_name)
  ) %>%
  # 3. Ungroup the data (good practice after summarise)
  ungroup()
R3_ToE_count$year <- R3_ToE_count$first_year # fix naming scheme of years for ToE/ToD
R3_ToE_count <- R3_ToE_count[,-1]

R3_ToD_count <- R3_summary %>%
  # 1. Group the data by the year the species was first observed
  group_by(last_year) %>%
  # 2. Count the number of unique accepted names in each group
  summarise(
    R3_disappearances = n_distinct(accepted_name)
  ) %>%
  ungroup()
R3_ToD_count$year <- R3_ToD_count$last_year # fix naming scheme of years for ToE/ToD
R3_ToD_count <- R3_ToD_count[,-1]

comp_sync_yearly_data <- data.frame(year = (1972:2018), trophic_synchrony = trophic_comp_sync_indices$TrophicSynchrony, rarity_synchrony = rare_sync_indices$SynchIndex, species_synchrony = spec_sync_indices$SynchIndex, trophic_compensation = trophic_comp_sync_indices$TrophicCompensation, rarity_compensation = rare_comp_indices$CompensationIndex, species_compensation = spec_comp_indices$CompensationIndex)

########### Envr yearly trends ##############

SCS_temp_volatility <- generate_temperature_volatility(scotian_shelf_summer_core)
SCS_temp_volatility_refined <- SCS_temp_volatility %>%
  select(c(year, CV_SBT, CV_SST))

envr_yearly_data <- SCS_haul_stats %>%
  group_by(year) %>%
  summarize(SBT = mean(sbt, na.rm = TRUE), SST = mean(sst, na.rm = TRUE), Depth = mean(depth, na.rm = TRUE), .groups = 'drop')

envr_yearly_data <- envr_yearly_data[1:51,]

envr_yearly_data <- left_join(envr_yearly_data, SCS_temp_volatility_refined, by = 'year')

mean_SST_yearly <- envr_yearly_data %>%
  select(c(year, SST))

mean_SBT_yearly <- envr_yearly_data %>%
  select(c(year, SBT))

CV_SST_yearly <- envr_yearly_data %>%
  select(c(year, CV_SST))

CV_SBT_yearly <- envr_yearly_data %>%
  select(c(year, CV_SBT))

evenness_yearly_data <- evenness_yearly_data[1:51,]

########## add NAFO landings ##############

Annual_catch_statistics <- read_csv("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/FISHGLOB Data Project/NAFO landings data/Annual catch statistics.csv")
names(Annual_catch_statistics)[ncol(Annual_catch_statistics)] <- "metric_tonnage"
names(Annual_catch_statistics)[1] <- "year"

total_NAFO_landings <- Annual_catch_statistics %>%
  group_by(year) %>%
  summarize(total_catch_tonnage = sum(metric_tonnage))

# install.packages("remotes")
#remotes::install_github("MarEcosystemApproaches/marea")
library(marea)
# search available Marea datasets
library(knitr)
kable(marea_metadata())

data("eco_indicators")

fishing_pressure_ea <- eco_indicators[eco_indicators@data$region %in% c("4VN","4X", "4W", "4VS")]

# Combine data
fishing_pressure_df <- data.frame(
  year = fishing_pressure_ea@data$year,
  nafo_div = fishing_pressure_ea@data$region,
  FP_all = fishing_pressure_ea@data$FP_ALL_value,
  FP_groundfish = fishing_pressure_ea@data$FP_GROUNDFISH_value
  )

annual_fishing_pressure_df <- fishing_pressure_df %>%
  group_by(year) %>%
  summarise(
    total_FP_all = mean(FP_all, na.rm = TRUE), 
    total_FP_groundfish = mean(FP_groundfish, na.rm = TRUE)
  ) 
#plot(annual_fishing_pressure_df$year, annual_fishing_pressure_df$total_FP_all, type = "l", col = "blue", xlab = "Year", ylab = "Fishing Pressure (All)", main = "Annual Fishing Pressure (All)")
#plot(annual_fishing_pressure_df$year, annual_fishing_pressure_df$total_FP_groundfish, type = "l", col = "blue", xlab = "Year", ylab = "Fishing Pressure (All)", main = "Annual Fishing Pressure (All)")
# note fishing pressure here calculated as biomass / landings according to Bundy, A., Gomez, C., & Cook, A. M. (n.d.). Guidance framework for the selection and evaluation of ecological indicators.
yearly_groundfish_fishing_pressure <- annual_fishing_pressure_df %>%
  filter(year >= 1970 & year <= 2020) %>%
  select(c(year, total_FP_groundfish))

################################################

data_list <- list(
  MTC_yearly_data,
  troph_yearly_data,
  cpua_yearly_data,
  bodywgt_yearly_data,
  size_spec_yearly_data,
  CV_biomass_yearly_data, 
  R4_ToE_count, 
  R4_ToD_count, 
  beta_div_yearly_data,
  evenness_yearly_data, 
  species_richness_yearly_data, 
  comp_sync_yearly_data, 
  total_NAFO_landings,
  yearly_groundfish_fishing_pressure,
  mean_SST_yearly, 
  mean_SBT_yearly, 
  CV_SST_yearly, 
  CV_SBT_yearly
)

# combine all yearly data to prepare for structural break analysis
structural_break_test <- data_list %>%
  reduce(left_join, by = "year")

structural_break_test[8:12][is.na(structural_break_test[8:12])] <- 0 # fill in ToE/ToD count with 0s when none occur

################# Print colnames of structural_break_test for variable selection #######################
print(colnames(structural_break_test))

######### Test for single change in relationship between two indicators ############
structural_break_df <- structural_break_test %>%
  mutate(time_index = 1:n()) # Add a simple sequential index

# Test for a break in the relationship: MTC_year = c0 + c1 * avg_cpua + error
bp_test <- Fstats(avg_weight ~ avg_cpua, data = structural_break_df)
plot(bp_test)

# Find the single most significant break point
bp_model <- breakpoints(bp_test) # Bai-Perron approach
summary(bp_model)

# The year of the break point:
year_of_break <- structural_break_df$year[bp_model$breakpoints]
cat("Most likely break year:", year_of_break, "\n")


############# top break points between all indicator pairs ################

df <- structural_break_test

# 1. Identify all indicator columns (excluding 'year')
indicators <- names(df)[names(df) != "year"]

# 2. Create all ordered pairs (Y ~ X)
all_pairs <- expand.grid(Dependent_Y = indicators, Independent_X = indicators, 
                         stringsAsFactors = FALSE) %>%
  filter(Dependent_Y != Independent_X)

results <- list()
for (i in 1:nrow(all_pairs)) {
  Y <- all_pairs$Dependent_Y[i]
  X <- all_pairs$Independent_X[i]
  
  formula_str <- paste(Y, "~", X)
  formula_obj <- as.formula(formula_str)
  
  # 1. Calculate F-statistics and find the single optimal break
  fstats_model <- Fstats(formula_obj, data = df)
  # Ensure the data has a 'year' column for the break point, as used later
  bp_model <- breakpoints(fstats_model, h = 0.10, break.max = 1) 
  
  break_index <- bp_model$breakpoints
  
  if (!is.na(break_index)) {
    # 2. Run the full linear model
    lm_model <- lm(formula_obj, data = df)
    
    # 3. Run a Chow test using the lm_model and the detected break index
    test_result <- sctest(lm_model, type = "Chow", point = break_index)
    
    results[[i]] <- data.frame(
      Dependent_Y = Y,
      Independent_X = X,
      Break_Year = df$year[break_index],
      F_Statistic = as.numeric(test_result$statistic),
      P_Value = as.numeric(test_result$p.value),
      Is_Significant_at_05 = as.numeric(test_result$p.value) < 0.05
    )
  } else {
    # ... (Handling for no break remains the same)
    results[[i]] <- data.frame(
      Dependent_Y = Y,
      Independent_X = X,
      Break_Year = NA,
      F_Statistic = NA,
      P_Value = NA,
      Is_Significant_at_05 = FALSE
    )
  }
}

# Combine all results into a single data frame
results_df <- do.call(rbind, results)

# ----------------------------------------------------------------------
## Print All Significant Relationships
# ----------------------------------------------------------------------

# Filter for all relationships that are significant at the 0.05 level and have a break
all_significant_breaks <- results_df %>%
  filter(!is.na(P_Value) & P_Value < 0.05) %>%
  arrange(P_Value) # Sort by P-Value for clear presentation

options(max.print = 10000)
print(all_significant_breaks)

# This histogram still uses all significant breaks
struc_break_hist <- ggplot(all_significant_breaks, aes(x = Break_Year)) +
  geom_histogram(binwidth = 1, fill = "lightgray", color = "black") +
  geom_vline(aes(xintercept = median(Break_Year)), linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(limits = c(1970, 2020)) +
  labs(
    title = "Distribution of significant structural breaks",
    subtitle = paste("Red line = Median break year"),
    x = "Year",
    y = "Frequency (number of indicator pairs)"
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
    panel.spacing.y = unit(1.5, "lines"),
    panel.spacing.x = unit(1.5, "lines")
  )
print(struc_break_hist)


# ----------------------------------------------------------------------
## Plot Segmented Regressions for the 30 Most Significant
# ----------------------------------------------------------------------

top_30_breaks <- all_significant_breaks %>% head(30)

# Create an empty list to store the plots
plot_list <- list()

for (j in 1:nrow(top_30_breaks)) {
  pair <- top_30_breaks[j, ]
  Y_sig <- pair$Dependent_Y
  X_sig <- pair$Independent_X
  BP_year <- pair$Break_Year
  
  # Prepare data for plotting the segmented model
  df_plot <- df %>%
    mutate(segment = factor(ifelse(year <= BP_year, "Before Break", "After Break")))
  
  # Create the plot
  p <- ggplot(df_plot, aes(x = .data[[X_sig]], y = .data[[Y_sig]])) +
    geom_point(aes(color = segment)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = segment), linewidth = 1.2) +
    labs(
      title = paste0(j, ". ", Y_sig, " vs ", X_sig),
      subtitle = paste0("Structural Break: ", BP_year, " (P-Value: ", format.pval(pair$P_Value, digits = 3), ")"),
      x = X_sig,
      y = Y_sig,
      color = "Time Period"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8, face = "bold"),
      legend.position = if(j == 1) "bottom" else "none",
      legend.text = element_text(size = 12)
    )
  
  plot_list[[j]] <- p
}

# Combine all plots
combined_plots <- wrap_plots(plot_list, nrow = 5, ncol = 6, guides = "collect") +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = 'Segmented Regression Models for Top 30 Significant Relationships',
    subtitle = 'Chow Test structural breaks (P < 0.05)',
    caption = 'Source: Your Data',
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold")
    )
  )

print(combined_plots)

# print examples of most significant structural breaks
most_significant_structural_break <- plot_list[[1]] + plot_list[[2]] + plot_annotation(
  title = 'Example Segmented Regression Models for Significant Relationships',
  subtitle = 'Chow Test structural breaks (P < 0.05)',
  theme = theme(
    plot.title = element_text(size = 32, face = "bold"),
    plot.subtitle = element_text(size = 28),
    legend.position = "bottom",
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 28, face = "bold")
  ))

print(most_significant_structural_break)

# 1. Ensure individual plots have their specific axis titles
# (Assuming plot_list[[1]] is SBT and plot_list[[15]] is SST for example)
struc_break_1 <- plot_list[[1]] + labs(y = "Yearly average haul evenness", x = "Yearly average community total density (kg/km^2)", title = "Structural break between community evenness and total density") + theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )
struc_break_2 <- plot_list[[2]] + labs(y = "Yearly Average Haul Evenness", x = "Annual MTC (°C)", title = "Structural break between community evenness and annual mean temperature of the catch") + theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )
'
# 2. Combine using patchwork
most_significant_structural_break <- p1 + p2 + 
  plot_layout(guides = "collect") + # This pulls both legends into one
  plot_annotation(
    title = Example Segmented Regression Models for Significant Relationships
    subtitle = Chow Test structural breaks (P < 0.05)
    theme_bw() + 
      theme(
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        # 1. Make the axis Titles bigger 
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        # 2. Make the axis Numbers (labels) bigger
        axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        # 3. Legend
        legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 16),
        legend.position = "bottom"
      )
  )

print(most_significant_structural_break)
'


# 1. Update the X-axis labels to include the proper superscript squared
x_label_density <- expression(paste("Yearly Average Community Total Density (kg ", km^-2, ")"))

# 2. Apply the theme and labels to the structural break plots INDIVIDUALLY
# This ensures consistency and avoids the patchwork & error
p_break1 <- struc_break_1 + labs(title = "Structural break point analysis: evenness & total density") +
  labs(x = x_label_density) + 
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )
print(p_break1)
p_break2 <- struc_break_2  + labs(title = "Structural break point analysis: evenness & MTC") +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )
print(p_break2)
# 3. Assemble the 4-panel plot
# Layout: A & B on top row, C & D on bottom row
combined_dynamic_relationships_plot <- (indicator_corr_plot + struc_break_hist) / (p_break1 + p_break2) +
  plot_annotation(tag_levels = 'A') 

print(combined_dynamic_relationships_plot)

####### structural break table ##########
colnames(all_significant_breaks)[colnames(all_significant_breaks) == "Dependent_Y"] <- "Dependent variable 1"
colnames(all_significant_breaks)[colnames(all_significant_breaks) == "Independent_X"] <- "Independent variable 2"
colnames(all_significant_breaks)[colnames(all_significant_breaks) == "Break_Year"] <- "Break year"
colnames(all_significant_breaks)[colnames(all_significant_breaks) == "F_Statistic"] <- "F statistic"
colnames(all_significant_breaks)[colnames(all_significant_breaks) == "P_Value"] <- "p value"
colnames(all_significant_breaks)[colnames(all_significant_breaks) == "Is_Significant_at_05"] <- "Significant if alpha = 0.05"

significant_breaks_table <- all_significant_breaks
print(significant_breaks_table)

##########################################
##### Explanatory variable table #########
##########################################

# make faceted figure of SST, SBT, landings

total_NAFO_landings

yearly_explanatory_variables <- left_join(total_NAFO_landings, envr_yearly_data, by = "year")
yearly_explanatory_variables <- left_join(yearly_explanatory_variables, yearly_groundfish_fishing_pressure, by = "year") 
# make into faceted plot of scatterplots


# 1. Create a mapping for your variable names to nice titles
var_names <- c(
  "total_catch_tonnage" = "Total Annual Landings (metric tonnes)",
  "SBT"    = "Sea Bottom Temperature (°C)",
  "CV_SBT" = "Coefficient of Variation in Sea Bottom Temp",
  "total_FP_groundfish" = "Groundfish Fishing Pressure (Biomass / Landings)"
)

# 2. Reshape the data
plot_data <- yearly_explanatory_variables %>%
  select(c(year, SBT, CV_SBT, total_FP_groundfish, total_catch_tonnage)) %>%
  pivot_longer(
    cols = c(SBT, CV_SBT, total_FP_groundfish, total_catch_tonnage), 
    names_to = "variable", 
    values_to = "value"
  ) %>%
  # THIS IS THE KEY STEP:
  # Define the order of the subplots by setting factor levels
  mutate(variable = factor(variable, levels = c(
    "SBT", 
    "CV_SBT", 
    "total_FP_groundfish", 
    "total_catch_tonnage"
  )))

# 3. Plot with individual facet labels
combined_explanatory_plot <- ggplot(plot_data, aes(x = year, y = value)) +
  geom_line(linetype = "dashed", color = "black", linewidth = 0.3, alpha = 0.6) +
  geom_point(color = "black", size = 1.5) +
  
  # Trend line and Confidence Interval
  geom_smooth(
    method = "loess", 
    se = TRUE, 
    color = "red", 
    fill = "grey", 
    alpha = 0.4, 
    linewidth = 1
  ) +
  
  # Facet using the labeller to change individual titles
  facet_wrap(
    ~variable, 
    scales = "free_y", 
    labeller = as_labeller(var_names)
  ) +
  
  # Styling
  labs(
    x = "Year",
    y = NULL # Remove the general Y title since facets have specific headers
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger 
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
    # --- THIS PART CONTROLS THE FACET TITLES ---
    strip.text = element_text(size = 18, face = "bold"), 
    # Optional: remove the grey box behind the titles for a cleaner look
    strip.background = element_blank()
    # -------------------------------------------
  )

##################################################
######## faceted plots for chapter 1 #############
##################################################

print(combined_biomass_plot)
print(combined_comm_struc_plot)
print(combined_biodiversity_plot)
print(combined_rarity_plot)
print(combined_ToE_ToD_plot)
print(combined_trophic_plot)
print(comp_sync_index_plot)
print(combined_dynamic_relationships_plot)
print(combined_explanatory_plot)

####################################
###### trophic comp/sync table ####
####################################

compared_clades <- c("Species", "Rarity groups", "Trophic groups")
# 1. Combine and format Compensation to 3 sig figs
compensation <- signif(
  c(SCS_interspecific_compensation_total$compensation_index, 
    SCS_rarity_compensation_index$compensation_index, 
    SCS_trophic_compensation_index$compensation_index), 
  digits = 3
)

# 2. Combine and format Synchrony to 3 sig figs
synchrony <- signif(
  c(SCS_interspecific_synchrony_total, 
    SCS_rarity_synchrony_index, 
    SCS_trophic_synchrony), 
  digits = 3
)
comp_sync_table <- data.frame(cbind(compared_clades, compensation, synchrony))

colnames(comp_sync_table)[colnames(comp_sync_table) == "compared_clades"] <- "Compared clades"
colnames(comp_sync_table)[colnames(comp_sync_table) == "compensation"] <- "Compensation index"
colnames(comp_sync_table)[colnames(comp_sync_table) == "synchrony"] <- "Synchrony index"

'
comp_sync_table <- comp_sync_table %>% 
  meps_theme()
'
print(comp_sync_table)


##############################################
######### Commmunity vignettes ###############
##############################################
# start with full community, find highest average yearly average cpua

species_total_avg_CPUA <- yearly_average_community %>%
  # 1. Pivot from wide to long (gather species columns)
  pivot_longer(
    cols = -Time, 
    names_to = "Species", 
    values_to = "CPUA"
  ) %>%
  # 2. Group by Species to calculate the overall average
  group_by(Species) %>%
  summarize(Average_CPUA = mean(CPUA, na.rm = TRUE)) %>%
  ungroup()

top_5_overall <- species_total_avg_CPUA %>%
  slice_max(Average_CPUA, n = 5)
print(top_5_overall)

# now do trophic groups
trophic_total_average <- trophic_source_table %>% 
  select(Species, trophic_group)
yearly_totals_to_match <- species_total_avg_CPUA %>%
  select(Species, Average_CPUA)
trophic_total_average <- left_join(trophic_total_average, yearly_totals_to_match, by = "Species")

# now species with trophic information have total average yearly average CPUA, find most spatially dominant members of each trophic group
top_5_per_trophic_group <- trophic_total_average %>%
  group_by(trophic_group) %>%
  slice_max(order_by = Average_CPUA, n = 5, with_ties = FALSE) %>%
  ungroup()
print(top_5_per_trophic_group)

#now do rarity groups
SCS_haul_prevalence$accepted_name <- as.character(SCS_haul_prevalence$accepted_name)
rarity_total_average <- SCS_haul_prevalence %>% 
  select(accepted_name, rarity_group)
names(rarity_total_average)[1] <- "Species"
rarity_total_average<-left_join(rarity_total_average, yearly_totals_to_match, by = "Species")

# now species with rarity information have total average yearly average CPUA, find most spatially dominant members of each trophic group
top_5_per_rarity_group <- rarity_total_average %>%
  group_by(rarity_group) %>%
  slice_max(order_by = Average_CPUA, n = 5, with_ties = FALSE) %>%
  ungroup()
print(top_5_per_rarity_group)

print(top_5_overall)
print(top_5_per_rarity_group)
print(top_5_per_trophic_group)

names(top_5_per_rarity_group)[2] <- "Group"
names(top_5_per_trophic_group)[2] <- "Group"
top_5_overall$Group <- "Overall"

names(top_5_per_rarity_group)[3] <- "Total Average CPUA"
names(top_5_per_trophic_group)[3] <- "Total Average CPUA"
top_5_overall$Group <- "Overall"
names(top_5_overall)[2] <- "Total Average CPUA"
top_5_overall <- top_5_overall[, c("Species", "Group", "Total Average CPUA")]

group_vignettes <- rbind(top_5_overall, top_5_per_rarity_group, top_5_per_trophic_group)
print(group_vignettes, n = 40)


# Create the gt table
species_table <- group_vignettes %>%
  gt() %>%
  # Basic formatting: making numbers look clean
  fmt_number(columns = `Total Average CPUA`, decimals = 1) %>%
  # Apply conditional coloring to the 'Group' column
  tab_style(
    style = cell_fill(color = "lightgreen"),
    locations = cells_body(rows = Group %in% c("Overall", "TL4", "R3"))
  ) %>%
  tab_style(
    style = cell_fill(color = "tomato"), # Red
    locations = cells_body(rows = Group %in% c("TL3", "R1"))
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_body(rows = Group %in% c("TL3.5", "R2"))
  ) %>%
  tab_style(
    style = cell_fill(color = "plum"), # Light Purple
    locations = cells_body(rows = Group == "R4")
  ) %>%
  # Adding a header and labels
  tab_header(
    title = md("Dominant members of investigated communities"),
    subtitle = "Highest survey mean yearly average CPUA"
  ) %>%
  cols_label(
    Species = md("**Assigned name**"),
    Group = md("**Community**"),
    `Total Average CPUA` = md("**CPUA (kg/km^2)**")
  )

# Display the table
community_table_species <- unique(group_vignettes$Species)


###########################
### add common names ######
###########################
alpha_id_list <- list()
for (species in community_table_species){
  records <- scotian_shelf_summer_core %>%
    filter(accepted_name == species) %>%
    pull(unique(aphia_id))
  alpha_id_list <- c(alpha_id_list, records)
}

# Create a unique lookup table from your dataframe
ref_lookup <- scotian_shelf_summer_core %>%
  select(accepted_name, aphia_id, family) %>%
  distinct() # This removes the duplicate rows

# Convert to dataframe
species_df <- data.frame(accepted_name = community_table_species)

# Perform a left join to pull in the aphia_id
results <- species_df %>%
  left_join(ref_lookup, by = "accepted_name")

my_ids <- results$aphia_id
common_names_list <- wm_common_id_(id = my_ids)
english_names <- common_names_list %>%
  bind_rows(.id = "aphia_id") %>%           # Flatten list to a dataframe
  filter(language == "English") %>%        # Filter for English
  group_by(aphia_id) 

unique_english_names <- english_names %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup()

aphia_list <- my_ids

'

# STEP 1: Get Scientific Names from WoRMS 
# rfishbase functions are optimized for scientific names
species_records <- wm_record(id = aphia_list)
sci_names <- species_records$scientificname

# STEP 2: Query FishBase for the "official" common names
# We use common_names and filter for the "Primary" name in FishBase
fb_names <- common_names(sci_names) %>%
  filter(Language == "English") %>%
  group_by(Species) %>%
  # FishBase often has many common names; 
  # we prioritize PreferredName if that column exists, or just the first entry
  slice(1) %>% 
  ungroup() %>%
  select(ScientificName = Species, FishBase_Common_Name = ComName)

# Result
print(fb_names)
'
unique_english_names <- unique_english_names[,2:3]
unique_english_names <- unique_english_names %>%
  mutate(id = as.numeric(id))

community_table_data <- results %>%
  left_join(unique_english_names, by = c("aphia_id" = "id"))

# deal with 2 missing values
community_table_data <- community_table_data %>%
  mutate(vernacular = case_when(
    accepted_name == "Sebastes" ~ "Redfish",
    accepted_name == "Myzopsetta ferruginea" ~ "Yellowtail flounder",
    TRUE ~ vernacular # This keeps the original value for everything else
  ))

# add population trends over time

# 1. Reshape the data to 'long' format
long_data <- relative_body_condition_df %>%
  pivot_longer(
    cols = -Year,                 # Keep year as is, pivot everything else
    names_to = "accepted_name", 
    values_to = "weight"
  ) %>%
  # 2. Filter for species in your community list
  filter(accepted_name %in% community_table_species) %>%
  # Remove NAs so they don't interfere with the regression
  drop_na(weight)

# 3. Calculate the slopes
coef_change_df <- long_data %>%
  group_by(accepted_name) %>%
  # Run regression (Weight ~ Year) for each species
  reframe(tidy(lm(weight ~ Year, data = cur_data()))) %>%
  # 4. Filter for only the slope (the 'year' term) and clean up
  filter(term == "Year") %>%
  select(accepted_name, coef_change = estimate)

# View the result
print(community_table_data)

community_table_data <- left_join(community_table_data, coef_change_df, by = "accepted_name")


# now combine community table data with right groupings
group_vignettes <- group_vignettes[-(1:5), ]
community_table_for_printing <- left_join(group_vignettes, community_table_data, by = c("Species" = "accepted_name"))
colnames(community_table_for_printing)[colnames(community_table_for_printing) == "Species"] <- "Taxon"
community_table_for_printing <- community_table_for_printing[, c("Taxon", "vernacular", "family", "aphia_id", "Group", "Total Average CPUA", "coef_change")]
colnames(community_table_for_printing)[colnames(community_table_for_printing) == "vernacular"] <- "Vernacular"
colnames(community_table_for_printing)[colnames(community_table_for_printing) == "family"] <- "Family"
colnames(community_table_for_printing)[colnames(community_table_for_printing) == "coef_change"] <- "Wgt slope"
colnames(community_table_for_printing)[colnames(community_table_for_printing) == "aphia_id"] <- "aphia ID"
colnames(community_table_for_printing)[colnames(community_table_for_printing) == "Total Average CPUA"] <- "Avg CPUA"

community_table_for_printing$`Avg CPUA` <- signif(community_table_for_printing$`Avg CPUA`, digits = 4)
community_table_for_printing$`Wgt slope` <- signif(community_table_for_printing$`Wgt slope`, digits = 3)

######## add abundance trends ######
unique_names <- unique(community_table_for_printing$Taxon)

comm_table_abundances <- SCS_relative_abundance_trends %>%
  filter(accepted_name %in% unique_names) %>%
  select(accepted_name, total_slope)

community_table_for_printing <- left_join(community_table_for_printing, comm_table_abundances, by = c("Taxon" = "accepted_name"))
colnames(community_table_for_printing)[colnames(community_table_for_printing) == "total_slope"] <- "Num slope"

# format for MEPS

community_vignettes_table <- community_table_for_printing
print(community_vignettes_table)


'
evenness_plot <- ggplot(yearly_evenness_summary, aes(x = year, y = avg_evenness)) +
  # Shaded confidence interval
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "cornflowerblue", alpha = 0.2) +
  # Yearly average line and points
  geom_line(color = "black", linewidth = 1, linetype = "dashed") +
  geom_point(color = "black", size = 2) +
  # 5-year sliding average line
  geom_line(aes(y = rolling_avg), color = "red", linewidth = 1.2) +
  labs(
    title = "Yearly Average Evenness of Hauls Over Time",
    x = "Year",
    y = "Average Evenness per Haul"
  ) +
  theme_bw()


combined_biodiversity_plot <- (barplot_richness + species_richness_plot) / (evenness_plot + beta_div_plot) +
  plot_annotation(tag_levels = A) +  # this ensures all panels have borders
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white")
  )

print(combined_biodiversity_plot)
'

#####################################################
######## how to calculate COVERAGE ##################
#####################################################

# Fishglob totals
nhauls_fishglob <- length(unique(FishGlob_public_clean$haul_id))
ntaxa_fishglob <- length(unique(FishGlob_public_clean$accepted_name))

# SCS totals
nhauls_SCS <- length(unique(scotian_shelf_year_round$haul_id))
ntaxa_SCS <- length(unique(scotian_shelf_year_round$accepted_name))

# SCS summer
nhauls_SCS_summer <- length(unique(scotian_shelf_summer$haul_id))
ntaxa_SCS_summer <- length(unique(scotian_shelf_summer$accepted_name))

# SCS Summer core
nhauls_SCS_summer_core <- length(unique(scotian_shelf_summer_core$haul_id))
ntaxa_SCS_summer_core <- length(unique(scotian_shelf_summer_core$accepted_name))

# SCS species rank
nhauls_SCS_summer_core_by_species <- length(unique(scotian_shelf_summer_core_by_species$haul_id))
ntaxa_SCS_summer_core_by_species <- length(unique(scotian_shelf_summer_core_by_species$accepted_name))

# SCS trophic
nhauls_SCS_summer_core_by_species_with_troph <- length(unique(scotian_with_troph$haul_id))
ntaxa_SCS_summer_core_by_species_with_troph <- length(unique(scotian_with_troph$accepted_name))
'
# SCS MTC1
nhauls_SCS_MTC1 <- length(unique(SCS_temp_hauls$haul_id))
ntaxa_SCS_MTC1 <- length(unique(SCS_temp_hauls$accepted_name))

# FISHGLOB MTC1
nhauls_Fish_MTC1 <- length(unique(Fish_temp_hauls$haul_id))
ntaxa_Fish_MTC1 <- length(unique(Fish_temp_hauls$accepted_name))

# FISHGLOB MTC2
nhauls_Fish_MTC2 <- length(unique(FishGlob_temp2_hauls$haul_id))
ntaxa_Fish_MTC2 <- length(unique(FishGlob_temp2_hauls$accepted_name))
'


# SCS MTC3
nhauls_SCS_MTC3 <- length(unique(SCS_temp3_hauls$haul_id))
ntaxa_SCS_MTC3 <- length(unique(SCS_temp3_hauls$accepted_name))


# construct MEPS format table for coverage
data_subsets <- c("FISHGLOB public dataset", "Scotian Shelf (SCS)", "SCS summer", "SCS summer (4Vn, 4Vs, 4W, 4X)*", "* with rarity information", "* by species", "* with thermal preference", "* with trophic information")
subset_number_of_hauls <- c(nhauls_fishglob, nhauls_SCS, nhauls_SCS_summer, nhauls_SCS_summer_core, nhauls_SCS_summer_core, nhauls_SCS_summer_core_by_species, nhauls_SCS_MTC3, nhauls_SCS_summer_core_by_species_with_troph)
subset_number_of_taxa <- c(ntaxa_fishglob, ntaxa_SCS, ntaxa_SCS_summer, ntaxa_SCS_summer_core, ntaxa_SCS_summer_core, ntaxa_SCS_summer_core_by_species, ntaxa_SCS_MTC3, ntaxa_SCS_summer_core_by_species_with_troph)
coverage_table <- data.frame(cbind(data_subsets, subset_number_of_hauls, subset_number_of_taxa))
colnames(coverage_table)[colnames(coverage_table) == "data_subsets"] <- "Data subset"
colnames(coverage_table)[colnames(coverage_table) == "subset_number_of_hauls"] <- "Number of hauls"
colnames(coverage_table)[colnames(coverage_table) == "subset_number_of_taxa"] <- "Number of taxa"

'coverage_table <- coverage_table %>% 
  meps_theme()'

print(coverage_table)
'
################# Check crossover of MTC1 and MTC2 on SCS ######################
MTC1_taxa <- unique(SCS_temp_hauls$accepted_name)
MTC2_taxa <- unique(SCS_temp2_hauls$accepted_name)
MTC_improvement <- setdiff(MTC2_taxa, MTC1_taxa) # finds taxa described in 2 not described in 1 so these would be added things, improvement
MTC_unexpected <- setdiff(MTC1_taxa, MTC2_taxa) # finds taxa described in 1 not described in 2 so these would be unexpected given expected improvement

################# Check crossover of MTC1 and MTC2 across FISHGLOB ######################
MTC1_taxa_Fish <- unique(Fish_temp_hauls$accepted_name)
MTC2_taxa_Fish <- unique(FishGlob_temp2_hauls$accepted_name)
MTC_improvement_Fish <- setdiff(MTC2_taxa_Fish, MTC1_taxa_Fish) # finds taxa described in 2 not described in 1 so these would be added things, improvement
MTC_unexpected_Fish <- setdiff(MTC1_taxa_Fish, MTC2_taxa_Fish) # finds taxa described in 1 not described in 2 so these would be unexpected given expected improvement

print(MTC_unexpected)
print(MTC_unexpected_Fish)

################# Check agreement of MTC1 and MTC2 across SCS ############################
MTC_shared_taxa <- intersect(MTC1_taxa, MTC2_taxa)

MTC1_shared_key <- simplified_thermal_dictionary %>%
  filter(accepted_name %in% MTC_shared_taxa) %>%
  distinct(accepted_name, estim_thermal_mean)

MTC2_shared_key <- SCS_temp2_hauls %>%
  filter(accepted_name %in% MTC_shared_taxa) %>%
  distinct(accepted_name, thetao_mean_mean)

MTC2_shared_key$estim_thermal_mean <- MTC2_shared_key$thetao_mean_mean

# Combine the two dataframes
# df1 and df2 should both have "accepted_name" and "estim_thermal_mean"
combined_MTC_strat_df <- inner_join(MTC1_shared_key, MTC2_shared_key, by = "accepted_name", suffix = c("_source1", "_source2"))

ggplot(combined_MTC_strat_df, aes(x = estim_thermal_mean_source1, y = estim_thermal_mean_source2)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + # 1:1 line
  labs(title = "Comparison of Thermal Mean Estimates",
       x = "Source 1 Temperature (°C)",
       y = "Source 2 Temperature (°C)") +
  theme_minimal()
'

############ make rank table for survey richness report ###############
rank_table <- as.data.frame(table(scotian_shelf_summer_core$rank))
colnames(rank_table) <- c("Rank", "Catch frequency")

rank_table <- rank_table %>%
  arrange(desc(`Catch frequency`)) %>%
  mutate(
    Proportion = `Catch frequency` / sum(`Catch frequency`),
    Percentage = Proportion * 100
  )

#######################################################################################
######### V3 Refined: Univariate Change & Gradient Forest Drivers #####################
#######################################################################################

# Load necessary libraries
library(bcp)
library(gradientForest)
library(tidyverse)
library(patchwork)
library(corrplot)
library(ggrepel)

# --- 1. Data Setup & Sanitization ---
univariate_indicator_df <- structural_break_test %>%
  select(-c(total_catch_tonnage, SST, SBT, CV_SST, CV_SBT, R4_disappearances, R4_appearances))

explanatory_variable_df <- structural_break_test %>%
  select(c(year, total_catch_tonnage, SST, SBT, CV_SST, CV_SBT))

years <- univariate_indicator_df$year
indicators <- univariate_indicator_df %>% select(-year) 

# --- 2. Refined BCP: Multiple Peak Detection ---
PROB_THRESHOLD <- 0.7 # Industry standard for "Strong" evidence
max_prob_years <- data.frame()

for(col_name in colnames(indicators)) {
  temp_df <- data.frame(year = years, val = indicators[[col_name]])
  clean_temp <- temp_df %>% filter(!is.na(val))
  
  if(nrow(clean_temp) < 6) next 
  
  bcp_fit <- bcp(clean_temp$val)
  prob_vector <- bcp_fit$posterior.prob
  
  # Identify all years where probability is a local maximum above threshold
  peak_indices <- which(
    prob_vector >= PROB_THRESHOLD & 
      prob_vector >= lag(prob_vector, default = 0) & 
      prob_vector >= lead(prob_vector, default = 0)
  )
  
  if(length(peak_indices) == 0) next
  
  for(peak_idx in peak_indices) {
    # Skip peaks at the very edges (need enough data to calculate means)
    if(peak_idx < 3 || peak_idx > (nrow(clean_temp) - 2)) next
    
    mean_before <- mean(clean_temp$val[1:peak_idx], na.rm = TRUE)
    mean_after  <- mean(clean_temp$val[(peak_idx+1):nrow(clean_temp)], na.rm = TRUE)
    
    percent_diff <- if(abs(mean_before) > 1e-5) {
      ((mean_after - mean_before) / mean_before) * 100
    } else { NA }
    
    max_prob_years <- rbind(max_prob_years, data.frame(
      Indicator = col_name, 
      Peak_Year = clean_temp$year[peak_idx], 
      Prob = prob_vector[peak_idx],
      Percent_Shift = percent_diff,
      Mean_Pre = mean_before,
      Mean_Post = mean_after
    ))
  }
}

# --- 2.5 Prepare Lagged Predictors & Sanitize for GF Updated: Include Year to test for temporal confounding ---
base_preds <- c("total_catch_tonnage", "SST", "SBT", "CV_SST", "CV_SBT")

# Add year into the explanatory pool
explanatory_with_time <- explanatory_variable_df 

explanatory_lagged <- explanatory_with_time %>%
  arrange(year) %>%
  mutate(across(all_of(base_preds), 
                list(lag_1 = ~lag(., 1), lag_2 = ~lag(., 2), lag_3 = ~lag(., 3)), 
                .names = "{.col}_lag_{.fn}"))

# Combine and force numeric to fix density error
combined_clean <- inner_join(univariate_indicator_df, explanatory_lagged, by = "year") %>%
  mutate(across(everything(), ~as.numeric(as.character(.)))) %>%
  drop_na()


# --- 3. Run Repeatable Gradient Forest ---

# Ensure year is in the predictor list
# We include year as a numeric column to compete with environmental variables
all_responses <- colnames(indicators)
valid_responses <- all_responses[sapply(combined_clean[all_responses], function(x) length(unique(x)) > 5)]
preds_with_year <- c("year", colnames(combined_clean)[!(colnames(combined_clean) %in% c("year", all_responses))])

set.seed(123)
PROB_CUTOFF <- 0.5 
IMP_CUTOFF <- 0.001 
gf <- gradientForest(data = as.data.frame(combined_clean), 
                     predictor.vars = preds_with_year, 
                     response.vars = valid_responses, 
                     ntree = 2000, 
                     transform = NULL, 
                     compact = TRUE)

# --- 4. Threshold & Direction Extraction ---
gf_summary <- data.frame()
actual_indicators <- names(gf$result)

for(ind in actual_indicators) {
  imp <- importance(gf, type = "Weighted", response = ind)
  if (length(imp) == 0 || max(imp) <= 0) next
  
  best_driver <- names(which.max(imp))
  all_cum_imp <- cumimp(gf, predictor = best_driver, type = "Species")
  dat <- all_cum_imp[[ind]]
  
  if (!is.null(dat)) {
    peak_threshold <- dat$x[which.max(diff(dat$y))]
    correlation <- cor(combined_clean[[ind]], combined_clean[[best_driver]], use="complete.obs")
    direction <- ifelse(correlation > 0, "Increasing", "Decreasing")
    
    gf_summary <- rbind(gf_summary, data.frame(
      Indicator = ind,
      Top_Driver = best_driver,
      Importance = as.numeric(max(imp)),
      Threshold = peak_threshold,
      Response_Direction = direction
    ))
  }
}

# --- 5. Final Classification & Threshold Mapping ---
find_threshold_year <- function(driver_name, threshold_val, data) {
  matching_years <- data %>%
    filter(!!sym(driver_name) <= threshold_val) %>% 
    pull(year)
  if(length(matching_years) > 0) return(min(matching_years)) else return(NA)
}

classification_df <- inner_join(max_prob_years, gf_summary, by = "Indicator") %>%
  mutate(
    Lag = ifelse(grepl("_lag_", Top_Driver), gsub(".*_lag_", "", Top_Driver), "0"),
    Driver_Base = gsub("_lag_.*", "", Top_Driver)
  ) %>%
  rowwise() %>%
  mutate(
    Threshold_Year = find_threshold_year(Top_Driver, Threshold, combined_clean),
    Change_Type = case_when(
      (Prob > PROB_CUTOFF & Importance > IMP_CUTOFF) ~ "Phase Shifted",
      (Prob <= PROB_CUTOFF & Importance > IMP_CUTOFF) ~ "Gradual (Driver-led)",
      (Prob > PROB_CUTOFF & Importance <= IMP_CUTOFF) ~ "Abrupt (Unexplained)",
      TRUE ~ "Gradual/Stochastic"
    )
  ) %>%
  ungroup()

# --- 6. Tables & Visualizations ---

# Table 1: BCP Results
Bayesian_change_point_table <- max_prob_years %>%
  rename("BCP year" = Peak_Year, `Probability of shift` = Prob, 
         `Mean before shift` = Mean_Pre, `Mean after shift` = Mean_Post, 
         `Percent shift (%)` = Percent_Shift) %>%
  arrange(`Indicator`, `BCP year`, desc(`Probability of shift`)) %>%  
  mutate(across(where(is.numeric), ~round(., 3))) 

# Table 2: Final Classification
classification_of_indicators_table <- classification_df %>%
  select(Indicator, Threshold, Threshold_Year, 
         Change_Type, Driver_Base, Lag, Response_Direction, Importance) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  arrange(desc(`Threshold`))%>%
  meps_theme()

# Plot 1: Correlation Matrix
# Identify all variables for the matrix ---
# We want: Year + Base Predictors + Indicators that made it into the GF model
base_preds_names <- c("year", "total_catch_tonnage", "SST", "SBT", "CV_SST", "CV_SBT")
actual_indicators <- names(gf$result) # Only those with enough importance/variance

# Combine them into one selection
corr_selection <- cbind(explanatory_with_time, indicators)

# Calculate Correlation Matrix
indicator_corr <- cor(corr_selection, use = "complete.obs")

# --- 2. Plotting with ggcorrplot ---
yearly_indicator_corr_plot <- ggcorrplot(
  indicator_corr,
  hc.order = TRUE,         # Groups variables with similar temporal patterns
  type = "upper",
  title = "Correlation between yearly shelf metrics", 
  lab = TRUE,  
  lab_size = 2.5,          # Slightly smaller to accommodate more variables
  colors = c("#2166AC", "white", "#B2182B"), 
  tl.srt = 45,  
  show.legend = TRUE) + 
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title.x = element_blank(), # Standard for corrplots
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 14)
  )

print(yearly_indicator_corr_plot)

###########################################################################
######## V3 corr plot edit ############
###########################################################################

# =================================================================
# Fisheries Sensitivity Analysis: SBT and total_FP_groundfish
# Framework: Generalized Additive Models (GAM) with Tensor Products
# =================================================================

# 1. Load Required Libraries
if (!require("mgcv")) install.packages("mgcv")
if (!require("tidyverse")) install.packages("tidyverse")

library(mgcv)      # For non-linear modeling
library(tidyverse) # For data manipulation and plotting

# 2. Define the Biological Indicators
bio_indicators <- c(
  "MTC_year", "overall_yearly_avg_trophic_level", "avg_cpua", 
  "avg_weight", "avg_kappa", "CV_biomass", "dissimilarity", 
  "ave_evenness", "ave_species_richness", "trophic_synchrony", 
  "rarity_synchrony", "species_synchrony", "trophic_compensation", 
  "rarity_compensation", "species_compensation"
)

# 3. Define the Model Fitting Function
# Uses te() to handle the high correlation between total_FP_groundfish and year
fit_fisheries_gam <- function(indicator, data) {
  # Note: Variable names must match your dataframe exactly (case-sensitive)
  form <- as.formula(paste(indicator, "~ s(SBT, k=5) + te(total_FP_groundfish, year, k=5)"))
  
  # We use REML for better smoothing parameter selection
  model <- gam(form, data = data, method = "REML", family = gaussian())
  return(model)
}

# 4. Define the Sensitivity Calculation Function
calculate_sensitivity <- function(model, data) {
  # Define the reference point (median values across the 51-year series)
  # na.rm = TRUE is used to handle any missing years in your dataset
  baseline_data <- data.frame(
    SBT = median(data$SBT, na.rm = TRUE),
    total_FP_groundfish = median(data$total_FP_groundfish, na.rm = TRUE),
    year = median(data$year, na.rm = TRUE)
  )
  
  base_pred <- predict(model, newdata = baseline_data)
  
  # Calculate ranges for 10% perturbation
  sbt_range <- diff(range(data$SBT, na.rm = TRUE))
  fp_range  <- diff(range(data$total_FP_groundfish, na.rm = TRUE))
  
  # Sensitivity to SBT: Predict change based on a 10% range increase
  sbt_perturb <- baseline_data %>% mutate(SBT = SBT + (0.1 * sbt_range))
  sens_sbt <- predict(model, newdata = sbt_perturb) - base_pred
  
  # Sensitivity to Fishing: Predict change based on a 10% range increase
  fp_perturb <- baseline_data %>% mutate(total_FP_groundfish = total_FP_groundfish + (0.1 * fp_range))
  sens_fp <- predict(model, newdata = fp_perturb) - base_pred
  
  return(data.frame(
    SBT_Sens = as.numeric(sens_sbt), 
    FP_Sens = as.numeric(sens_fp)
  ))
}

# 5. Execute the Loop Across All Indicators
results_list <- list()

# Ensure your dataframe is named yearly_change
yearly_change <- corr_selection 
# Standardize all biological indicators (Z-score: mean=0, sd=1)
# This allows you to compare sensitivity across different units
yearly_change_logged <- yearly_change %>%
  mutate(across(all_of(bio_indicators), ~ log(. + 1))) # Log-transform to handle skewness, add 1 to avoid log(0)
yearly_change_std <- yearly_change_logged %>%
  mutate(across(all_of(bio_indicators), ~ as.numeric(scale(.))))

# Now run the rest of the script using 'yearly_change_std'
for (ind in bio_indicators) {
  message(paste("Analyzing:", ind))
  
  # Wrap in try() to catch errors (e.g., if an indicator has zero variance)
  fit <- try(fit_fisheries_gam(ind, yearly_change), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    sens <- calculate_sensitivity(fit, yearly_change)
    sens$indicator <- ind
    results_list[[ind]] <- sens
  } else {
    message(paste("Warning: Model for", ind, "failed to converge."))
  }
}

# 6. Final Data Formatting and Visualization
sensitivity_summary <- bind_rows(results_list)

# Reshape for ggplot
plot_data <- sensitivity_summary %>%
  pivot_longer(cols = c(SBT_Sens, FP_Sens), names_to = "Driver", values_to = "Impact")

ggplot(plot_data, aes(x = indicator, y = Impact, fill = Driver)) +
  geom_col(position = "dodge") +
  coord_flip() + # Better for long indicator names
  theme_minimal() +
  scale_fill_manual(
    values = c("SBT_Sens" = "firebrick3", "FP_Sens" = "dodgerblue4"),
    labels = c("Fishing Pressure", "Sea Bottom Temp")
  ) +
  labs(
    title = "Indicator Sensitivity Analysis",
    subtitle = "Predicted response to a 10% increase in environmental/fishing drivers",
    x = "Biological Indicator",
    y = "Magnitude of Predicted Change"
  )


#########################################################################
# find maximal correlation for each indicator with other bio indicators #
#########################################################################
bio_corr_df <- corr_selection %>%
  select(-c(year, total_catch_tonnage, total_FP_groundfish, SST, SBT, CV_SST, CV_SBT)) # Keep only the biological indicators

# Calculate Correlation Matrix
bio_indicator_corr <- cor(bio_corr_df, use = "complete.obs")

# 1. Create a copy so we don't mess up the original matrix
mat <- as.matrix(bio_indicator_corr)

# 2. Set diagonal to 0 to ignore self-correlation
diag(mat) <- 0

# 3. Find the index of the max absolute value for each row
max_idx <- apply(abs(mat), 1, which.max)

# 4. Extract the actual values (with signs) using matrix indexing
max_values <- mat[cbind(1:nrow(mat), max_idx)]

# 5. Build the final dataframe
maximal_correlation_table <- data.frame(
  indicator = rownames(mat),
  most_correlated_with = colnames(mat)[max_idx],
  correlation = max_values
)
#########################################################################################
# USING HAUL VARS find maximal correlation for each indicator with other bio indicators #
#########################################################################################
haul_bio_corr_df <- SCS_haul_stats %>%
  select(c("species_richness", "evenness", "haul_total_cpua", "haul_avg_weight", "body_weight_coeff", "haul_avg_troph", "MTC_haul")) # Keep only the biological indicators

# Calculate Correlation Matrix
haul_bio_indicator_corr <- cor(haul_bio_corr_df, use = "pairwise.complete.obs")

# 1. Create a copy so we don't mess up the original matrix
mat <- as.matrix(haul_bio_indicator_corr)

# 2. Set diagonal to 0 to ignore self-correlation
diag(mat) <- 0

# 3. Find the index of the max absolute value for each row
max_idx <- apply(abs(mat), 1, which.max)

# 4. Extract the actual values (with signs) using matrix indexing
max_values <- mat[cbind(1:nrow(mat), max_idx)]

# 5. Build the final dataframe
haul_maximal_correlation_table <- data.frame(
  indicator = rownames(mat),
  most_correlated_with = colnames(mat)[max_idx],
  correlation = max_values
)

#########################################################################
# find maximal correlation for each indicator with predictors ###########
#########################################################################
bio_corr_df <- corr_selection

# Calculate Correlation Matrix
predictor_indicator_corr <- cor(bio_corr_df, use = "complete.obs")

# 1. Create a copy so we don't mess up the original matrix
mat <- as.matrix(predictor_indicator_corr)

# Define predictors
predictors <- c("SBT", "total_FP_groundfish", "CV_SBT", "CV_SST", "SST" )

# Subset to predictor columns
sub_mat <- mat[, predictors]

removes <- c(predictors, "year", "total_catch_tonnage")
# Remove predictors from rows
sub_mat <- sub_mat[!rownames(sub_mat) %in% removes, ]

# Max correlation per row
max_idx <- apply(abs(sub_mat), 1, which.max)
max_values <- sub_mat[cbind(1:nrow(sub_mat), max_idx)]

# Final table
maximal_predictor_correlation_table <- data.frame(
  indicator = rownames(sub_mat),
  most_correlated_predictor = colnames(sub_mat)[max_idx],
  correlation = max_values
)
#####################################################################################
# USING HAULS find maximal correlation for each indicator with predictors ###########
#####################################################################################
haul_predictor_corr_df <- SCS_haul_stats %>%
  select(c("species_richness", "evenness", "haul_total_cpua", "haul_avg_weight", "body_weight_coeff", "haul_avg_troph", "MTC_haul", "sst", "sbt", "latitude", "longitude", "depth", "area_swept")) # Keep only the biological indicators

# Calculate Correlation Matrix
haul_predictor_corr <- cor(haul_predictor_corr_df, use = "complete.obs")

# 1. Create a copy so we don't mess up the original matrix
mat <- haul_predictor_corr

# Define predictors
predictors <- c("sst", "sbt", "latitude", "longitude", "depth", "area_swept")

# Subset to predictor columns
sub_mat <- mat[, predictors]

removes <- c(predictors)
# Remove predictors from rows
sub_mat <- sub_mat[!rownames(sub_mat) %in% removes, ]

# Max correlation per row
max_idx <- apply(abs(sub_mat), 1, which.max)
max_values <- sub_mat[cbind(1:nrow(sub_mat), max_idx)]

# Final table
maximal_haul_predictor_correlation_table <- data.frame(
  indicator = rownames(sub_mat),
  most_correlated_predictor = colnames(sub_mat)[max_idx],
  correlation = max_values
)


######################################################################
# Plot 2: BCP Histogram
BCP_histogram <- ggplot(max_prob_years, aes(x = Peak_Year)) +
  geom_bar(fill = "lightgray", color = "black", alpha = 0.8) +
  geom_vline(aes(xintercept = median(Peak_Year)), linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(limits = c(1970, 2020)) +
  labs(title = "Distribution of significant Bayesian change points",
       subtitle = paste("Red line = Median shift year"),
       x = "Year", y = "Frequency (number of indicators)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )



# Plot 3: Community Sensitivity
sensitivity_data <- classification_df %>%
  filter(Driver_Base == "total_catch_tonnage")

community_sensitivity_plot <- ggplot(sensitivity_data, aes(x = Threshold, y = reorder(Indicator, Threshold))) +
  geom_segment(aes(x = min(Threshold)*0.8, xmax = max(Threshold)*1.1, y = Indicator, yend = Indicator), 
               color = "grey90", size = 0.5) +
  geom_point(aes(color = Change_Type, size = Importance)) +
  geom_text_repel(aes(label = paste0(round(Threshold, 0), " t")), size = 3.5) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Community Sensitivity to Fishing Pressure", x = "Total Catch Tonnage Threshold", y = "") +
  theme_bw()

# --- 7. Community Threshold Shape & Trend Validation (ORDERED) ---

# 1. Identify the top 8 indicators based on Importance
top_indicators <- gf_summary %>%
  arrange(desc(Importance)) %>%
  slice(1:8) %>%
  pull(Indicator)

# 2. Pre-calculate thresholds to allow sorting
threshold_info <- data.frame()

for(ind in top_indicators) {
  best_d <- gf_summary %>% filter(Indicator == ind) %>% pull(Top_Driver)
  all_cum_imp <- cumimp(gf, predictor = best_d, type = "Species")
  
  if (ind %in% names(all_cum_imp)) {
    raw_dat <- all_cum_imp[[ind]]
    # Calculate the threshold (steepest part of the curve)
    rates_of_change <- diff(raw_dat$y) / diff(raw_dat$x)
    peak_val <- raw_dat$x[which.max(rates_of_change)]
    
    threshold_info <- rbind(threshold_info, data.frame(
      Indicator = ind,
      Threshold = peak_val,
      Driver = best_d
    ))
  }
}

# 3. Sort indicators by descending Threshold value
ordered_indicators <- threshold_info %>%
  arrange(desc(Threshold))

plot_list <- list()

# 4. Loop through the ORDERED list
for(i in 1:nrow(ordered_indicators)) {
  ind <- ordered_indicators$Indicator[i]
  best_d <- ordered_indicators$Driver[i]
  peak_val <- ordered_indicators$Threshold[i]
  
  all_cum_imp <- cumimp(gf, predictor = best_d, type = "Species")
  raw_dat <- all_cum_imp[[ind]]
  plot_df <- data.frame(x = raw_dat$x, y = raw_dat$y)
  
  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_line(linewidth = 1.2, color = "steelblue") +
    geom_vline(xintercept = peak_val, linetype = "dashed", color = "red") +
    annotate("text", x = peak_val, y = max(plot_df$y)*0.1, 
             label = paste("Threshold:", round(peak_val, 1)), 
             angle = 90, vjust = -0.5, size = 3, color = "red") +
    labs(title = paste(ind),
         subtitle = paste("Top Driver:", best_d),
         x = "Driver Value", 
         y = "Cum. Importance") +
    theme_bw() +
    theme(plot.title = element_text(size = 11, face = "bold"),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8))
  
  plot_list[[ind]] <- p
}

# 5. Use patchwork to wrap the ordered list into a grid
if(length(plot_list) > 0) {
  # Increased ncol to 3 or 4 might be better for 8 plots, keeping 2 as per your code
  threshold_shapes_plot <- wrap_plots(plot_list, ncol = 2) + 
    plot_annotation(title = "Biological Sensitivity Curves (Ordered by Descending Threshold)",
                    subtitle = "Higher threshold values (e.g., higher Catch or SST) appear first")
  print(threshold_shapes_plot)
}

# --- 8. Landing vs. Year Importance Comparison ---
# This visually demonstrates your finding that Landings > Year
imp_all <- importance(gf)
imp_df <- data.frame(variable = names(imp_all), importance = as.numeric(imp_all))

comp_plot <- ggplot(imp_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity", aes(fill = ifelse(variable == "year", "Time", 
                                                ifelse(grepl("catch", variable), "Fishery", "Environment")))) +
  coord_flip() +
  scale_fill_manual(values = c("Fishery" = "firebrick", "Time" = "black", "Environment" = "grey70")) +
  labs(title = "Overall Driver Importance across Community",
       subtitle = "Landings (Red) underperforming Year (Black) suggesting spurious relationship",
       x = "Predictor Variable", y = "Weighted Mean R-squared", fill = "Category") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )

# Print outputs

print(threshold_shapes_plot)
print(comp_plot)

print(Bayesian_change_point_table)
print(BCP_histogram)

print(classification_of_indicators_table)

#print(community_sensitivity_plot)

print(yearly_indicator_corr_plot)

#######################################
###### BCP Visualization for Fig 10 ###
#######################################

# 1. Select the Indicator (e.g., 'avg_kappa' from your image)
target_ind <- colnames(indicators)[4]

# 1. Select the Indicator and the "Primary" shift to highlight in P1

all_shifts <- max_prob_years %>% filter(Indicator == target_ind)

# Define the primary shift (e.g., the one with the highest probability)
primary_shift <- all_shifts %>% filter(Prob == max(Prob)) %>% slice(1)

# 2. Prepare the data
plot_data <- data.frame(
  year = years, 
  val = indicators[[target_ind]]
) %>% filter(!is.na(val))

# Re-run BCP to get the continuous probability vector for P2
bcp_fit_single <- bcp(plot_data$val)
plot_data$prob <- bcp_fit_single$posterior.prob

# --- 3. Upper Panel (P1): Focus on ONLY the Primary Shift ---
BCP_primary <- ggplot(plot_data, aes(x = year, y = val)) +
  geom_line(color = "grey80", alpha = 0.8) +
  geom_point(size = 1.5) +
  # Single vertical line for the primary shift
  geom_vline(xintercept = primary_shift$Peak_Year, 
             color = "darkred", linetype = "solid", size = 1) +
  # Mean segments for ONLY this primary shift
  geom_segment(aes(x = min(year), xend = primary_shift$Peak_Year, 
                   y = primary_shift$Mean_Pre, yend = primary_shift$Mean_Pre), 
               color = "blue", linetype = "dashed", size = 1.2) +
  geom_segment(aes(x = primary_shift$Peak_Year + 1, xend = max(year), 
                   y = primary_shift$Mean_Post, yend = primary_shift$Mean_Post), 
               color = "red", linetype = "dashed", size = 1.2) +
  annotate("text", x = primary_shift$Peak_Year, y = max(plot_data$val, na.rm=T), 
           label = paste("Primary shift:", primary_shift$Peak_Year), 
           angle = 90, vjust = 1.5, hjust = 1.5, color = "darkred", size = 3.5) +
  labs(title = "Bayesian change point analysis: average body weight",
       subtitle = paste0("Bayesian change point: ", primary_shift$Peak_Year, 
                         " (Posterior Probability = ", round(primary_shift$Prob, 3), ")"),
       y = "Haul average body weight (kg)", x = "Year") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )

# --- 4. Lower Panel (P2): Keep ALL Detected Shifts ---
BCP_prob <- ggplot(plot_data, aes(x = year, y = prob)) +
  geom_area(fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue", size = 0.8) +
  # Reference threshold
  geom_hline(yintercept = 0.7, linetype = "dotted", color = "red", size = 1) +
  # Mark ALL peaks identified in the multi-peak loop
  geom_point(data = all_shifts, aes(x = Peak_Year, y = Prob), 
             color = "red", size = 2) +
  geom_text(data = all_shifts, aes(x = Peak_Year, y = Prob, label = Peak_Year), 
            vjust = -1, size = 3) +
  labs(title = "Posterior probability of shift: average body weight", y = "Posterior Probability", x = "Year") +
  ylim(0, 1.1) + # Leave room for labels
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )

########## Add univariate stuff to dynamic relationships plot
combined_corr_plot<- indicator_corr_plot + yearly_indicator_corr_plot + plot_annotation(tag_levels = 'A')
print(combined_corr_plot)

combined_dynamic_relationships_plot <- (BCP_primary + BCP_histogram) /  (p_break1 + struc_break_hist) + 
  plot_annotation(tag_levels = 'A')

print(combined_dynamic_relationships_plot)


##################################################
############ TABLES for chapter 1 ################
##################################################
print(rank_table)

##################################################
############ TABLES formatted for MEPS ###########
##################################################
# New MEPS Table Function
new_meps_theme <- function(data, title = NULL, subtitle = NULL) {
  data %>%
    gt() %>%
    # Typography (MEPS accepts Helvetica/Arial or Times New Roman)
    tab_options(
      table.font.names = "Times New Roman", 
      table.font.size = px(14),
      column_labels.font.weight = "bold",
      table.align = "center",
      table.border.top.color = "black",
      table.border.top.width = px(2),
      table.border.bottom.color = "black",
      table.border.bottom.width = px(2),
      column_labels.border.top.color = "black",
      column_labels.border.top.width = px(2),
      column_labels.border.bottom.color = "black",
      column_labels.border.bottom.width = px(1),
      table_body.border.bottom.color = "black",
      table_body.border.bottom.width = px(2),
      heading.border.bottom.style = "none"
    ) %>%
    # Alignment & Spacing
    cols_align(align = "center", columns = everything()) %>%
    cols_align(align = "left", columns = 1) %>% 
    tab_style(
      style = cell_text(indent = pct(2)),
      locations = cells_body(columns = 1)
    )
}

# Load data
T1_quantities <- read_csv("C:/Users/besch/OneDrive - Dalhousie University/Masters Thesis - Worm Lab/Lit Review/Manuscript elements/Edit 2 for committee/Figures save/Tables/T1 quantities.csv")

# Generate table
indicator_list_table <- T1_quantities %>%
  new_meps_theme()
# --- MEPS COMPLIANT EXPORT STEP ---
# Save directly as an editable Word document
gtsave(
  data = indicator_list_table,
  filename = "C:/Users/besch/Desktop/FishGlob project/MEPS formatted tables/Table_1.docx"
)

T2_comp_sync_table <- comp_sync_table %>%
  new_meps_theme()
gtsave(data = T2_comp_sync_table,
  filename = "C:/Users/besch/Desktop/FishGlob project/MEPS formatted tables/Table_2.docx"
)

TS2_coverage_table <- coverage_table %>% 
  new_meps_theme()
gtsave(data = TS2_coverage_table,
       filename = "C:/Users/besch/Desktop/FishGlob project/MEPS formatted tables/TS2.docx")

TS3_comm_vignettes <- community_vignettes_table %>%
  new_meps_theme()
gtsave(data = TS3_comm_vignettes,
       filename = "C:/Users/besch/Desktop/FishGlob project/MEPS formatted tables/TS3.docx")

TS4_cluster_table <- cluster_community_table %>%
  new_meps_theme()
gtsave(data = TS4_cluster_table,
       filename = "C:/Users/besch/Desktop/FishGlob project/MEPS formatted tables/TS4.docx")

TS5_BCP <- Bayesian_change_point_table %>%
  new_meps_theme()
gtsave(data = TS5_BCP,
       filename = "C:/Users/besch/Desktop/FishGlob project/MEPS formatted tables/TS5.docx")

TS6_struc_breaks <- significant_breaks_table %>%
  new_meps_theme()
gtsave(data = TS6_struc_breaks,
       filename = "C:/Users/besch/Desktop/FishGlob project/MEPS formatted tables/TS6.docx")

#################################################################################
########## Check whats driving late stage body weight gains among T4 ############
#################################################################################
library(dplyr)
library(broom)
T4_SCS <- scotian_with_troph %>%
  filter(trophic_group == "TL4" & body_weight > 0 & body_weight < 1000) # filter to TL4 and remove outliers in body weight
# 1. Calculate the threshold: Global TL4 average body weight (post-2000)
global_tl4_mean <- T4_SCS %>%
  filter(year >= 2000) %>%
  summarise(mu = mean(body_weight, na.rm = TRUE)) %>%
  pull(mu)

# 2. Identify all positive contributors
all_positive_drivers <- T4_SCS %>%
  filter(year >= 2000) %>%
  group_by(accepted_name) %>%
  filter(n_distinct(year) > 3) %>%
  summarise(
    mean_bw = mean(body_weight, na.rm = TRUE),
    weight_slope = tidy(lm(body_weight ~ year))$estimate[2],
    biomass_slope = tidy(lm(wgt_cpua ~ year))$estimate[2],
    .groups = "drop"
  ) %>%
  # Filter: Must be increasing in weight OR (Larger than avg AND increasing in biomass)
  filter(
    weight_slope > 0 | (mean_bw > global_tl4_mean & biomass_slope > 0)
  ) %>%
  mutate(
    is_naturally_heavier = mean_bw > global_tl4_mean,
    primary_mechanism = case_when(
      weight_slope > 0 & biomass_slope > 0 ~ "Both Growth & Biomass",
      weight_slope > 0 ~ "Individual Weight Gain",
      TRUE ~ "Increasing Biomass of Large Individuals"
    )
  ) %>%
  arrange(desc(weight_slope))

# Export the full list

T4_weight_gain_drivers <- all_positive_drivers %>%
  rename(Species = accepted_name, 
         `Mean body weight (kg)` = mean_bw,
         `Wgt slope` = weight_slope,
         `Biomass slope` = biomass_slope,
         `Naturally heavier than avg` = is_naturally_heavier,
         `Recovery mechanism` = primary_mechanism
  )

T4_weight_gain_drivers <- T4_weight_gain_drivers %>%
  mutate(Vernacular = case_when(
    Species == "Tetronarce nobiliana"      ~ "Atlantic torpedo ray",
    Species == "Macrourus berglax"         ~ "roughhead grenadier",
    Species == "Hippoglossus hippoglossus" ~ "Atlantic halibut",
    Species == "Coryphaenoides rupestris"  ~ "roundnose grenadier",
    Species == "Cottunculus microps"       ~ "polar sculpin",
    Species == "Merluccius albidus"        ~ "offshore hake",
    Species == "Centroscymnus coelolepis"  ~ "Portuguese dogfish",
    Species == "Cyttopsis rosea"           ~ "rosy dory",
    Species == "Petromyzon marinus"        ~ "sea lamprey",
    Species == "Polymixia nobilis"         ~ "stout beardfish",
    Species == "Lophius americanus"        ~ "monkfish",
    Species == "Anoplogaster cornuta"      ~ "fangtooth",
    Species == "Merluccius bilinearis"     ~ "silver hake",
    Species == "Chauliodus sloani"         ~ "Sloanes viperfish",
    Species == "Paraliparis calidus"       ~ "lowfin snailfish",
    TRUE                                   ~ "unknown"
  ))
T4_weight_gain_drivers <- T4_weight_gain_drivers %>%
  select(c("Species", "Vernacular", `Naturally heavier than avg`, `Mean body weight (kg)`, `Recovery mechanism`, `Wgt slope`, `Biomass slope`))

# Apply to the dataframe
T4_weight_gain_drivers <- T4_weight_gain_drivers %>%
  mutate(across(where(is.numeric), ~ {
    # 1. Round to 4 significant figures
    s <- signif(., 4)
    # 2. Convert to character, stripping trailing zeros/dots
    char_s <- format(s, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
    # 3. Clean up any lingering decimal points (e.g., "10.")
    sub("\\.$", "", char_s)
  }))

T4_weight_gain_drivers <- T4_weight_gain_drivers %>%
  new_meps_theme()

gtsave(
  data = T4_weight_gain_drivers,
  filename = "C:/Users/besch/Desktop/FishGlob project/MEPS formatted tables/S7_T4_drivers.docx"
)

########################################################
######## faceted plots formatted for MEPS ##############
########################################################
theme_ir <- function() {
  theme_minimal(base_size = 9, base_family = "Helvetica") %+replace% 
    theme(
      # Axis titles and text set to the recommended 8-9 pt range
      axis.title      = element_text(size = 9, face = "plain"),
      axis.text       = element_text(size = 8),
      
      # Rotates the X-axis text 45 degrees and aligns it perfectly to the ticks
      axis.text.x     = element_text(hjust = 1, vjust = 1),
      
      # Subpanel tags (A, B, C...) should be clear
      plot.tag        = element_text(size = 10, face = "bold"), 
      
      # Legend text alignment
      legend.title    = element_text(size = 9),
      legend.text     = element_text(size = 8),
      
      # Clean borders typical for scientific journals
      axis.line       = element_line(color = "black", linewidth = 0.5),
      axis.ticks      = element_line(color = "black", linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Apply the updated theme to your plots
indicator_corr_plot <- indicator_corr_plot + 
  theme_ir() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

yearly_indicator_corr_plot <- yearly_indicator_corr_plot + 
  theme_ir() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Combine them horizontally with patchwork
combined_corr_plot <- indicator_corr_plot + yearly_indicator_corr_plot + plot_annotation(tag_levels = 'A')

# Print the final result
print(combined_corr_plot)

###########################################
######## More MEPS themed figures #########
###########################################
library(cowplot)

# --- 1. DEFINE COMPLIANT MEPS THEME ---
theme_meps <- function(base_size = 9, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace% {
    theme(
      text = element_text(family = base_family, size = base_size),
      axis.title = element_text(size = base_size, face = "plain"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      
      # MEPS style prefers distinct black axis boundary lines
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      
      # Strip out background grids for a clean presentation
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Plot margins (ensures nothing clips when packed together)
      plot.margin = margin(t = 5, r = 10, b = 5, l = 10, unit = "pt")
    )
  }
}

# --- 2. BASE PLOTTING FUNCTION FOR SUBPLOTS ---
# This ensures lines, points, and loess curves are identical across all panels
make_meps_sub_plot <- function(data, y_var, y_label) {
  ggplot(data, aes(x = year, y = .data[[y_var]])) +
    geom_line(linetype = "dashed", color = "black", linewidth = 0.3, alpha = 0.6) +
    geom_point(color = "black", size = 1.2) +
    geom_smooth(
      method = "loess", 
      se = TRUE, 
      color = "red", 
      fill = "grey", 
      alpha = 0.4, 
      linewidth = 0.8
    ) +
    labs(x = "Year", y = y_label) +
    theme_meps()
}

# --- 3. GENERATE THE 4 INDIVIDUAL PANELS ---
# Units are formatted into standard SI parentheses layout required by MEPS
p_sbt <- make_meps_sub_plot(
  yearly_explanatory_variables, 
  "SBT", 
  "Sea bottom temperature (°C)"
)

p_cv  <- make_meps_sub_plot(
  yearly_explanatory_variables, 
  "CV_SBT", 
  "CV of sea bottom temperature"
)

p_fp  <- make_meps_sub_plot(
  yearly_explanatory_variables, 
  "total_FP_groundfish", 
  "Groundfish fishing pressure (landings / biomass)"
)

p_cat <- make_meps_sub_plot(
  yearly_explanatory_variables, 
  "total_catch_tonnage", 
  "Total annual landings (metric tonnes)"
)

# --- 4. COMBINE AND AUTOMATICALLY TAG PANELS (A, B, C, D) ---
# Four panels are grouped into a full-width page block (169 mm width)
F3_combined_explanatory_plot <- plot_grid(
  p_sbt, p_cv, p_fp, p_cat,
  ncol = 2, 
  nrow = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 10,
  label_fontfamily = "Helvetica",
  label_fontface = "bold",
  align = "v" # Keeps axes vertically aligned nicely
)

# Preview the figure in your R session
print(F3_combined_explanatory_plot)

#####################################################
f4_top <- community_cpua_plot + 
  theme_meps() + 
  labs(x = NULL)

f4_bottom <- CV_community_cpua_plot + 
  theme_meps()

F4_combined_biomass_plot <- plot_grid(
  f4_top, f4_bottom,
  ncol = 1, 
  nrow = 2,
  labels = c("A", "B"),
  label_size = 10,
  label_fontfamily = "Helvetica",
  label_fontface = "bold",
  align = "v"
)

library(cowplot)
library(ggplot2)

# Set up your target directory pathway
export_dir <- "C:/Users/besch/Desktop/FishGlob project/MEPS formatted figures/"
if(!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)


# =========================================================================
# FIGURE 5: COMBINED COMMUNITY STRUCTURE PLOT (Custom Grid Spacing)
# =========================================================================

# 1. Format individual panels with MEPS theme and prune redundant top-row x-labels
p5_a <- community_body_weight_plot_haul_avg + theme_meps() + labs(x = NULL)
p5_b <- community_size_spectra_plot         + theme_meps() + labs(x = NULL)
p5_c <- proportional_change_plot            + theme_meps() + labs(x = NULL)
p5_d <- MTC_plot                            + theme_meps()
p5_e <- PCoA_plot                           + theme_meps()

# 2. Assemble Top Row (3 columns)
f5_top_row <- plot_grid(
  p5_a, p5_b, p5_c,
  ncol = 3,
  labels = c("A", "B", "C"),
  label_size = 10,
  label_fontfamily = "Helvetica",
  label_fontface = "bold"
)

# 3. Assemble Bottom Row (2 columns, matching outer horizontal boundaries)
f5_bottom_row <- plot_grid(
  p5_d, p5_e,
  ncol = 2,
  labels = c("D", "E"),
  label_size = 10,
  label_fontfamily = "Helvetica",
  label_fontface = "bold"
)

# 4. Combine rows vertically and align axes
F5_combined_comm_struc_plot <- plot_grid(
  f5_top_row,
  f5_bottom_row,
  ncol = 1,
  nrow = 2,
  rel_heights = c(1, 1), # Ensures rows take up equal vertical space
  align = "v"            # Aligns left and right margins perfectly
)

# 5
print(F5_combined_comm_struc_plot)


# =========================================================================
# FIGURE 6: COMBINED BIODIVERSITY PLOT (2x2 Grid)
# =========================================================================
F6_combined_biodiversity_plot <- plot_grid(
  barplot_richness     + theme_meps() + labs(x = NULL),
  species_richness_plot + theme_meps() + labs(x = NULL),
  evenness_plot        + theme_meps(),
  beta_div_plot        + theme_meps(),
  ncol = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 10,
  label_fontfamily = "Helvetica",
  label_fontface = "bold",
  align = "hv"
)

print(F6_combined_biodiversity_plot)

# =========================================================================
# FIGURE 7: COMBINED RARITY PLOT (Horizontal Split)
# =========================================================================
F7_combined_rarity_plot <- plot_grid(
  boxplot_bwgt_R + theme_meps(),
  boxplot_MTC_R  + theme_meps() + labs(y = "Estimated thermal mean (°C)"),
  ncol = 2,
  labels = c("A", "B"),
  label_size = 10,
  label_fontfamily = "Helvetica",
  label_fontface = "bold",
  align = "h"
)
print(F7_combined_rarity_plot)


# =========================================================================
# FIGURE 8: COMBINED TIME OF EMERGENCE / DISAPPEARANCE (Horizontal Split)
# =========================================================================
F8_combined_ToE_ToD_plot <- plot_grid(
  p3 + theme_meps(),
  p4 + theme_meps() + labs(y = NULL), # Strip out redundant left axis labels
  ncol = 2,
  labels = c("A", "B"),
  label_size = 10,
  label_fontfamily = "Helvetica",
  label_fontface = "bold",
  align = "h"
)
print(F8_combined_ToE_ToD_plot)


# =========================================================================
# FIGURE 10: COMP/SYNC INDEX PLOT (Vertical Split)
# =========================================================================
F10_comp_sync_index_plot <- plot_grid(
  p_comp + theme_meps() + labs(x = NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)),
  p_sync + theme_meps() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)),
  ncol = 1,
  nrow = 2,
  labels = c("A", "B"),
  label_size = 10,
  label_fontfamily = "Helvetica",
  label_fontface = "bold",
  align = "v"
)
print(F10_comp_sync_index_plot)


# =========================================================================
# FIGURE 12: COMBINED DYNAMIC RELATIONSHIPS PLOT (2x2 Grid)
# =========================================================================
F12_combined_dynamic_relationships_plot <- plot_grid(
  BCP_primary      + theme_meps(),
  BCP_histogram    + theme_meps() + labs(x = NULL),
  p_break1         + theme_meps(),
  struc_break_hist + theme_meps(),
  ncol = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 10,
  label_fontfamily = "Helvetica",
  label_fontface = "bold",
  align = "hv"
)

print(F12_combined_dynamic_relationships_plot)

# =========================================================================
# FIGURE 9: Trophic dynamics (1x2 Grid)
# =========================================================================
F9_trophic_plot <- plot_grid(
  T_level_hist + theme_meps(),
  yearly_trophic_level_plot + theme_meps(),
  ncol = 2,
  nrow = 1,
  labels = c("A", "B"),
  label_size = 10,
  label_fontfamily = "Helvetica",
  label_fontface = "bold",
  align = "v"
)
print(F9_trophic_plot)

#####################################################
##### spatial autocorrelation check V2 ##############
#####################################################

autocorr_data <- SCS_haul_stats %>%
  select(year, longitude, latitude, haul_avg_troph, species_richness, evenness, haul_total_cpua, haul_avg_weight, body_weight_coeff)
library(tidyverse)
library(spdep)
library(sf)

# 1. Convert dataframe to an sf object (assuming WGS84)
autocorr_sf <- st_as_sf(autocorr_data, coords = c("longitude", "latitude"), crs = 4326)

# 2. Define the variables of interest
bio_vars <- c("haul_avg_troph", "species_richness", "evenness", 
              "haul_total_cpua", "haul_avg_weight", "body_weight_coeff")

yearly_moran_results <- autocorr_data %>%
  group_by(year) %>%
  group_modify(~ {
    sub_sf <- st_as_sf(.x, coords = c("longitude", "latitude"), crs = 4326)
    
    map_df(bio_vars, function(v) {
      # 1. Subset data to remove NAs for this specific variable
      valid_data <- sub_sf %>% filter(!is.na(.data[[v]]))
      
      # 2. Check if we have enough points left to actually calculate neighbors
      # (Need at least k+1 points, where k=5)
      if (nrow(valid_data) <= 5) {
        return(tibble(variable = v, moran_i = NA, p_value = NA))
      }
      
      # 3. Re-calculate neighbors for the subset of valid points
      nb <- knn2nb(knearneigh(st_coordinates(valid_data), k = 5))
      listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
      
      # 4. Run the test
      test <- moran.test(valid_data[[v]], listw, 
                         alternative = "greater", 
                         zero.policy = TRUE)
      
      tibble(
        variable = v,
        moran_i = test$estimate[1],
        p_value = test$p.value
      )
    })
  })

print(yearly_moran_results)

# 1. Clean and prepare data for plotting
plot_data <- yearly_moran_results %>%
  filter(!is.na(moran_i)) %>%
  mutate(significant = ifelse(p_value < 0.05, "p < 0.05", "NS"))

# 2. Create the plot
ggplot(plot_data, aes(x = year, y = moran_i, group = variable)) +
  # Add a reference line at 0 (Random distribution)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  # Lines for trends
  geom_line(color = "steelblue", alpha = 0.5) +
  # Points colored by significance
  geom_point(aes(color = significant), size = 2.5) +
  # Separate panels for each variable
  facet_wrap(~variable, scales = "free_y") +
  # Formatting
  scale_color_manual(values = c("p < 0.05" = "#D55E00", "NS" = "#56B4E9")) +
  theme_minimal() +
  labs(
    title = "Yearly spatial autocorrelation trends via nearest neighbor approach",
    subtitle = "Positive values indicate clustering; Orange points represent statistical significance",
    x = "Year",
    y = "Moran's I Index",
    color = "Significance"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )

###########################
## now whole time series ##
###########################

library(sf)
library(spdep)
library(purrr)

# 1. Jitter coordinates slightly to handle identical haul locations
# st_jitter moves points by a tiny fraction of a degree (~1m)
autocorr_sf_jittered <- st_jitter(autocorr_sf, factor = 0.00001)
coords_all <- st_coordinates(autocorr_sf_jittered)

# 2. Global Analysis with NA handling
whole_series_results <- map_df(bio_vars, function(v) {
  
  # Remove NAs for this specific variable
  valid_sf <- autocorr_sf_jittered %>% filter(!is.na(.data[[v]]))
  valid_coords <- st_coordinates(valid_sf)
  
  # Re-calculate weights for the subset of valid points
  nb_all <- knn2nb(knearneigh(valid_coords, k = 5))
  listw_all <- nb2listw(nb_all, style = "W", zero.policy = TRUE)
  
  # Run the test
  test <- moran.test(valid_sf[[v]], listw_all, 
                     alternative = "greater", 
                     zero.policy = TRUE)
  
  tibble(
    variable = v,
    moran_i = test$estimate[1],
    p_value = test$p.value
  )
})

print(whole_series_results)

ggplot(whole_series_results, aes(x = reorder(variable, moran_i), y = moran_i)) +
  geom_col(aes(fill = p_value < 0.05)) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey70"), 
                    labels = c("TRUE" = "p < 0.05", "FALSE" = "NS"),
                    name = "Significance") +
  labs(
    title = "Whole series spatial autocorrelation via nearest neighbor approach",
    x = "Biological Variable",
    y = "Moran's I Index"
  ) + 
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # 1. Make the axis Titles bigger
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    # 2. Make the axis Numbers (labels) bigger
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    # 3. Legend
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16)
  )

#####################################################
##### Spatial Autocorrelation V3: Distance-Based ####
#####################################################
library(tidyverse)
library(spdep)
library(sf)


'
# Function to find the "Distance to Independence" (Correlogram)
find_autocorr_dist <- function(sf_data, var_name, max_dist = 100, step = 10) {
  # Extract coordinates once
  coords <- st_coordinates(sf_data)
  distances <- seq(step, max_dist, by = step)
  
  results <- map_df(distances, function(d) {
    # 1. Define neighbors within distance d (d in km if using longlat=TRUE)
    # Note: dnearneigh returns a  object
    nb <- dnearneigh(coords, 0, d, longlat = TRUE)
    
    # 2. Safety check: Are there ANY connections at this distance?
    # If a distance is too small, no hauls may touch each other.
    if(all(card(nb) == 0)) {
      return(tibble(dist_km = d, moran_i = NA, p_value = NA))
    }
    
    # 3. Create listw with zero.policy=TRUE to handle isolated points (islands)
    listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
    
    # 4. Calculate Morans I
# Using tryCatch because some distance bins might have insufficient variance
test_res <- tryCatch({
  test <- moran.test(sf_data[[var_name]], listw, 
                     alternative = "greater", 
                     zero.policy = TRUE, 
                     randomization = TRUE)
  tibble(dist_km = d, moran_i = test$estimate[1], p_value = test$p.value)
}, error = function(e) {
  tibble(dist_km = d, moran_i = NA, p_value = NA)
})

return(test_res)
})
return(results)
}

#------------------------------------------------------------------------------
# 1. Yearly Spatial Range Analysis
#------------------------------------------------------------------------------
yearly_dist_results <- autocorr_data %>%
  group_by(year) %>%
  group_modify(~ {
    sub_sf <- st_as_sf(.x, coords = c("longitude", "latitude"), crs = 4326)
    
    map_df(bio_vars, function(v) {
      valid_sf <- sub_sf %>% filter(!is.na(.data[[v]]))
      
      # Use a reasonable max_dist for yearly hauls (e.g., 200km)
      dist_stats <- find_autocorr_dist(valid_sf, v, max_dist = 200, step = 20)
      
      # Identify the maximum distance where spatial clustering is still significant
      max_sig_dist <- dist_stats %>% 
        filter(p_value < 0.05) %>% 
        summarise(limit = max(dist_km, na.rm = TRUE)) %>% 
        pull(limit)
      
      # Handle cases where no distance was significant
      spatial_range <- if(is.infinite(max_sig_dist) || is.na(max_sig_dist)) 0 else max_sig_dist
      
      tibble(variable = v, spatial_range_km = spatial_range)
    })
  })

#------------------------------------------------------------------------------
# 2. Global Correlogram Analysis
#------------------------------------------------------------------------------
# Warning: Running the whole series (thousands of points) with dnearneigh 
# at large distances can be memory-intensive. Using 10km steps up to 300km.
whole_series_dist <- map_df(bio_vars, function(v) {
  valid_sf <- autocorr_sf_jittered %>% filter(!is.na(.data[[v]]))
  
  cat("Processing Correlogram for:", v, "\n")
  
  correlogram <- find_autocorr_dist(valid_sf, v, max_dist = 200, step = 10)
  correlogram %>% mutate(variable = v)
})

#------------------------------------------------------------------------------
# 3. Visualization: The Correlogram
#------------------------------------------------------------------------------
ggplot(whole_series_dist %>% filter(!is.na(moran_i)), 
       aes(x = dist_km, y = moran_i, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(linewidth = 1) +
  # Use shapes to denote significance (easier to see than alpha)
  geom_point(aes(shape = p_value < 0.05), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), 
                     name = "Significant (p < 0.05)") +
  theme_bw() +
  facet_wrap(~variable, scales = "free_y") +
  labs(title = "Spatial Correlogram: Patch Size of Ecological Indicators",
       subtitle = "Distance where Morans I crosses zero indicates the spatial limit of autocorrelation",
       x = Distance (km), 
       y = Morans I Index) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))
'


# 1. New Robust Global Function
find_autocorr_dist_robust <- function(sf_data, var_name, max_dist = 300, step = 10) {
  
  # Subsample if the dataset is huge (e.g., > 2000 points) to prevent memory crashes
  if(nrow(sf_data) > 2000) {
    sf_data <- sf_data[sample(1:nrow(sf_data), 2000), ]
  }
  
  coords <- st_coordinates(sf_data)
  distances <- seq(step, max_dist, by = step)
  
  results <- map_df(distances, function(d) {
    nb <- dnearneigh(coords, 0, d, longlat = TRUE)
    
    # Check for empty neighbor sets
    if(all(card(nb) == 0)) return(NULL)
    
    listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
    
    # Run test with numeric enforcement
    test_res <- tryCatch({
      test <- moran.test(as.numeric(sf_data[[var_name]]), listw, 
                         alternative = "greater", zero.policy = TRUE)
      
      tibble(dist_km = as.numeric(d), 
             moran_i = as.numeric(test$estimate[1]), 
             p_value = as.numeric(test$p.value))
    }, error = function(e) {
      NULL # Return NULL so map_df just ignores failed bins
    })
    
    return(test_res)
  })
  return(results)
}

# 2. Re-run Global Analysis
whole_series_dist <- map_df(bio_vars, function(v) {
  # Drop NAs and check for variance
  valid_sf <- autocorr_sf_jittered %>% 
    filter(!is.na(.data[[v]]))
  
  if(var(valid_sf[[v]], na.rm=TRUE) == 0) return(NULL)
  
  cat("Processing Correlogram for:", v, "...\n")
  find_autocorr_dist_robust(valid_sf, v) %>% mutate(variable = v)
})

# 3. Visualization (Now protected against empty data)
if (nrow(whole_series_dist) > 0) {
  ggplot(whole_series_dist, aes(x = dist_km, y = moran_i, color = variable)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(aes(shape = p_value < 0.05), size = 3) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), name = "Significant") +
    facet_wrap(~variable, scales = "free_y") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14),
      # 1. Make the axis Titles bigger
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold"),
      # 2. Make the axis Numbers (labels) bigger
      axis.text.y = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 14, color = "black"),
      # 3. Legend
      legend.position = "right",
      legend.text = element_text(size = 14),
      legend.title = element_text(face = "bold", size = 16)
    ) + 
    labs(title = "Spatial correlograms of haul level biological indicators", x = "Patch size - Distance (km)", y = "Moran's I")
} else {
  print("Error: No significant autocorrelation found or data points too sparse.")
}


#############################################################################
#### V2 Facet plot of biological haul variables ########
##############################################################################


# Required libraries
library(tidyverse)
library(viridis)
library(patchwork)
library(sf)

# 1. Define variables and Threshold
vars_to_map <- c("haul_avg_troph", "species_richness", "evenness", 
                 "haul_total_cpua", "haul_avg_weight", "body_weight_coeff")
sd_threshold <- 3  # Remove data points beyond 3 Standard Deviations

# --- STEP A: HISTOGRAMS OF RAW DATA ---
hist_list <- map(vars_to_map, function(var_name) {
  ggplot(SCS_haul_stats, aes(x = .data[[var_name]])) +
    geom_histogram(fill = "steelblue", color = "white", bins = 30) +
    labs(title = str_replace_all(var_name, "_", " ") %>% str_to_title(),
         x = NULL, y = "Count") +
    theme_minimal()
})

print(wrap_plots(hist_list, ncol = 3) + 
        plot_annotation(title = "Raw Biological Variables: Distribution Before Cleaning"))


# --- STEP B: SYSTEMATIC OUTLIER REMOVAL ---
# We create a 'cleaned' version of the dataframe
SCS_cleaned <- SCS_haul_stats

for(var in vars_to_map) {
  # Calculate mean and SD, ignoring NAs
  m <- mean(SCS_cleaned[[var]], na.rm = TRUE)
  s <- sd(SCS_cleaned[[var]], na.rm = TRUE)
  
  # Identify the bounds
  upper <- m + (sd_threshold * s)
  lower <- m - (sd_threshold * s)
  
  # Filter the data
  SCS_cleaned <- SCS_cleaned %>%
    filter(.data[[var]] >= lower & .data[[var]] <= upper | is.na(.data[[var]]))
  
  message(paste("Variable:", var, "| Bounds:", round(lower, 2), "to", round(upper, 2)))
}


# --- STEP C: SMOOTH MAPS (CLEANED DATA) ---
generate_smooth_map <- function(var_name, df) {
  
  data_subset <- df %>% filter(!is.na(.data[[var_name]]))
  
  ggplot() +
    # Layer 1: Landmass (Ensure na_map is pre-loaded)
    geom_polygon(data = na_map, aes(x = long, y = lat, group = group), 
                 fill = "grey85", color = "grey70") +
    
    # Layer 2: Diffused Biological Surface
    stat_summary_2d(data = data_subset, 
                    aes(x = longitude, y = latitude, z = .data[[var_name]]),
                    fun = mean, bins = 500) + 
    scale_fill_viridis(option = "magma", name = "Mean",
                       guide = guide_colorbar(barheight = unit(3, "cm"))) +
    
    # Layer 3: NAFO Divisions (Ensure nafo_sf is pre-loaded)
    geom_sf(data = nafo_sf, fill = NA, color = "black", 
            linetype = "dashed", linewidth = 0.2) +
    
    coord_sf(xlim = c(-69.0, -56.0), ylim = c(41.0, 48.0), expand = FALSE) +
    labs(title = str_replace_all(var_name, "_", " ") %>% str_to_title()) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title = element_blank(),
      panel.background = element_rect(fill = "aliceblue"),
      panel.grid = element_blank()
    )
}

# Generate and print final map
plot_list_smooth <- map(vars_to_map, ~generate_smooth_map(.x, SCS_cleaned))

final_smooth_plot <- wrap_plots(plot_list_smooth, ncol = 2) + 
  plot_annotation(
    title = "Spatially diffused haul level biological indicators (outliers removed)",
    subtitle = paste("Outliers defined as >", sd_threshold, "SD from the mean"),
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

print(final_smooth_plot)

#############################################################################
######### Kriging - spatial coverage of haul averages on SCS ################
#############################################################################

test_year_for_kriging <- SCS_haul_stats %>%
  filter(year == 1985) 

# Example: your dataframe should be named df and have columns:
# longitude, latitude, avg_troph, species_richness, evenness, total_biomass_cpua

# Convert to sf and then sp object
df_sf <- st_as_sf(test_year_for_kriging, coords = c("longitude", "latitude"), crs = 4326)
df_sp <- as(df_sf, "Spatial")

# Project to UTM zone 20N (suitable for Nova Scotia)
df_sp_utm <- spTransform(df_sp, CRS("+proj=utm +zone=20 +datum=WGS84 +units=km"))

# Create a regular grid over the region for prediction
bbox <- bbox(df_sp_utm)
x_range <- seq(bbox[1,1], bbox[1,2], by = 5) # 5 km grid spacing
y_range <- seq(bbox[2,1], bbox[2,2], by = 5)
grid <- expand.grid(x = x_range, y = y_range)
coordinates(grid) <- ~x+y
proj4string(grid) <- CRS(proj4string(df_sp_utm))

# Automated variogram model selection function
auto_fit_variogram <- function(vgm_emp, models = c("Sph", "Exp", "Gau", "Lin")) {
  best_model <- NULL
  best_sse <- Inf
  best_fit <- NULL
  for (m in models) {
    fit <- try(fit.variogram(vgm_emp, model = vgm(m)), silent = TRUE)
    if (!inherits(fit, "try-error") && !any(is.na(fit))) {
      # Check for negative or zero range/psill
      if (any(fit$range <= 0) || any(fit$psill < 0)) next
      fitted_vals <- variogramLine(fit, dist_vector = vgm_emp$dist)$gamma
      sse <- sum((vgm_emp$gamma - fitted_vals)^2, na.rm = TRUE)
      if (sse < best_sse) {
        best_sse <- sse
        best_model <- m
        best_fit <- fit
      }
    }
  }
  if (is.null(best_fit)) {
    warning("No suitable variogram model found, using arithmetic mean.")
    return(list(fit = NULL, model = "arithmetic mean"))
  }
  return(list(fit = best_fit, model = best_model))
}

# Function to compute spatially weighted mean and record model
spatially_weighted_mean_and_model <- function(varname, sp_points, grid) {
  vgm_model <- variogram(as.formula(paste(varname, "~ 1")), data = sp_points)
  fit_result <- auto_fit_variogram(vgm_model)
  fit_vgm <- fit_result$fit
  model_used <- fit_result$model
  if (is.null(fit_vgm)) {
    mean_val <- mean(sp_points@data[[varname]], na.rm = TRUE)
  } else {
    kriged <- krige(as.formula(paste(varname, "~ 1")), sp_points, grid, model = fit_vgm)
    mean_val <- mean(kriged$var1.pred, na.rm = TRUE)
  }
  return(list(mean = mean_val, model = model_used))
}

# Variables to process
vars <- c("haul_avg_troph", "haul_avg_weight", "species_richness", "genus_richness", "evenness", "haul_total_cpua", "body_weight_coeff", "MTC_haul")

# Compute spatially weighted means and models
results <- lapply(vars, function(v) spatially_weighted_mean_and_model(v, df_sp_utm, grid))
spatial_means_df <- data.frame(
  variable = vars,
  spatially_weighted_mean = sapply(results, function(x) x$mean),
  variogram_model = sapply(results, function(x) x$model)
)

print(spatial_means_df)



vgm_model <- variogram(haul_avg_troph ~ 1, data = df_sp_utm) 
plot(vgm_model)
vgm_model <- variogram(species_richness ~ 1, data = df_sp_utm) 
plot(vgm_model)
vgm_model <- variogram(evenness ~ 1, data = df_sp_utm) 
plot(vgm_model)
vgm_model <- variogram(haul_total_cpua ~ 1, data = df_sp_utm) 
plot(vgm_model)

#############################################################################################
####### try across entire dataset calculating yearly spatially weighted averages ############
#############################################################################################

vars <- c("haul_avg_troph", "haul_avg_weight", "species_richness", "genus_richness", "evenness", "haul_total_cpua", "body_weight_coeff", "MTC_haul")
years <- sort(unique(SCS_haul_stats$year))

# Function to fit variogram and extract range, with specified scale
extract_variogram_range <- function(varname, sp_points) {
  # Remove NAs
  valid_idx <- !is.na(sp_points@data[[varname]])
  sp_points_valid <- sp_points[valid_idx, ]
  if (nrow(sp_points_valid) < 5) return(NA)
  vgm_emp <- variogram(
    as.formula(paste(varname, "~ 1")),
    data = sp_points_valid,
    cutoff = 100,   # maximum distance (km)
    width = 10      # bin width (km)
  )
  # Try several models
  models <- c("Sph", "Exp", "Gau", "Lin")
  for (m in models) {
    fit <- try(fit.variogram(vgm_emp, model = vgm(m)), silent = TRUE)
    if (!inherits(fit, "try-error") && !any(is.na(fit))) {
      if (any(fit$range <= 0) || any(fit$psill < 0)) next
      return(fit$range[2]) # The second row is the structural component
    }
  }
  return(NA)
}

# Main loop over years and variables
range_results <- data.frame()
for (yr in years) {
  df_year <- SCS_haul_stats %>%
    filter(year == yr)
  if (nrow(df_year) < 5) next
  df_sf <- st_as_sf(df_year, coords = c("longitude", "latitude"), crs = 4326)
  df_sp <- as(df_sf, "Spatial")
  df_sp_utm <- spTransform(df_sp, CRS("+proj=utm +zone=20 +datum=WGS84 +units=km"))
  for (v in vars) {
    range_val <- extract_variogram_range(v, df_sp_utm)
    range_results <- rbind(range_results, data.frame(
      year = yr,
      variable = v,
      range_km = range_val
    ))
  }
}

print(range_results)

# Optional: Plot the range over years for each variable
ggplot(range_results, aes(x = year, y = range_km, color = variable)) +
  geom_point() + geom_line() +
  labs(title = "Estimated Range of Spatial Autocorrelation (up to 100km, 10km bins)",
       x = "Year", y = "Range (km)", color = "Variable") +
  theme_minimal()

############# this would suggest no spatial autocorrelation for all variables for all time within a year so check for spatial autocorrelation another way #############

vars <- c("haul_avg_troph", "haul_avg_weight", "species_richness", "genus_richness", "evenness", "haul_total_cpua", "body_weight_coeff", "MTC_haul")
years <- sort(unique(SCS_haul_stats$year))

moran_results <- data.frame()

for (yr in years) {
  df_year <- SCS_haul_stats %>%
    filter(year == yr)
  if (nrow(df_year) < 5) next  # skip years with too few points
  coords <- as.matrix(df_year[, c("longitude", "latitude")])
  # Create spatial weights (k-nearest neighbors, k=4)
  nb <- knn2nb(knearneigh(coords, k = min(4, nrow(df_year)-1)))
  lw <- nb2listw(nb, style = "W")
  for (v in vars) {
    if (all(is.na(df_year[[v]]))) next
    x <- df_year[[v]]
    if (sum(!is.na(x)) < 5) next
    mi <- moran.test(x, lw, na.action = na.omit)
    moran_results <- rbind(moran_results, data.frame(
      year = yr,
      variable = v,
      moran_I = mi$estimate[1],
      expected_I = mi$estimate[2],
      sd_I = mi$estimate[3],
      p_value = mi$p.value
    ))
  }
}

print(moran_results)

significant_moran <- moran_results %>%
  filter(p_value < 0.05)

print(significant_moran)

######### try to reconcile kriging and morins I #############################

# Function to fit variogram and extract range
extract_variogram_range <- function(varname, sp_points) {
  # Remove NAs
  valid_idx <- !is.na(sp_points@data[[varname]])
  sp_points_valid <- sp_points[valid_idx, ]
  if (nrow(sp_points_valid) < 5) return(NA)
  vgm_emp <- variogram(as.formula(paste(varname, "~ 1")), data = sp_points_valid)
  # Try several models
  models <- c("Sph", "Exp", "Gau", "Lin")
  for (m in models) {
    fit <- try(fit.variogram(vgm_emp, model = vgm(m)), silent = TRUE)
    if (!inherits(fit, "try-error") && !any(is.na(fit))) {
      if (any(fit$range <= 0) || any(fit$psill < 0)) next
      return(fit$range[2]) # The second row is the structural component
    }
  }
  return(NA)
}

# Main loop over years and variables
range_results <- data.frame()
for (yr in years) {
  df_year <- SCS_haul_stats %>%
    filter(year == yr)
  if (nrow(df_year) < 5) next
  df_sf <- st_as_sf(df_year, coords = c("longitude", "latitude"), crs = 4326)
  df_sp <- as(df_sf, "Spatial")
  df_sp_utm <- spTransform(df_sp, CRS("+proj=utm +zone=20 +datum=WGS84 +units=km"))
  for (v in vars) {
    range_val <- extract_variogram_range(v, df_sp_utm)
    range_results <- rbind(range_results, data.frame(
      year = yr,
      variable = v,
      range_km = range_val
    ))
  }
}

print(range_results)

# Optional: Plot the range over years for each variable
ggplot(range_results, aes(x = year, y = range_km, color = variable)) +
  geom_point() + geom_line() +
  labs(title = "Estimated Range of Spatial Autocorrelation by Year",
       x = "Year", y = "Range (km)", color = "Variable") +
  theme_minimal()

#############################################################################
library(tidyverse)
library(sf)
library(gstat)
library(sp)

# 1. Define variables and a container for results
vars_to_analyze <- c("haul_avg_troph", "species_richness", "evenness", "haul_total_cpua")
all_years_results <- list()

# 2. Start the Loop
unique_years <- sort(unique(SCS_haul_stats$year))

for (yr in unique_years) {
  cat("--- Processing Year:", yr, "---\n")
  
  # Filter data for the specific year
  yearly_data <- SCS_haul_stats %>% filter(year == yr)
  
  # Convert to SF and Project to UTM 20N (Crucial for distance in km)
  # Ensure your coordinates are named 'longitude' and 'latitude'
  yearly_sf <- st_as_sf(yearly_data, coords = c("longitude", "latitude"), crs = 4326)
  yearly_sp_utm <- as(st_transform(yearly_sf, 32620), "Spatial") 
  
  # Define prediction grid (5km resolution)
  grid_res <- 5 
  bbox_yr <- bbox(yearly_sp_utm)
  grid_yr <- expand.grid(
    x = seq(bbox_yr[1,1], bbox_yr[1,2], by = grid_res),
    y = seq(bbox_yr[2,1], bbox_yr[2,2], by = grid_res)
  )
  coordinates(grid_yr) <- ~x+y
  proj4string(grid_yr) <- CRS(proj4string(yearly_sp_utm))
  
  # Process each variable for this year
  year_stats <- map_df(vars_to_analyze, function(v) {
    
    # Calculate Empirical Variogram
    vgm_emp <- variogram(as.formula(paste(v, "~ 1")), data = yearly_sp_utm)
    
    # Auto-fit (Trying Spherical and Exponential models)
    fit_res <- try(fit.variogram(vgm_emp, model = vgm(c("Sph", "Exp"))), silent = TRUE)
    
    if (inherits(fit_res, "try-error")) {
      return(tibble(year = yr, variable = v, spatial_range_km = NA, weighted_mean = NA))
    }
    
    # Perform Kriging
    kriged <- krige(as.formula(paste(v, "~ 1")), yearly_sp_utm, grid_yr, model = fit_res)
    
    # Extract Range and Spatially Weighted Mean
    # Note: 'range' in gstat is the distance parameter of the model
    tibble(
      year = yr,
      variable = v,
      spatial_range_km = fit_res$range[2], 
      weighted_mean = mean(kriged$var1.pred, na.rm = TRUE)
    )
  })
  
  all_years_results[[as.character(yr)]] <- year_stats
}

# 3. Combine results into a single data frame
kriging_time_series <- bind_rows(all_years_results)














