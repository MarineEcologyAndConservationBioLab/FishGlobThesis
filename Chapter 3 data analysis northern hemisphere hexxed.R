#____________________________________#
######## workspace setup  ########
#____________________________________#

library(readr)
library(dplyr)
library(ggplot2)
library(grid)
library(lubridate) # to add year, month columns to datasets with time / date character columns
library(purrr)
library(rnaturalearth) # Used to get a cleaner map background
library(sf)
library(stringr)
library(tidyr)
library(vegan)
library(worrms)

setwd("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/Additional survey datasets")

#____________________________________#
##### import primary datasets and standardize metadata  #######
#____________________________________#


#____________________________________#
######## Initialize FISHGLOB  ########
#_____________________________________#

load("C:/Users/besch/Desktop/flagged FISHGLOB/FishGlob_public_std_clean.RData")
FishGlob_hexxed <- data
# 1. Identify the species names used in Nor-BTS that lack AphiaIDs
nor_bts_species <- FishGlob_hexxed %>%
  filter(survey_unit == "Nor-BTS", is.na(aphia_id) | aphia_id == 0 | aphia_id == "") %>%
  pull(accepted_name) %>%
  unique()
# Unique Nor-BTS names
nor_bts_names <- FishGlob_hexxed %>%
  filter(survey_unit == "Nor-BTS") %>%
  distinct(accepted_name) %>%
  filter(!is.na(accepted_name))

# Query WoRMS

nor_bts_lookup <- map_dfr(
  nor_bts_names$accepted_name,
  function(sp){
    
    tryCatch({
      
      # exact scientific-name search
      res <- wm_records_name(sp)
      
      if(length(res) == 0){
        
        # fallback fuzzy search
        res <- wm_records_taxamatch(sp)
        
      }
      
      if(length(res) == 0){
        
        return(
          tibble(
            accepted_name = sp,
            aphia_id_new = NA_character_
          )
        )
        
      }
      
      res <- bind_rows(res)
      
      if(nrow(res) == 0){
        
        return(
          tibble(
            accepted_name = sp,
            aphia_id_new = NA_character_
          )
        )
        
      }
      
      accepted_row <- res %>%
        filter(status == "accepted")
      
      if(nrow(accepted_row) == 0){
        accepted_row <- res[1, ]
      }
      
      tibble(
        accepted_name = sp,
        aphia_id_new =
          as.character(
            accepted_row$AphiaID[1]
          )
      )
      
    }, error = function(e){
      
      tibble(
        accepted_name = sp,
        aphia_id_new = NA_character_
      )
      
    })
    
  }
)
# Update the main dataset
FishGlob_hexxed <- FishGlob_hexxed %>%
  mutate(
    aphia_id = as.character(aphia_id)
  ) %>%
  left_join(
    nor_bts_lookup,
    by = "accepted_name"
  ) %>%
  mutate(
    aphia_id = case_when(
      is.na(aphia_id) ~ aphia_id_new,
      aphia_id == "" ~ aphia_id_new,
      aphia_id == "0" ~ aphia_id_new,
      TRUE ~ aphia_id
    )
  ) %>%
  select(-aphia_id_new)

missing_aphia_nor <- FishGlob_hexxed%>%
  filter(survey_unit == "Nor-BTS")%>%
  summarize(
    n_species = length(unique(accepted_name)),
    n_aphia = length(unique(aphia_id))
  )
#____________________________________#
# ____ Initialize EU Med and Black Sea ____#
#____________________________________#

EU_Med_Black_Sea <- read_csv("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/Additional survey datasets/Med and Black Sea EU/2024_MEDBSsurvey/Demersal/TA.csv")

TB <- read_csv("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/Additional survey datasets/Med and Black Sea EU/2024_MEDBSsurvey/Demersal/TB.csv")
TA <- read_csv("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/Additional survey datasets/Med and Black Sea EU/2024_MEDBSsurvey/Demersal/TA.csv")
identifiers <- intersect(colnames(TB), colnames(TA))

combined_EU_data <- left_join(TB, TA, by = identifiers)

dim(combined_EU_data)

combined_EU_data <- combined_EU_data %>%
  filter(validity == "V" & catfau %in% c("A", "Aa", "Ae","Ao")) %>% # removes invalid hauls and keeps only Fish, Fish Agnatha, Fish Chondrichthyes, Fish Osteichthyes
  mutate (
    haul_id = paste0(name_of_survey, "_", country, "_", area, "_", vessel, "_", year, "_", month, "_", day, "_", haul_number),
  ) 

combined_EU_data <- combined_EU_data %>%
  rename (
    num = nbtot,
    wgt = ptot,
    survey_unit = name_of_survey,
    haul_dur = distance, 
    timestamp = shooting_time, 
    sub_area = area
  ) %>%
  mutate(
    latitude = (shooting_latitude + hauling_latitude) / 200,
    longitude = (shooting_longitude + hauling_longitude) / 200,
    depth = (shooting_depth + hauling_depth) / 2,
    sbt = (bottom_temperature_beginning + bottom_temperature_end) / 2,
    wgt_cpue = wgt / haul_dur,
    num_cpue = num / haul_dur,
    source = "EU_JRC", 
    survey = "EU_Med_Black_Sea", 
    continent = "Europe",
    combined_MEDITS_code = paste0(genus,species)
  ) 

MEDITS_code_guide <- read_csv("Med and Black Sea EU/MEDITS_code_guide.txt") #import naming convention of species to ID fish taxa
full_MEDITS_code_guide <- read_table("Med and Black Sea EU/full_MEDITS_code_guide.txt")
MEDITS_code_key <- full_MEDITS_code_guide %>%
  rename(
    combined_MEDITS_code = MEDITS
  ) %>%
  mutate(
    accepted_name = paste0(Code, " ", Scientific)
  ) %>%
  select(
    combined_MEDITS_code, accepted_name)

combined_EU_data <- left_join(combined_EU_data, MEDITS_code_key, by = "combined_MEDITS_code")

combined_EU_data_FISHGLOB_formatted <- combined_EU_data %>%
  select(
    survey, source, timestamp, haul_id, country, sub_area, continent, year, month, day, latitude, longitude, haul_dur, gear, depth, sbt, num, num_cpue, wgt, wgt_cpue, survey_unit, accepted_name 
  )


unique_species <- unique(combined_EU_data_FISHGLOB_formatted$accepted_name) %>%
  unique() %>%
  na.omit() 

build_worms_taxonomy <- function(species_names){
  
  species_names <- species_names %>%
    unique() %>%
    na.omit() %>%
    trimws()
  
  purrr::map_df(species_names, function(sp){
    
    tryCatch({
      
      cat("\nLOOKING UP:", sp, "\n")
      
      res <- wm_records_taxamatch(sp)[[1]]
      
      print(res)
      
      if(is.null(res) || nrow(res) == 0){
        
        cat("NO MATCH\n")
        
        return(
          tibble(
            accepted_name = sp,
            aphia_id = NA_character_,
            kingdom = NA_character_,
            phylum = NA_character_,
            class = NA_character_,
            order = NA_character_,
            family = NA_character_,
            genus = NA_character_,
            species = NA_character_,
            rank = NA_character_
          )
        )
        
      }
      
      match <- res[1, ]
      
      print(match)
      
      tibble(
        accepted_name = sp,
        aphia_id = as.character(match$AphiaID),
        kingdom = match$kingdom,
        phylum = match$phylum,
        class = match$class,
        order = match$order,
        family = match$family,
        genus = match$genus,
        species = match$scientificname,
        rank = match$rank
      )
      
    }, error = function(e){
      
      cat("ERROR:", e$message, "\n")
      
      tibble(
        accepted_name = sp,
        aphia_id = NA_character_,
        kingdom = NA_character_,
        phylum = NA_character_,
        class = NA_character_,
        order = NA_character_,
        family = NA_character_,
        genus = NA_character_,
        species = NA_character_,
        rank = NA_character_
      )
      
    })
    
  })
  
}
taxonomy_lookup <- build_worms_taxonomy(unique_species)
combined_EU_data_FISHGLOB_formatted <- combined_EU_data_FISHGLOB_formatted %>%
  left_join(taxonomy_lookup, by = "accepted_name")

MEDITS <- combined_EU_data_FISHGLOB_formatted %>%
  filter(survey_unit == "MEDITS" & month %in% c(5, 6, 7, 8, 9)) # MEDITS is primarily conducted in May-September)
SOLEMON <- combined_EU_data_FISHGLOB_formatted %>%
  filter(survey_unit == "SOLEMON")

'
ROMANIA_SPR <- combined_EU_data_FISHGLOB_formatted %>%
  filter(survey_unit == "spring")
ROMANIA_AUT <- combined_EU_data_FISHGLOB_formatted %>%
  filter(survey_unit == "autumn")
BTSBS_AUT <- combined_EU_data_FISHGLOB_formatted %>%
  filter(survey_unit == "BTSBS-AUT")
BTSBS_SPR <- combined_EU_data_FISHGLOB_formatted %>%
  filter(survey_unit == "BTSBS-SPR")
  
# these surveys only have 1 species in them Scophthalmus maximus	 so no community level thermal affinity
'
#____________________________________#
##### Initialize Faroese trawl survey   #######
#____________________________________#

FaroeGroundFishData_schatteman <- read_csv("~/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/Additional survey datasets/Faroe islands/final data receipt/FaroeGroundFishData_schatteman.csv")
'
colnames(FaroeGroundFishData_schatteman)
colnames(FishGlob_public_clean)
'

FaroeseSurveys <- FaroeGroundFishData_schatteman %>%
  mutate(
    year = year(timestamp),
    month = month(timestamp),
    day = day(timestamp),
    num_cpue = num / haul_dur, 
    wgt_cpue = wgt / haul_dur
  )

# 1. Clean and prepare unique names
unique_species <- unique(FaroeseSurveys$accepted_name) %>%
  unique() %>%
  na.omit()

build_worms_taxonomy <- function(species_names){
  
  species_names <- species_names %>%
    unique() %>%
    na.omit() %>%
    trimws()
  
  purrr::map_df(species_names, function(sp){
    
    tryCatch({
      
      cat("\nLOOKING UP:", sp, "\n")
      
      res <- wm_records_taxamatch(sp)[[1]]
      
      print(res)
      
      if(is.null(res) || nrow(res) == 0){
        
        cat("NO MATCH\n")
        
        return(
          tibble(
            accepted_name = sp,
            aphia_id = NA_character_,
            kingdom = NA_character_,
            phylum = NA_character_,
            class = NA_character_,
            order = NA_character_,
            family = NA_character_,
            genus = NA_character_,
            species = NA_character_,
            rank = NA_character_
          )
        )
        
      }
      
      match <- res[1, ]
      
      print(match)
      
      tibble(
        accepted_name = sp,
        kingdom = match$kingdom,
        phylum = match$phylum,
        class = match$class,
        order = match$order,
        family = match$family,
        genus = match$genus,
        species = match$scientificname,
        rank = match$rank
      )
      
    }, error = function(e){
      
      cat("ERROR:", e$message, "\n")
      
      tibble(
        accepted_name = sp,
        aphia_id = NA_character_,
        kingdom = NA_character_,
        phylum = NA_character_,
        class = NA_character_,
        order = NA_character_,
        family = NA_character_,
        genus = NA_character_,
        species = NA_character_,
        rank = NA_character_
      )
      
    })
    
  })
  
}


taxonomy_lookup <- build_worms_taxonomy(unique_species)

FaroeseSurveys_FISHGLOB_formatted <- FaroeseSurveys %>%
  left_join(taxonomy_lookup, by = "accepted_name") %>%
  mutate(accepted_name = species)

"
SURVEY_UNIT (4 unique):	
yvirlitstrol_lg_var = overviewtrawl Faroe Plateau spring
yvirlitstrol_lg_summar = overviewtrawl Faroe Plateau summer
yvirlitstrol_fb_var = overviewtrawl Faroe Bank spring
yvirlitstrol_fb_summar = overviewtrawl Faroe Bank summer

var = spring (which is a bit misleading, since the survey is conducted in February-March)
"

FaroePlSpr <- FaroeseSurveys_FISHGLOB_formatted %>%
  filter(survey_unit == "yvirlitstrol_lg_var")
FaroePlSum <- FaroeseSurveys_FISHGLOB_formatted %>%
  filter(survey_unit == "yvirlitstrol_lg_summar")
FaroeBaSpr <- FaroeseSurveys_FISHGLOB_formatted %>%
  filter(survey_unit == "yvirlitstrol_fb_var")
FaroeBaSum <- FaroeseSurveys_FISHGLOB_formatted %>%
  filter(survey_unit == "yvirlitstrol_fb_summar")

#____________________________________#
##### Initialize GINR Greenland Fish-Shrimp West 2007-2023  #######
#____________________________________#

GINR_fish_shrimp_west <- read_delim("Greenland Greenland Institute of Natural Resources/df_SF.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)
GINR_fish_shrimp_west <- GINR_fish_shrimp_west %>%
  rename(
    year = Year,
    station = Station,
    depth = DEPTH,
    haul_dur = `Haul Duration`,
    area_swept = `SWEPTAREA_KM2`, 
    sbt = TEMP_BOT,
    sst = `TEMP SUR`,
    accepted_name = `NameLatin`, 
    wgt = Weight_Kg, 
    num = Number, 
    latitude = LAT, 
    longitude = LON, 
    gear = GearID, 
    haul_id = KeySt,
  ) %>%
  mutate(
    wgt_cpua = wgt/area_swept, 
    num_cpua = num/area_swept, 
    wgt_cpue = wgt/haul_dur,
    num_cpue = num/haul_dur,
    Survey = "GINR fish shrimp", # Assign survey name, 
    survey_unit = "GINR fish shrimp west",
    source = "GINR", 
    country = "Greenland", 
    month = month(DATE), # This gives numbers (e.g., 7),
    day = day(DATE)
  )

# 1. Clean and prepare unique names from the Bering Sea data
unique_species <- GINR_fish_shrimp_west$accepted_name %>%
  unique() %>%
  na.omit()

build_worms_taxonomy <- function(species_names){
  
  species_names <- species_names %>%
    unique() %>%
    na.omit() %>%
    trimws()
  
  purrr::map_df(species_names, function(sp){
    
    tryCatch({
      
      cat("\nLOOKING UP:", sp, "\n")
      
      res <- wm_records_taxamatch(sp)[[1]]
      
      print(res)
      
      if(is.null(res) || nrow(res) == 0){
        
        cat("NO MATCH\n")
        
        return(
          tibble(
            accepted_name = sp,
            aphia_id = NA_character_,
            kingdom = NA_character_,
            phylum = NA_character_,
            class = NA_character_,
            order = NA_character_,
            family = NA_character_,
            genus = NA_character_,
            species = NA_character_,
            rank = NA_character_
          )
        )
        
      }
      
      match <- res[1, ]
      
      print(match)
      
      tibble(
        accepted_name = sp,
        aphia_id = match$AphiaID,
        kingdom = match$kingdom,
        phylum = match$phylum,
        class = match$class,
        order = match$order,
        family = match$family,
        genus = match$genus,
        species = match$scientificname
      )
      
    }, error = function(e){
      
      cat("ERROR:", e$message, "\n")
      
      tibble(
        accepted_name = sp,
        aphia_id = NA_character_,
        kingdom = NA_character_,
        phylum = NA_character_,
        class = NA_character_,
        order = NA_character_,
        family = NA_character_,
        genus = NA_character_,
        species = NA_character_
      )
      
    })
    
  })
  
}


taxonomy_lookup <- build_worms_taxonomy(unique_species)

GINR_fish_shrimp_west <- GINR_fish_shrimp_west %>%
  left_join(taxonomy_lookup, by = "accepted_name")

#____________________________________#
##### Initialize GINR Greenland Fish-Shrimp East  #######
#____________________________________#

GINR_fish_shrimp_east <- read_delim("GINR fish shrimp east/df_SFE.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)

GINR_fish_shrimp_east <- GINR_fish_shrimp_east %>%
  rename(
    year = Year,
    station = Station,
    depth = DEPTH,
    haul_dur = `Haul Duration`,
    area_swept = `SWEPTAREA_KM2`, 
    sbt = TEMP_BOT,
    sst = `TEMP SUR`,
    accepted_name = `NameLatin`, 
    wgt = Weight_Kg, 
    num = Number, 
    latitude = LAT, 
    longitude = LON, 
    gear = GearID, 
    haul_id = KeySt,
  ) %>%
  mutate(
    wgt_cpua = wgt/area_swept, 
    num_cpua = num/area_swept, 
    wgt_cpue = wgt/haul_dur,
    num_cpue = num/haul_dur,
    Survey = "GINR fish shrimp", # Assign survey name, 
    survey_unit = "GINR fish shrimp east",
    source = "GINR", 
    country = "Greenland", 
    month = month(DATE), # This gives numbers (e.g., 7),
    day = day(DATE)
  )

# 1. Clean and prepare unique names from the Bering Sea data
unique_species <- GINR_fish_shrimp_east$accepted_name %>%
  unique() %>%
  na.omit()

build_worms_taxonomy <- function(species_names){
  
  species_names <- species_names %>%
    unique() %>%
    na.omit() %>%
    trimws()
  
  purrr::map_df(species_names, function(sp){
    
    tryCatch({
      
      cat("\nLOOKING UP:", sp, "\n")
      
      res <- wm_records_taxamatch(sp)[[1]]
      
      print(res)
      
      if(is.null(res) || nrow(res) == 0){
        
        cat("NO MATCH\n")
        
        return(
          tibble(
            accepted_name = sp,
            aphia_id = NA_character_,
            kingdom = NA_character_,
            phylum = NA_character_,
            class = NA_character_,
            order = NA_character_,
            family = NA_character_,
            genus = NA_character_,
            species = NA_character_,
            rank = NA_character_
          )
        )
        
      }
      
      match <- res[1, ]
      
      print(match)
      
      tibble(
        accepted_name = sp,
        aphia_id = match$AphiaID,
        kingdom = match$kingdom,
        phylum = match$phylum,
        class = match$class,
        order = match$order,
        family = match$family,
        genus = match$genus,
        species = match$scientificname
      )
      
    }, error = function(e){
      
      cat("ERROR:", e$message, "\n")
      
      tibble(
        accepted_name = sp,
        aphia_id = NA_character_,
        kingdom = NA_character_,
        phylum = NA_character_,
        class = NA_character_,
        order = NA_character_,
        family = NA_character_,
        genus = NA_character_,
        species = NA_character_
      )
      
    })
    
  })
  
}


taxonomy_lookup <- build_worms_taxonomy(unique_species)

GINR_fish_shrimp_east <- GINR_fish_shrimp_east %>%
  left_join(taxonomy_lookup, by = "accepted_name")

#____________________________________#
##### Initialize NOAA Northeastern Bering Sea (Extension)  #######
#____________________________________#

NE_Bering_Sea <- read_csv("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/Additional survey datasets/NOAA Northeastern Bering Sea extension/CATCH AND HAUL DATA .csv")

NE_Bering_Sea <- NE_Bering_Sea %>%
  rename(
    year = `Survey year`,
    haul_id = `Haul ID`,
    station = `Station ID`,
    stratum = `Stratum ID`,
    depth = `Depth (m)`,
    haul_dur = `Tow duration (decimal hr)`,
    area_swept = `Area swept (km)`, 
    sbt = `Bottom temperature (degrees Celsius)`,
    sst = `Surface temperature (degrees Celsius)`,
    verbatim_name = `Taxon common name`,
    accepted_name = `Taxon scientific name`, 
    aphia_id = `World register of marine species (WoRMS) taxonomic serial number`,
    wgt = `Taxon weight (kg)`, 
    num = `Taxon count`, 
    wgt_cpua = `Weight CPUE (kg/km2)`, 
    num_cpua = `Number CPUE (no/km2)`, 
  ) %>%
  mutate(
    Survey = "Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension", # Assign survey name, 
    survey_unit = "Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension",
    source = "NOAA", 
    country = "United States", 
    latitude = (`Start latitude (decimal degrees)` + `End latitude (decimal degrees)`)/2,
    longitude = (`Start longitude (decimal degrees)` + `End longitude (decimal degrees)`)/2,
    # 1. Convert the text to a real Date object first
    date_formatted = dmy_hm(`Date and time`),
    # 2. Extract the Year and Month
    year  = year(date_formatted),
    month = month(date_formatted), # This gives numbers (e.g., 7),
  )

# 1. Clean and prepare unique names from the Bering Sea data
unique_species <- NE_Bering_Sea$accepted_name %>%
  unique() %>%
  na.omit()

build_worms_taxonomy <- function(species_names){
  
  species_names <- species_names %>%
    unique() %>%
    na.omit() %>%
    trimws()
  
  purrr::map_df(species_names, function(sp){
    
    tryCatch({
      
      cat("\nLOOKING UP:", sp, "\n")
      
      res <- wm_records_taxamatch(sp)[[1]]
      
      print(res)
      
      if(is.null(res) || nrow(res) == 0){
        
        cat("NO MATCH\n")
        
        return(
          tibble(
            accepted_name = sp,
            aphia_id = NA_character_,
            kingdom = NA_character_,
            phylum = NA_character_,
            class = NA_character_,
            order = NA_character_,
            family = NA_character_,
            genus = NA_character_,
            species = NA_character_,
            rank = NA_character_
          )
        )
        
      }
      
      match <- res[1, ]
      
      print(match)
      
      tibble(
        accepted_name = sp,
        kingdom = match$kingdom,
        phylum = match$phylum,
        class = match$class,
        order = match$order,
        family = match$family,
        genus = match$genus,
        species = match$scientificname
      )
      
    }, error = function(e){
      
      cat("ERROR:", e$message, "\n")
      
      tibble(
        accepted_name = sp,
        kingdom = NA_character_,
        phylum = NA_character_,
        class = NA_character_,
        order = NA_character_,
        family = NA_character_,
        genus = NA_character_,
        species = NA_character_
      )
      
    })
    
  })
  
}


taxonomy_lookup <- build_worms_taxonomy(unique_species)

NE_Bering_Sea <- NE_Bering_Sea %>%
  left_join(taxonomy_lookup, by = "accepted_name")

#____________________________________#
##### Initialize Eastern Bering Slope  #######
#____________________________________#

Bering_Slope <- read_csv("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/Additional survey datasets/NOAA Eastern Bering Sea Upper Continental Slope/CATCH AND HAUL DATA .csv")

Bering_Slope <- Bering_Slope %>%
  rename(
    year = `Survey year`,
    haul_id = `Haul ID`,
    station = `Station ID`,
    stratum = `Stratum ID`,
    depth = `Depth (m)`,
    haul_dur = `Tow duration (decimal hr)`,
    area_swept = `Area swept (km)`, 
    sbt = `Bottom temperature (degrees Celsius)`,
    sst = `Surface temperature (degrees Celsius)`,
    verbatim_name = `Taxon common name`,
    accepted_name = `Taxon scientific name`, 
    aphia_id = `World register of marine species (WoRMS) taxonomic serial number`,
    wgt = `Taxon weight (kg)`, 
    num = `Taxon count`, 
    wgt_cpua = `Weight CPUE (kg/km2)`, 
    num_cpua = `Number CPUE (no/km2)`, 
  ) %>%
  mutate(
    Survey = "Eastern Bering Sea Slope Bottom Trawl Survey", # Assign survey name, 
    survey_unit = "Eastern Bering Sea Slope Bottom Trawl Survey",
    source = "NOAA", 
    country = "United States", 
    latitude = (`Start latitude (decimal degrees)` + `End latitude (decimal degrees)`)/2,
    longitude = (`Start longitude (decimal degrees)` + `End longitude (decimal degrees)`)/2,
    # 1. Convert the text to a real Date object first
    date_formatted = dmy_hm(`Date and time`),
    # 2. Extract the Year and Month
    year  = year(date_formatted),
    month = month(date_formatted), # This gives numbers (e.g., 7),
  )

# 1. Clean and prepare unique names from the Bering Sea data
unique_species <- Bering_Slope$accepted_name %>%
  unique() %>%
  na.omit()

build_worms_taxonomy <- function(species_names){
  
  species_names <- species_names %>%
    unique() %>%
    na.omit() %>%
    trimws()
  
  purrr::map_df(species_names, function(sp){
    
    tryCatch({
      
      cat("\nLOOKING UP:", sp, "\n")
      
      res <- wm_records_taxamatch(sp)[[1]]
      
      print(res)
      
      if(is.null(res) || nrow(res) == 0){
        
        cat("NO MATCH\n")
        
        return(
          tibble(
            accepted_name = sp,
            aphia_id = NA_character_,
            kingdom = NA_character_,
            phylum = NA_character_,
            class = NA_character_,
            order = NA_character_,
            family = NA_character_,
            genus = NA_character_,
            species = NA_character_,
            rank = NA_character_
          )
        )
        
      }
      
      match <- res[1, ]
      
      print(match)
      
      tibble(
        accepted_name = sp,
        kingdom = match$kingdom,
        phylum = match$phylum,
        class = match$class,
        order = match$order,
        family = match$family,
        genus = match$genus,
        species = match$scientificname
      )
      
    }, error = function(e){
      
      cat("ERROR:", e$message, "\n")
      
      tibble(
        accepted_name = sp,
        kingdom = NA_character_,
        phylum = NA_character_,
        class = NA_character_,
        order = NA_character_,
        family = NA_character_,
        genus = NA_character_,
        species = NA_character_
      )
      
    })
    
  })
  
}


taxonomy_lookup <- build_worms_taxonomy(unique_species)

Bering_Slope <- Bering_Slope %>%
  left_join(taxonomy_lookup, by = "accepted_name")

#____________________________________#
##### Initialize DFO Newfoundland (Newfoundland)  #######
#____________________________________#

# Load DFO Newfoundland Data
EDR_1996_2019_OldVessels <- read_csv("1996-2022 EDR DFO Newfoundland/1996-2022 EDR/EDR_1996-2019_OldVessels.csv")
EDR_1996_2019_OldVessels <- EDR_1996_2019_OldVessels %>%
  filter(EndLatitude != 0 & EndLongitude != 0) # Remove rows with zero coordinates)
EDR_1996_2019_OldVessels$Latitude <- (EDR_1996_2019_OldVessels$EndLatitude + EDR_1996_2019_OldVessels$StartLatitude)/2000
EDR_1996_2019_OldVessels$Longitude <- (EDR_1996_2019_OldVessels$EndLongitude + EDR_1996_2019_OldVessels$StartLongitude)/-2000

EDR_2020_2022_NewVessels <- read_csv("1996-2022 EDR DFO Newfoundland/1996-2022 EDR/EDR_2020-2022_NewVessels.csv")
EDR_2020_2022_NewVessels$Latitude <- (EDR_2020_2022_NewVessels$EndLatitude + EDR_2020_2022_NewVessels$StartLatitude)/200
EDR_2020_2022_NewVessels$Longitude <- (EDR_2020_2022_NewVessels$EndLongitude + EDR_2020_2022_NewVessels$StartLongitude)/-200

EDR_2020_2022_OldVessels <- read_csv("1996-2022 EDR DFO Newfoundland/1996-2022 EDR/EDR_2020-2022_OldVessels.csv")
EDR_2020_2022_OldVessels$Latitude <- (EDR_2020_2022_OldVessels$EndLatitude + EDR_2020_2022_OldVessels$StartLatitude)/200
EDR_2020_2022_OldVessels$Longitude <- (EDR_2020_2022_OldVessels$EndLongitude + EDR_2020_2022_OldVessels$StartLongitude)/-200

class(EDR_1996_2019_OldVessels$MidPointDateTime)
class(EDR_2020_2022_NewVessels$MidPointDateTime)
class(EDR_2020_2022_OldVessels$MidPointDateTime)

# Convert the character columns to POSIXct
# Convert the New Vessels
EDR_2020_2022_NewVessels$MidPointDateTime <- as.POSIXct(EDR_2020_2022_NewVessels$MidPointDateTime, 
                                                        format = "%m/%d/%Y %H:%M")

# Convert the Old Vessels (2020-2022)
EDR_2020_2022_OldVessels$MidPointDateTime <- as.POSIXct(EDR_2020_2022_OldVessels$MidPointDateTime, 
                                                        format = "%m/%d/%Y %H:%M")

# Now combine them!
NewfoundlandSurvey <- rbind(EDR_1996_2019_OldVessels, 
                            EDR_2020_2022_NewVessels, 
                            EDR_2020_2022_OldVessels)

NewfoundlandSurvey_FISHGLOB_formatted <- NewfoundlandSurvey %>%
  rename_with(tolower) 
#colnames(NewfoundlandSurvey_FISHGLOB_formatted)

library(stringr)
NewfoundlandSurvey_FISHGLOB_formatted <- NewfoundlandSurvey_FISHGLOB_formatted %>%
  # 1. Rename: new_name = existing_name
  rename(
    sbt = bottomtemp,
    haul_dur = distance,
    depth = meandepth,
    num = totalcount,
    wgt = totalweight,
    gear = gearcodeid,
    common_name = speciesname,
    accepted_name = scientificname, # this will need to be fixed later 
    order = orders
  ) %>%
  # 2. Mutate: Create new columns or modify renamed ones
  mutate(
    survey  = "DFO-NL", 
    survey_unit = "DFO-NL",
    source  = "DFO", 
    month = month(midpointdatetime), # Extract month as a number (e.g., 4),
    num_cpue = num / haul_dur, 
    wgt_cpue = wgt / haul_dur, 
    haul_id = str_remove(uniqueid, "-[0-9]+$")
  )

# taxonomic repair
library(dplyr)
library(stringr)
library(purrr)
library(worrms)

# Build Newfoundland AphiaID lookup using all available taxonomy columns
# Priority:
# 1. ScientificName species match
# 2. ScientificName with comments removed
# 3. ScientificName "SP." -> Genus
# 4. Genus
# 5. Family
# 6. Orders
# 7. Class
# 8. Phylum

# Create unique taxonomic records
nl_taxa <- NewfoundlandSurvey_FISHGLOB_formatted %>%
  select(
    accepted_name,
    genus,
    family,
    order,
    class,
    phylum
  ) %>%
  distinct()

# Helper cleaning function
clean_taxon <- function(x){
  
  if(is.na(x) || x == "") return(NA_character_)
  
  x %>%
    toupper() %>%
    str_remove_all("\\([^\\)]*\\)") %>%   # remove brackets/explanations
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

# Robust WoRMS lookup
lookup_aphia <- function(accepted_name,
                         genus,
                         family,
                         order,
                         class,
                         phylum){
  
  candidates <- c()
  
  sci_clean <- clean_taxon(accepted_name)
  
  if(!is.na(sci_clean))
    candidates <- c(candidates, sci_clean)
  
  if(!is.na(sci_clean) &&
     str_detect(sci_clean, "\\bSP\\.?\\b")){
    candidates <- c(candidates, clean_taxon(genus))
  }
  
  candidates <- c(
    candidates,
    clean_taxon(genus),
    clean_taxon(family),
    clean_taxon(order),
    clean_taxon(class),
    clean_taxon(phylum)
  )
  
  candidates <- unique(
    candidates[
      !is.na(candidates) &
        candidates != "" &
        candidates != "INVERTEBRATA" &
        candidates != "MULTIPLE FAMILIES"
    ]
  )
  
  for(candidate in candidates){
    
    res <- tryCatch(
      wm_records_taxamatch(candidate)[[1]],
      error = function(e) NULL
    )
    
    if(!is.null(res) && nrow(res) > 0){
      
      accepted <- res %>%
        filter(tolower(status) == "accepted")
      
      if(nrow(accepted) == 0){
        accepted <- res[1, ]
      }
      
      return(
        tibble(
          accepted_name = accepted_name,
          aphia_id = as.character(accepted$AphiaID[1]),
          worms_accepted_name = accepted$scientificname[1],
          rank = accepted$rank[1],
          kingdom = accepted$kingdom[1],
          phylum_worms = accepted$phylum[1],
          class_worms = accepted$class[1],
          order_worms = accepted$order[1],
          family_worms = accepted$family[1],
          genus_worms = accepted$genus[1]
        )
      )
    }
  }
  
  tibble(
    accepted_name = accepted_name,
    aphia_id = NA_character_,
    worms_accepted_name = NA_character_,
    rank = NA_character_,
    kingdom = NA_character_,
    phylum_worms = NA_character_,
    class_worms = NA_character_,
    order_worms = NA_character_,
    family_worms = NA_character_,
    genus_worms = NA_character_
  )
}

nl_aphia_lookup <- purrr::pmap_df(
  nl_taxa,
  lookup_aphia
)

'
##############################################################################
# Run lookup once for each unique ScientificName
##############################################################################


##############################################################################
# Diagnostics
##############################################################################

cat(
  "\nUnique Newfoundland taxa:",
  nrow(nl_taxa),
  "\n"
)

cat(
  "Matched AphiaIDs:",
  sum(!is.na(nl_aphia_lookup$aphia_id)),
  "\n"
)

cat(
  "Unmatched taxa:",
  sum(is.na(nl_aphia_lookup$aphia_id)),
  "\n"
)

# Review failures if desired
unmatched_taxa <- nl_aphia_lookup %>%
  filter(is.na(aphia_id))

print(unmatched_taxa)
'


# Join back to Newfoundland dataset
# Adds ONLY aphia_id

# Join WoRMS taxonomy back to Newfoundland dataset

NewfoundlandSurvey_FISHGLOB_formatted <-
  NewfoundlandSurvey_FISHGLOB_formatted %>%
  left_join(
    nl_aphia_lookup,
    by = "accepted_name"
  ) %>%
  mutate(
    
    # Preserve original survey taxonomy
    original_name = accepted_name,
    
    # Replace accepted_name with WoRMS accepted name
    accepted_name = coalesce(
      worms_accepted_name,
      accepted_name
    ),
    
    # Optionally replace higher taxonomy with WoRMS taxonomy
    genus = coalesce(genus_worms, genus),
    family = coalesce(family_worms, family),
    order = coalesce(order_worms, order),
    class = coalesce(class_worms, class),
    phylum = coalesce(phylum_worms, phylum)
  ) %>%
  select(
    -worms_accepted_name,
    -genus_worms,
    -family_worms,
    -order_worms,
    -class_worms,
    -phylum_worms
  )
'
# Coverage check

NewfoundlandSurvey_FISHGLOB_formatted %>%
  distinct(original_name, accepted_name, rank) %>%
  filter(original_name != accepted_name) %>%
  print(n = 100)
NewfoundlandSurvey_FISHGLOB_formatted %>%
  count(rank, sort = TRUE)
'

# seasonal split
'
haul_month <- NewfoundlandSurvey_FISHGLOB_formatted %>%
  group_by(month) %>%
  summarise(n_haul = n_distinct(haul_id), .groups = "drop")

ggplot(haul_month, aes(x = factor(month, levels = 1:12), y = n_haul)) +
  geom_col(fill = "steelblue") +
  labs(
    x = "Month",
    y = "Number of Unique Haul IDs",
    title = "Unique Haul IDs per Month"
  ) +
  theme_bw()
'
NewfoundlandSurvey_FISHGLOB_formatted_spring <- NewfoundlandSurvey_FISHGLOB_formatted %>%
  filter(month %in% c(3, 4, 5, 6, 7)) %>% # March, April, May, June, July 
  mutate(survey_unit = "DFO-NL-Spring")

NewfoundlandSurvey_FISHGLOB_formatted_winter <- NewfoundlandSurvey_FISHGLOB_formatted %>%
  filter(month %in% c(9, 10, 11, 12, 1)) %>% # September, October, November, December, January
  mutate(survey_unit = "DFO-NL-Winter")

#____________________________________#
##### Initialize DFO Arctic (Greenland Halibut & Northern Shrimp) #######
#____________________________________#

## load trawl, species, catch, and length weight files to merge into final dataset and length weight
# USE Trawl_ID column to match hauls, OS_ID to match taxa

DFO_Arctic_trawl_data <- read_csv("DFO Arctic Greenland Halibut and Northern Shrimp/greenland halibut and northern shrimp trawl data.csv")

DFO_Arctic_trawl_data <- DFO_Arctic_trawl_data %>%
  mutate(`SweptArea(km2)` = `SweptArea(m2)` / 1000000)

# now do catch data 
DFO_Arctic_catch_data <- read_csv("DFO Arctic Greenland Halibut and Northern Shrimp/greenland halibut and northern shrimp catch data.csv")

DFO_Arctic_combined_catch_trawl <- left_join(DFO_Arctic_catch_data, DFO_Arctic_trawl_data, by = "Trawl_ID")

# now do species taxonomy data 
DFO_Arctic_species_info_data <- read_csv("DFO Arctic Greenland Halibut and Northern Shrimp/greenland halibut and northern shrimp species info data.csv")

DFO_Arctic_combined_catch_trawl_species <- left_join(DFO_Arctic_combined_catch_trawl, DFO_Arctic_species_info_data, by = "OS_ID")

unique(is.finite(c(DFO_Arctic_combined_catch_trawl_species$...16, DFO_Arctic_combined_catch_trawl_species$...15, DFO_Arctic_combined_catch_trawl_species$...14, DFO_Arctic_combined_catch_trawl_species$...13)))
DFO_Arctic_combined_catch_trawl_species <- DFO_Arctic_combined_catch_trawl_species %>%
  select(-...13, -...14, -...15, -...16)

# now do catch data 
DFO_Arctic_length_weight_data <- read_csv("DFO Arctic Greenland Halibut and Northern Shrimp/greenland halibut and northern shrimp length weight data.csv")

# Prepare Greenland Halibut Survey Data
Greenland_Halibut_Survey_Data <- DFO_Arctic_combined_catch_trawl_species %>%
  filter(Gear %in% c("Alfredo-03", 'Bacalao')) 

# Prepare Northern Shrimp Survey Data
Northern_Shrimp_Survey_Data <- DFO_Arctic_combined_catch_trawl_species %>%
  filter(Gear %in% c("Cosmos 2000", '21" Campelen', 'STD Campelen')) 

Greenland_Halibut_Survey_data_FISHGLOB_formatted <- Greenland_Halibut_Survey_Data %>%
  rename_with(tolower) 
#colnames(Greenland_Halibut_Survey_data_FISHGLOB_formatted)

Greenland_Halibut_Survey_data_FISHGLOB_formatted <- Greenland_Halibut_Survey_data_FISHGLOB_formatted %>%
  # 1. Rename: new_name = existing_name
  rename(
    haul_id = trawl_id,
    depth  = mean_depth,
    num = `total set count`,
    wgt = `total set wt(kg)`,
    aphia_id = aphiaid, 
    sbt = temp_ctd_degc
  ) %>%
  # 2. Mutate: Create new columns or modify renamed ones
  mutate(
    survey  = "DFO-GH", 
    survey_unit = "DFO-GH",
    source  = "DFO", 
    area_swept = `sweptarea(km2)`, 
    wgt_cpua = wgt / `sweptarea(km2)`, 
    num_cpua = num / `sweptarea(km2)`,
    # Taking the accepted name from the WoRMS column
    accepted_name = wormsnamecurr,
    # 1. Convert the text to a real Date object first
    date_formatted = mdy(`start date`),
    # 2. Extract the Year and Month
    year  = year(date_formatted),
    month = month(date_formatted) # This gives numbers (e.g., 7)
  )

Northern_Shrimp_Survey_data_FISHGLOB_formatted <- Northern_Shrimp_Survey_Data %>%
  rename_with(tolower) 
#colnames(Northern_Shrimp_Survey_data_FISHGLOB_formatted)

Northern_Shrimp_Survey_data_FISHGLOB_formatted <- Northern_Shrimp_Survey_data_FISHGLOB_formatted %>%
  # 1. Rename: new_name = existing_name
  rename(
    haul_id = trawl_id,
    depth  = mean_depth,
    num = `total set count`,
    wgt = `total set wt(kg)`,
    aphia_id = aphiaid, 
    sbt = temp_ctd_degc
  ) %>%
  # 2. Mutate: Create new columns or modify renamed ones
  mutate(
    survey  = "DFO-NS", 
    survey_unit = "DFO-NS",
    source  = "DFO", 
    area_swept = `sweptarea(km2)`, 
    wgt_cpua = wgt / `sweptarea(km2)`, 
    num_cpua = num / `sweptarea(km2)`,
    # Taking the accepted name from the WoRMS column
    accepted_name = wormsnamecurr, 
    # 1. Convert the text to a real Date object first
    date_formatted = mdy(`start date`),
    # 2. Extract the Year and Month
    year  = year(date_formatted),
    month = month(date_formatted) # This gives numbers (e.g., 7)
  )

#____________________________________#
##### Initialize German Greenlandic Survey #######
#____________________________________#

NetzFi <- read_csv("Thunen German Greenland Survey/NetzFi_Brian16.01.26.csv")
NetzFi <- NetzFi %>%
  mutate(Latitude = (GBENDE + GBSTART)/2)
NetzFi <- NetzFi %>%
  mutate(Longitude = (GLENDE + GLSTART)/2)
NetzFi <- NetzFi %>%
  mutate(Depth = (FTMIN + FTMAX)/2)
NetzFi <- NetzFi %>%
  unite("haul_ID", JAHR, REISENR, STATION, sep = "_", remove = FALSE)

print(colnames(NetzFi))

NetzFi <- NetzFi %>%
  rename(
    year = JAHR,
    ship_number = REISENR,
    'distance (nautical miles)' = DISTANZ,
    duration = SCHLDAUER,
    speed = SCHLGESCHW,
    sbt = TBODEN,
    sst = TOBERFL, 
    bottom_salinity = SBODEN, 
    surface_salinity = SOBERFL
  )

# Prepare Greenland Survey Data
ThunenGreendland_for_plotting <- NetzFi %>%
  mutate(Survey = "Thunen Greenland") %>% # Assign survey name
  select(
    Latitude,
    Longitude,
    Survey
  )

FishFi <- read_csv("Thunen German Greenland Survey/FishFi_Brian20.01.2026.csv")
FishFi <- FishFi %>%
  unite("haul_ID", JAHR, REISENR, STATION, sep = "_", remove = FALSE)
FishFi <- FishFi %>%
  rename(
    accepted_name = FARTNAME, 
    weight = GESAMTKG, 
    num = GESAMTSTCK
  )

LengFi <- read_csv("Thunen German Greenland Survey/LengFi_Brian06.03.2026.csv")
LengFi <- LengFi %>%
  unite("haul_ID", JAHR, REISENR, STATION, sep = "_", remove = FALSE)
LengFi <- LengFi %>%
  rename(
    length_class = LAENGE,
    class_num = LANZAHL
  )

# use haul_ID to match catch and trawl data

# Attempt splitting German Greenland survey 

# Load NAFO divisions for north WEST atlantic rfmo divisions
nafo_sf <- sf::st_read("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/DFO advice on SCS subregioning/NAFO divisions/Divisions.shp", quiet = TRUE)

# Load FAO divisions for north EAST atlantic rfmo divisions
FAO27_sf <- sf::st_read("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/DFO advice on SCS subregioning/FAO 27 subdivisions/ICES_Areas_20160601_cut_dense_3857.shp", quiet = TRUE)

# Load FAO37 subdivisions called GSAs for Med / Black Sea
FAO37_GSAs_sf <- sf::st_read("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/DFO advice on SCS subregioning/General Fisheries Commission for the Mediterranean subregions/gfcm_gsa.shp", quiet = TRUE)


library(sf)
library(dplyr)
sf_use_s2(FALSE)

# 1. Convert NetzFi to a spatial object (sf)
# We use CRS 4326 (WGS84) which is standard for Lat/Lon
NetzFi_sf <- NetzFi %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

# 2. Ensure shapefiles are in the same CRS (4326)
# Your FAO file mentions 3857 in the name, so we must transform them to match
nafo_sf <- st_transform(nafo_sf, 4326)
FAO27_sf <- st_transform(FAO27_sf, 4326)
FAO37_GSAs_sf <- st_transform(FAO37_GSAs_sf, 4326)

# 3. Spatial Join: This "tags" each haul with the NAFO 'ZONE' and FAO 'Area_27'
NetzFi_sf <- NetzFi_sf %>%
  st_join(nafo_sf %>% select(ZONE), join = st_intersects) %>%
  st_join(FAO27_sf %>% select(Area_27), join = st_intersects)

# 4. Create the Greenland_Region column based on your logic
NetzFi <- NetzFi_sf %>%
  mutate(Greenland_Region = case_when(
    # West Greenland Logic (NAFO 1B-1F)
    ZONE %in% c("1B", "1C", "1D", "1E", "1F") ~ "West Greenland",
    
    # East Greenland Logic (FAO 14.b.1 or 14.b.2)
    # Note: Check your shapefile attributes to see if they use "14.b.1" or "XIVb1"
    Area_27 %in% c("14.b.1", "14.b.2", "XIVb1", "XIVb2") ~ "East Greenland",
    
    TRUE ~ "Other"
  )) %>%
  # Convert back to a regular dataframe if you don't want to keep it as an sf object
  st_drop_geometry()

# Check the results
table(NetzFi$Greenland_Region)

# Join FishFi and NetzFi 
Thunen_combined_trawl_species <- left_join(FishFi, NetzFi, by = "haul_ID")


# Prepare Greenland Survey Data
ThunenGreendlandWest <- Thunen_combined_trawl_species %>%
  filter(Greenland_Region == "West Greenland") %>% # Filter for West Greenland
  mutate(Survey = "Thunen Greenland West")  # Assign survey name

ThunenGreendlandEast <- Thunen_combined_trawl_species %>%
  filter(Greenland_Region == "East Greenland") %>% # Filter for West Greenland
  mutate(Survey = "Thunen Greenland East") # Assign survey name


ThunenGreendlandEast_FISHGLOB_formatted <- ThunenGreendlandEast %>%
  rename_with(tolower) 
colnames(ThunenGreendlandEast_FISHGLOB_formatted)
ThunenGreendlandWest_FISHGLOB_formatted <- ThunenGreendlandWest %>%
  rename_with(tolower) 

library(worrms)
library(dplyr)
library(stringr)
library(purrr)

# 1. Clean and prepare unique names
unique_species <- c(
  ThunenGreendlandEast_FISHGLOB_formatted$accepted_name, 
  ThunenGreendlandWest_FISHGLOB_formatted$accepted_name
) %>%
  unique() %>%
  na.omit() 


build_worms_taxonomy <- function(species_names){
  
  species_names <- species_names %>%
    unique() %>%
    na.omit() %>%
    trimws()
  
  purrr::map_df(species_names, function(sp){
    
    tryCatch({
      
      cat("\nLOOKING UP:", sp, "\n")
      
      res <- wm_records_taxamatch(sp)[[1]]
      
      print(res)
      
      if(is.null(res) || nrow(res) == 0){
        
        cat("NO MATCH\n")
        
        return(
          tibble(
            accepted_name = sp,
            aphia_id = NA_character_,
            kingdom = NA_character_,
            phylum = NA_character_,
            class = NA_character_,
            order = NA_character_,
            family = NA_character_,
            genus = NA_character_,
            species = NA_character_,
            rank = NA_character_
          )
        )
        
      }
      
      match <- res[1, ]
      
      print(match)
      
      tibble(
        accepted_name = sp,
        aphia_id = as.character(match$AphiaID),
        kingdom = match$kingdom,
        phylum = match$phylum,
        class = match$class,
        order = match$order,
        family = match$family,
        genus = match$genus,
        species = match$scientificname,
        rank = match$rank
      )
      
    }, error = function(e){
      
      cat("ERROR:", e$message, "\n")
      
      tibble(
        accepted_name = sp,
        aphia_id = NA_character_,
        kingdom = NA_character_,
        phylum = NA_character_,
        class = NA_character_,
        order = NA_character_,
        family = NA_character_,
        genus = NA_character_,
        species = NA_character_,
        rank = NA_character_
      )
      
    })
    
  })
  
}

taxonomy_lookup <- build_worms_taxonomy(unique_species)

ThunenGreendlandEast_FISHGLOB_formatted <- ThunenGreendlandEast_FISHGLOB_formatted %>%
  left_join(taxonomy_lookup, by = "accepted_name") %>%
  rename(
    station = station.x,
    month = monat.x,     # Using . rather than monat.x if renamed earlier
    wgt = weight,
    gear = netz.x,
    quarter = quartal.x,
    haul_dur = `distance (nautical miles)`
  ) %>%
  mutate(
    survey  = "Thunen Greenland", 
    survey_unit = "East Greenland",
    source  = "Thunen", 
    num_cpue = num / haul_dur,
    wgt_cpue = wgt / haul_dur,
    accepted_name = species
  ) 

colnames(ThunenGreendlandWest)


ThunenGreendlandWest_FISHGLOB_formatted <- ThunenGreendlandWest_FISHGLOB_formatted %>%
  left_join(taxonomy_lookup, by = "accepted_name") %>%
# 1. Rename: new_name = existing_name
  rename(
    station = station.x,
    month = monat.x,
    wgt = weight,
    gear = netz.x,
    quarter = quartal.x,
    haul_dur = `distance (nautical miles)`
  ) %>%
  # 2. Mutate: Create new columns or modify renamed ones
  mutate(
    survey  = "Thunen Greenland", 
    survey_unit = "West Greenland",
    source  = "Thunen", 
    num_cpue = num / haul_dur,
    wgt_cpue = wgt / haul_dur, 
    accepted_name = species
  )

#____________________________________#
##### Initialize Icelandic Survey #######
#____________________________________#

Icelandic_groundfish_survey <- read_csv("Iceland/Icelandic groundfish survey_jan2026.csv")

Icelandic_groundfish_survey_FISHGLOB_formatted <- Icelandic_groundfish_survey %>%
  rename_with(tolower) 

Icelandic_groundfish_survey_FISHGLOB_formatted <- Icelandic_groundfish_survey_FISHGLOB_formatted %>%
  # 1. Rename: new_name = existing_name
  rename(
    haul_id = tow_id,
    station = station_nr,
    month = mon,
    latitude = lat, 
    longitude = lon,
    sbt = bottom_temp,
    sst = surface_temp,
    haul_dur = towlength,
    num = n,
    wgt = estimate.kg,
    wgt_cpua = estimate.kg.km2, 
    aphia_id = worms_id
  ) %>%
  # 2. Mutate: Create new columns or modify renamed ones
  mutate(
    survey  = "IS", 
    survey_unit = "IS",
    source  = "MFRI", 
    accepted_name = species,
    phylum = "Chordata", 
    num_cpue = num / haul_dur
  )

#____________________________________#
##### Initialize Norway Survey #######
#____________________________________#

load("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/Additional survey datasets/Norway Institute of Marine Research/Data_Kysttokt_article_Moving_North.RData")
Norway_catch <- Kys_species_clean
Norway_hauls <- Kys_stations_clean
Norway_combined <- left_join(Norway_catch, Norway_hauls, by = "year_serialno")
Norway_combined_survey <- Norway_combined %>%
  mutate(
    trawled_distance_km = as.numeric(trawled_distance_km),
    trawl_opening_horizontal_m = as.numeric(trawl_opening_horizontal_m), 
    Latitude = as.numeric(lat), 
    Longitude = as.numeric(long), 
    Survey = "IMR Norwegian Sea", 
    area_swept = trawled_distance_km * trawl_opening_horizontal_m / 1000
  ) %>%
  select(-c(lat, long))

Norway_combined_survey_FISHGLOB_formatted <- Norway_combined_survey %>%
  rename_with(tolower) %>%
  rename(
    haul_id = year_serialno,
    haul_dur = trawled_distance_km,
    num = count,
    wgt = weight_kg,
    accepted_name = species
  ) %>%
  mutate(
    # 1. Fix the Year: Extract first 4 chars from haul_id and convert to numeric
    year = as.numeric(substr(haul_id, 1, 4)),
    # If it's already numeric, as.numeric() does nothing. 
    # If it's a character, it converts it.
    wgt = as.numeric(wgt), 
    area_swept = as.numeric(area_swept),
    num = as.numeric(num),
    
    survey_unit = "IMR Norwegian Sea",
    source  = "IMR", 
    
    # Calculations
    wgt_cpua = wgt / area_swept, 
    num_cpua = num / area_swept, 
    
    country = "Norway",
    
    # Clean encoding and whitespace
    accepted_name = iconv(accepted_name, from = "latin1", to = "UTF-8", sub = ""),
    accepted_name = trimws(accepted_name)
  )


# 1. Clean and prepare unique names
unique_species <- c(Norway_combined_survey_FISHGLOB_formatted$accepted_name) %>%
  unique() %>%
  na.omit() 

build_worms_taxonomy <- function(species_names){
  
  species_names <- species_names %>%
    unique() %>%
    na.omit() %>%
    trimws()
  
  purrr::map_df(species_names, function(sp){
    
    tryCatch({
      
      cat("\nLOOKING UP:", sp, "\n")
      
      res <- wm_records_taxamatch(sp)[[1]]
      
      print(res)
      
      if(is.null(res) || nrow(res) == 0){
        
        cat("NO MATCH\n")
        
        return(
          tibble(
            accepted_name = sp,
            aphia_id = NA_character_,
            kingdom = NA_character_,
            phylum = NA_character_,
            class = NA_character_,
            order = NA_character_,
            family = NA_character_,
            genus = NA_character_,
            species = NA_character_,
            rank = NA_character_
          )
        )
        
      }
      
      match <- res[1, ]
      
      print(match)
      
      tibble(
        accepted_name = sp,
        aphia_id = as.character(match$AphiaID),
        kingdom = match$kingdom,
        phylum = match$phylum,
        class = match$class,
        order = match$order,
        family = match$family,
        genus = match$genus,
        species = match$scientificname,
        rank = match$rank
      )
      
    }, error = function(e){
      
      cat("ERROR:", e$message, "\n")
      
      tibble(
        accepted_name = sp,
        aphia_id = NA_character_,
        kingdom = NA_character_,
        phylum = NA_character_,
        class = NA_character_,
        order = NA_character_,
        family = NA_character_,
        genus = NA_character_,
        species = NA_character_,
        rank = NA_character_
      )
      
    })
    
  })
  
}

taxonomy_lookup <- build_worms_taxonomy(unique_species)

Norway_combined_survey_FISHGLOB_formatted <- left_join(Norway_combined_survey_FISHGLOB_formatted, taxonomy_lookup, by = "accepted_name")

#____________________________________#
##### Standardize datasets and check completeness #######
#____________________________________#

external_dataset_list <- c(Norway_combined_survey_FISHGLOB_formatted, Icelandic_groundfish_survey_FISHGLOB_formatted, NewfoundlandSurvey, Greenland_Halibut_Survey_data_FISHGLOB_formatted, Northern_Shrimp_Survey_data_FISHGLOB_formatted, ThunenGreendlandWest_FISHGLOB_formatted, ThunenGreendlandEast_FISHGLOB_formatted, NE_Bering_Sea, Bering_Slope, FaroePlSpr, FaroePlSum, FaroeBaSpr, FaroeBaSum, GINR_fish_shrimp_west, GINR_fish_shrimp_east)

# Harmonize all external datasets! 
library(tibble)
library(dplyr)

template_df <- tibble(
  # --- Character (Text) Columns ---
  survey            = character(),
  source            = character(),
  haul_id           = character(),
  country           = character(),
  sub_area          = character(),
  continent         = character(),
  stat_rec          = character(),
  station           = character(),
  stratum           = character(),
  season            = character(),
  gear              = character(),
  verbatim_name     = character(),
  accepted_name     = character(),
  kingdom           = character(),
  phylum            = character(),
  class             = character(),
  order             = character(),
  family            = character(),
  genus             = character(),
  rank              = character(),
  survey_unit       = character(),
  verbatim_aphia_id = character(), 
  aphia_id          = character(),
  SpecCode          = character(),
  
  # --- Numeric (Double) Columns ---
  year              = numeric(),
  month             = numeric(),
  day               = numeric(),
  quarter           = numeric(),
  latitude          = numeric(),
  longitude         = numeric(),
  haul_dur          = numeric(),
  area_swept        = numeric(),
  depth             = numeric(),
  sbt               = numeric(),
  sst               = numeric(),
  num               = numeric(),
  num_cpue          = numeric(),
  num_cpua          = numeric(),
  wgt               = numeric(),
  wgt_cpue          = numeric(),
  wgt_cpua          = numeric(),
  
  # --- Date/Time Columns ---
  timestamp         = as.POSIXct(character())
  
)

# 2. Your function stays mostly the same, 
# but we'll use 'as.character' on haul_id to be safe
standardize_to_template <- function(df, template) {
  # 1. Identify common columns
  common_cols <- intersect(names(df), names(template))
  
  # 2. Loop through common columns and match types to template
  for (col in common_cols) {
    template_type <- class(template[[col]])[1]
    
    if (template_type == "character") {
      df[[col]] <- as.character(df[[col]])
    } else if (template_type == "numeric" || template_type == "double") {
      df[[col]] <- as.numeric(df[[col]])
    } else if (template_type == "integer") {
      df[[col]] <- as.integer(df[[col]])
    } else if (template_type == "POSIXct" || template_type == "POSIXt") {
      # Try to preserve date-time if possible, otherwise keep as is
      if (!inherits(df[[col]], "POSIXct")) {
        df[[col]] <- as.POSIXct(df[[col]])
      }
    }
  }
  
  # 3. Finalize the structure
  df %>%
    select(any_of(names(template))) %>%
    bind_rows(template, .) %>%
    select(all_of(names(template)))
}
# Put your dataframes into a named list
external_survey_list <- list(
  GreenlandHalibut = Greenland_Halibut_Survey_data_FISHGLOB_formatted,
  NorthernShrimp = Northern_Shrimp_Survey_data_FISHGLOB_formatted,
  Newfoundland_spring = NewfoundlandSurvey_FISHGLOB_formatted_spring,
  Newfoundland_winter = NewfoundlandSurvey_FISHGLOB_formatted_winter,
  Iceland = Icelandic_groundfish_survey_FISHGLOB_formatted,
  GreenlandEast = ThunenGreendlandEast_FISHGLOB_formatted, 
  GreenlandWest = ThunenGreendlandWest_FISHGLOB_formatted, 
  NorwegianSea = Norway_combined_survey_FISHGLOB_formatted, 
  NEBeringSea = NE_Bering_Sea,
  Bering_Slope = Bering_Slope,
  Medits = MEDITS,
  Solemon = SOLEMON,
  FaroePlateauSpring = FaroePlSpr, 
  FaroePlateauSummer = FaroePlSum, 
  FaroeBankSpring = FaroeBaSpr, 
  FaroeBankSummer = FaroeBaSum, 
  GINR_west = GINR_fish_shrimp_west, 
  GINR_east = GINR_fish_shrimp_east
)

# Apply the "Auto-Typer" to your list
standardized_list <- purrr::map(external_survey_list, ~standardize_to_template(.x, template_df))

# Combine everything
External_surveys_FISHGLOB_compatible <- bind_rows(standardized_list)

library(dplyr)
library(tidyr)
library(purrr)

analyze_seasonal_drift <- function(data) {
  
  # 1. Summarize monthly effort distribution per survey unit over time
  effort_summary <- data %>%
    group_by(survey_unit, year, month) %>%
    summarise(num_hauls = n(), .groups = "drop")
  
  # 2. Extract unique survey units to evaluate
  units <- unique(effort_summary$survey_unit)
  
  report <- map_df(units, function(unit) {
    unit_data <- effort_summary %>% filter(survey_unit == unit)
    
    # Calculate global profile for this unit
    total_months_active <- unique(unit_data$month)
    num_unique_months <- length(total_months_active)
    
    # FIX: Changed fill_defaults to values_fill
    matrix_view <- unit_data %>%
      pivot_wider(id_cols = year, names_from = month, values_from = num_hauls, values_fill = 0)
    
    # Check for multi-quarter sampling footprints (potential distinct seasons)
    quarters_active <- data %>% 
      filter(survey_unit == unit) %>% 
      pull(quarter) %>% 
      unique()
    
    # Metric A: Check for wide within-year seasonal spreads (e.g., Spring AND Fall)
    max_month_span <- if(num_unique_months > 1) max(total_months_active) - min(total_months_active) else 0
    
    # Metric B: Test if month distribution changes over time (correlation between year and month)
    correlation_test <- tryCatch({
      raw_vectors <- data %>% filter(survey_unit == unit)
      cor(raw_vectors$year, raw_vectors$month, method = "spearman", use = "complete.obs")
    }, error = function(e) 0)
    
    # Decision Matrix Logic
    action_required <- "Keep Combined"
    reason = "Consistent seasonal window."
    
    if (length(quarters_active) > 1 && max_month_span >= 4) {
      action_required <- "DIVIDE"
      reason <- paste("Multi-seasonal profile detected across quarters:", paste(sort(quarters_active), collapse=", "))
    } else if (!is.na(correlation_test) && abs(correlation_test) > 0.3) {
      action_required <- "DIVIDE"
      reason <- paste0("Temporal drift detected (r = ", round(correlation_test, 2), "): months sampled have shifted across years.")
    } else if (num_unique_months > 3) {
      action_required <- "REVIEW"
      reason <- "Wide monthly footprint within years; inspect for macro-seasonal contamination."
    }
    
    tibble(
      survey_unit = unit,
      unique_months_total = num_unique_months,
      months_list = paste(sort(total_months_active), collapse = ","),
      quarters_active = paste(sort(quarters_active), collapse = ","),
      temporal_drift_index = round(correlation_test, 3),
      recommendation = action_required,
      justification = reason
    )
  })
  
  return(report)
}

# Run the analyzer on your dataset
seasonal_report <- analyze_seasonal_drift(External_surveys_FISHGLOB_compatible)

# View units flagged for separation
print(seasonal_report %>% filter(recommendation %in% c("DIVIDE", "REVIEW")))

library(ggplot2)
library(dplyr)
library(grid)


External_surveys_FISHGLOB_compatible %>%
  count(survey_unit, year, month) %>%
  ggplot(aes(year, factor(month), fill = n)) +
  geom_tile() +
  facet_wrap(~survey_unit, scales = "free_x")

# 1. Isolate the specific survey units flagged for DIVIDE or REVIEW
flagged_units <- c("DFO-GH", "DFO-NL", "DFO-NS", "East Greenland", "IMR Norwegian Sea")

audit_data <- External_surveys_FISHGLOB_compatible %>% 
  filter(survey_unit %in% flagged_units)

'
# 2. Build a faceted heatmap to inspect the seasonal profiles simultaneously
ggplot(audit_data, aes(x = factor(month), y = factor(year))) +
  geom_bin2d(aes(fill = ..count..)) +
  scale_fill_viridis_c(name = "Haul Count", option = "plasma") +
  facet_wrap(~survey_unit, scales = "free_y", ncol = 2) +
  theme_minimal() +
  labs(
    title = "Seasonal Drift & Footprint Audit",
    subtitle = "Vertical shifts = temporal drift over years | Wide horizontal blocks = seasonal spread",
    x = "Month of Year",
    y = "Year"
  ) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 6), # Shrink year labels to prevent overlapping
    panel.background = element_rect(fill = "gray98", color = NA)
  )
'


'
###### month fix #############################
na_diagnosis_report <- External_surveys_FISHGLOB_compatible %>%
  # Create a flag that catches both true logical NAs and string "NA" text artifacts
  mutate(is_month_na = is.na(month) | month == "NA" | month == "") %>%
  group_by(survey, survey_unit) %>%
  summarise(
    total_records = n(),
    na_month_count = sum(is_month_na),
    pct_missing    = round((na_month_count / total_records) * 100, 2),
    .groups = "drop"
  ) %>%
  # Filter to show only the survey units that actually contain missing months
  filter(na_month_count > 0) %>%
  arrange(desc(na_month_count))

# 2. Print the final diagnosis report to the console
print(na_diagnosis_report)
library(dplyr)
library(tidyr)

# 1. Standardize literal "NA" strings to true logical NAs
External_surveys_standardized <- External_surveys_FISHGLOB_compatible %>%
  mutate(
    month = case_when(
      month == "NA" ~ NA_real_,
      month == ""   ~ NA_real_,
      TRUE          ~ as.numeric(month)
    )
  )

# 2. Build a composite context dictionary
# We group by haul_id AND the structural metadata to prevent cross-survey collisions
haul_month_dictionary <- External_surveys_standardized %>%
  filter(!is.na(month)) %>%
  distinct(survey, survey_unit, year, haul_id, valid_month = month)

# 3. Join back using the multi-key composite match
# This guarantees a clean 1-to-many join, completely eliminating the warning
External_surveys_patched <- External_surveys_standardized %>%
  left_join(
    haul_month_dictionary, 
    by = c("survey", "survey_unit", "year", "haul_id")
  ) %>%
  mutate(
    # Patch the missing months
    month = coalesce(month, valid_month)
  ) %>%
  select(-valid_month)

# Final validation audit on the patched dataset
patched_na_audit <- External_surveys_patched %>%
  mutate(is_month_na = is.na(month)) %>%
  group_by(survey, survey_unit) %>%
  summarise(
    total_records = n(),
    remaining_nas = sum(is_month_na),
    .groups = "drop"
  ) %>%
  filter(survey_unit %in% c("IMR Norwegian Sea", "DFO-NL"))

print(patched_na_audit)

# final conclusion: missing months of newfoundland and imr cannot be salvaged using other catches in same hauls

'

# Calculate the fraction of finite (non-NA) values
survey_unit_completeness_report <- External_surveys_FISHGLOB_compatible %>%
  group_by(survey_unit) %>%
  summarise(across(everything(), ~ mean(!is.na(.x))))

completeness_long <- survey_unit_completeness_report %>%
  pivot_longer(cols = -survey_unit, names_to = "column_name", values_to = "fraction_complete") %>%
  filter(fraction_complete <= 1.0) %>%  # Show only columns that have at least some missing data
  arrange(survey_unit, fraction_complete)

# This gives you a targeted list of where data is missing
print(completeness_long)
incomplete_data_fields <- completeness_long %>%
  filter(fraction_complete < 0.5)


worst_offenders <- incomplete_data_fields %>%
  group_by(column_name) %>%
  summarise(num_incomplete = n()) %>%
  arrange(num_incomplete)

completeness_of_datafields <- completeness_long %>%
  group_by(column_name) %>%
  summarise(frac_complete = mean(fraction_complete)) %>%
  arrange(frac_complete)


FishGlob_hexxed <- FishGlob_hexxed [, 1:42]

Combined_global_SBTS_data <- rbind(FishGlob_hexxed, External_surveys_FISHGLOB_compatible)
Combined_FISH_only <- Combined_global_SBTS_data %>%
  filter(phylum == 'Chordata', 
         !is.na(latitude), 
         !is.na(longitude)) %>%
  mutate(year = as.integer(year), 
         month = as.integer(month), 
         day = as.integer(day), 
  )

# Fix for Nor-BTS and other FishGlob units where AphiaID might be missing or character
Combined_FISH_only <- Combined_FISH_only %>%
  mutate(
    # Ensure aphia_id is character for matching
    aphia_id = as.character(aphia_id),
    # If aphia_id is NA, try to fill it using the accepted_name and WoRMS
    # (Note: You already have a map_df/wm_records_taxamatch logic for external surveys,
    # apply a similar bridge for the FishGlob units if aphia_id is missing)
  )

#_________________________________________________#
##### Implement survey unit name changes ######
#_________________________________________________#

Naming_convention_survey_units_Chap_2 <- read_csv("C:/Users/besch/Desktop/FishGlob project/Naming convention survey units Chap 2.csv")
library(gt)

Naming_table_gt <- Naming_convention_survey_units_Chap_2 %>%
  gt() %>%
  tab_header(
    title = "Survey Unit Naming Convention"
  )

gtsave(
  Naming_table_gt,
  "Survey_Unit_Naming_Convention.html"
)

old_names <- unique(Naming_convention_survey_units_Chap_2$`Previous Name`)
current_names <- unique(Combined_FISH_only$survey_unit)
matches <- current_names[current_names %in% old_names]
fix <- current_names[!current_names %in% old_names]

# Create named lookup vector: names = old names, values = new names
survey_lookup <- Naming_convention_survey_units_Chap_2 %>%
  select(`Previous Name`, `Name`) %>%
  tibble::deframe()

# Replace survey_unit values with new names where a match exists
Combined_FISH_only <- Combined_FISH_only %>%
  mutate(
    survey_unit = ifelse(
      survey_unit %in% names(survey_lookup),
      survey_lookup[survey_unit],
      survey_unit
    )
  )

# Check results
Combined_FISH_only %>%
  count(survey_unit, sort = TRUE)

#_________________________________________________#
##### Identify temporal crossover of surveys ######
#_________________________________________________#

Combined_FISH_only$survey_unit <- as.factor(Combined_FISH_only$survey_unit)
Combined_FISH_only$year <- as.numeric(Combined_FISH_only$year)
# 1. Summarize the data to get the range for each unit
plot_data <- Combined_FISH_only %>%
  group_by(survey_unit) %>%
  summarize(
    start_year = min(year),
    end_year = max(year)
  ) %>%
  # Optional: Sort by start year so the chart looks organized
  arrange(start_year) %>%
  mutate(survey_unit = factor(survey_unit, levels = survey_unit))

# 2. Create the visualization
timeline_plot <- ggplot(plot_data, aes(x = start_year, xend = end_year, y = survey_unit, yend = survey_unit)) +
  # Draw the horizontal lines (the "duration")
  geom_segment(size = 3, color = "skyblue", lineend = "round") +
  # Add points at the start and end for clarity
  geom_point(aes(x = start_year), size = 3) +
  geom_point(aes(x = end_year), size = 3) +
  # Formatting
  labs(
    x = "Year",
    y = "Survey Unit"
  )  +
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

ggsave(
  plot = timeline_plot,
  filename = "Survey_Unit_Timespans.png",
  width = 9,
  height = 12,
  dpi = 200
)

#_________________________________________________#
##### Plot haul locations across all surveys ######
#_________________________________________________#

library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)

global_haul_plot_data <- Combined_FISH_only %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  group_by(survey_unit, year, haul_id) %>%
  summarize(
    Latitude = first(latitude), 
    Longitude = first(longitude), 
    .groups = "drop"
  )

global_haul_plot_data <- global_haul_plot_data %>%
  mutate(
    Longitude_plot = ifelse(Longitude > 60, Longitude - 360, Longitude)
  )
library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
library(ggrepel)

# Generate survey unit IDs
survey_units <- sort(unique(global_haul_plot_data$survey_unit))

survey_lookup <- data.frame(
  survey_unit = survey_units,
  survey_id = seq_along(survey_units)
)

# Add IDs to plotting data
global_haul_plot_data <- global_haul_plot_data %>%
  left_join(survey_lookup, by = "survey_unit")

# Create legend labels with ID numbers
survey_labels <- setNames(
  paste0(
    survey_lookup$survey_id,
    ": ",
    survey_lookup$survey_unit
  ),
  survey_lookup$survey_unit
)

# More saturated colors
base_colors <- c(
  brewer.pal(8, "Dark2"),
  brewer.pal(9, "Set1"),
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Accent")
)

survey_colors <- setNames(
  colorRampPalette(base_colors)(length(survey_units)),
  survey_units
)

# Create centroid locations for survey labels
survey_centroids <- global_haul_plot_data %>%
  group_by(survey_unit, survey_id) %>%
  summarise(
    Longitude = median(Longitude),
    Latitude = median(Latitude),
    .groups = "drop"
  )

# World map
world_map_wrapped <- ne_countries(
  scale = "medium",
  returnclass = "sf"
) %>%
  st_make_valid() %>%
  st_wrap_dateline(
    options = c(
      "WRAPDATELINE=YES",
      "DATELINEOFFSET=180"
    )
  )

# Plot limits
lon_lims <- c(-195, 50)

lat_lims <- c(
  min(global_haul_plot_data$Latitude, na.rm = TRUE) - 2,
  max(global_haul_plot_data$Latitude, na.rm = TRUE) + 2
)
highlight_countries <- c(
  "Canada",
  "Germany",
  "Iceland",
  "Norway",
  "United States of America",
  "Greenland",
  "France",
  "Italy",
  "Spain",
  "Malta",
  "Greece",
  "Slovenia",
  "Cyprus",
  "Croatia",
  "Faeroe Is.",
  "Ireland",
  "United Kingdom",
  "Portugal",
  "Denmark",
  "Netherlands",
  "Belgium",
  "Sweden",
  "Poland",
  "Lithuania",
  "Latvia",
  "Estonia"
)
library(grid)
world_map_plot <- world_map_wrapped %>%
  filter(
    continent %in% c(
      "North America",
      "Europe",
      "Africa",
      "Asia"
    ) |
      name %in% c(
        "Greenland",
        "Iceland"
      )
  ) %>%
  mutate(
    country_fill = ifelse(
      name %in% highlight_countries,
      "Survey Country",
      "Other Country"
    )
  )
library(grid)
global_haul_map <- ggplot() +
  
  # Land masses (keep countries blue, remove legend)
  geom_sf(
    data = world_map_plot,
    aes(fill = country_fill),
    color = "black",
    linewidth = 0.2,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(
    values = c(
      "Survey Country" = "lightblue",
      "Other Country" = "grey85"
    )
  ) +
  
  # Haul locations
  geom_point(
    data = global_haul_plot_data,
    aes(
      x = Longitude,
      y = Latitude,
      color = survey_unit
    ),
    alpha = 0.8,
    size = 1.0
  ) +
  
  # Larger survey ID labels with overlap avoidance
  geom_text_repel(
    data = survey_centroids,
    aes(
      x = Longitude,
      y = Latitude,
      label = survey_id
    ),
    size = 5,
    fontface = "bold",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey40",
    segment.size = 0.4,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = survey_colors,
    labels = survey_labels,
    name = "Survey Unit"
  ) +
  
  coord_sf(
    xlim = lon_lims,
    ylim = lat_lims,
    crs = 4326,
    expand = FALSE,
    label_graticule = "SW"
  ) +
  
  scale_x_continuous(
    breaks = seq(-190, 50, by = 15)
  ) +
  
  scale_y_continuous(
    breaks = seq(
      floor(min(lat_lims)),
      ceiling(max(lat_lims)),
      by = 5
    )
  ) +
  
  labs(
    title = "All Unique Haul Locations by Survey Unit",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  guides(
    color = guide_legend(
      ncol = 8,
      byrow = TRUE,
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold"
    ),
    
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    
    panel.grid.major = element_line(
      color = "grey70",
      linewidth = 0.4
    ),
    
    panel.grid.minor = element_line(
      color = "grey90",
      linewidth = 0.25
    ),
    
    legend.position = "bottom",
    legend.box = "vertical",
    
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    
    legend.text = element_text(
      size = 10
    )
  )

print(global_haul_map)
ggsave(
  "combined_global_haul_map.png",
  plot = global_haul_map,
  width = 24,
  height = 12,
  dpi = 300
)

'
this version plots survey units by westerliness instead of alphabetical order
#_________________________________________________#
##### Plot haul locations across all surveys ######
#_________________________________________________#

library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
library(ggrepel)
library(lwgeom)
library(grid)

# Create unique haul-location dataset
global_haul_plot_data <- Combined_FISH_only %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  group_by(survey_unit, year, haul_id) %>%
  summarize(
    Latitude = first(latitude),
    Longitude = first(longitude),
    .groups = "drop"
  )

# Shift longitudes so Pacific surveys plot correctly on the left
global_haul_plot_data <- global_haul_plot_data %>%
  mutate(
    Longitude_plot = ifelse(Longitude > 60, Longitude - 360, Longitude)
  )

#--------------------------------------------------#
# Assign survey IDs from west --> east
#--------------------------------------------------#

survey_lookup <- global_haul_plot_data %>%
  group_by(survey_unit) %>%
  summarise(
    representative_longitude = median(Longitude_plot, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(representative_longitude) %>%
  mutate(
    survey_id = row_number()
  ) %>%
  select(survey_unit, survey_id)

survey_units <- survey_lookup$survey_unit

# Add IDs to plotting data
global_haul_plot_data <- global_haul_plot_data %>%
  left_join(survey_lookup, by = "survey_unit")

# Force west-to-east ordering everywhere
global_haul_plot_data$survey_unit <- factor(
  global_haul_plot_data$survey_unit,
  levels = survey_lookup$survey_unit
)

# Legend labels
survey_labels <- setNames(
  paste0(
    survey_lookup$survey_id,
    ": ",
    survey_lookup$survey_unit
  ),
  survey_lookup$survey_unit
)

#--------------------------------------------------#
# Colors
#--------------------------------------------------#

base_colors <- c(
  brewer.pal(8, "Dark2"),
  brewer.pal(9, "Set1"),
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Accent")
)

survey_colors <- setNames(
  colorRampPalette(base_colors)(length(survey_units)),
  survey_units
)

#--------------------------------------------------#
# Survey label locations
#--------------------------------------------------#

survey_centroids <- global_haul_plot_data %>%
  group_by(survey_unit, survey_id) %>%
  summarise(
    Longitude = median(Longitude_plot, na.rm = TRUE),
    Latitude = median(Latitude, na.rm = TRUE),
    .groups = "drop"
  )

#--------------------------------------------------#
# World map
#--------------------------------------------------#

world_map_wrapped <- ne_countries(
  scale = "medium",
  returnclass = "sf"
) %>%
  st_make_valid() %>%
  st_wrap_dateline(
    options = c(
      "WRAPDATELINE=YES",
      "DATELINEOFFSET=180"
    )
  )

#--------------------------------------------------#
# Map extent
#--------------------------------------------------#

lon_lims <- c(-195, 50)

lat_lims <- c(
  min(global_haul_plot_data$Latitude, na.rm = TRUE) - 2,
  max(global_haul_plot_data$Latitude, na.rm = TRUE) + 2
)

highlight_countries <- c(
  "Canada",
  "Germany",
  "Iceland",
  "Norway",
  "United States of America",
  "Greenland",
  "France",
  "Italy",
  "Spain",
  "Malta",
  "Greece",
  "Slovenia",
  "Cyprus",
  "Croatia",
  "Faeroe Is.",
  "Ireland",
  "United Kingdom",
  "Portugal",
  "Denmark",
  "Netherlands",
  "Belgium",
  "Sweden",
  "Poland",
  "Lithuania",
  "Latvia",
  "Estonia"
)

world_map_plot <- world_map_wrapped %>%
  filter(
    continent %in% c(
      "North America",
      "Europe",
      "Africa",
      "Asia"
    ) |
      name %in% c(
        "Greenland",
        "Iceland"
      )
  ) %>%
  mutate(
    country_fill = ifelse(
      name %in% highlight_countries,
      "Survey Country",
      "Other Country"
    )
  )

#--------------------------------------------------#
# Plot
#--------------------------------------------------#

global_haul_map <- ggplot() +
  
  geom_sf(
    data = world_map_plot,
    aes(fill = country_fill),
    color = "black",
    linewidth = 0.2,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(
    values = c(
      "Survey Country" = "lightblue",
      "Other Country" = "grey85"
    )
  ) +
  
  geom_point(
    data = global_haul_plot_data,
    aes(
      x = Longitude_plot,
      y = Latitude,
      color = survey_unit
    ),
    alpha = 0.8,
    size = 1.0
  ) +
  
  geom_text_repel(
    data = survey_centroids,
    aes(
      x = Longitude,
      y = Latitude,
      label = survey_id
    ),
    size = 5,
    fontface = "bold",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey40",
    segment.size = 0.4,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = survey_colors,
    labels = survey_labels,
    name = "Survey Unit",
    breaks = survey_lookup$survey_unit
  ) +
  
  coord_sf(
    xlim = lon_lims,
    ylim = lat_lims,
    crs = 4326,
    expand = FALSE,
    label_graticule = "SW"
  ) +
  
  scale_x_continuous(
    breaks = seq(-190, 50, by = 10)
  ) +
  
  scale_y_continuous(
    breaks = seq(
      floor(min(lat_lims)),
      ceiling(max(lat_lims)),
      by = 5
    )
  ) +
  
  labs(
    title = "All Unique Haul Locations by Survey Unit",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  guides(
    color = guide_legend(
      ncol = 4,
      byrow = TRUE,
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold"
    ),
    
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    
    panel.grid.major = element_line(
      color = "grey70",
      linewidth = 0.4
    ),
    
    panel.grid.minor = element_line(
      color = "grey90",
      linewidth = 0.25
    ),
    
    legend.position = "bottom",
    legend.box = "vertical",
    
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.6, "cm")
  )

print(global_haul_map)

ggsave(
  "combined_global_haul_map.png",
  plot = global_haul_map,
  width = 16,
  height = 10,
  dpi = 300
)
'

#_________________________________________________#
##### Plot selected external survey units ########
#_________________________________________________#

library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(ggrepel)
library(RColorBrewer)
library(lwgeom)
library(grid)

# Survey units to include

selected_surveys <- c(
  "E Bering Slope",
  "NE Bering Sea",
  "Newfoundland Labrador Q1",
  "Newfoundland Labrador Q2",
  "Baffin - Davis",
  "Hudson Strait",
  "W Greenland GINR",
  "E Greenland GINR",
  "E Greenland Thunen",
  "W Greenland Thunen",
  "Iceland",
  "Norwegian Sea",
  "Faroe Plateau Q1",
  "Faroe Plateau Q3",
  "Faroe Bank Q1",
  "Faroe Bank Q3",
  "Mediterranean Sea",
  "Adriatic Sea"
)

# Prepare haul locations

plot_data <- Combined_FISH_only %>%
  filter(
    survey_unit %in% selected_surveys,
    !is.na(latitude),
    !is.na(longitude)
  ) %>%
  group_by(survey_unit, haul_id, year) %>%
  summarise(
    Latitude = first(latitude),
    Longitude = first(longitude),
    .groups = "drop"
  ) %>%
  mutate(
    Longitude_plot = ifelse(Longitude > 60, Longitude - 360, Longitude)
  )

# Assign survey IDs west -> east

survey_lookup <- plot_data %>%
  group_by(survey_unit) %>%
  summarise(
    representative_longitude =
      median(Longitude_plot, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(representative_longitude) %>%
  mutate(
    survey_id = row_number()
  )

plot_data <- plot_data %>%
  left_join(
    survey_lookup %>%
      select(survey_unit, survey_id),
    by = "survey_unit"
  )

plot_data$survey_unit <- factor(
  plot_data$survey_unit,
  levels = survey_lookup$survey_unit
)

# Survey labels and colors

survey_labels <- setNames(
  paste0(
    survey_lookup$survey_id,
    ": ",
    survey_lookup$survey_unit
  ),
  survey_lookup$survey_unit
)

base_colors <- c(
  brewer.pal(8, "Dark2"),
  brewer.pal(9, "Set1"),
  brewer.pal(12, "Paired")
)

survey_colors <- setNames(
  colorRampPalette(base_colors)(
    length(survey_lookup$survey_unit)
  ),
  survey_lookup$survey_unit
)

# Label locations

survey_centroids <- plot_data %>%
  group_by(
    survey_unit,
    survey_id
  ) %>%
  summarise(
    Longitude = median(
      Longitude_plot,
      na.rm = TRUE
    ),
    Latitude = median(
      Latitude,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# World map

world_map_wrapped <- ne_countries(
  scale = "medium",
  returnclass = "sf"
) %>%
  st_make_valid() %>%
  st_wrap_dateline(
    options = c(
      "WRAPDATELINE=YES",
      "DATELINEOFFSET=180"
    )
  )

highlight_countries <- c(
  "Canada",
  "Germany",
  "Iceland",
  "Norway",
  "United States of America",
  "Greenland",
  "France",
  "Italy",
  "Spain",
  "Malta",
  "Greece",
  "Slovenia",
  "Cyprus",
  "Croatia",
  "Faeroe Is."
)

world_map_plot <- world_map_wrapped %>%
  filter(
    continent %in% c(
      "North America",
      "Europe",
      "Africa",
      "Asia"
    ) |
      name %in% c(
        "Greenland",
        "Iceland"
      )
  ) %>%
  mutate(
    country_fill = ifelse(
      name %in% highlight_countries,
      "Survey Country",
      "Other Country"
    )
  )

# Plot limits

lon_lims <- c(-195, 40)

lat_lims <- c(
  min(plot_data$Latitude, na.rm = TRUE) - 2,
  max(plot_data$Latitude, na.rm = TRUE) + 2
)

# Plot

external_survey_map <- ggplot() +
  
  geom_sf(
    data = world_map_plot,
    aes(fill = country_fill),
    color = "black",
    linewidth = 0.2,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(
    values = c(
      "Survey Country" = "lightblue",
      "Other Country" = "grey85"
    )
  ) +
  
  geom_point(
    data = plot_data,
    aes(
      x = Longitude_plot,
      y = Latitude,
      color = survey_unit
    ),
    alpha = 0.75,
    size = 1.2
  ) +
  
  geom_text_repel(
    data = survey_centroids,
    aes(
      x = Longitude,
      y = Latitude,
      label = survey_id
    ),
    size = 5,
    fontface = "bold",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey40",
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = survey_colors,
    labels = survey_labels,
    breaks = survey_lookup$survey_unit,
    name = "Survey Unit"
  ) +
  
  coord_sf(
    xlim = lon_lims,
    ylim = lat_lims,
    crs = 4326,
    expand = FALSE
  ) +
  
  scale_x_continuous(
    breaks = seq(-190, 40, by = 10)
  ) +
  
  labs(
    title = "Selected External Survey Units",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  guides(
    color = guide_legend(
      ncol = 4,
      byrow = TRUE,
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold"
    ),
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      size = 10
    )
  )

print(external_survey_map)

ggsave(
  "selected_external_surveys_map.png",
  plot = external_survey_map,
  width = 16,
  height = 10,
  dpi = 300
)

#_________________________________________________#
##### optimize spatial flagging at scale of individual survey units ######
#_________________________________________________#

library(dggridR)
library(dplyr)
library(purrr)
library(tidyr)
library(CFtime)
optimize_survey_hexes <- function(data, 
                                  consistency_threshold = 0.7, 
                                  retention_goal = 0.95) {
  
  dggs_temp <- dggridR::dgconstruct(res = 12, metric = TRUE)
  res_table <- dggridR::dggetres(dggs_temp) %>%
    filter(res >= 6 & res <= 13) 
  
  survey_list <- split(data, data$survey_unit)
  
  # Change 1: purrr::map
  results <- purrr::map(survey_list, function(sub_data) {
    
    sub_data <- sub_data %>%
      group_by(haul_id) %>%
      mutate(
        latitude = mean(latitude, na.rm = TRUE),
        longitude = mean(longitude, na.rm = TRUE)
      ) %>%
      ungroup()
    
    total_hauls_initial <- nrow(sub_data)
    total_years <- n_distinct(sub_data$year)
    
    # Change 2: purrr::map_df
    res_testing <- purrr::map_df(res_table$res, function(r) {
      dggs_res <- dggridR::dgconstruct(res = r, metric = TRUE)
      actual_area <- res_table$area_km[res_table$res == r]
      
      temp <- sub_data %>%
        mutate(hex_id = dggridR::dgGEO_to_SEQNUM(dggs_res, longitude, latitude)$seqnum) %>%
        group_by(hex_id) %>%
        mutate(years_sampled = n_distinct(year),
               is_valid = (years_sampled / total_years) >= consistency_threshold) %>%
        ungroup()
      
      data.frame(
        res = r,
        area_km2 = actual_area,
        rows_kept = sum(temp$is_valid),
        hexes_kept = n_distinct(temp$hex_id[temp$is_valid]),
        total_hexes = n_distinct(temp$hex_id)
      )
    })
    
    max_rows_achievable <- max(res_testing$rows_kept)
    
    
    optimal_res_row <- res_testing %>%
      filter(rows_kept >= (max_rows_achievable * retention_goal)) %>%
      filter(area_km2 == min(area_km2)) %>%
      dplyr::slice(1)
    
    if(nrow(optimal_res_row) == 0) {
      optimal_res_row <- res_testing %>%
        filter(area_km2 == max(area_km2)) %>%
        dplyr::slice(1)
    }
    
    
    final_dggs <- dggridR::dgconstruct(res = optimal_res_row$res, metric = TRUE)
    
    final_data <- sub_data %>%
      mutate(
        hexagon_id = dggridR::dgGEO_to_SEQNUM(final_dggs, longitude, latitude)$seqnum,
        hex_area_km2 = optimal_res_row$area_km2
      ) %>%
      group_by(hexagon_id) %>%
      mutate(
        years_sampled = n_distinct(year),
        spatial_valid = (years_sampled / total_years) >= consistency_threshold
      ) %>%
      ungroup()
    
    summary_tab <- data.frame(
      survey_unit = unique(sub_data$survey_unit),
      opt_res = optimal_res_row$res,
      hex_size_km2 = round(optimal_res_row$area_km2, 2),
      n_unique_hexes = optimal_res_row$hexes_kept,
      prop_rows_kept = round(optimal_res_row$rows_kept / total_hauls_initial, 3),
      prop_hexes_kept = round(optimal_res_row$hexes_kept / optimal_res_row$total_hexes, 3),
      total_area_before = optimal_res_row$total_hexes * optimal_res_row$area_km2,
      total_area_kept = optimal_res_row$hexes_kept * optimal_res_row$area_km2
    )
    
    return(list(data = final_data, summary = summary_tab))
  })
  
  # Change 3 & 4: purrr::map_df
  return(list(
    full_data = purrr::map_df(results, ~ .x$data),
    summary_table = purrr::map_df(results, ~ .x$summary)
  ))
}

# --- Execute the function ---
optimized_output <- optimize_survey_hexes(Combined_FISH_only, consistency_threshold = 0.70)

# The updated data with hex assignments
Combined_FISH_optimized <- optimized_output$full_data

# The summary table you requested
optimization_summary <- optimized_output$summary_table
library(gt)
library(dplyr)
meps_theme <- function(data, title = NULL, subtitle = NULL) {
  data %>%
    gt() %>%
    # Typography
    tab_options(
      table.font.names = "Times New Roman", # Or Arial, depending on final layout
      table.font.size = px(12),
      column_labels.font.weight = "bold",
      table.align = "left"
    ) %>%
    # Borders (The MEPS "Sandwich" Look)
    tab_options(
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
    cols_align(align = "left", columns = 1) %>% # Keep first column left-aligned
    tab_style(
      style = cell_text(indent = pct(2)),
      locations = cells_body(columns = 1)
    )
}

optimization_summary <- optimization_summary %>%
  arrange(hex_size_km2)%>%
  meps_theme()
print(optimization_summary)

library(dplyr)
library(gt)
library(scales)

optimization_summary_gt <-
  optimized_output$summary_table %>%
  arrange(hex_size_km2) %>%
  transmute(
    `Survey Unit` = survey_unit,
    
    `Optimal Resolution` = opt_res,
    
    `Selected Hexagon Area (km²)` =
      comma(round(hex_size_km2, 0)),
    
    `Number of Hexagons` =
      comma(n_unique_hexes),
    
    `Hauls Retained (%)` =
      percent(prop_rows_kept, accuracy = 0.1),
    
    `Hexagons Retained (%)` =
      percent(prop_hexes_kept, accuracy = 0.1),
    
    `Survey Area Retained (km²)` =
      comma(round(total_area_kept, 0)),
    
    `Area Retained (%)` =
      percent(
        total_area_kept / total_area_before,
        accuracy = 0.1
      )
  ) %>%
  gt() %>%
  tab_header(
    title = md("**Survey-Specific Spatial Optimization Summary**"),
    subtitle =
      "Optimal DGGS resolution selected to maximize spatial consistency while retaining ≥95% of achievable haul observations"
  ) %>%
  cols_align(
    align = "center",
    columns = -`Survey Unit`
  ) %>%
  cols_width(
    `Survey Unit` ~ px(220)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  opt_row_striping() %>%
  fmt_missing(
    columns = everything(),
    missing_text = "—"
  ) %>%
  tab_source_note(
    source_note =
      "Valid hexagons were defined as cells sampled during ≥70% of available survey years. Area retained (%) represents the fraction of the original survey footprint remaining after spatial consistency filtering."
  )

print(optimization_summary_gt)

gtsave(
  optimization_summary_gt,
  "Survey_Spatial_Optimization_Summary.html"
)

# note high variation in proportion of hexes kept but high proportion of catches kept; super high variation in optimal hex size; gonnna have to choose a different threshold, hex strategy or fix these across all survey units and report cuts
# fishglob used hexes of size 7 = 23,320km^2 or size 8 = 7,770km^2 which do actually make sense but for two groups of survey units

#_________________________________________________#
##### Generate spatially valid dataset with identified hexes ######
#_________________________________________________#

library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
library(ggrepel)
library(lwgeom)
library(grid)
library(dggridR)

Combined_FISH_optimized_spatially_valid <- Combined_FISH_optimized %>%
  filter(spatial_valid == TRUE) %>%
  select(-spatial_valid)

Combined_FISH_optimized_spatially_valid <- Combined_FISH_optimized_spatially_valid %>%
  mutate(
    wgt = as.numeric(wgt),
    num = as.numeric(num),
    wgt_cpue = as.numeric(wgt_cpue),
    wgt_cpua = as.numeric(wgt_cpua),
    num_cpue = as.numeric(num_cpue),
    num_cpua = as.numeric(num_cpua),
    year = as.numeric(year),
    month = as.numeric(month),
    day = as.numeric(day),
    sst = as.numeric(sst),
    sbt = as.numeric(sbt),
    depth = as.numeric(depth),
    area_swept = as.numeric(area_swept),
    haul_dur = as.numeric(haul_dur),
    wgt_combined = ifelse(
      !is.na(wgt_cpua) & wgt_cpua >= 0,
      wgt_cpua,
      wgt_cpue
    )
  )
Combined_FISH_optimized_spatially_valid <-
  Combined_FISH_optimized_spatially_valid %>%
  mutate(
    haul_id = paste(survey_unit, year, haul_id, sep = "_")
  )

##############################################################################
######## Q2 prep: Fix missing MTC models ##########################
# Aleutian islands, e bering sea, and gulf of alaska lack S combined, bit of T and fish as well
# rockall lacks fishing pressure
##############################################################################
##############################################################################
# APHIA-ID COVERAGE AUDIT FOR SURVEYS DROPPED FROM HEX_DAT
##############################################################################

Combined_FISH_optimized_spatially_valid %>%
  
  filter(
    survey_unit %in% c(
      "Aleutian Islands",
      "E Bering Sea",
      "Gulf of Alaska",
      "Rockall"
    )
  ) %>%
  
  summarise(
    
    n_rows = n(),
    
    n_wgt =
      sum(!is.na(wgt)),
    
    pct_wgt =
      round(
        100 * mean(!is.na(wgt),
                   na.rm = TRUE),
      ),
    
    n_wgt_cpua =
      sum(!is.na(wgt_cpua)),
    
    pct_wgt_cpua =
      round(
        100 * mean(!is.na(wgt_cpua),
                   na.rm = TRUE),
      ),
    
    n_wgt_cpue =
      sum(!is.na(wgt_cpue)),
    
    pct_wgt_cpue =
      round(
        100 * mean(!is.na(wgt_cpue),
                   na.rm = TRUE),
      ),
    
    n_haul_dur =
      sum(!is.na(haul_dur)),
    
    pct_haul_dur =
      round(
        100 * mean(!is.na(haul_dur),
                   na.rm = TRUE),
      ), 
    n_area =
      sum(!is.na(area_swept)),
    
    pct_area =
      round(
        100 * mean(!is.na(area_swept),
                   na.rm = TRUE),
      )
    
  )

#_________________________________________________#
##### Rebuild optimized hexagon geometries ########
#_________________________________________________#

res_lookup <- optimized_output$summary_table %>%
  select(
    survey_unit,
    opt_res
  ) %>%
  distinct()

hex_polygons_list <- list()

for (i in seq_len(nrow(res_lookup))) {
  
  s_unit <- res_lookup$survey_unit[i]
  s_res <- res_lookup$opt_res[i]
  
  ids <- Combined_FISH_optimized_spatially_valid %>%
    filter(survey_unit == s_unit) %>%
    pull(hexagon_id) %>%
    unique()
  
  if (length(ids) > 0) {
    
    dggs <- dgconstruct(res = s_res)
    
    polys <- dgcellstogrid(
      dggs,
      ids
    ) %>%
      st_as_sf() %>%
      st_make_valid() %>%
      st_wrap_dateline(
        options = c(
          "WRAPDATELINE=YES",
          "DATELINEOFFSET=180"
        )
      ) %>%
      mutate(
        survey_unit = s_unit,
        hexagon_id = ids,
        resolution = s_res
      )
    
    hex_polygons_list[[s_unit]] <- polys
  }
}

all_hex_geometries <- do.call(
  rbind,
  hex_polygons_list
) %>%
  st_as_sf() %>%
  st_cast("MULTIPOLYGON") %>%
  st_make_valid()

print(all_hex_geometries)

all_hex_geometries %>%
  st_drop_geometry() %>%
  count(survey_unit, sort = TRUE)

n_distinct(all_hex_geometries$hexagon_id)

#_________________________________________________#
##### Prepare plotting data #######################
#_________________________________________________#

hex_plot_data <- all_hex_geometries

# Use actual haul locations to determine label placement
# This avoids dateline-centroid problems (e.g. Aleutian Islands)

survey_centroids <- Combined_FISH_optimized_spatially_valid %>%
  group_by(survey_unit) %>%
  summarise(
    Longitude = median(longitude, na.rm = TRUE),
    Latitude = median(latitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Longitude_plot = ifelse(
      Longitude > 60,
      Longitude - 360,
      Longitude
    )
  )

#--------------------------------------------------#
# Survey IDs west --> east
#--------------------------------------------------#

survey_lookup <- survey_centroids %>%
  arrange(Longitude_plot) %>%
  mutate(
    survey_id = row_number()
  ) %>%
  select(
    survey_unit,
    survey_id
  )

hex_plot_data <- hex_plot_data %>%
  left_join(
    survey_lookup,
    by = "survey_unit"
  )

survey_centroids <- survey_centroids %>%
  left_join(
    survey_lookup,
    by = "survey_unit"
  )

ordered_surveys <- survey_lookup$survey_unit

hex_plot_data$survey_unit <- factor(
  hex_plot_data$survey_unit,
  levels = ordered_surveys
)

survey_centroids$survey_unit <- factor(
  survey_centroids$survey_unit,
  levels = ordered_surveys
)

survey_labels <- setNames(
  paste0(
    survey_lookup$survey_id,
    ": ",
    survey_lookup$survey_unit
  ),
  survey_lookup$survey_unit
)

#_________________________________________________#
##### Survey colours #############################
#_________________________________________________#

base_colors <- c(
  brewer.pal(8, "Dark2"),
  brewer.pal(12, "Paired"),
  brewer.pal(9, "Set1"),
  brewer.pal(8, "Accent")
)

survey_colors <- setNames(
  colorRampPalette(base_colors)(
    length(ordered_surveys)
  ),
  ordered_surveys
)

#_________________________________________________#
##### World map ##################################
#_________________________________________________#

world_map_wrapped <- ne_countries(
  scale = "medium",
  returnclass = "sf"
) %>%
  st_make_valid() %>%
  st_wrap_dateline(
    options = c(
      "WRAPDATELINE=YES",
      "DATELINEOFFSET=180"
    )
  )

highlight_countries <- c(
  "Canada",
  "Germany",
  "Iceland",
  "Norway",
  "United States of America",
  "Greenland",
  "France",
  "Italy",
  "Spain",
  "Malta",
  "Greece",
  "Slovenia",
  "Cyprus",
  "Croatia",
  "Faeroe Is.",
  "Ireland",
  "United Kingdom",
  "Portugal",
  "Denmark",
  "Netherlands",
  "Belgium",
  "Sweden",
  "Poland",
  "Lithuania",
  "Latvia",
  "Estonia"
)

world_map_plot <- world_map_wrapped %>%
  filter(
    continent %in% c(
      "North America",
      "Europe",
      "Africa",
      "Asia"
    ) |
      name %in% c(
        "Greenland",
        "Iceland"
      )
  )

#_________________________________________________#
##### Map limits #################################
#_________________________________________________#

lon_lims <- c(-195, 50)

lat_lims <- c(
  st_bbox(hex_plot_data)$ymin - 2,
  st_bbox(hex_plot_data)$ymax + 2
)

#_________________________________________________#
##### Plot #######################################
#_________________________________________________#

global_hex_map <- ggplot() +
  
  # Non-data-provider countries
  geom_sf(
    data = world_map_plot %>%
      filter(!(name %in% highlight_countries)),
    fill = "grey85",
    color = "black",
    linewidth = 0.25,
    show.legend = FALSE
  ) +
  
  # Data-provider countries only
  geom_sf(
    data = world_map_plot %>%
      filter(name %in% highlight_countries),
    fill = "lightblue",
    color = "black",
    linewidth = 0.25,
    show.legend = FALSE
  ) +
  
  # Hexagons
  geom_sf(
    data = hex_plot_data,
    aes(fill = survey_unit),
    color = "black",
    linewidth = 0.18,
    alpha = 0.95
  ) +
  
  # Survey numbers
  geom_text_repel(
    data = survey_centroids,
    aes(
      x = Longitude_plot,
      y = Latitude,
      label = survey_id
    ),
    size = 5,
    fontface = "bold",
    color = "black",
    box.padding = 0.4,
    point.padding = 0.2,
    segment.color = "grey30",
    segment.size = 0.3,
    max.overlaps = Inf
  ) +
  
  # Survey colours
  scale_fill_manual(
    values = survey_colors,
    labels = survey_labels,
    breaks = ordered_surveys,
    name = "Survey Unit"
  ) +
  
  coord_sf(
    xlim = lon_lims,
    ylim = lat_lims,
    crs = 4326,
    expand = FALSE
  ) +
  
  scale_x_continuous(
    breaks = seq(-180, 40, by = 15)
  ) +
  
  scale_y_continuous(
    breaks = seq(
      floor(lat_lims[1] / 5) * 5,
      ceiling(lat_lims[2] / 5) * 5,
      by = 5
    )
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  
  guides(
    fill = guide_legend(
      ncol = 10,
      byrow = TRUE,
      override.aes = list(
        alpha = 1,
        linewidth = 0.5
      )
    )
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    
    panel.grid.major = element_line(
      color = "grey75",
      linewidth = 0.35
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    
    legend.title = element_text(
      size = 12,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.key.height = grid::unit(0.20, "cm"),
    legend.key.width = grid::unit(0.30, "cm"),
    
    plot.margin = margin(
      3, 3, 3, 3
    )
  )

print(global_hex_map)

ggsave(
  "Combined global hex map.png",
  plot = global_hex_map,
  width = 19,
  height = 8,
  dpi = 400
)

##_______________________________________##
### facet data integration figures ###
##_______________________________________##


library(dplyr)
library(ggplot2)
library(patchwork)

#--------------------------------------------------#
# Prepare timeline data using same survey order
# as the map (west -> east)
#--------------------------------------------------#

timeline_data <- Combined_FISH_only %>%
  group_by(survey_unit) %>%
  summarise(
    start_year = min(year, na.rm = TRUE),
    end_year   = max(year, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    survey_lookup,
    by = "survey_unit"
  ) %>%
  arrange(survey_id)

#--------------------------------------------------#
# Timeline plot
#--------------------------------------------------#

timeline_plot <- ggplot(
  timeline_data,
  aes(
    x = start_year,
    xend = end_year,
    y = factor(survey_id),
    yend = factor(survey_id)
  )
) +
  
  geom_segment(
    linewidth = 2.5,
    color = "skyblue3",
    lineend = "round"
  ) +
  
  geom_point(
    aes(x = start_year),
    size = 2.5
  ) +
  
  geom_point(
    aes(x = end_year),
    size = 2.5
  ) +
  
  labs(
    x = "Year",
    y = "Survey Unit ID"
  ) +
  scale_x_continuous(
    breaks = seq(
      floor(min(timeline_data$start_year) / 10) * 10,
      ceiling(max(timeline_data$end_year) / 10) * 10,
      by = 10
    )
  ) +
  scale_y_discrete(
    limits = rev
  ) +
  
  theme_bw() +
  
  theme(
    axis.title.x = element_text(
      size = 14,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 14,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 10,
      color = "black"
    ),
    
    panel.grid.minor = element_blank(),
    
    plot.margin = margin(
      t = 5,
      r = 10,
      b = 35,
      l = 5
    ),
    
    legend.position = "none"
  )

#--------------------------------------------------#
# Remove duplicate legend from timeline
# Keep legend only from map
#--------------------------------------------------#

timeline_plot <- timeline_plot +
  theme(
    legend.position = "none"
  )


#--------------------------------------------------#
# Combined figure
#--------------------------------------------------#

combined_data_integration_plot <-
  timeline_plot +
  global_hex_map +
  plot_layout(
    widths = c(1.25, 3),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "bottom",
    
    # Larger legend
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    # Panel labels
    plot.tag = element_text(
      size = 20,
      face = "bold"
    ),
    
    plot.tag.position = c(
      0.02,
      0.98
    )
  )

#--------------------------------------------------#
# Export
#--------------------------------------------------#

print(combined_data_integration_plot)

ggsave(
  filename = "Survey_Timeline_and_Map.png",
  plot = combined_data_integration_plot,
  width = 24,
  height = 10,
  dpi = 400,
  bg = "white"
)
##_______________________________________##
### add secondary data to hexed dataset ###
##_______________________________________##

## LME shapefiles from NAFO

shapefile_path <- "C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/FISHGLOB Data Project/LME data/lme66.shp"
library(sf)
library(dplyr)

assign_lme <- function(df,
                       lat_col,
                       lon_col,
                       shapefile_path,
                       rockall_name = "Rockall") {
  
  # Read LME polygons
  lmes <- st_read(shapefile_path, quiet = TRUE)
  
  # Assign CRS if missing
  if (is.na(st_crs(lmes))) {
    st_crs(lmes) <- 4326
  }
  
  # Create sf points
  pts <- st_as_sf(
    df,
    coords = c(lon_col, lat_col),
    crs = 4326,
    remove = FALSE
  )
  
  # Match CRS
  lmes <- st_transform(lmes, st_crs(pts))
  
  ############################################################
  # FIRST PASS: assign polygons containing each haul
  ############################################################
  
  joined <- st_join(
    pts,
    lmes,
    join = st_within,
    left = TRUE
  )
  
  ############################################################
  # SECOND PASS: assign nearest LME to unmatched hauls
  # except Rockall
  ############################################################
  
  missing_idx <- which(
    is.na(joined$LME_NAME) &
      joined$survey_unit != rockall_name
  )
  
  if(length(missing_idx) > 0) {
    
    nearest_idx <- st_nearest_feature(
      joined[missing_idx, ],
      lmes
    )
    
    joined$LME_NAME[missing_idx] <-
      lmes$LME_NAME[nearest_idx]
  }
  
  joined %>%
    st_drop_geometry() %>%
    as.data.frame()
}

# Replace 'LME66.shp' with your actual filename
Combined_FISH_optimized_spatially_valid <- assign_lme(Combined_FISH_optimized_spatially_valid, 
                                                      lat_col = "latitude", 
                                                      lon_col = "longitude", 
                                                      shapefile_path = shapefile_path)


# Trophic level derivation

weighted_average_troph <- function(target_species, trophic_table){
  target_trophic_table <- trophic_table %>% filter(Species == target_species)
  observations <- target_trophic_table$Troph
  num_observations <- length(observations)
  se2 <- target_trophic_table$seTroph^2
  weights <- rep(1, num_observations)
  weights <- weights / se2
  weight_sum <- sum(weights)
  weights_obs <- mapply(`*`, observations, weights) # `*` is the element wise multiplication function
  sum_weights_obs <- sum(weights_obs)
  weighted_average <- sum_weights_obs / weight_sum
  se_weighted_average <- sqrt(1 / weight_sum)
  avg_trophic_info <- data.frame(Species = target_species, WE_trophic_estimate = weighted_average, SE_WE_trophic_estimate = se_weighted_average)
  return(avg_trophic_info)
}

assign_troph_fn <- function(df, my_dict) {
  key_list <- df$accepted_name
  value_list <- my_dict[key_list]
  df$Troph <- value_list
  return(df)
}

assign_SEtroph_fn <- function(df, my_dict) {
  key_list <- df$accepted_name
  value_list <- my_dict[key_list]
  df$SE_Troph <- value_list
  return(df)
}

species_list <- unique(Combined_FISH_optimized_spatially_valid$accepted_name) 
library(rfishbase)
species_list <- validate_names(species_list)

species_info <- species(species_list)

diet_table <- diet(species_list) # diet breakdown of global species including by lifestage

trophic_key <- data.frame(Species = diet_table$Species,Troph = diet_table$Troph, seTroph = diet_table$seTroph,  Remark = diet_table$Remark, Source = "Diet table") # trophic breakdown of global species with only relevant columns kept from diet info

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
    x = "Weighted Average Trophic Level Estimate from Fishbase Ecology and Diet Tables",
    y = "Frequency (Number of Species)"
  ) +
  scale_x_continuous(limits = c(2.0, 5.0)) +
  scale_y_continuous(limits = c(0, 175)) +
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


weighted_average_troph <- function(target_species, trophic_table){
  target_trophic_table <- trophic_table %>% filter(Species == target_species)
  observations <- target_trophic_table$Troph
  num_observations <- length(observations)
  se2 <- target_trophic_table$seTroph^2
  weights <- rep(1, num_observations)
  weights <- weights / se2
  weight_sum <- sum(weights)
  weights_obs <- mapply(`*`, observations, weights) # `*` is the element wise multiplication function
  sum_weights_obs <- sum(weights_obs)
  weighted_average <- sum_weights_obs / weight_sum
  se_weighted_average <- sqrt(1 / weight_sum)
  avg_trophic_info <- data.frame(Species = target_species, WE_trophic_estimate = weighted_average, SE_WE_trophic_estimate = se_weighted_average)
  return(avg_trophic_info)
}

assign_troph_fn <- function(df, my_dict) {
  key_list <- df$accepted_name
  value_list <- my_dict[key_list]
  df$Troph <- value_list
  return(df)
}

assign_SEtroph_fn <- function(df, my_dict) {
  key_list <- df$accepted_name
  value_list <- my_dict[key_list]
  df$SE_Troph <- value_list
  return(df)
}

Combined_FISH_optimized_spatially_valid <- assign_SEtroph_fn(Combined_FISH_optimized_spatially_valid, SE_troph_dict) # all species get trophic level estimate
Combined_FISH_optimized_spatially_valid <- assign_troph_fn(Combined_FISH_optimized_spatially_valid, troph_dict)

Combined_FISH_optimized_spatially_valid <- Combined_FISH_optimized_spatially_valid %>%
  mutate(
    trophic_group = case_when(
      Troph < 3 ~ "TL3",
      Troph >= 3.0 & Troph < 3.5 ~ "TL3.5",
      Troph >= 3.5 & Troph < 4.0 ~ "TL4",
      Troph >= 4.0 ~ "TL5",
      is.na(Troph) ~ NA )
  )

trophic_source_table <- avg_trophic_key
trophic_source_table <- trophic_source_table %>%
  mutate(
    trophic_group = case_when(
      WE_trophic_estimate < 3 ~ "TL3",
      WE_trophic_estimate >= 3.0 & WE_trophic_estimate < 3.5 ~ "TL3.5",
      WE_trophic_estimate >= 3.5 & WE_trophic_estimate < 4.0 ~ "TL4",
      WE_trophic_estimate >= 3.0 & WE_trophic_estimate < 4.0 ~ "TL5"
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

'
# Add missing bottom depth to haul data
#
# Rules:
#   1) Match to nearest etopo 2022 ice surface 15 arc seconds observatio
'

library(terra)
library(dplyr)

# Paths to all ETOPO tiles
etopo_files <- c(
  "~/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/Additional survey datasets/exportImage.tiff",
  "~/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/Additional survey datasets/exportImage (1).tiff",
  "~/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/Additional survey datasets/exportImage (2).tiff",
  "~/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/Additional survey datasets/exportImage (3).tiff",
  "~/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/Additional survey datasets/exportImage (4).tiff"
)

# Load rastersf
etopo_list <- lapply(etopo_files, rast)

# Unique hauls with missing depth
missing_depth_hauls <- Combined_FISH_optimized_spatially_valid %>%
  filter(is.na(depth)) %>%
  distinct(haul_id, latitude, longitude)

# Spatial points
pts <- vect(
  missing_depth_hauls,
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)

# Initialize depth vector
missing_depth_hauls$etopo_depth <- NA_real_

# Extract from each raster sequentially
for(i in seq_along(etopo_list)) {
  
  idx <- which(is.na(missing_depth_hauls$etopo_depth))
  
  if(length(idx) == 0) break
  
  vals <- extract(
    etopo_list[[i]],
    pts[idx, ],
    method = "bilinear"
  )[,2]
  
  missing_depth_hauls$etopo_depth[idx] <- vals
}

# Convert bathymetry to positive depth
missing_depth_hauls <- missing_depth_hauls %>%
  mutate(
    etopo_depth = abs(etopo_depth)
  )

# Join back and fill missing values
Combined_FISH_optimized_spatially_valid <-
  Combined_FISH_optimized_spatially_valid %>%
  left_join(
    missing_depth_hauls %>%
      select(haul_id, etopo_depth),
    by = "haul_id"
  ) %>%
  mutate(
    depth_imputed = is.na(depth) & !is.na(etopo_depth)
  ) 

# Join back and fill missing values
Combined_FISH_optimized_spatially_valid <-
  Combined_FISH_optimized_spatially_valid %>%
  mutate(
    depth = coalesce(depth, etopo_depth)
  ) 

#####_________________________________________####
### Begin analysis by haul and hex communities ###
#####_________________________________________####

unique_hexes <- unique(Combined_FISH_optimized_spatially_valid$hexagon_id)
unique_hauls <- unique(Combined_FISH_optimized_spatially_valid$haul_id)
missing_aphias <- Combined_FISH_optimized_spatially_valid %>%
  filter(is.na(aphia_id))


# data cleaning of valid set

spatial_valid_cleaned <- Combined_FISH_optimized_spatially_valid %>%
  filter(
    # 1. ESSENTIAL IDENTIFIERS
    !is.na(haul_id),
    !is.na(aphia_id) | !is.na(accepted_name),
    
    # 2. THE INCLUSIVE GATE
    # Keep the row if it has valid CPUA data OR valid CPUE data
    !is.na(wgt_combined) & wgt_combined >= 0,
    
    # 3. GENERAL SANITY
    wgt >= 0 | is.na(wgt),
    num >= 0 | is.na(num)
  )

test_Combined_FISH_optimized_spatially_valid <- Combined_FISH_optimized_spatially_valid[1:1000,]
# first logistical variables

global_spatially_valid_haul_stats <- spatial_valid_cleaned %>%
  group_by(haul_id) %>%
  summarise(
    # --- Metadata ---
    hex_id = first(hexagon_id), 
    survey_unit = first(survey_unit),
    year = first(year),
    month = first(month),
    day = first(day),
    lat = first(latitude), 
    lon = first(longitude),
    lme_name = first(LME_NAME),
    lme_number = first(LME_NUMBER),
    depth = first(depth),
    gear = first(gear),
    sst = first(sst), 
    sbt = if_else(first(sbt) > 40, NA_real_, first(sbt)),
    
    # --- Aggregated Biomass ---
    cpue = sum(wgt_cpue, na.rm = TRUE), 
    cpua = sum(wgt_cpua, na.rm = TRUE),
    num_cpue = sum(num_cpue, na.rm = TRUE),
    num_cpua = sum(num_cpua, na.rm = TRUE),
    wgt = sum(wgt, na.rm = TRUE), 
    num = sum(num, na.rm = TRUE), 
    
    # --- Not effort standardized Richness (S) ---
    S_aphia = n_distinct(aphia_id),
    
    # --- Diversity & Evenness ---
    # Evenness for CPUA
    evenness_cpua = if(S_aphia > 1) {
      diversity(wgt_cpua[!is.na(wgt_cpua)]) / log(S_aphia)
    } else { 0 },
    
    # Evenness for CPUE
    evenness_cpue = if(S_aphia > 1) {
      diversity(wgt_cpue[!is.na(wgt_cpue)]) / log(S_aphia)
    } else { 0 },
    
    .groups = "drop"
  ) %>%
  # --- The "Combined" Priority Logic ---
  mutate(
    # Prefer area-standardized richness density
    S_combined = S_aphia,
    
    evenness_combined = ifelse(
      evenness_cpua > 0,
      evenness_cpua,
      evenness_cpue
    ),
    
    density_combined = ifelse(
      cpua > 0,
      cpua,
      cpue
    ),
    
    body_weight_combined = case_when(
      wgt > 0 & num > 0       ~ wgt / num,
      cpua > 0 & num_cpua > 0 ~ cpua / num_cpua,
      cpue > 0 & num_cpue > 0 ~ cpue / num_cpue,
      TRUE                    ~ NA_real_
    )
  )

'
# Add GLORYS bottom temperature to haul community data
#
# Rules:
#   1) Match to nearest GLORYS grid cell
#   2) Use exact YYYY-MM-DD when day is known
#   3) Use YYYY-MM-15 when day is missing
#   4) Only assign temperatures for dates within the
#      GLORYS temporal coverage period
#   5) If an in-range lookup date is not found exactly,
#      use the nearest available GLORYS date
#   6) For dates outside GLORYS coverage, leave NA
'
library(ncdf4)
library(data.table)
library(FNN)

############################################################
# Open GLORYS file
############################################################

nc <- nc_open(
  "C:/Users/besch/cmems_mod_glo_phy_my_0.083deg_P1D-m_bottomT_179.92W-36.75E_23.08N-88.83N_1993-01-01-2026-05-26.nc"
)

lon  <- ncvar_get(nc, "longitude")
lat  <- ncvar_get(nc, "latitude")
time <- ncvar_get(nc, "time")

dates <- as.Date(
  as.POSIXct(
    time * 3600,
    origin = "1950-01-01",
    tz = "UTC"
  )
)

glorys_min_date <- min(dates)
glorys_max_date <- max(dates)

############################################################
# Working copy
############################################################

hauls <- copy(global_spatially_valid_haul_stats)

############################################################
# Spatial matching
############################################################

grid <- expand.grid(
  lon = lon,
  lat = lat
)

nn <- get.knnx(
  data  = as.matrix(grid),
  query = as.matrix(hauls[, c("lon", "lat")]),
  k = 1
)

hauls$grid_index <- nn$nn.index[, 1]

hauls$grid_lon <- grid$lon[hauls$grid_index]
hauls$grid_lat <- grid$lat[hauls$grid_index]

nlon <- length(lon)

hauls$lon_idx <- ((hauls$grid_index - 1) %% nlon) + 1
hauls$lat_idx <- ((hauls$grid_index - 1) %/% nlon) + 1

############################################################
# Temporal matching
############################################################

# Exact date when available, otherwise use 15th
hauls$lookup_date <- as.Date(
  ifelse(
    is.na(hauls$day),
    sprintf(
      "%04d-%02d-15",
      as.integer(hauls$year),
      as.integer(hauls$month)
    ),
    sprintf(
      "%04d-%02d-%02d",
      as.integer(hauls$year),
      as.integer(hauls$month),
      as.integer(hauls$day)
    )
  )
)

# Flag records within GLORYS temporal coverage
within_glorys_period <-
  !is.na(hauls$lookup_date) &
  hauls$lookup_date >= glorys_min_date &
  hauls$lookup_date <= glorys_max_date

hauls$time_idx <- NA_integer_

# Exact date lookup
date_lookup <- data.table(
  lookup_date = dates,
  time_idx = seq_along(dates)
)

date_matches <- merge(
  data.table(
    row_id = seq_len(nrow(hauls)),
    lookup_date = hauls$lookup_date
  ),
  date_lookup,
  by = "lookup_date",
  all.x = TRUE,
  sort = FALSE
)

hauls$time_idx[date_matches$row_id] <- date_matches$time_idx

# Use nearest GLORYS day only for dates inside coverage
need_nearest <-
  within_glorys_period &
  is.na(hauls$time_idx)

hauls$time_idx[need_nearest] <- vapply(
  hauls$lookup_date[need_nearest],
  function(d) {
    which.min(abs(as.numeric(dates - d)))
  },
  integer(1)
)

############################################################
# Extract GLORYS temperatures
############################################################

hauls$glorys_bottomT <- NA_real_

valid <-
  !is.na(hauls$lon_idx) &
  !is.na(hauls$lat_idx) &
  !is.na(hauls$time_idx)

for (i in which(valid)) {
  
  hauls$glorys_bottomT[i] <- ncvar_get(
    nc,
    "bottomT",
    start = c(
      as.integer(hauls$lon_idx[i]),
      as.integer(hauls$lat_idx[i]),
      as.integer(hauls$time_idx[i])
    ),
    count = c(1L, 1L, 1L)
  )
  
}

nc_close(nc)

############################################################
# Save back to master dataset
############################################################

global_spatially_valid_haul_stats$glorys_bottomT <-
  hauls$glorys_bottomT

############################################################
# Diagnostics
############################################################

cat("\n================ GLORYS COVERAGE ================\n")

cat(
  "\nGLORYS start:",
  as.character(glorys_min_date),
  "\n"
)

cat(
  "GLORYS end:  ",
  as.character(glorys_max_date),
  "\n"
)

cat(
  "\nOutside temporal coverage:",
  sum(!within_glorys_period),
  "\n"
)

cat(
  "Inside temporal coverage:",
  sum(within_glorys_period),
  "\n"
)

table(
  used_midmonth = is.na(global_spatially_valid_haul_stats$day),
  glorys_found = !is.na(global_spatially_valid_haul_stats$glorys_bottomT)
)

library(ggplot2)

glorys_sbt_agreement <- ggplot(
  global_spatially_valid_haul_stats,
  aes(x = sbt, y = glorys_bottomT)
) +
  geom_point(size = 0.4) +
  geom_abline(intercept = 0, slope = 1, colour = "red", linewidth = 1) +
  labs(
    title = "Agreement between haul level SBT measured by CTD and GLORYS model hindcasting",
    x = "CTD SBT",
    y = "GLORYS Bottom Temperature"
  ) +
  theme_bw()

ggsave(
  "Glorys vs sbt CTD agreement plot.png",
  plot = glorys_sbt_agreement,
  width = 9,
  height = 9,
  dpi = 300
)
############################################################
# Diagnose missing values
############################################################

temporal_failure <-
  !within_glorys_period

spatial_failure <-
  is.na(hauls$lon_idx) |
  is.na(hauls$lat_idx)

temp_failure <-
  is.na(hauls$glorys_bottomT)

cat("\n==================== SUMMARY ====================\n")

cat(
  "\nTotal hauls: ",
  nrow(hauls),
  "\n",
  sep = ""
)

cat(
  "Missing GLORYS bottomT: ",
  sum(temp_failure),
  "\n",
  sep = ""
)

cat(
  "\nOutside GLORYS temporal coverage: ",
  sum(temporal_failure),
  "\n",
  sep = ""
)

cat(
  "Spatial matching failures: ",
  sum(spatial_failure),
  "\n",
  sep = ""
)

cat(
  "Missing bottomT after successful spatial and temporal matching: ",
  sum(
    temp_failure &
      !temporal_failure &
      !spatial_failure
  ),
  "\n",
  sep = ""
)

############################################################
# Breakdown table
############################################################

diagnostics <- data.frame(
  temporal_failure = temporal_failure,
  spatial_failure = spatial_failure,
  missing_bottomT = temp_failure
)

cat("\n================ BREAKDOWN TABLE ================\n")

print(
  with(
    diagnostics,
    table(
      temporal_failure,
      spatial_failure,
      missing_bottomT
    )
  )
)

############################################################
# Examine failures that occurred despite valid matching
############################################################

failed_after_match <- hauls[
  temp_failure &
    !temporal_failure &
    !spatial_failure,
]

cat("\n=========== FAILED DESPITE VALID MATCH ===========\n")

cat(
  "Count: ",
  nrow(failed_after_match),
  "\n",
  sep = ""
)

if (nrow(failed_after_match) > 0) {
  
  cols_to_show <- intersect(
    c(
      "lon",
      "lat",
      "grid_lon",
      "grid_lat",
      "year",
      "month",
      "day",
      "lookup_date",
      "time_idx"
    ),
    names(failed_after_match)
  )
  
  print(
    head(
      failed_after_match[
        ,
        cols_to_show,
        drop = FALSE
      ]
    )
  )
  
  cat("\nFailures by year:\n")
  print(
    sort(
      table(failed_after_match$year),
      decreasing = TRUE
    )
  )
  
  cat("\nRange of failed longitudes:\n")
  print(summary(failed_after_match$lon))
  
  cat("\nRange of failed latitudes:\n")
  print(summary(failed_after_match$lat))
}

############################################################
# Final interpretation
############################################################

cat("\n================== CONCLUSION ==================\n")

n_temporal <- sum(temporal_failure)
n_spatial  <- sum(spatial_failure)

n_data_gap <- sum(
  temp_failure &
    !temporal_failure &
    !spatial_failure
)

if (n_temporal > 0) {
  cat(
    "Temporal coverage limitation:",
    n_temporal,
    "hauls occurred outside the GLORYS period and were intentionally left as NA.\n"
  )
}

if (n_spatial > 0) {
  cat(
    "Spatial coverage limitation:",
    n_spatial,
    "hauls could not be assigned a GLORYS grid cell.\n"
  )
}

if (n_data_gap > 0) {
  cat(
    "Masked-cell issue:",
    n_data_gap,
    "hauls had valid dates and valid grid cells but GLORYS returned NA.\n"
  )
}

if (
  n_temporal == 0 &&
  n_spatial == 0 &&
  n_data_gap > 0
) {
  cat(
    "\nInterpretation:\n",
    "Missing values are not caused by temporal coverage or domain coverage.\n",
    "They are likely due to coastal/land/masked GLORYS cells where bottomT is unavailable.\n",
    sep = ""
  )
}




#####################################################
##### add average trophic level to haul dataset #####
#####################################################
# attempt to assign trophic levels to species with good trophic information
keep_values <- avg_trophic_key$Species
# Filter the global dataframe to analyze only species with trophic information
troph_filtered_Combined_FISH_optimized_spatially_valid <- Combined_FISH_optimized_spatially_valid %>% # generates new df of SCS year round only including species whose trophic information could be found
  filter(accepted_name %in% keep_values)

haul_species_data <- troph_filtered_Combined_FISH_optimized_spatially_valid %>%
  mutate(wgt_cpua = replace_na(wgt_cpua, 0))

# 3. Propagate uncertainty to the Haul Level
library(dplyr)

haul_troph_summary <- haul_species_data %>%
  # Filter out NAs to prevent the whole calculation from returning NA
  filter(!is.na(Troph), !is.na(SE_Troph), !is.na(wgt_combined)) %>%
  group_by(haul_id, year) %>%
  summarize(
    # 1. Calculate the weighted average Trophic Level
    # Note: Using wgt_combined as an additional weighting factor 
    # if you want to account for the biomass of that species in the haul
    avg_troph = sum(Troph * wgt_combined) / sum(wgt_combined),
    
    # 2. Propagate the Standard Error
    # We use the inverse-variance formula: SE_total = sqrt(1 / sum(weights))
    # where weight = 1 / (SE^2)
    total_se_troph = sqrt(1 / sum(1 / (SE_Troph^2), na.rm = TRUE)),
    
    # Optional: Keep track of how many species were used for this haul
    n_species_troph = n()
  ) %>%
  ungroup()


global_spatially_valid_haul_stats <- global_spatially_valid_haul_stats %>%
  left_join(haul_troph_summary, by = c("haul_id", "year"))

#####################################################
###### thermal preferences from AquaX ###############
#####################################################

'Citation: Reygondeau G, Egorova Y, Boerder K, Tittensor DP, Kaschner K, Kesner-Reyes K, et al. (2026) AquaX: An enhanced and revised AquaMaps framework to model marine species distributions and biodiversity. PLoS One 21(2): e0335823. https://doi.org/10.1371/journal.pone.0335823
'
# V3 MTC
load("C:/Users/besch/Desktop/New Niche Files/META_02122025.Rdata")

# Strategy: use seafloor occ niche files when available, then use pelagic occ niche files if unavailable

id_list_for_MTC <- unique(Combined_FISH_optimized_spatially_valid$aphia_id) # 2303 things
META_FISHGLOB <- META %>%
  filter(AphiaID %in% id_list_for_MTC)

length(unique(META_FISHGLOB$AphiaID)) # 1795 things with niche data


############# add thermal data from individual r data files ########

# 1. Setup paths
data_dir <- "C:/Users/besch/Desktop/New Niche Files/NICHE_SEAFLOOR_OCC/"
AphiaID_list <- sort(unique(META_FISHGLOB$AphiaID))

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
FISHGLOB_META_and_THERMALS <- META_FISHGLOB %>%
  mutate(AphiaID = as.double(AphiaID)) %>% # Ensure types match
  left_join(seafloor_occ_thermals_wide_MEANS_ONLY, by = "AphiaID")

# Check how many species now have thermal data
sum(!is.na(FISHGLOB_META_and_THERMALS$thetao_mean_mean)) # 1401 aphiaIDs are benthic

reduced_FISHGLOB_META_and_THERMALS <- FISHGLOB_META_and_THERMALS %>%
  select(c(AphiaID, thetao_mean_mean, thetao_mean_min,thetao_mean_max,thetao_mean_q01,thetao_mean_q05,thetao_mean_q25,thetao_mean_q50,thetao_mean_q75,thetao_mean_q95, thetao_mean_q99))

# Now need thermal niches for species without benthic models
# 1. Setup paths
data_dir <- "C:/Users/besch/Desktop/New Niche Files/NICHE_PELAGIC_OCC/"
AphiaID_list <- sort(unique(META_FISHGLOB$AphiaID))

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
FISHGLOB_META_and_THERMALS <- META_FISHGLOB %>%
  mutate(AphiaID = as.double(AphiaID)) %>% # Ensure types match
  left_join(combined_pelagic_seafloor_occ_data, by = "AphiaID")

# Check how many species now have thermal data
sum(!is.na(FISHGLOB_META_and_THERMALS$thetao_mean_mean)) # 1749

missing_AphiaID <- setdiff(as.character(unique(META_FISHGLOB$AphiaID)), 
                           unique(combined_pelagic_seafloor_occ_data$AphiaID))
print(missing_AphiaID) 

# NO SDM FOR 
'
[1] "126367"  "126438"  "126447"  "126454"  "126499"  "126658"  "126790"  "126859"  "126862"  "127229"  "159277" 
[12] "272334"  "272616"  "272858"  "273991"  "274361"  "274402"  "274422"  "274454"  "274483"  "274546"  "274796" 
[23] "274835"  "275632"  "275694"  "279409"  "279415"  "279557"  "279680"  "280313"  "280948"  "281139"  "281630" 
[34] "282359"  "298167"  "306617"  "306769"  "578660"  "1021330" "1042864" "1384967" "1493316" "103700"  "103815" 
[45] "103937"  "1307332" "159019" 
'

colnames(FISHGLOB_META_and_THERMALS)

reduced_FISHGLOB_META_and_THERMALS <- FISHGLOB_META_and_THERMALS %>%
  select(c(AphiaID, thetao_mean_mean, thetao_mean_min,thetao_mean_max,thetao_mean_q01,thetao_mean_q05,thetao_mean_q25,thetao_mean_q50,thetao_mean_q75,thetao_mean_q95, thetao_mean_q99))
reduced_FISHGLOB_META_and_THERMALS$AphiaID <- as.character(reduced_FISHGLOB_META_and_THERMALS$AphiaID)

##############################################################
#### calculate coverage of MTC3 for global #########
##############################################################
global_temp3_hauls <- left_join(
  Combined_FISH_optimized_spatially_valid, 
  reduced_FISHGLOB_META_and_THERMALS, 
  by = c("aphia_id" = "AphiaID")
)


global_temp3_hauls <- global_temp3_hauls %>%
  filter(!is.na(thetao_mean_mean))

nhaul_global_MTC3 <- length(unique(global_temp3_hauls$haul_id)) # 276,550
ntaxa_global_MTC3 <- length(unique(global_temp3_hauls$aphia_id)) # 1749

######## MTC Derivation (Haul Level) ################

# 1. Enhanced Thermal Dictionary Preparation
# This uses multiple points of the distribution to reach a consensus SD
simplified_thermal_dictionary <- reduced_FISHGLOB_META_and_THERMALS %>%
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

########################### Combine thermal data into main dataset #################################

Combined_FISH_optimized_spatially_valid <- left_join( # note NA's present thrughout!
  Combined_FISH_optimized_spatially_valid, 
  simplified_thermal_dictionary, 
  by = c("aphia_id" = "AphiaID")
)

# 2. Join and MTC Calculation (Logic remains consistent)
haul_species_data <- global_temp3_hauls %>%
  mutate(wgt_combined = replace_na(wgt_combined, 0)) %>%
  left_join(simplified_thermal_dictionary, by = c("aphia_id" = "AphiaID")) %>%
  filter(!is.na(estim_thermal_mean))

# 1. Join and Calculate Composite Weights
haul_mtc_results <- global_temp3_hauls %>%
  mutate(wgt_combined = replace_na(wgt_combined, 0)) %>%
  left_join(simplified_thermal_dictionary, by = c("aphia_id" = "AphiaID")) %>%
  filter(!is.na(estim_thermal_mean), !is.na(Variance_Ti)) %>%
  mutate(
    # Create a composite weight: Biomass divided by Uncertainty
    # We add a tiny epsilon (1e-10) to avoid division by zero if Var is 0
    comp_weight = wgt_combined / (Variance_Ti + 1e-10)
  )

# 2. Propagate uncertainty to the Haul Level
haul_mtc_summary <- haul_mtc_results %>%
  group_by(haul_id, year) %>%
  summarise(
    # MTC now influenced by both biomass and precision
    MTC_haul = sum(estim_thermal_mean * comp_weight) / sum(comp_weight),
    
    # Propagated Variance for a Weighted Mean
    # Formula: Sum(Wi^2 * Var_i) / (Sum Wi)^2
    MTC_variance_Ti = sum(comp_weight^2 * Variance_Ti) / (sum(comp_weight)^2),
    
    total_wgt = sum(wgt_combined),
    n_species = n(),
    .groups = "drop"
  ) %>%
  filter(total_wgt > 0)

# 3. Final Join (Using both haul_id and year for safety)
global_spatially_valid_haul_stats <- global_spatially_valid_haul_stats %>%
  left_join(
    haul_mtc_summary %>% select(haul_id, year, MTC_haul, MTC_variance_Ti), 
    by = c("haul_id", "year")
  )

##############################################
##### assemble trait dataset ################
##############################################

###### Trait database from Beukhof et al 2019 #######
Atlantic_Pacific_trait_database <- read.csv("C:/Users/besch/OneDrive - Dalhousie University/Documents/Dalhousie Biology/Post acceptance/Masters Thesis/Data files/FISHGLOB Data Project/trait data for North Atlantic and Northeast Pacific groundfish/TraitDatabase.csv")
colnames(Atlantic_Pacific_trait_database)

unique_fish_taxa_combined <- unique(Combined_FISH_optimized_spatially_valid$accepted_name) # 2,315
traited_fish <- unique(Atlantic_Pacific_trait_database$taxon) # 2,010 taxa
traited_fish_in_study <- intersect(traited_fish, unique_fish_taxa_combined) #1,740 taxa

global_fish_trait_database <- Atlantic_Pacific_trait_database %>%
  filter(taxon %in% traited_fish_in_study) 

"
Continue trait investigation by matching trait data to taxa dependent on spatial scale
"
possible_traits_to_test <- colnames(global_fish_trait_database)

'
possible_traits_to_test
 [1] "family"                   "genus"                    "species"                  "taxonomic.rank"           "taxon"                    "LME"                     
 [7] "FAO"                      "habitat"                  "reference.habitat"        "feeding.mode"             "reference.feeding.mode"   "tl"                      
[13] "level.tl"                 "reference.tl"             "body.shape"               "fin.shape"                "AR"                       "level.AR"                
[19] "reference.AR"             "offspring.size"           "level.offspring.size"     "reference.offspring.size" "spawning.type"            "reference.spawning.type" 
[25] "age.maturity"             "level.age.maturity"       "reference.age.maturity"   "fecundity"                "level.fecundity"          "reference.fecundity"     
[31] "length.infinity"          "growth.coefficient"       "level.growth"             "reference.growth"         "length.max"               "level.length.max"        
[37] "reference.length.max"     "age.max"                  "level.age.max"            "reference.age.max"       
'
library(dplyr)
library(purrr)

# Define the columns you want to pull from the trait data
trait_cols <- c("habitat", "feeding.mode", "tl", "body.shape", "fin.shape", 
                "AR", "offspring.size", "spawning.type", "age.maturity", 
                "fecundity", "length.infinity", "growth.coefficient", 
                "length.max", "age.max")

# Helper function to summarize trait data at different levels
summarize_traits <- function(df, group_cols) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(across(all_of(trait_cols), 
                     ~if(is.numeric(.x)) mean(.x, na.rm = TRUE) else first(.x)),
              .groups = "drop")
}

# Create lookup levels
# 1. Species + LME (Highest accuracy)
traits_sp_lme <- summarize_traits(global_fish_trait_database, c("taxon", "LME"))

# 2. Genus + LME
traits_gen_lme <- summarize_traits(global_fish_trait_database, c("genus", "LME"))

# 3. Family + LME
traits_fam_lme <- summarize_traits(global_fish_trait_database, c("family", "LME"))

# 4. Global averages (Species/Genus/Family) for when LME data is missing
traits_sp_global  <- summarize_traits(global_fish_trait_database, "taxon")
traits_gen_global <- summarize_traits(global_fish_trait_database, "genus")
traits_fam_global <- summarize_traits(global_fish_trait_database, "family")

Combined_FISH_optimized_spatially_valid <- Combined_FISH_optimized_spatially_valid %>%
  # Ensure LME_NUMBER is an integer to match the trait data
  mutate(LME_NUMBER = as.integer(LME_NUMBER)) %>%
  
  # Join 1: Species + LME (Use LME_NUMBER now)
  left_join(traits_sp_lme, by = c("accepted_name" = "taxon", "LME_NUMBER" = "LME"), suffix = c("", ".sp_lme")) %>%
  
  # Join 2: Genus + LME
  left_join(traits_gen_lme, by = c("genus" = "genus", "LME_NUMBER" = "LME"), suffix = c("", ".gen_lme")) %>%
  
  # Join 3: Family + LME
  left_join(traits_fam_lme, by = c("family" = "family", "LME_NUMBER" = "LME"), suffix = c("", ".fam_lme")) %>%
  
  # Join 4: Global Species (No LME needed here)
  left_join(traits_sp_global, by = c("accepted_name" = "taxon"), suffix = c("", ".sp_glob")) %>%
  
  # Join 5: Global Genus
  left_join(traits_gen_global, by = c("genus" = "genus"), suffix = c("", ".gen_glob"))

# 1. Define the traits you are working with (from your earlier list)
trait_cols <- c("habitat", "feeding.mode", "tl", "body.shape", "fin.shape", 
                "AR", "offspring.size", "spawning.type", "age.maturity", 
                "fecundity", "length.infinity", "growth.coefficient", 
                "length.max", "age.max")

# 2. Collapse the columns
# This loop takes the 'best' version and fills the original column
for (trait in trait_cols) {
  
  # Identify the extra columns created by the joins
  # Note: The first join likely didn't have a suffix, or had one you specified
  sp_lme   <- trait 
  gen_lme  <- paste0(trait, ".gen_lme")
  fam_lme  <- paste0(trait, ".fam_lme")
  sp_glob  <- paste0(trait, ".sp_glob")
  gen_glob <- paste0(trait, ".gen_glob")
  
  Combined_FISH_optimized_spatially_valid <- Combined_FISH_optimized_spatially_valid %>%
    mutate(!!trait := coalesce(!!sym(sp_lme),   # 1st Priority: Sp + LME
                               !!sym(gen_lme),  # 2nd Priority: Gen + LME
                               !!sym(fam_lme),  # 3rd Priority: Fam + LME
                               !!sym(sp_glob),  # 4th Priority: Sp Global
                               !!sym(gen_glob)  # 5th Priority: Gen Global
    ))
}

Combined_FISH_optimized_spatially_valid <- Combined_FISH_optimized_spatially_valid %>%
  mutate(trait_resolution = case_when(
    !is.na(tl) ~ "Species_LME",
    !is.na(tl.gen_lme) ~ "Genus_LME",
    !is.na(tl.fam_lme) ~ "Family_LME",
    !is.na(tl.sp_glob) ~ "Species_Global",
    TRUE ~ "None"
  ))

# 3. Drop the 'ghost' columns to clean up the workspace
Combined_FISH_optimized_spatially_valid <- Combined_FISH_optimized_spatially_valid %>%
  select(-ends_with(".gen_lme"), 
         -ends_with(".fam_lme"), 
         -ends_with(".sp_glob"), 
         -ends_with(".gen_glob"))

# now with catch level traits, build haul community stats weighting by combined weight
# Function to find the category with the highest total biomass
get_weighted_dominant <- function(x, w) {
  # Clean out NAs from both trait and weight
  valid_indices <- !is.na(x) & !is.na(w)
  x <- x[valid_indices]
  w <- w[valid_indices]
  
  if (length(x) == 0) return(NA_character_)
  
  # Aggregate weight by category
  weights_by_cat <- aggregate(w ~ x, FUN = sum)
  # Return the category with the max weight
  weights_by_cat$x[which.max(weights_by_cat$w)]
}

# Identify which columns are which
num_traits <- trait_cols[sapply(Combined_FISH_optimized_spatially_valid[trait_cols], is.numeric)]
cat_traits <- trait_cols[!sapply(Combined_FISH_optimized_spatially_valid[trait_cols], is.numeric)]

Haul_Weighted_Trait_Summary <- Combined_FISH_optimized_spatially_valid %>%
  group_by(haul_id, year) %>% 
  summarise(
    # 1. Weighted Mean for Numeric Traits
    across(
      all_of(num_traits), 
      ~ weighted.mean(.x, w = wgt_combined, na.rm = TRUE)
    ),
    # 2. Weighted Dominance for Categorical Traits
    across(
      all_of(cat_traits), 
      ~ get_weighted_dominant(.x, w = wgt_combined)
    ),
    .groups = "drop"
  )

global_spatially_valid_haul_stats <- left_join(global_spatially_valid_haul_stats, Haul_Weighted_Trait_Summary, by = c("haul_id", "year")) # now hauls have dominant or mean trait data

get_dominant <- function(x) {
  x <- x[!is.na(x)] 
  if (length(x) == 0) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

global_spatially_valid_hex_stats <- global_spatially_valid_haul_stats %>%
  group_by(survey_unit, year, hex_id) %>%
  summarise(
    lme_number = first(lme_number),
    lme_name = first(lme_name),
    lat = mean(lat, na.rm = TRUE),
    lon = mean(lon, na.rm = TRUE),
    depth = mean(depth, na.rm = TRUE),
    CV_sst = sd(sst, na.rm = TRUE) / mean(sst, na.rm = TRUE),
    sst = mean(sst, na.rm = TRUE),
    CV_sbt = sd(sbt, na.rm = TRUE) / mean(sbt, na.rm = TRUE),
    
    ########################################################
    # TEMPERATURES
    ########################################################
    
    sbt = if (
      all(is.na(sbt))
    ) {
      NA_real_
    } else {
      mean(sbt, na.rm = TRUE)
    },
    
    glorys_bottomT = if (
      all(is.na(glorys_bottomT))
    ) {
      NA_real_
    } else {
      mean(glorys_bottomT, na.rm = TRUE)
    },
    
    ########################################################
    
    evenness_combined = mean(evenness_combined, na.rm = TRUE),
    CV_density_combined = sd(density_combined, na.rm = TRUE) / mean(density_combined, na.rm = TRUE),
    density_combined = mean(density_combined, na.rm = TRUE), 
    S_combined = mean(S_combined, na.rm = TRUE),
    body_weight_combined = mean(body_weight_combined, na.rm = TRUE),
    
    
    tl = mean(tl, na.rm = TRUE), 
    AR = mean(AR, na.rm = TRUE),
    offspring.size = mean(offspring.size, na.rm = TRUE),
    age.maturity = mean(age.maturity, na.rm = TRUE), 
    fecundity = mean(age.maturity, na.rm = TRUE), 
    length.infinity = mean(length.infinity, na.rm = TRUE), 
    growth.coefficient = mean(growth.coefficient, na.rm = TRUE),
    length_max = mean(length.max, na.rm = TRUE),
    age.max = mean(age.max, na.rm = TRUE),
    
    # Categorical traits (Dominant value across hauls)
    dom_habitat    = get_dominant(habitat),
    dom_feeding    = get_dominant(feeding.mode),
    dom_spawning   = get_dominant(spawning.type),
    dom_body.shape = get_dominant(body.shape), 
    dom_fin.shape = get_dominant(fin.shape),
    gear = get_dominant(gear),
    .groups = "drop"
  )

################################################
##### Begin MTC change test ####################
################################################

calculate_yearly_hex_mtc <- function(df) {
  
  yearly_hex_stats <- df %>%
    # 1. Filter out rows where MTC or Variance is missing
    filter(!is.na(MTC_haul), !is.na(MTC_variance_Ti), MTC_variance_Ti > 0) %>%
    
    group_by(survey_unit, year, hex_id) %>%
    summarise(
      # Calculate weights (Inverse Variance)
      # We sum the weights first for the denominator
      total_weight = sum(1 / MTC_variance_Ti, na.rm = TRUE),
      
      # 2. Calculate Weighted Average MTC
      # Sum of (Value * Weight) / Total Weight
      MTC_hex_avg = sum(MTC_haul * (1 / MTC_variance_Ti), na.rm = TRUE) / total_weight,
      
      # 3. Calculate the Variance of the weighted mean
      # This represents the uncertainty of the hexagon's yearly estimate
      MTC_hex_variance = 1 / total_weight,
      
      # 4. Metadata: Keep track of haul density per hex
      n_hauls = n(),
      
      # Optional: Carry over LME info (assuming one hex doesn't straddle LMEs)
      lme_name = first(lme_name),
      
      .groups = "drop"
    ) %>%
    # Clean up helper columns
    select(-total_weight)
  
  return(yearly_hex_stats)
}

# Apply the function
global_MTC_by_hex <- calculate_yearly_hex_mtc(global_spatially_valid_haul_stats)


calculate_survey_yearly_mtc <- function(hex_df) {
  
  survey_yearly_stats <- hex_df %>%
    # 1. Filter out invalid variances to avoid infinite weights
    # Note: We are now checking 'MTC_hex_variance'
    filter(!is.na(MTC_hex_avg), !is.na(MTC_hex_variance), MTC_hex_variance > 0) %>%
    
    group_by(survey_unit, year) %>%
    summarise(
      # FIX: Use MTC_hex_variance instead of MTC_variance_Ti
      total_survey_weight = sum(1 / MTC_hex_variance, na.rm = TRUE),
      
      # 2. Survey Unit Weighted Average MTC
      MTC_survey_avg = sum(MTC_hex_avg * (1 / MTC_hex_variance), na.rm = TRUE) / total_survey_weight,
      
      # 3. Survey Unit Variance (Uncertainty of the survey-wide mean)
      MTC_survey_variance = 1 / total_survey_weight,
      
      # 4. Metadata for context
      n_hexagons_sampled = n_distinct(hex_id),
      total_hauls_in_survey = sum(n_hauls, na.rm = TRUE),
      
      # Optional: List LMEs covered by this survey unit
      lmes_included = paste(sort(unique(lme_name)), collapse = "; "),
      
      .groups = "drop"
    ) %>%
    select(-total_survey_weight)
  
  return(survey_yearly_stats)
}

# Now run it again
global_MTC_by_year <- calculate_survey_yearly_mtc(global_MTC_by_hex)

library(dplyr)
library(ggplot2)

# Define the region mapping
library(dplyr)
library(ggplot2)

assign_ocean_region_comprehensive <- function(df) {
  df %>%
    mutate(ocean_region = case_when(
      # --- PACIFIC ---
      survey_unit %in% c("AI", "EBS", "GOA", "WCANN", "WCTRI", 
                         "DFO-HS", "DFO-SOG", "DFO-QCS", "DFO-WCHG", "DFO-WCVI") ~ "Pacific",
      
      # --- NORTHWEST ATLANTIC / ARCTIC ---
      # Includes DFO (Canada), NEUS/SEUS/GMEX (USA), and Greenland
      survey_unit %in% c("DFO-GH", "DFO-NS", "GSL-N", "GSL-S", 
                         "SCS-FALL", "SCS-SPRING", "SCS-SUMMER", 
                         "NEUS-Fall", "NEUS-Spring", "SEUS-fall", "SEUS-spring", "SEUS-summer",
                         "GMEX-Fall", "GMEX-Summer", 
                         "West Greenland", "East Greenland", "DF0-NL-Winter", "DFO-NL-Spring") ~ "Northwest Atlantic / Arctic",
      
      # --- NORTHEAST ATLANTIC / ARCTIC ---
      # Includes ICES areas, Norway, Iceland, and Iberia
      survey_unit %in% c("BITS-1", "BITS-4", "EVHOE", "FR-CGFS", "IE-IGFS", 
                         "NIGFS", "NS-IBTS-1", "NS-IBTS-3", "PT-IBTS", 
                         "ROCKALL", "SWC-IBTS-1", "SWC-IBTS-4", 
                         "IMR Norwegian Sea", "IS", "SP-ARSA-1", 
                         "SP-ARSA-4", "SP-NORTH", "SP-PORC", 
                         "yvirlitstrol_fb_summar", "yvirlitstrol_fb_var", "yvirlitstrol_lg_summar", "yvirlitstrol_lg_var") ~ "Northeast Atlantic / Arctic",
      
      TRUE ~ "Other"
    ))
}

# Apply to your aggregated data
global_MTC_by_year <- assign_ocean_region_comprehensive(global_MTC_by_year)

library(RColorBrewer)

# 1. Ensure the mapping is applied
global_MTC_by_year <- assign_ocean_region_comprehensive(global_MTC_by_year)

# 2. Define the regions to loop through
target_regions <- c("Pacific", "Northwest Atlantic / Arctic", "Northeast Atlantic / Arctic")

for (reg in target_regions) {
  
  # Filter data for the specific region and calculate CI
  region_data <- global_MTC_by_year %>% 
    filter(ocean_region == reg) %>%
    mutate(
      se = sqrt(MTC_survey_variance),
      ci_l = MTC_survey_avg - (1.96 * se),
      ci_u = MTC_survey_avg + (1.96 * se)
    )
  
  # Create the plot
  p <- ggplot(region_data, aes(x = year, y = MTC_survey_avg, color = survey_unit, fill = survey_unit)) +
    # Confidence intervals
    geom_ribbon(aes(ymin = ci_l, ymax = ci_u), alpha = 0.1, color = NA) +
    # Trend lines
    geom_line() +
    # Constant size points
    geom_point(size = 1) +
    # Standard BW theme
    theme_bw() +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      # 1. Make the axis Titles bigger
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold"),
      # 2. Make the axis Numbers (labels) bigger
      axis.text.y = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 14, color = "black"),
      legend.text = element_text(size = 14),
      legend.title = element_text(face = "bold", size = 16)
    ) +
    labs(
      title = reg,
      x = "Year",
      y = "MTC (°C)",
      color = "Survey Unit"
    ) +
    # Ensure fill legend is removed to keep it clean
    guides(fill = "none")
  
  print(p)
}

################################################################################
##### NORTH AMERICA-CENTERED PERSISTENT MTC ANIMATION ##########################
################################################################################

library(dggridR)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rnaturalearth)
library(magick)

# 1. DEFINE PROJECTION ---------------------------------------------------------
# We use a Lambert Azimuthal Equal Area or a simple longlat centered at -100
# but the key is how we handle the 'cut' line to avoid horizontal streaks.
target_crs <- "+proj=longlat +datum=WGS84 +no_defs"

# 2. RECONSTRUCT & WRAP DATA ---------------------------------------------------
res_lookup <- optimized_output$summary_table %>% 
  select(survey_unit, opt_res) %>% 
  distinct()

hex_polygons_list <- list()

for(i in 1:nrow(res_lookup)) {
  s_unit <- res_lookup$survey_unit[i]
  s_res  <- res_lookup$opt_res[i]
  
  ids <- global_MTC_by_hex %>% 
    filter(survey_unit == s_unit) %>% 
    pull(hex_id) %>% 
    unique()
  
  if(length(ids) > 0) {
    dggs <- dgconstruct(res = s_res)
    polys <- dgcellstogrid(dggs, ids) %>% 
      st_as_sf() %>%
      # This handles the internal polygons crossing the dateline 
      # before we project/center them.
      st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")) %>%
      mutate(survey_unit = s_unit, hex_id = ids)
    
    hex_polygons_list[[s_unit]] <- polys
  }
}

all_hex_geometries <- do.call(rbind, hex_polygons_list) %>%
  st_cast("MULTIPOLYGON")

# 3. PERSISTENT DATA (Carry forward MTC) ---------------------------------------
full_year_range <- min(global_MTC_by_hex$year):max(global_MTC_by_hex$year)

global_MTC_persistent <- global_MTC_by_hex %>%
  group_by(survey_unit, hex_id) %>%
  complete(year = full_year_range) %>%
  fill(MTC_hex_avg, .direction = "down") %>%
  filter(!is.na(MTC_hex_avg)) %>%
  ungroup()

# 4. JOIN & RECENTERING --------------------------------------------------------
animation_ready_data <- all_hex_geometries %>%
  inner_join(global_MTC_persistent, by = c("survey_unit", "hex_id")) %>%
  mutate(year = as.integer(year))

# Calculate the bounding box in standard coordinates first
bbox <- st_bbox(animation_ready_data)

# Load and prepare world map
# We use st_break_antimeridian to ensure Russia/Alaska don't stretch across the map
world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid()

# 5. RENDER FRAMES -------------------------------------------------------------
dir.create("mtc_na_centered_frames", showWarnings = FALSE)
all_years <- sort(unique(animation_ready_data$year))

# Set the limits manually based on your survey units to ensure North America is center
# This covers roughly Russia (East) to Europe (West) centered on the Americas
lon_lims <- c(-195, 50)

lat_lims <- c(
  st_bbox(hex_plot_data)$ymin - 2,
  st_bbox(hex_plot_data)$ymax + 2
)

mtc_limits <- c(min(animation_ready_data$MTC_hex_avg), max(animation_ready_data$MTC_hex_avg))

for (yr in all_years) {
  yr_data <- animation_ready_data %>% filter(year == yr)
  
  p <- ggplot() +
    # The world map
    #geom_sf(data = world_map, fill = "gray95", color = "gray80", linewidth = 0.1) +
    # Non-data-provider countries
    geom_sf(
      data = world_map_plot %>%
        filter(!(name %in% highlight_countries)),
      fill = "grey95",
      color = "black",
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    
    # Data-provider countries only
    geom_sf(
      data = world_map_plot %>%
        filter(name %in% highlight_countries),
      fill = "grey55",
      color = "black",
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    # The hex data
    geom_sf(data = yr_data, aes(fill = MTC_hex_avg), color = NA) +
    scale_fill_viridis_c(
      option = "turbo",
      limits = mtc_limits,
      oob = scales::squish,
      na.value = "transparent",
      name = "CTI (°C)"
    ) +
    # This coord_sf call with a specific CRS and limits centers the view
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE, crs = 4326) +

        scale_x_continuous(
      breaks = seq(-180, 40, by = 15)
    ) +
    
    scale_y_continuous(
      breaks = seq(
        floor(lat_lims[1] / 5) * 5,
        ceiling(lat_lims[2] / 5) * 5,
        by = 5
      )
    ) +
    labs(title = paste("Year:", yr)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 28, face = "bold"),
      axis.text = element_text(size = 12),
      legend.key.height = grid::unit(2, "cm"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    )
  
  ggsave(filename = sprintf("mtc_na_centered_frames/frame_%04d.png", yr),
         plot = p, width = 14, height = 9, dpi = 400)
}

# 6. STITCH --------------------------------------------------------------------
img_sequence <- image_read(list.files("mtc_na_centered_frames", full.names = TRUE))
animation <- image_animate(img_sequence, fps = 2, dispose = "background")
image_write(animation, "MTC_NA_Centered_HighRes.gif")

###############################################################################
##### HEX TEMPERATURE ANIMATION (SBT PRIORITY, GLORYS FALLBACK) ################
###############################################################################

library(dplyr)
library(tidyr)
library(dggridR)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(magick)
library(viridis)

############################################################
# 1. CREATE HEX TEMPERATURE VARIABLE
############################################################

global_temperature_by_hex <- global_spatially_valid_hex_stats %>%
  mutate(
    hex_temperature = case_when(
      !is.na(sbt) ~ sbt,
      is.na(sbt) & !is.na(glorys_bottomT) ~ glorys_bottomT,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(hex_temperature)) %>%
  select(
    survey_unit,
    year,
    hex_id,
    lme_name,
    hex_temperature
  ) %>%
  mutate(
    hex_id = as.character(hex_id)
  )

############################################################
# 2. REBUILD HEX GEOMETRIES
############################################################

res_lookup <- optimized_output$summary_table %>%
  select(survey_unit, opt_res) %>%
  distinct()

hex_polygons_list <- list()

for(i in seq_len(nrow(res_lookup))) {
  
  s_unit <- res_lookup$survey_unit[i]
  s_res  <- res_lookup$opt_res[i]
  
  survey_hexes <- global_temperature_by_hex %>%
    filter(survey_unit == s_unit) %>%
    distinct(hex_id)
  
  ids_char <- survey_hexes$hex_id
  
  # DGGRID requires numeric IDs
  ids_num <- suppressWarnings(as.numeric(ids_char))
  
  if(length(ids_num) == 0) next
  
  if(any(is.na(ids_num))) {
    
    warning(
      paste(
        "Survey",
        s_unit,
        "contains non-numeric hex IDs and was skipped."
      )
    )
    
    next
  }
  
  dggs <- dgconstruct(res = s_res)
  
  polys <- dgcellstogrid(
    dggs,
    ids_num
  ) %>%
    st_as_sf()
  
  # Verify geometry count matches ID count
  if(nrow(polys) != length(ids_char)) {
    
    warning(
      paste(
        "Geometry count mismatch for",
        s_unit,
        "- polygons:",
        nrow(polys),
        "IDs:",
        length(ids_char)
      )
    )
  }
  
  # Preserve original character hex IDs
  polys$hex_id <- ids_char[seq_len(min(
    nrow(polys),
    length(ids_char)
  ))]
  
  polys <- polys %>%
    st_wrap_dateline(
      options = c(
        "WRAPDATELINE=YES",
        "DATELINEOFFSET=180"
      )
    ) %>%
    st_make_valid() %>%
    mutate(
      survey_unit = s_unit
    )
  
  hex_polygons_list[[as.character(s_unit)]] <- polys
}

all_hex_geometries <- do.call(
  rbind,
  hex_polygons_list
) %>%
  st_make_valid() %>%
  st_cast("MULTIPOLYGON") %>%
  mutate(
    hex_id = as.character(hex_id)
  )

############################################################
# 3. PERSIST TEMPERATURE THROUGH TIME
############################################################

full_year_range <-
  min(global_temperature_by_hex$year):
  max(global_temperature_by_hex$year)

global_temperature_persistent <-
  global_temperature_by_hex %>%
  group_by(survey_unit, hex_id) %>%
  complete(year = full_year_range) %>%
  arrange(year, .by_group = TRUE) %>%
  fill(hex_temperature, .direction = "down") %>%
  filter(!is.na(hex_temperature)) %>%
  ungroup()

############################################################
# 4. JOIN GEOMETRY
############################################################

missing_hexes <- anti_join(
  global_temperature_persistent,
  st_drop_geometry(all_hex_geometries),
  by = c("survey_unit", "hex_id")
)

cat(
  "\nMissing geometry matches:",
  nrow(missing_hexes),
  "\n"
)

if(nrow(missing_hexes) > 0) {
  
  cat("\nSurveys with missing geometry:\n")
  
  print(
    missing_hexes %>%
      count(survey_unit, sort = TRUE),
    n = Inf
  )
}

animation_ready_data <-
  all_hex_geometries %>%
  inner_join(
    global_temperature_persistent,
    by = c("survey_unit", "hex_id")
  ) %>%
  mutate(
    year = as.integer(year)
  )

cat(
  "\nRows in animation dataset:",
  nrow(animation_ready_data),
  "\n"
)

############################################################
# 4B. DIAGNOSTICS
############################################################

cat("\n========================================\n")
cat("ANIMATION DIAGNOSTICS\n")
cat("========================================\n")

survey_diagnostics <-
  animation_ready_data %>%
  st_drop_geometry() %>%
  group_by(survey_unit) %>%
  summarise(
    first_anim_year = min(year),
    last_anim_year  = max(year),
    n_years         = n_distinct(year),
    n_hex           = n_distinct(hex_id),
    .groups = "drop"
  ) %>%
  arrange(first_anim_year)

print(
  survey_diagnostics,
  n = nrow(survey_diagnostics)
)

late_surveys <-
  survey_diagnostics %>%
  filter(
    first_anim_year >= max(animation_ready_data$year) - 2
  )

if(nrow(late_surveys) > 0) {
  
  cat(
    "\nWARNING: Survey units only appearing near final years:\n"
  )
  
  print(
    late_surveys,
    n = nrow(late_surveys)
  )
  
}

n_invalid <- sum(
  !st_is_valid(animation_ready_data)
)

cat(
  "\nInvalid geometries:",
  n_invalid,
  "\n"
)

cat(
  "Missing geometry matches:",
  nrow(missing_hexes),
  "\n"
)

if(nrow(missing_hexes) > 0) {
  
  print(
    missing_hexes %>%
      count(survey_unit, sort = TRUE),
    n = Inf
  )
  
}

final_year <- max(animation_ready_data$year)

cat(
  "\nFinal animation year:",
  final_year,
  "\n"
)

final_surveys <-
  animation_ready_data %>%
  st_drop_geometry() %>%
  filter(year == final_year) %>%
  count(survey_unit, sort = TRUE)


cat("\nSurvey units present in final frame:\n")

'
animation_ready_data %>%
  st_drop_geometry() %>%
  filter(year == final_year) %>%
  count(survey_unit, sort = TRUE) %>%
  print(n = Inf)

cat("\n========================================\n")
cat("END DIAGNOSTICS\n")
cat("========================================\n\n")
'


############################################################
# 5. WORLD MAP
############################################################

world_map <- ne_countries(
  scale = "medium",
  returnclass = "sf"
) %>%
  st_make_valid()

bbox <- st_bbox(animation_ready_data)

lon_lims <- c(-190, 40)

lat_lims <- c(
  bbox["ymin"] - 5,
  bbox["ymax"] + 5
)

############################################################
# 6. FIXED TEMPERATURE SCALE
############################################################

animation_ready_data <- animation_ready_data %>%
  filter(
    is.finite(hex_temperature)
  )

temp_limits <- range(c(-5, 40),
  na.rm = TRUE
)

print(temp_limits)

############################################################
# 7. CREATE FRAMES
############################################################

dir.create(
  "temperature_na_centered_frames",
  showWarnings = FALSE
)

all_years <- sort(
  unique(animation_ready_data$year)
)

for (yr in all_years) {
  
  yr_data <-
    animation_ready_data %>%
    filter(year == yr)
  
  p <- ggplot() +
    
    geom_sf(
      data = world_map,
      fill = "gray95",
      color = "gray80",
      linewidth = 0.1
    ) +
    
    geom_sf(
      data = yr_data,
      aes(fill = hex_temperature),
      color = NA
    ) +
    
    scale_fill_viridis_c(
      option = "turbo",
      limits = temp_limits,
      oob = scales::squish,
      na.value = "transparent",
      name = "Bottom\nTemp (°C)"
    ) +
    
    coord_sf(
      xlim = lon_lims,
      ylim = lat_lims,
      expand = FALSE,
      crs = 4326
    ) +
    
    theme_bw() +
    
    labs(
      title = paste("Bottom Temperature:", yr)
    ) +
    
    theme(
      plot.title = element_text(
        size = 28,
        face = "bold"
      ),
      axis.text = element_text(size = 12),
      legend.key.height = grid::unit(2, "cm"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    )
  
  ggsave(
    filename = sprintf(
      "temperature_na_centered_frames/frame_%04d.png",
      yr
    ),
    plot = p,
    width = 14,
    height = 9,
    dpi = 300
  )
}

############################################################
# 8. CREATE GIF
############################################################

frame_files <- list.files(
  "temperature_na_centered_frames",
  pattern = "\\.png$",
  full.names = TRUE
)

frame_files <- sort(frame_files)

img_sequence <- image_read(frame_files)

animation <- image_animate(
  img_sequence,
  fps = 1,
  dispose = "background"
)

image_write(
  animation,
  "Hex_Bottom_Temperature_NA_Centered.gif"
)


##_________________________________##
### Coverage estimates for pre-Q1 ###
##_________________________________##
print(range(Combined_FISH_optimized_spatially_valid$year))

haul_level <- Combined_FISH_optimized_spatially_valid %>%
  distinct(survey_unit, year, haul_id, area_swept, haul_dur)


area_swept_summary <- haul_level %>%
  filter(!is.na(area_swept)) %>%
  group_by(survey_unit, year) %>%
  summarise(
    n_hauls_with_area = n(),
    total_area_swept = sum(area_swept, na.rm = TRUE),
    mean_area_swept = mean(area_swept, na.rm = TRUE),
    median_area_swept = median(area_swept, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(survey_unit, year)

# Total across all survey units and years
total_area_swept_all <- sum(
  area_swept_summary$total_area_swept,
  na.rm = TRUE
)

missing_area_swept_hauls <- haul_level %>%
  filter(is.na(area_swept))

survey_effort_completeness <- haul_level %>%
  group_by(survey_unit) %>%
  summarise(
    total_unique_hauls = n(),
    
    hauls_with_area_swept =
      sum(!is.na(area_swept)),
    
    pct_area_swept_available =
      100 * hauls_with_area_swept / total_unique_hauls,
    
    hauls_with_haul_dur =
      sum(!is.na(haul_dur)),
    
    pct_haul_dur_available =
      100 * hauls_with_haul_dur / total_unique_hauls,
    
    hauls_with_both =
      sum(!is.na(area_swept) & !is.na(haul_dur)),
    
    pct_both_available =
      100 * hauls_with_both / total_unique_hauls,
    
    .groups = "drop"
  ) %>%
  arrange(desc(pct_area_swept_available))

survey_year_completeness <- haul_level %>%
  group_by(survey_unit, year) %>%
  summarise(
    total_unique_hauls = n(),
    pct_area_swept_available =
      100 * mean(!is.na(area_swept)),
    pct_haul_dur_available =
      100 * mean(!is.na(haul_dur)),
    pct_both_available =
      100 * mean(!is.na(area_swept) & !is.na(haul_dur)),
    .groups = "drop"
  )

area_swept_summary
survey_effort_completeness
missing_area_swept_hauls

cat(
  "\nEstimated total area swept across all survey-unit/year combinations:",
  format(round(total_area_swept_all, 2), big.mark = ","),
  "\n"
)

# try to estimate total area of all spatial hexagons
library(units)

# Calculate area (km2) for each hexagon
hex_area_summary <- all_hex_geometries %>%
  st_transform(6933) %>%  # equal-area projection
  mutate(area_km2 = set_units(st_area(.), km^2)) %>%
  st_drop_geometry() %>%
  distinct(survey_unit, hex_id, .keep_all = TRUE) %>%
  group_by(survey_unit) %>%
  summarise(
    n_hexes = n(),
    total_area_km2 = sum(area_km2),
    .groups = "drop"
  ) %>%
  mutate(total_area_km2 = as.numeric(total_area_km2))

hex_area_summary

all_hex_geometries %>% #shows that hexes are being counted only once per survey unit so its total sampling area of each survey unit
  st_drop_geometry() %>%
  count(survey_unit, hex_id) %>%
  filter(n > 1)

grand_total_area_km2 <- all_hex_geometries %>%
  st_transform(6933) %>%
  mutate(area_km2 = set_units(st_area(.), km^2)) %>%
  st_drop_geometry() %>%
  distinct(survey_unit, hex_id, .keep_all = TRUE) %>%
  summarise(total_area_km2 = sum(area_km2)) %>%
  pull(total_area_km2) %>%
  as.numeric()

grand_total_area_km2

# generate percentages of coverage for each secondary dataset
library(dplyr)
library(tidyr)

#------------------------------------------------------------
# One record per Aphia ID
#------------------------------------------------------------
aphia_traits <- Combined_FISH_optimized_spatially_valid %>%
  select(
    aphia_id,
    accepted_name,
    rank,
    Troph,
    estim_thermal_mean,
    trait_resolution
  ) %>%
  filter(!is.na(aphia_id)) %>%
  distinct()

#------------------------------------------------------------
# Species-level summary
#------------------------------------------------------------
n_total_aphia <- n_distinct(aphia_traits$aphia_id)

n_species_aphia <- aphia_traits %>%
  filter(tolower(rank) == "species") %>%
  summarise(n = n_distinct(aphia_id)) %>%
  pull(n)

rank_summary <- tibble(
  total_unique_aphia_ids = n_total_aphia,
  species_level_ids = n_species_aphia,
  non_species_ids = n_total_aphia - n_species_aphia,
  pct_species_level =
    round(100 * n_species_aphia / n_total_aphia, 2),
  pct_non_species =
    round(100 * (n_total_aphia - n_species_aphia) / n_total_aphia, 2)
)

#------------------------------------------------------------
# Trait coverage
#------------------------------------------------------------
trait_coverage <- tibble(
  variable = c(
    "Troph",
    "estim_thermal_mean",
    "trait_resolution"
  )
) %>%
  rowwise() %>%
  mutate(
    
    # coverage across all Aphia IDs
    n_aphia_with_trait =
      n_distinct(
        aphia_traits$aphia_id[
          !is.na(aphia_traits[[variable]])
        ]
      ),
    
    pct_all_aphia =
      round(
        100 * n_aphia_with_trait / n_total_aphia,
        2
      ),
    
    # coverage among species-level Aphia IDs only
    n_species_with_trait =
      n_distinct(
        aphia_traits$aphia_id[
          tolower(aphia_traits$rank) == "species" &
            !is.na(aphia_traits[[variable]])
        ]
      ),
    
    pct_species_aphia =
      round(
        100 * n_species_with_trait / n_species_aphia,
        2
      )
  ) %>%
  ungroup() %>%
  select(
    variable,
    n_aphia_with_trait,
    pct_all_aphia,
    n_species_with_trait,
    pct_species_aphia
  )

#------------------------------------------------------------
# Results
#------------------------------------------------------------
rank_summary
trait_coverage

###_________________________________________________####
##### Q1: Estimate rates of thermophilization ##########
###_________________________________________________####
# begin with global and survey unit scale hex level linear regressions of MTC change

global_spatially_valid_hex_stats <- left_join(
  global_spatially_valid_hex_stats,
  global_MTC_by_hex,
  by = c("survey_unit", "year", "hex_id", "lme_name")
)
############################################################
#### THERMOPHILIZATION MODEL FRAMEWORK #####################
############################################################

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(mgcv)
library(gratia)
library(AICcmodavg)
library(sf)
library(ggplot2)
library(gt)
library(lwgeom)
library(rnaturalearth)
library(rnaturalearthdata)

############################################################
# COMMON PALETTE
############################################################

trend_palette <- c(
  "Significant Increase" = "#F06292",
  "Significant Decrease" = "#FFD54F",
  "No Significant Trend" = "#64B5F6"
)

############################################################
# COMMON TREND CLASSIFIER
############################################################

classify_trend <- function(rate, p_value){
  
  case_when(
    p_value < 0.05 & rate > 0 ~
      "Significant Increase",
    
    p_value < 0.05 & rate < 0 ~
      "Significant Decrease",
    
    TRUE ~
      "No Significant Trend"
  )
}

############################################################
# STANDARDIZED SURVEY RESULT OBJECT
############################################################

format_survey_output <- function(
    survey_unit,
    method,
    rate,
    conf_low,
    conf_high,
    p_value,
    model_type = NA,
    r_squared = NA,
    adj_r2 = NA,
    deviance_explained = NA,
    edf = NA,
    n_obs = NA,
    n_hex = NA,
    year_range = NA
){
  
  tibble(
    
    survey_unit = survey_unit,
    
    method = method,
    
    annual_rate = rate,
    
    conf_low = conf_low,
    
    conf_high = conf_high,
    
    p_value = p_value,
    
    trend_cat =
      classify_trend(
        rate,
        p_value
      ),
    
    model_type = model_type,
    
    r_squared = r_squared,
    
    adj_r2 = adj_r2,
    
    deviance_explained =
      deviance_explained,
    
    edf = edf,
    
    n_obs = n_obs,
    
    n_hex = n_hex,
    
    year_range = year_range
  )
}

############################################################
# STANDARDIZED HEX OUTPUT OBJECT
############################################################

format_hex_output <- function(
    survey_unit,
    hex_id,
    method,
    rate,
    conf_low,
    conf_high,
    p_value
){
  
  tibble(
    
    survey_unit,
    
    hex_id,
    
    method,
    
    annual_rate = rate,
    
    conf_low,
    
    conf_high,
    
    p_value,
    
    trend_class =
      classify_trend(
        rate,
        p_value
      )
  )
}

############################################################
#### SECTION 2: MODEL FITTING FUNCTIONS (CORRECTED) ########
############################################################

############################################################
# 2A. WEIGHTED LINEAR TREND
############################################################

fit_linear_trend <- function(dat){
  
  dat <- dat %>%
    dplyr::filter(
      !is.na(MTC_hex_avg),
      !is.na(MTC_hex_variance),
      MTC_hex_variance > 0,
      !is.na(year)
    ) %>%
    dplyr::arrange(year)
  
  if(dplyr::n_distinct(dat$year) < 4){
    return(NULL)
  }
  
  model <- lm(
    MTC_hex_avg ~ year,
    data = dat,
    weights = 1 / MTC_hex_variance
  )
  
  td <- broom::tidy(
    model,
    conf.int = TRUE
  )
  
  gl <- broom::glance(model)
  
  yr <- td %>%
    dplyr::filter(term == "year")
  
  trend <- tibble::tibble(
    
    annual_rate = yr$estimate,
    
    conf_low = yr$conf.low,
    
    conf_high = yr$conf.high,
    
    std_error = yr$std.error,
    
    p_value = yr$p.value,
    
    model_type = "Weighted Linear",
    
    r_squared = gl$r.squared,
    
    adj_r2 = gl$adj.r.squared,
    
    deviance_explained = NA_real_,
    
    edf = 1,
    
    n_obs = nrow(dat),
    
    n_hex =
      if ("hex_id" %in% names(dat)) {
        dplyr::n_distinct(dat$hex_id)
      } else {
        1L
      },
    
    year_range = paste0(
      min(dat$year),
      "-",
      max(dat$year)
    )
    
  )
  
  return(
    list(
      model = model,
      trend = trend
    )
  )
}

############################################################
# 2B. HARMONIC MODEL WITH AICc SELECTION
############################################################

fit_harmonic_trend <- function(dat){
  
  dat <- dat %>%
    dplyr::filter(
      !is.na(MTC_hex_avg),
      !is.na(MTC_hex_variance),
      MTC_hex_variance > 0,
      !is.na(year)
    ) %>%
    dplyr::arrange(year)
  
  if(dplyr::n_distinct(dat$year) < 8){
    return(NULL)
  }
  
  ##########################################################
  # HELPER FOR AICc
  ##########################################################
  
  calc_AICc <- function(model){
    
    k <- length(stats::coef(model))
    
    n <- stats::nobs(model)
    
    aic <- stats::AIC(model)
    
    if((n - k - 1) <= 0){
      return(Inf)
    }
    
    aic +
      (
        2 * k * (k + 1)
      ) /
      (n - k - 1)
  }
  
  ##########################################################
  # PERIODS
  ##########################################################
  
  duration <- max(dat$year) -
    min(dat$year)
  
  periods <- c(
    5,
    10,
    20,
    40
  )
  
  periods <- periods[
    periods <= duration
  ]
  
  t <- dat$year -
    min(dat$year)
  
  ##########################################################
  # BUILD HARMONIC TERMS
  ##########################################################
  
  if(length(periods) > 0){
    
    for(p in periods){
      
      dat[[paste0("sin_", p)]] <-
        sin(
          2 * pi * t / p
        )
      
      dat[[paste0("cos_", p)]] <-
        cos(
          2 * pi * t / p
        )
    }
  }
  
  ##########################################################
  # MODEL LIST
  ##########################################################
  
  wts <- 1 /
    dat$MTC_hex_variance
  
  model_list <- list()
  
  ##########################################################
  # BASE LINEAR
  ##########################################################
  
  model_list[["linear"]] <-
    lm(
      MTC_hex_avg ~ year,
      data = dat,
      weights = wts
    )
  
  ##########################################################
  # SINGLE HARMONICS
  ##########################################################
  
  if(length(periods) > 0){
    
    for(p in periods){
      
      ff <- as.formula(
        paste0(
          "MTC_hex_avg ~ year + ",
          "sin_", p,
          " + cos_", p
        )
      )
      
      model_list[[paste0(
        "harmonic_",
        p
      )]] <-
        lm(
          ff,
          data = dat,
          weights = wts
        )
    }
  }
  
  ##########################################################
  # DOUBLE HARMONICS
  ##########################################################
  
  if(length(periods) >= 2){
    
    combos <- combn(
      periods,
      2,
      simplify = FALSE
    )
    
    for(cc in combos){
      
      p1 <- cc[1]
      p2 <- cc[2]
      
      ff <- as.formula(
        paste0(
          "MTC_hex_avg ~ year + ",
          "sin_", p1,
          " + cos_", p1,
          " + ",
          "sin_", p2,
          " + cos_", p2
        )
      )
      
      model_list[[paste0(
        "harmonic_",
        p1,
        "_",
        p2
      )]] <-
        lm(
          ff,
          data = dat,
          weights = wts
        )
    }
  }
  
  ##########################################################
  # AICc MODEL SELECTION
  ##########################################################
  
  selection <- tibble::tibble(
    model_name = names(model_list),
    AICc = sapply(
      model_list,
      calc_AICc
    )
  ) %>%
    dplyr::arrange(AICc) %>%
    dplyr::mutate(
      delta_AICc =
        AICc - min(AICc)
    )
  
  best_name <-
    selection$model_name[1]
  
  best_model <-
    model_list[[best_name]]
  
  ##########################################################
  # EXTRACT YEAR COEFFICIENT
  ##########################################################
  
  td <- broom::tidy(
    best_model,
    conf.int = TRUE
  )
  
  gl <- broom::glance(
    best_model
  )
  
  yr <- td %>%
    dplyr::filter(
      term == "year"
    )
  
  trend <- tibble::tibble(
    
    annual_rate = yr$estimate,
    
    conf_low = yr$conf.low,
    
    conf_high = yr$conf.high,
    
    std_error = yr$std.error,
    
    p_value = yr$p.value,
    
    model_type = best_name,
    
    r_squared = gl$r.squared,
    
    adj_r2 = gl$adj.r.squared,
    
    deviance_explained = NA_real_,
    
    edf = NA_real_,
    
    n_obs = nrow(dat),
    
    n_hex =
      if ("hex_id" %in% names(dat)) {
        dplyr::n_distinct(dat$hex_id)
      } else {
        1L
      },
    
    year_range =
      paste0(
        min(dat$year),
        "-",
        max(dat$year)
      )
  )
  
  return(
    list(
      model = best_model,
      selection = selection,
      trend = trend
    )
  )
}

############################################################
# 2C. WGAMM MODEL
# RATE = MEAN FIRST DERIVATIVE
############################################################

fit_wgamm_trend <- function(dat){
  
  dat <- dat %>%
    dplyr::filter(
      !is.na(MTC_hex_avg),
      !is.na(MTC_hex_variance),
      MTC_hex_variance > 0,
      !is.na(year)
    ) %>%
    dplyr::arrange(year)
  
  if(dplyr::n_distinct(dat$year) < 10){
    return(NULL)
  }
  
  ##########################################################
  # FIT GAM
  ##########################################################
  
  model <- mgcv::gam(
    
    MTC_hex_avg ~
      
      s(
        year,
        k = min(
          5,
          dplyr::n_distinct(dat$year) - 1
        )
      ),
    
    data = dat,
    
    weights =
      1 /
      MTC_hex_variance,
    
    method = "REML"
  )
  
  ##########################################################
  # DERIVATIVES
  ##########################################################
  
  derivs <- gratia::derivatives(
    model,
    term = "s(year)",
    interval = "confidence"
  )
  
  annual_rate <-
    mean(
      derivs$.derivative,
      na.rm = TRUE
    )
  
  conf_low <-
    mean(
      derivs$.lower_ci,
      na.rm = TRUE
    )
  
  conf_high <-
    mean(
      derivs$.upper_ci,
      na.rm = TRUE
    )
  
  std_error <-
    (
      conf_high -
        conf_low
    ) /
    (2 * 1.96)
  
  ##########################################################
  # MODEL SUMMARY
  ##########################################################
  
  sm <- summary(model)
  
  trend <- tibble::tibble(
    
    annual_rate =
      annual_rate,
    
    conf_low =
      conf_low,
    
    conf_high =
      conf_high,
    
    std_error =
      std_error,
    
    p_value =
      sm$s.table[
        1,
        "p-value"
      ],
    
    model_type =
      "wGAMM",
    
    r_squared =
      sm$r.sq,
    
    adj_r2 =
      NA_real_,
    
    deviance_explained =
      sm$dev.expl,
    
    edf =
      sm$s.table[
        1,
        "edf"
      ],
    
    n_obs =
      nrow(dat),
    
    n_hex =
      if ("hex_id" %in% names(dat)) {
        dplyr::n_distinct(dat$hex_id)
      } else {
        1L
      },
    
    year_range =
      paste0(
        min(dat$year),
        "-",
        max(dat$year)
      )
  )
  
  return(
    list(
      model = model,
      derivatives = derivs,
      trend = trend
    )
  )
}

############################################################
#### SECTION 3: SURVEY-LEVEL MODELS ########################
############################################################

############################################################
# HELPER FUNCTION
############################################################

run_survey_models <- function(
    data,
    fit_function,
    method_name
){
  
  data %>%
    group_by(survey_unit) %>%
    nest() %>%
    mutate(
      fit = map(
        data,
        fit_function
      )
    ) %>%
    filter(
      !map_lgl(
        fit,
        is.null
      )
    ) %>%
    mutate(
      
      model =
        map(
          fit,
          "model"
        ),
      
      trend =
        map(
          fit,
          "trend"
        )
      
    ) %>%
    select(
      survey_unit,
      model,
      trend
    ) %>%
    unnest(trend) %>%
    mutate(
      method = method_name,
      
      trend_cat =
        classify_trend(
          annual_rate,
          p_value
        )
    )
}

############################################################
# LINEAR SURVEY MODELS
############################################################

linear_survey_results <-
  run_survey_models(
    data =
      global_spatially_valid_hex_stats,
    
    fit_function =
      fit_linear_trend,
    
    method_name =
      "Linear"
  )

############################################################
# HARMONIC SURVEY MODELS
############################################################

harmonic_survey_results <-
  global_spatially_valid_hex_stats %>%
  group_by(survey_unit) %>%
  nest() %>%
  mutate(
    fit =
      map(
        data,
        fit_harmonic_trend
      )
  ) %>%
  filter(
    !map_lgl(
      fit,
      is.null
    )
  ) %>%
  mutate(
    
    model =
      map(
        fit,
        "model"
      ),
    
    trend =
      map(
        fit,
        "trend"
      ),
    
    model_selection =
      map(
        fit,
        "selection"
      )
    
  ) %>%
  select(
    survey_unit,
    model,
    trend
  ) %>%
  unnest(trend) %>%
  mutate(
    
    method = "Harmonic",
    
    trend_cat =
      classify_trend(
        annual_rate,
        p_value
      )
    
  )

############################################################
# WGAMM SURVEY MODELS
############################################################

wgamm_survey_results <-
  run_survey_models(
    data =
      global_spatially_valid_hex_stats,
    
    fit_function =
      fit_wgamm_trend,
    
    method_name =
      "wGAMM"
  )

############################################################
# COMBINED SURVEY RESULTS
############################################################

survey_results_all <-
  bind_rows(
    
    linear_survey_results,
    
    harmonic_survey_results,
    
    wgamm_survey_results
    
  ) %>%
  arrange(
    survey_unit,
    method
  )

############################################################
# VERIFY STRUCTURE
############################################################

glimpse(
  survey_results_all
)

############################################################
# MASTER SURVEY TABLE
############################################################

survey_model_comparison_table <-
  survey_results_all %>%
  mutate(
    
    estimate_ci =
      sprintf(
        "%.4f (%.4f, %.4f)",
        annual_rate,
        conf_low,
        conf_high
      )
    
  ) %>%
  select(
    
    survey_unit,
    
    method,
    
    annual_rate,
    
    estimate_ci,
    
    p_value,
    
    trend_cat,
    
    model_type,
    
    r_squared,
    
    adj_r2,
    
    deviance_explained,
    
    edf,
    
    n_obs,
    
    n_hex,
    
    year_range
    
  ) %>%
  arrange(
    survey_unit,
    method
  )

############################################################
# WIDE COMPARISON TABLE (MANUSCRIPT READY)
############################################################

survey_model_wide <-
  survey_results_all %>%
  mutate(
    
    estimate_ci =
      sprintf(
        "[%.4f, %.4f]",
        conf_low,
        conf_high
      )
    
  ) %>%
  select(
    
    survey_unit,
    
    method,
    
    annual_rate,
    
    estimate_ci,
    
    p_value,
    
    trend_cat,
    
    model_type,
    
    r_squared,
    
    adj_r2,
    
    deviance_explained,
    
    edf
    
  ) %>%
  pivot_wider(
    
    names_from = method,
    
    values_from = c(
      annual_rate,
      estimate_ci,
      p_value,
      trend_cat,
      model_type,
      r_squared,
      adj_r2,
      deviance_explained,
      edf
    )
    
  )

############################################################
# VIEW TABLES
############################################################

print(
  survey_model_comparison_table
)

print(
  survey_model_wide
)

############################################################
# OPTIONAL GT TABLE
############################################################

survey_model_gt <-
  survey_model_comparison_table %>%
  gt() %>%
  tab_header(
    title =
      md(
        "**Survey Unit Thermophilization Comparison Across Models**"
      )
  ) %>%
  fmt_number(
    columns =
      c(
        annual_rate,
        r_squared,
        adj_r2,
        deviance_explained,
        edf
      ),
    decimals = 4
  )

survey_model_gt

############################################################
#### SECTION 4: HEX-LEVEL MODELS ###########################
############################################################

############################################################
# HELPER FUNCTION
############################################################

run_hex_models <- function(
    data,
    fit_function,
    method_name,
    min_years
){
  
  data %>%
    group_by(
      survey_unit,
      hex_id
    ) %>%
    filter(
      n_distinct(year) >= min_years
    ) %>%
    nest() %>%
    mutate(
      fit =
        map(
          data,
          fit_function
        )
    ) %>%
    filter(
      !map_lgl(
        fit,
        is.null
      )
    ) %>%
    mutate(
      trend =
        map(
          fit,
          "trend"
        )
    ) %>%
    select(
      survey_unit,
      hex_id,
      trend
    ) %>%
    unnest(trend) %>%
    mutate(
      
      method =
        method_name,
      
      trend_class =
        classify_trend(
          annual_rate,
          p_value
        )
      
    )
}

############################################################
# LINEAR HEX MODELS
############################################################

linear_hex_results <-
  run_hex_models(
    data =
      global_spatially_valid_hex_stats,
    
    fit_function =
      fit_linear_trend,
    
    method_name =
      "Linear",
    
    min_years = 4
  )

############################################################
# HARMONIC HEX MODELS
############################################################

harmonic_hex_results <-
  global_spatially_valid_hex_stats %>%
  group_by(
    survey_unit,
    hex_id
  ) %>%
  filter(
    n_distinct(year) >= 8
  ) %>%
  nest() %>%
  mutate(
    fit =
      map(
        data,
        fit_harmonic_trend
      )
  ) %>%
  filter(
    !map_lgl(
      fit,
      is.null
    )
  ) %>%
  mutate(
    trend =
      map(
        fit,
        "trend"
      )
  ) %>%
  select(
    survey_unit,
    hex_id,
    trend
  ) %>%
  unnest(trend) %>%
  mutate(
    
    method =
      "Harmonic",
    
    trend_class =
      classify_trend(
        annual_rate,
        p_value
      )
    
  )

############################################################
# WGAMM HEX MODELS
############################################################

wgamm_hex_results <-
  run_hex_models(
    data =
      global_spatially_valid_hex_stats,
    
    fit_function =
      fit_wgamm_trend,
    
    method_name =
      "wGAMM",
    
    min_years = 10
  )

############################################################
# COMBINED HEX RESULTS
############################################################

hex_results_all <-
  bind_rows(
    
    linear_hex_results,
    
    harmonic_hex_results,
    
    wgamm_hex_results
    
  ) %>%
  mutate(
    hex_id =
      as.character(hex_id)
  )

############################################################
# VERIFY OUTPUT
############################################################

glimpse(
  hex_results_all
)

############################################################
# SUMMARY BY METHOD
############################################################

hex_method_summary <-
  hex_results_all %>%
  group_by(method) %>%
  summarise(
    
    total_hexes =
      n(),
    
    significant_hexes =
      sum(
        p_value < 0.05,
        na.rm = TRUE
      ),
    
    increasing_hexes =
      sum(
        p_value < 0.05 &
          annual_rate > 0,
        na.rm = TRUE
      ),
    
    decreasing_hexes =
      sum(
        p_value < 0.05 &
          annual_rate < 0,
        na.rm = TRUE
      ),
    
    mean_rate =
      mean(
        annual_rate,
        na.rm = TRUE
      ),
    
    median_rate =
      median(
        annual_rate,
        na.rm = TRUE
      ),
    
    min_rate =
      min(
        annual_rate,
        na.rm = TRUE
      ),
    
    max_rate =
      max(
        annual_rate,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

print(
  hex_method_summary
)

############################################################
# GLOBAL COLOR SCALE FOR ALL HEX MAPS
############################################################

global_max_hex_rate <-
  max(
    abs(
      hex_results_all$annual_rate
    ),
    na.rm = TRUE
  )

print(
  global_max_hex_rate
)

############################################################
# JOIN TO GEOMETRIES
############################################################

all_hex_geometries_join <-
  all_hex_geometries %>%
  mutate(
    hex_id =
      as.character(hex_id)
  )

hex_map_data <-
  all_hex_geometries_join %>%
  left_join(
    hex_results_all,
    by = c(
      "survey_unit",
      "hex_id"
    )
  ) %>%
  mutate(
    
    slope_sig =
      ifelse(
        p_value < 0.05,
        annual_rate,
        NA_real_
      )
    
  )

############################################################
# WRAP ACROSS DATELINE
############################################################

hex_map_data <-
  hex_map_data %>%
  st_make_valid() %>%
  st_wrap_dateline(
    options = c(
      "WRAPDATELINE=YES",
      "DATELINEOFFSET=180"
    )
  ) %>%
  st_cast(
    "MULTIPOLYGON"
  )

############################################################
# VERIFY FACET MAP OBJECT
############################################################

hex_map_data

############################################################
#### SECTION 5: FACETED LOLLIPOP COMPARISON ###############
############################################################

############################################################
# COMMON SURVEY ORDER ACROSS ALL METHODS
############################################################

survey_order <-
  
  survey_results_all %>%
  
  group_by(
    survey_unit
  ) %>%
  
  summarise(
    
    mean_rate =
      mean(
        annual_rate,
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(
    mean_rate
  ) %>%
  
  pull(
    survey_unit
  )

############################################################
# PREPARE LOLLIPOP DATA
############################################################

lollipop_data <-
  survey_results_all %>%
  mutate(
    
    survey_unit =
      factor(
        survey_unit,
        levels = survey_order
      ),
    
    model_shape =
      case_when(
        
        method == "Linear"
        ~ "Linear",
        
        method == "Harmonic" &
          grepl(
            "^harmonic_[0-9]+$",
            model_type
          )
        ~ "Single Harmonic",
        
        method == "Harmonic" &
          grepl(
            "^harmonic_[0-9]+_[0-9]+$",
            model_type
          )
        ~ "Double Harmonic",
        
        method == "Harmonic" &
          model_type == "linear"
        ~ "Linear",
        
        method == "wGAMM" &
          edf < 1.5
        ~ "Near Linear",
        
        method == "wGAMM" &
          edf < 3
        ~ "Moderately Nonlinear",
        
        method == "wGAMM"
        ~ "Strongly Nonlinear",
        
        TRUE
        ~ "Other"
      )
  )

############################################################
# HARMONIZED SURVEY MODEL TABLE
############################################################

survey_model_comparison_table <-
  lollipop_data %>%
  mutate(
    
    estimate_ci =
      sprintf(
        "[%.4f, %.4f]",
        conf_low,
        conf_high
      )
    
  ) %>%
  select(
    
    survey_unit,
    
    method,
    
    annual_rate,
    
    estimate_ci,
    
    p_value,
    
    trend_cat,
    
    model_shape,
    
    model_type,
    
    r_squared,
    
    adj_r2,
    
    deviance_explained,
    
    edf,
    
    n_obs,
    
    n_hex,
    
    year_range
    
  ) %>%
  arrange(
    method,
    desc(annual_rate)
  )

############################################################
# GT TABLE
############################################################

survey_model_gt <-
  survey_model_comparison_table %>%
  gt(
    groupname_col = "method"
  ) %>%
  
  cols_label(
    survey_unit = "Survey Unit",
    
    annual_rate =
      html("Thermophilization Rate (°C yr<sup>-1</sup>)"),
    
    estimate_ci =
      "95% Confidence Interval",
    
    p_value =
      "P-value",
    
    trend_cat =
      "Trend Classification",
    
    model_shape =
      "Model Type",
    
    model_type =
      "AICC Selected Model",
    
    r_squared =
      html("R<sup>2</sup>"),
    
    adj_r2 =
      html("Adjusted R<sup>2</sup>"),
    
    deviance_explained =
      "Deviance Explained",
    edf =
      "Effective Degrees of Freedom",
    
    n_obs =
      "Number of Hex-Years",
    
    n_hex =
      "Spatial Hexagons",
    
    year_range =
      "Survey Period"
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(
      columns = everything()
    )
  ) %>%
  
  tab_style(
    style = list(
      cell_text(
        weight = "bold",
        size = px(14)
      )
    ),
    locations = cells_row_groups()
  ) %>%
  
  tab_header(
    title = md(
      "**Thermophilization Rate Comparison across Models**"
    )
  ) %>%
  
  fmt_number(
    columns = c(
      annual_rate,
      p_value
    ),
    decimals = 4
  )
survey_model_gt # T2 of full model comparison!

############################################################
# COMMON X-AXIS LIMITS (IDENTICAL FOR ALL PANELS)
############################################################

max_abs_rate <- max(
  abs(
    c(
      lollipop_data$conf_low,
      lollipop_data$conf_high
    )
  ),
  na.rm = TRUE
)

common_xlim <- c(
  -0.25,0.25
)

lollipop_data <- lollipop_data %>%
  mutate(
    survey_unit = factor(
      survey_unit,
      levels = survey_order
    )
  )
############################################################
# ORDER SURVEYS BY HARMONIC THERMOPHILIZATION RATE
############################################################

harmonic_order <-
  lollipop_data %>%
  filter(method == "Harmonic") %>%
  arrange(desc(annual_rate)) %>%
  pull(survey_unit) %>%
  as.character() %>%
  unique()

survey_order <-
  c(
    harmonic_order,
    setdiff(
      unique(as.character(lollipop_data$survey_unit)),
      harmonic_order
    )
  )

lollipop_data <-
  lollipop_data %>%
  mutate(
    survey_unit = factor(
      as.character(survey_unit),
      levels = rev(survey_order)
    )
  )

############################################################
# COMMON THEME
############################################################

library(patchwork)

lollipop_theme <-
  theme_classic(base_size = 20) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 22,
      hjust = 0.5
    ),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position = "right"
  )

############################################################
# LINEAR PANEL
############################################################

p_linear <-
  lollipop_data %>%
  filter(method == "Linear") %>%
  ggplot(
    aes(
      x = annual_rate,
      y = survey_unit,
      colour = trend_cat
    )
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  
  geom_segment(
    aes(
      x = 0,
      xend = annual_rate,
      yend = survey_unit
    ),
    linewidth = 0.6,
    alpha = 0.6
  ) +
  
  geom_errorbar(
    aes(
      xmin = conf_low,
      xmax = conf_high
    ),
    width = 0,
    linewidth = 0.6
  ) +
  
  geom_point(
    aes(
      size = n_hex,
      shape = model_shape
    ),
    alpha = 0.95,
    stroke = 1
  ) +
  
  scale_colour_manual(
    values = trend_palette
  ) +
  
  scale_size_area(
    max_size = 10,
    name = "Hexes"
  ) +
  
  coord_cartesian(
    xlim = common_xlim
  ) +
  
  labs(
    title = "Linear",
    x = NULL,
    y = "Survey unit",
    colour = "Trend",
    shape = "Model Type"
  ) +
  
  lollipop_theme

############################################################
# HARMONIC PANEL
############################################################

p_harmonic <-
  lollipop_data %>%
  filter(method == "Harmonic") %>%
  ggplot(
    aes(
      x = annual_rate,
      y = survey_unit,
      colour = trend_cat
    )
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  
  geom_segment(
    aes(
      x = 0,
      xend = annual_rate,
      yend = survey_unit
    ),
    linewidth = 0.6,
    alpha = 0.6
  ) +
  
  geom_errorbar(
    aes(
      xmin = conf_low,
      xmax = conf_high
    ),
    width = 0,
    linewidth = 0.6
  ) +
  
  geom_point(
    aes(
      size = n_hex,
      shape = model_shape
    ),
    alpha = 0.95,
    stroke = 1
  ) +
  
  scale_colour_manual(
    values = trend_palette
  ) +
  
  scale_size_area(
    max_size = 10,
    name = "Hexes"
  ) +
  
  coord_cartesian(
    xlim = common_xlim
  ) +
  
  labs(
    title = "Harmonic",
    x = NULL,
    y = NULL,
    colour = "Trend",
    shape = "Model Type"
  ) +
  
  lollipop_theme +
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

############################################################
# WGAMM PANEL
############################################################

p_wgamm <-
  lollipop_data %>%
  filter(method == "wGAMM") %>%
  ggplot(
    aes(
      x = annual_rate,
      y = survey_unit,
      colour = trend_cat
    )
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  
  geom_segment(
    aes(
      x = 0,
      xend = annual_rate,
      yend = survey_unit
    ),
    linewidth = 0.6,
    alpha = 0.6
  ) +
  
  geom_errorbar(
    aes(
      xmin = conf_low,
      xmax = conf_high
    ),
    width = 0,
    linewidth = 0.6
  ) +
  
  geom_point(
    aes(
      size = n_hex,
      shape = model_shape
    ),
    alpha = 0.95,
    stroke = 1
  ) +
  
  scale_colour_manual(
    values = trend_palette
  ) +
  scale_size_area(
    max_size = 10,
    name = "Hexes"
  ) +
  
  coord_cartesian(
    xlim = common_xlim
  ) +
  
  labs(
    title = "wGAMM",
    x = NULL,
    y = NULL,
    colour = "Trend",
    shape = "Model Type"
  ) +
  
  lollipop_theme +
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

############################################################
# COMBINE PANELS WITH SHARED LEGEND
############################################################

faceted_lollipop_plot <-
  (
    p_linear |
      p_harmonic |
      p_wgamm
  ) +
  plot_layout(
    guides = "collect",
    widths = c(1.4, 1, 1)
  ) +
  plot_annotation(
    tag_levels = "A"
  ) +
  plot_annotation(
    caption = expression(
      paste(
        "Annual CTI Change (",
        degree,
        "C yr"^-1,
        ")"
      )
    ),
    theme = theme(
      plot.caption = element_text(
        size = 20,
        hjust = 0.5
      ),
      plot.tag = element_text(
        size = 24,
        face = "bold"
    )
  )
)
############################################################
# EXPORT
############################################################

ggsave(
  "Faceted_Thermophilization_Lollipop_Horizontal.png",
  faceted_lollipop_plot,
  width = 30,
  height = 15,
  dpi = 300
)

############################################################
# QUICK SUMMARY
############################################################

model_summary_counts <-
  survey_results_all %>%
  group_by(
    method,
    trend_cat
  ) %>%
  summarise(
    n_surveys = n(),
    .groups = "drop"
  )

print(
  model_summary_counts
)

############################################################
#### SECTION 6: FACETED HEX THERMOPHILIZATION MAPS #########
############################################################

library(ggpattern)
library(grid)

############################################################
# REMOVE NA METHOD PANEL + SET METHOD ORDER
############################################################

hex_map_data <-
  hex_map_data %>%
  filter(!is.na(method)) %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Linear",
        "Harmonic",
        "wGAMM"
      )
    ),
    
    significant = p_value < 0.05
  )

############################################################
# COMMON SURVEY ORDER
############################################################

hex_map_data <-
  hex_map_data %>%
  mutate(
    survey_unit =
      factor(
        survey_unit,
        levels = survey_order
      )
  )

############################################################
# COLOUR SCALE USING 95th PERCENTILE
############################################################

global_max_hex_rate <-
  quantile(
    abs(
      hex_map_data$annual_rate
    ),
    probs = 0.95,
    na.rm = TRUE
  )

print(global_max_hex_rate)
'
library(grid)
############################################################
# FACETED HEX MAP
############################################################

faceted_hex_map <-
  
  ggplot() +
  
  ##########################################################
# BASEMAP: NON-PROVIDER COUNTRIES
##########################################################

geom_sf(
  data =
    world_map_plot %>%
    filter(
      !(name %in% highlight_countries)
    ),
  fill = "grey85",
  colour = "black",
  linewidth = 0.20
) +
  
  ##########################################################
# BASEMAP: PROVIDER COUNTRIES
##########################################################

geom_sf(
  data =
    world_map_plot %>%
    filter(
      name %in% highlight_countries
    ),
  fill = "lightblue",
  colour = "black",
  linewidth = 0.20
) +
  
  ##########################################################
# NON-SIGNIFICANT HEXES
##########################################################

ggpattern::geom_sf_pattern(
  data =
    hex_map_data %>%
    filter(!significant),
  
  fill = "white",
  
  colour = NA,
  
  pattern = "stripe",
  
  pattern_fill = "grey35",
  
  pattern_colour = "grey35",
  
  pattern_angle = 45,
  
  pattern_density = 0.45,
  
  pattern_spacing = 0.025
) +
  
  ##########################################################
# SIGNIFICANT HEXES
##########################################################

geom_sf(
  data =
    hex_map_data %>%
    filter(significant),
  
  aes(
    fill = annual_rate
  ),
  
  colour = NA
) +
  
  ##########################################################
# FACETS
##########################################################

facet_wrap(
  ~method,
  ncol = 1
) +
  
  ##########################################################
# HIGH-CONTRAST COLOUR SCALE
##########################################################

scale_fill_gradient2(
  
  low = "#053061",
  
  mid = "#F7F7F7",
  
  high = "#67001F",
  
  midpoint = 0,
  
  limits = c(
    -global_max_hex_rate,
    global_max_hex_rate
  ),
  
  oob = scales::squish,
  
  name =
    expression(
      paste(
        "MTC Change (",
        degree,
        "C yr"^-1,
        ")"
      )
    )
) +
  
  ##########################################################
# COORDINATES
##########################################################

coord_sf(
  xlim = lon_lims,
  ylim = lat_lims,
  crs = 4326,
  expand = FALSE
) +
  
  scale_x_continuous(
    breaks =
      seq(
        -180,
        40,
        by = 20
      )
  ) +
  
  scale_y_continuous(
    breaks =
      seq(
        floor(lat_lims[1] / 5) * 5,
        ceiling(lat_lims[2] / 5) * 5,
        by = 5
      )
  ) +
  
  ##########################################################
# LABELS
##########################################################

labs(
  
  x = "Longitude",
  
  y = "Latitude",
  
  caption =
    "Cross-hatched hexagons indicate non-significant trends. Colour intensity reflects thermophilization rate."
) +
  
  ##########################################################
# THEME
##########################################################

theme_bw() +
  
  theme(
    
    strip.text =
      element_text(
        face = "bold",
        size = 22
      ),
    
    axis.title =
      element_text(
        face = "bold",
        size = 18
      ),
    
    axis.text =
      element_text(
        size = 14,
        colour = "black"
      ),
    
    legend.title =
      element_text(
        face = "bold",
        size = 18
      ),
    
    legend.text =
      element_text(
        size = 14
      ),
    
    plot.caption =
      element_text(
        size = 12,
        hjust = 0.5
      ),
    
    panel.grid.major =
      element_line(
        colour = "grey75",
        linewidth = 0.30
      ),
    
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "right"
  )

############################################################
# VIEW FIGURE
############################################################

print(
  faceted_hex_map
)

############################################################
# EXPORT
############################################################

ggsave(
  "Faceted_Hex_Thermophilization_Map.png",
  faceted_hex_map,
  width = 18,
  height = 24,
  dpi = 200
)
'
library(ggplot2)
library(sf)
library(dplyr)
library(grid)

############################################################
# CREATE THERMOPHILIZATION CLASSES
############################################################

positive_cut <- median(
  hex_map_data$annual_rate[
    hex_map_data$annual_rate > 0
  ],
  na.rm = TRUE
)

negative_cut <- median(
  hex_map_data$annual_rate[
    hex_map_data$annual_rate < 0
  ],
  na.rm = TRUE
)

hex_map_data <-
  
  hex_map_data %>%
  
  mutate(
    
    thermo_class = case_when(
      
      annual_rate <= negative_cut ~ "High negative",
      
      annual_rate < 0 ~ "Low negative",
      
      annual_rate < positive_cut ~ "Low positive",
      
      annual_rate >= positive_cut ~ "High positive",
      
      TRUE ~ NA_character_
      
    )
    
  )

############################################################
# FACETED HEX MAP
############################################################

faceted_hex_map <-
  
  ggplot() +
  
  ##########################################################
# BASEMAP: NON-PROVIDER COUNTRIES
##########################################################

geom_sf(
  
  data =
    world_map_plot %>%
    filter(
      !(name %in% highlight_countries)
    ),
  
  fill = "grey85",
  
  colour = "black",
  
  linewidth = 0.20
  
) +
  
  ##########################################################
# BASEMAP: PROVIDER COUNTRIES
##########################################################

geom_sf(
  
  data =
    world_map_plot %>%
    filter(
      name %in% highlight_countries
    ),
  
  fill = "lightblue",
  
  colour = "black",
  
  linewidth = 0.20
  
) +
  
  ##########################################################
# NON-SIGNIFICANT HEXES
##########################################################

geom_sf(
  
  data =
    hex_map_data %>%
    filter(!significant),
  
  fill = NA,
  
  colour = "grey85",
  
  linewidth = 0.15
  
) +
  
  ##########################################################
# SIGNIFICANT HEXES
##########################################################

geom_sf(
  
  data =
    hex_map_data %>%
    filter(significant),
  
  aes(
    fill = thermo_class
  ),
  
  colour = NA
  
) +
  
  ##########################################################
# THERMOPHILIZATION COLOUR SCALE
##########################################################

scale_fill_manual(
  
  values = c(
    
    "High negative" = "#08519C",
    
    "Low negative" = "#CBC3E3",
    
    "Low positive" = "#FCAE91",
    
    "High positive" = "#CB181D"
    
  ),
  
  breaks = c(
    "High negative",
    "Low negative",
    "Low positive",
    "High positive"
  ),
  
  drop = FALSE,
  
  name = "Thermophilization\nRate"
  
) +
  
  ##########################################################
# FACETS
##########################################################

facet_wrap(
  ~method,
  ncol = 1
) +
  
  ##########################################################
# COORDINATES
##########################################################

coord_sf(
  xlim = lon_lims,
  ylim = lat_lims,
  crs = 4326,
  expand = FALSE
) +
  
  scale_x_continuous(
    breaks =
      seq(
        -180,
        40,
        by = 20
      )
  ) +
  
  scale_y_continuous(
    breaks =
      seq(
        floor(lat_lims[1] / 5) * 5,
        ceiling(lat_lims[2] / 5) * 5,
        by = 5
      )
  ) +
  
  ##########################################################
# LABELS
##########################################################

labs(
  
  x = "Longitude",
  
  y = "Latitude"
  
) +
  
  ##########################################################
# THEME
##########################################################

theme_bw() +
  
  theme(
    
    strip.text =
      element_text(
        face = "bold",
        size = 22
      ),
    
    axis.title =
      element_text(
        face = "bold",
        size = 18
      ),
    
    axis.text =
      element_text(
        size = 14,
        colour = "black"
      ),
    
    legend.title =
      element_text(
        face = "bold",
        size = 18
      ),
    
    legend.text =
      element_text(
        size = 14
      ),
    
    plot.caption =
      element_text(
        size = 12,
        hjust = 0.5
      ),
    
    panel.grid.major =
      element_line(
        colour = "grey75",
        linewidth = 0.30
      ),
    
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "right"
    
  )

############################################################
# VIEW FIGURE
############################################################

print(
  faceted_hex_map
)

############################################################
# EXPORT
############################################################

ggsave(
  filename = "Faceted_Hex_Thermophilization_Map.png",
  plot = faceted_hex_map,
  width = 18,
  height = 24,
  dpi = 300,
  bg = "white"
)
############################################################
# HEX SUMMARY TABLE
############################################################

hex_summary_table <-
  
  hex_results_all %>%
  
  group_by(
    method
  ) %>%
  
  summarise(
    
    total_hexes =
      n(),
    
    significant_hexes =
      sum(
        p_value < 0.05,
        na.rm = TRUE
      ),
    
    increasing_hexes =
      sum(
        p_value < 0.05 &
          annual_rate > 0,
        na.rm = TRUE
      ),
    
    decreasing_hexes =
      sum(
        p_value < 0.05 &
          annual_rate < 0,
        na.rm = TRUE
      ),
    
    percent_significant =
      100 *
      significant_hexes /
      total_hexes,
    percent_increasing = 
      100*increasing_hexes / total_hexes, 
    percent_decreasing = 
      100*decreasing_hexes / total_hexes, 
    percent_insignificant = 
      1-percent_significant, 
    mean_rate =
      mean(
        annual_rate,
        na.rm = TRUE
      ),
    
    median_rate =
      median(
        annual_rate,
        na.rm = TRUE
      ),
    
    sd_rate =
      sd(
        annual_rate,
        na.rm = TRUE
      ),
    
    min_rate =
      min(
        annual_rate,
        na.rm = TRUE
      ),
    
    max_rate =
      max(
        annual_rate,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

print(
  hex_summary_table
)

############################################################
# OPTIONAL GT TABLE
############################################################

hex_summary_gt <-
  
  hex_summary_table %>%
  
  gt() %>%
  
  tab_header(
    title =
      md(
        "**Hex-Level Thermophilization Summary by Model Type**"
      )
  ) %>%
  
  fmt_number(
    columns = c(
      percent_significant,
      mean_rate,
      median_rate,
      sd_rate,
      min_rate,
      max_rate
    ),
    decimals = 2
  )

hex_summary_gt

############################################################
#### SECTION 7: REGRESSION FIT VISUALIZATIONS (CORRECTED) ##
############################################################

############################################################
# SURVEY-LEVEL OBSERVED TIME SERIES
############################################################

annual_survey_mtc <-
  
  global_spatially_valid_hex_stats %>%
  
  filter(
    !is.na(MTC_hex_avg),
    !is.na(MTC_hex_variance),
    MTC_hex_variance > 0
  ) %>%
  
  mutate(
    weight = 1 / MTC_hex_variance
  ) %>%
  
  group_by(
    survey_unit,
    year
  ) %>%
  
  summarise(
    
    observed =
      weighted.mean(
        MTC_hex_avg,
        weight
      ),
    
    variance =
      1 / sum(weight),
    
    n_hex = n(),
    
    .groups = "drop"
    
  )

############################################################
# REFIT SURVEY MODELS FOR VISUALIZATION
############################################################

linear_models_plot <-
  
  global_spatially_valid_hex_stats %>%
  
  group_by(survey_unit) %>%
  
  nest() %>%
  
  mutate(
    fit = map(
      data,
      fit_linear_trend
    )
  ) %>%
  
  filter(
    !map_lgl(
      fit,
      is.null
    )
  ) %>%
  
  mutate(
    model = map(
      fit,
      "model"
    )
  )

harmonic_models_plot <-
  
  global_spatially_valid_hex_stats %>%
  
  group_by(survey_unit) %>%
  
  nest() %>%
  
  mutate(
    fit = map(
      data,
      fit_harmonic_trend
    )
  ) %>%
  
  filter(
    !map_lgl(
      fit,
      is.null
    )
  ) %>%
  
  mutate(
    model = map(
      fit,
      "model"
    )
  )

wgamm_models_plot <-
  
  global_spatially_valid_hex_stats %>%
  
  group_by(survey_unit) %>%
  
  nest() %>%
  
  mutate(
    fit = map(
      data,
      fit_wgamm_trend
    )
  ) %>%
  
  filter(
    !map_lgl(
      fit,
      is.null
    )
  ) %>%
  
  mutate(
    model = map(
      fit,
      "model"
    )
  )

############################################################
# LINEAR PREDICTIONS
############################################################

linear_prediction_data <-
  
  linear_models_plot %>%
  
  mutate(
    
    prediction_data =
      
      map2(
        model,
        survey_unit,
        
        function(mod, su){
          
          obs <-
            
            annual_survey_mtc %>%
            
            filter(
              survey_unit == su
            )
          
          pred <-
            
            predict(
              mod,
              newdata =
                data.frame(
                  year = obs$year
                ),
              interval = "confidence"
            )
          
          tibble(
            
            year = obs$year,
            
            observed =
              obs$observed,
            
            fitted =
              pred[,1],
            
            lower =
              pred[,2],
            
            upper =
              pred[,3]
            
          )
        }
      )
    
  ) %>%
  
  select(
    survey_unit,
    prediction_data
  ) %>%
  
  tidyr::unnest(
    prediction_data
  )

############################################################
# HARMONIC PREDICTIONS
############################################################

harmonic_prediction_data <-
  
  harmonic_models_plot %>%
  
  mutate(
    
    prediction_data =
      
      map2(
        model,
        survey_unit,
        
        function(mod, su){
          
          obs <-
            
            annual_survey_mtc %>%
            
            filter(
              survey_unit == su
            )
          
          pred_dat <- obs
          
          duration <-
            max(obs$year) -
            min(obs$year)
          
          tt <-
            obs$year -
            min(obs$year)
          
          periods <- c(
            5,
            10,
            20,
            40
          )
          
          periods <-
            periods[
              periods <= duration
            ]
          
          if(length(periods) > 0){
            
            for(p in periods){
              
              pred_dat[[paste0(
                "sin_",
                p
              )]] <-
                sin(
                  2*pi*tt/p
                )
              
              pred_dat[[paste0(
                "cos_",
                p
              )]] <-
                cos(
                  2*pi*tt/p
                )
            }
          }
          
          pred <-
            
            predict(
              mod,
              newdata =
                pred_dat,
              interval =
                "confidence"
            )
          
          tibble(
            
            year =
              obs$year,
            
            observed =
              obs$observed,
            
            fitted =
              pred[,1],
            
            lower =
              pred[,2],
            
            upper =
              pred[,3]
            
          )
        }
      )
    
  ) %>%
  
  select(
    survey_unit,
    prediction_data
  ) %>%
  
  tidyr::unnest(
    prediction_data
  )

############################################################
# WGAMM PREDICTIONS
############################################################

wgamm_prediction_data <-
  
  wgamm_models_plot %>%
  
  mutate(
    
    prediction_data =
      
      map2(
        model,
        survey_unit,
        
        function(mod, su){
          
          obs <-
            
            annual_survey_mtc %>%
            
            filter(
              survey_unit == su
            )
          
          pred <-
            
            predict(
              mod,
              newdata =
                data.frame(
                  year =
                    obs$year
                ),
              se.fit = TRUE
            )
          
          tibble(
            
            year =
              obs$year,
            
            observed =
              obs$observed,
            
            fitted =
              pred$fit,
            
            lower =
              pred$fit -
              1.96 *
              pred$se.fit,
            
            upper =
              pred$fit +
              1.96 *
              pred$se.fit
            
          )
        }
      )
    
  ) %>%
  
  select(
    survey_unit,
    prediction_data
  ) %>%
  
  tidyr::unnest(
    prediction_data
  )

############################################################
# COMMON PLOTTING FUNCTION
############################################################

make_regression_facets <- function(
    dat,
    title_text
){
  
  ggplot(
    dat,
    aes(
      x = year,
      y = observed
    )
  ) +
    
    geom_ribbon(
      aes(
        ymin = lower,
        ymax = upper
      ),
      fill = "lightblue",
      alpha = 0.35
    ) +
    
    geom_line(
      aes(
        y = fitted
      ),
      colour = "red",
      linewidth = 0.8
    ) +
    
    geom_point(
      size = 1.2,
      alpha = 0.75
    ) +
    
    facet_wrap(
      ~ survey_unit,
      scales = "free_y"
    ) +
    
    theme_bw() +
    
    labs(
      title = title_text,
      x = "Year",
      y = "Survey Mean CTI (°C)"
    ) +
    
    theme(
      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5
        ),
      strip.text =
        element_text(
          face = "bold"
        )
    )
}

############################################################
# CREATE FIGURES
############################################################

linear_regression_facets <-
  make_regression_facets(
    linear_prediction_data,
    "Linear Thermophilization Models"
  )

harmonic_regression_facets <-
  make_regression_facets(
    harmonic_prediction_data,
    "Harmonic Thermophilization Models"
  )

wgamm_regression_facets <-
  make_regression_facets(
    wgamm_prediction_data,
    "Weighted GAM Thermophilization Models"
  )

print(linear_regression_facets)
print(harmonic_regression_facets)
print(wgamm_regression_facets)

############################################################
# EXPORT REGRESSION FIGURES
############################################################

ggsave(
  "Linear_Regression_Facets.png",
  linear_regression_facets,
  width = 18,
  height = 14,
  dpi = 400
)

ggsave(
  "Harmonic_Regression_Facets.png",
  harmonic_regression_facets,
  width = 18,
  height = 14,
  dpi = 400
)

ggsave(
  "WGAMM_Regression_Facets.png",
  wgamm_regression_facets,
  width = 18,
  height = 14,
  dpi = 400
)
############################################################
# MODEL FIT SUMMARY
############################################################

fit_summary_table <-
  
  survey_results_all %>%
  
  group_by(method) %>%
  
  summarise(
    
    mean_rate =
      mean(
        annual_rate,
        na.rm = TRUE
      ),
    
    median_rate =
      median(
        annual_rate,
        na.rm = TRUE
      ),
    
    sd_rate =
      sd(
        annual_rate,
        na.rm = TRUE
      ),
    
    significant_surveys =
      sum(
        p_value < 0.05,
        na.rm = TRUE
      ),
    
    total_surveys =
      n(),
    
    percent_significant =
      100 *
      significant_surveys /
      total_surveys,
    
    .groups = "drop"
    
  )

print(fit_summary_table)
############################################################
#### SECTION 8: METHOD AGREEMENT DIAGNOSTICS ###############
############################################################

############################################################
# CREATE SURVEY-LEVEL WIDE COMPARISON TABLE
############################################################

thermophilization_method_comparison <-
  
  survey_results_all %>%
  
  mutate(
    
    estimate_ci =
      sprintf(
        "%.4f [%.4f, %.4f]",
        annual_rate,
        conf_low,
        conf_high
      )
    
  ) %>%
  
  select(
    
    survey_unit,
    
    method,
    
    annual_rate,
    
    estimate_ci,
    
    p_value,
    
    trend_cat,
    
    model_type,
    
    r_squared,
    
    adj_r2,
    
    deviance_explained,
    
    edf,
    
    n_hex,
    
    year_range
    
  ) %>%
  
  pivot_wider(
    
    names_from = method,
    
    values_from = c(
      
      annual_rate,
      
      estimate_ci,
      
      p_value,
      
      trend_cat,
      
      model_type,
      
      r_squared,
      
      adj_r2,
      
      deviance_explained,
      
      edf
      
    )
    
  ) %>%
  
  arrange(
    survey_unit
  )

print(
  thermophilization_method_comparison
)

gtsave(
  thermophilization_method_comparison %>%
    gt() %>%
    tab_header(
      title =
        md(
          "**Thermophilization Model Comparison by Survey Unit**"
        )
    ),
  "Thermophilization_Model_Comparison_by_Survey_Unit.html"
)

###############################################################################
# Compact publication-ready comparison table
###############################################################################

thermophilization_method_comparison <-
  
  survey_results_all %>%
  
  mutate(
    
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    
    fit_value = case_when(
      !is.na(adj_r2) ~ paste0("Adj R² = ", round(adj_r2, 3)),
      !is.na(r_squared) ~ paste0("R² = ", round(r_squared, 3)),
      !is.na(deviance_explained) ~ paste0("Dev = ", round(deviance_explained, 3)),
      TRUE ~ NA_character_
    ),
    
    summary_cell = paste0(
      sprintf("%.4f [%.4f, %.4f]",
              annual_rate,
              conf_low,
              conf_high),
      "<br>",
      trend_cat,
      " ",
      significance,
      "<br>",
      fit_value
    )
  ) %>%
  
  select(
    survey_unit,
    method,
    summary_cell,
    n_hex,
    year_range
  ) %>%
  
  distinct() %>%
  
  pivot_wider(
    names_from = method,
    values_from = summary_cell
  ) %>%
  
  arrange(survey_unit)

###############################################################################
# Export
###############################################################################

thermophilization_method_comparison %>%
  gt() %>%
  tab_header(
    title = md(
      "**Thermophilization Model Comparison by Survey Unit**"
    ),
    subtitle = md(
      "Annual trend (95% CI), significance, and model fit"
    )
  ) %>%
  cols_label(
    survey_unit = "Survey Unit"
  ) %>%
  fmt_markdown(
    columns = everything()
  ) %>%
  gtsave(
    "Thermophilization_Model_Comparison_by_Survey_Unit.html"
  )

############################################################
# SURVEY-LEVEL RATE CORRELATIONS
############################################################

method_correlations <-
  
  survey_results_all %>%
  
  select(
    survey_unit,
    method,
    annual_rate
  ) %>%
  
  pivot_wider(
    names_from = method,
    values_from = annual_rate
  )

############################################################
# PAIRWISE CORRELATIONS
############################################################

cor_linear_harmonic <-
  cor(
    method_correlations$Linear,
    method_correlations$Harmonic,
    use = "complete.obs"
  )

cor_linear_wgamm <-
  cor(
    method_correlations$Linear,
    method_correlations$wGAMM,
    use = "complete.obs"
  )

cor_harmonic_wgamm <-
  cor(
    method_correlations$Harmonic,
    method_correlations$wGAMM,
    use = "complete.obs"
  )

method_correlation_summary <-
  tibble(
    
    comparison =
      c(
        "Linear vs Harmonic",
        "Linear vs wGAMM",
        "Harmonic vs wGAMM"
      ),
    
    correlation =
      c(
        cor_linear_harmonic,
        cor_linear_wgamm,
        cor_harmonic_wgamm
      )
    
  )

print(
  method_correlation_summary
)

############################################################
# RATE SCATTERPLOTS
############################################################

p_linear_harmonic <-
  
  ggplot(
    method_correlations,
    aes(
      Linear,
      Harmonic
    )
  ) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "red"
  ) +
  
  theme_bw() +
  
  labs(
    title = "Linear vs Harmonic",
    x = expression(
      paste(
        "Linear Rate (",
        degree,
        "C yr"^-1,
        ")"
      )
    ),
    y = expression(
      paste(
        "Harmonic Rate (",
        degree,
        "C yr"^-1,
        ")"
      )
    )
  )

p_linear_wgamm <-
  
  ggplot(
    method_correlations,
    aes(
      Linear,
      wGAMM
    )
  ) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "red"
  ) +
  
  theme_bw() +
  
  labs(
    title = "Linear vs wGAMM",
    x = expression(
      paste(
        "Linear Rate (",
        degree,
        "C yr"^-1,
        ")"
      )
    ),
    y = expression(
      paste(
        "wGAMM Rate (",
        degree,
        "C yr"^-1,
        ")"
      )
    )
  )

p_harmonic_wgamm <-
  
  ggplot(
    method_correlations,
    aes(
      Harmonic,
      wGAMM
    )
  ) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "red"
  ) +
  
  theme_bw() +
  
  labs(
    title = "Harmonic vs wGAMM",
    x = expression(
      paste(
        "Harmonic Rate (",
        degree,
        "C yr"^-1,
        ")"
      )
    ),
    y = expression(
      paste(
        "wGAMM Rate (",
        degree,
        "C yr"^-1,
        ")"
      )
    )
  )

############################################################
# TREND-DIRECTION AGREEMENT
############################################################

trend_agreement <-
  
  survey_results_all %>%
  
  mutate(
    
    direction =
      case_when(
        
        annual_rate > 0 ~
          "Positive",
        
        annual_rate < 0 ~
          "Negative",
        
        TRUE ~
          "Zero"
      )
    
  ) %>%
  
  select(
    survey_unit,
    method,
    direction
  ) %>%
  
  pivot_wider(
    names_from = method,
    values_from = direction
  ) %>%
  
  mutate(
    
    all_agree =
      Linear ==
      Harmonic &
      Harmonic ==
      wGAMM
    
  )

print(
  trend_agreement
)

############################################################
# TREND AGREEMENT SUMMARY
############################################################

method_agreement_summary <-
  
  trend_agreement %>%
  
  summarise(
    
    total_surveys =
      n(),
    
    surveys_all_agree =
      sum(
        all_agree,
        na.rm = TRUE
      ),
    
    percent_agreement =
      100 *
      surveys_all_agree /
      total_surveys
    
  )

print(
  method_agreement_summary
)

############################################################
# SIGNIFICANCE AGREEMENT
############################################################

significance_agreement <-
  
  survey_results_all %>%
  
  mutate(
    
    significant =
      p_value < 0.05
    
  ) %>%
  
  select(
    survey_unit,
    method,
    significant
  ) %>%
  
  pivot_wider(
    names_from = method,
    values_from = significant
  ) %>%
  
  mutate(
    
    all_significance_agree =
      Linear ==
      Harmonic &
      Harmonic ==
      wGAMM
    
  )

print(
  significance_agreement
)

############################################################
# MODEL-SPECIFIC SUMMARY
############################################################

model_performance_summary <-
  
  survey_results_all %>%
  
  group_by(
    method
  ) %>%
  
  summarise(
    
    n_surveys =
      n(),
    
    mean_rate =
      mean(
        annual_rate,
        na.rm = TRUE
      ),
    
    median_rate =
      median(
        annual_rate,
        na.rm = TRUE
      ),
    
    sd_rate =
      sd(
        annual_rate,
        na.rm = TRUE
      ),
    
    significant_surveys =
      sum(
        p_value < 0.05,
        na.rm = TRUE
      ),
    
    percent_significant =
      100 *
      significant_surveys /
      n_surveys,
    
    mean_abs_rate =
      mean(
        abs(
          annual_rate
        ),
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  )

print(
  model_performance_summary
)

############################################################
# MANUSCRIPT-READY GT TABLE
############################################################

comparison_gt <-
  
  thermophilization_method_comparison %>%
  
  gt() %>%
  
  tab_header(
    title =
      md(
        "**Comparison of Thermophilization Estimation Methods**"
      )
  )

comparison_gt

############################################################
# SAVE TABLES
############################################################

write.csv(
  thermophilization_method_comparison,
  "Thermophilization_Method_Comparison.csv",
  row.names = FALSE
)

write.csv(
  model_performance_summary,
  "Thermophilization_Model_Performance.csv",
  row.names = FALSE
)

write.csv(
  method_correlation_summary,
  "Thermophilization_Method_Correlations.csv",
  row.names = FALSE
)

write.csv(
  trend_agreement,
  "Thermophilization_Trend_Agreement.csv",
  row.names = FALSE
)

write.csv(
  significance_agreement,
  "Thermophilization_Significance_Agreement.csv",
  row.names = FALSE
)

############################################################
# FINAL SUMMARY PRINTS
############################################################

cat(
  "\n==================================================\n",
  "THERMOPHILIZATION METHOD COMPARISON SUMMARY\n",
  "==================================================\n"
)

print(
  method_correlation_summary
)

print(
  method_agreement_summary
)

print(
  model_performance_summary
)

print(
  hex_summary_table
)

############################################################
# FINAL OBJECTS CREATED
############################################################

# SURVEY-LEVEL RESULTS
# --------------------
# linear_survey_results
# harmonic_survey_results
# wgamm_survey_results
# survey_results_all

# SURVEY TABLES
# -------------
# survey_model_comparison_table
# survey_model_wide
# survey_model_gt
# thermophilization_method_comparison
# comparison_gt

# LOLLIPOP FIGURES
# ----------------
# faceted_lollipop_plot

# HEX RESULTS
# -----------
# linear_hex_results
# harmonic_hex_results
# wgamm_hex_results
# hex_results_all

# HEX MAPS
# --------
# faceted_hex_map

# REGRESSION PLOTS
# ----------------
# linear_regression_facets
# harmonic_regression_facets
# wgamm_regression_facets

# PREDICTION DATA
# ---------------
# linear_prediction_data
# harmonic_prediction_data
# wgamm_prediction_data

# MODEL AGREEMENT ANALYSES
# ------------------------
# method_correlations
# method_correlation_summary
# trend_agreement
# significance_agreement
# method_agreement_summary
# fit_summary_table
# model_performance_summary
# hex_summary_table

############################################################
# RECOMMENDED MANUSCRIPT OUTPUTS
############################################################

# Figure 1:
# Faceted_Thermophilization_Lollipop.png

# Figure 2:
# Faceted_Hex_Thermophilization_Map.png

# Figure 3:
# Linear_Regression_Facets.png

# Figure 4:
# Harmonic_Regression_Facets.png

# Figure 5:
# WGAMM_Regression_Facets.png

# Table 1:
# Thermophilization_Method_Comparison.csv

# Supplementary Table:
# Thermophilization_Model_Performance.csv

# Supplementary Table:
# Thermophilization_Method_Correlations.csv

##########################################################
##### FACETED Q1 THERMOPHILIZATION DIAGNOSTIC FIGURE #####
##########################################################
library(dplyr)
library(ggplot2)
library(patchwork)

target_survey_unit <- "Faroe Bank Q3"

############################################################
# SIGNIFICANCE FORMATTER
############################################################

fmt_p <- function(p){
  
  case_when(
    p < 0.001 ~ "<0.001***",
    p < 0.01  ~ paste0(sprintf("%.3f", p), "**"),
    p < 0.05  ~ paste0(sprintf("%.3f", p), "*"),
    TRUE      ~ sprintf("%.3f", p)
  )
  
}

############################################################
# MODEL RESULTS
############################################################

title_info <-
  
  survey_results_all %>%
  
  filter(
    survey_unit == target_survey_unit
  )

linear_rate <-
  
  title_info %>%
  
  filter(method == "Linear") %>%
  
  pull(annual_rate)

linear_p <-
  
  title_info %>%
  
  filter(method == "Linear") %>%
  
  pull(p_value)

harmonic_rate <-
  
  title_info %>%
  
  filter(method == "Harmonic") %>%
  
  pull(annual_rate)

harmonic_p <-
  
  title_info %>%
  
  filter(method == "Harmonic") %>%
  
  pull(p_value)

wgamm_rate <-
  
  title_info %>%
  
  filter(method == "wGAMM") %>%
  
  pull(annual_rate)

wgamm_p <-
  
  title_info %>%
  
  filter(method == "wGAMM") %>%
  
  pull(p_value)

############################################################
# REGRESSION TITLES
############################################################

linear_title <-
  
  paste0(
    "Linear (",
    sprintf("%.4f", linear_rate),
    " °C yr⁻¹; P ",
    fmt_p(linear_p),
    ")"
  )

harmonic_title <-
  
  paste0(
    "Harmonic (",
    sprintf("%.4f", harmonic_rate),
    " °C yr⁻¹; P ",
    fmt_p(harmonic_p),
    ")"
  )

wgamm_title <-
  
  paste0(
    "wGAMM (",
    sprintf("%.4f", wgamm_rate),
    " °C yr⁻¹; P ",
    fmt_p(wgamm_p),
    ")"
  )

############################################################
# UNWEIGHTED DENSITY DATA
############################################################

density_unweighted <-
  
  survey_results_all %>%
  
  mutate(
    
    method = factor(
      method,
      levels = c(
        "Linear",
        "Harmonic",
        "wGAMM"
      )
    )
    
  ) %>%
  
  filter(
    !is.na(annual_rate)
  )
survey_effort <-
  
  Combined_FISH_optimized_spatially_valid %>%
  
  group_by(survey_unit) %>%
  
  summarise(
    total_hauls = n_distinct(haul_id),
    .groups = "drop"
  )

density_dat <-
  
  survey_results_all %>%
  
  left_join(
    survey_effort,
    by = "survey_unit"
  ) %>%
  
  mutate(
    
    weight = sqrt(total_hauls),
    
    weight =
      weight /
      mean(weight, na.rm = TRUE)
    
  )
'
density_dat <-

  survey_results_all %>%

  left_join(
    survey_effort,
    by = "survey_unit"
  ) %>%

  mutate(

    se =
      case_when(
        !is.na(std_error) & std_error > 0 ~ std_error,
        TRUE ~ (conf_high - conf_low)/(2*1.96)
      ),

    weight =
      (1/se) *
      sqrt(total_hauls),

    weight =
      weight /
      mean(weight, na.rm = TRUE)

  )
'
############################################################
# COMMON X SCALE ACROSS ALL DENSITY PANELS
############################################################

density_xlims <-
  
  range(
    c(
      density_unweighted$annual_rate,
      density_dat$annual_rate
    ),
    na.rm = TRUE
  )

############################################################
# COMMON Y SCALE ACROSS ALL DENSITY PANELS
# (COLUMN A + COLUMN B USE IDENTICAL SCALE)
############################################################

unweighted_ymax <-
  
  max(
    
    sapply(
      
      c(
        "Linear",
        "Harmonic",
        "wGAMM"
      ),
      
      function(m){
        
        max(
          
          density(
            density_unweighted$annual_rate[
              density_unweighted$method == m
            ]
          )$y,
          
          na.rm = TRUE
          
        )
        
      }
      
    )
    
  )

weighted_ymax <-
  
  max(
    
    sapply(
      
      c(
        "Linear",
        "Harmonic",
        "wGAMM"
      ),
      
      function(m){
        
        sub <-
          
          density_dat %>%
          
          filter(
            method == m
          )
        
        max(
          
          density(
            sub$annual_rate,
            weights =
              sub$weight /
              sum(sub$weight)
          )$y,
          
          na.rm = TRUE
          
        )
        
      }
      
    )
    
  )

density_ymax <-
  
  max(
    unweighted_ymax,
    weighted_ymax
  ) * 1.10

############################################################
# DENSITY PANEL FUNCTION
############################################################

make_density_panel <- function(
    dat,
    method_name,
    weighted = FALSE,
    panel_letter = NULL,
    panel_subtitle,
    fill_colour
){
  
  density_sub <-
    
    dat %>%
    
    filter(
      method == method_name
    )
  
  mean_rate <-
    
    if(weighted){
      
      weighted.mean(
        density_sub$annual_rate,
        density_sub$weight,
        na.rm = TRUE
      )
      
    } else {
      
      mean(
        density_sub$annual_rate,
        na.rm = TRUE
      )
      
    }
  
  p <-
    
    ggplot(
      density_sub,
      aes(
        x = annual_rate
      )
    )
  
  if(weighted){
    
    p <-
      
      p +
      
      geom_density(
        aes(
          weight = weight
        ),
        fill = fill_colour,
        colour = "black",
        alpha = 0.80,
        linewidth = 0.8,
        adjust = 1
      )
    
  } else {
    
    p <-
      
      p +
      
      geom_density(
        fill = fill_colour,
        colour = "black",
        alpha = 0.80,
        linewidth = 0.8,
        adjust = 1
      )
    
  }
  
  p +
    
    geom_vline(
      xintercept = mean_rate,
      colour = "red",
      linetype = "dashed",
      linewidth = 0.9
    ) +
    
    annotate(
      "label",
      x = mean_rate,
      y = density_ymax * 0.91,
      label = paste0(
        "Mean = ",
        sprintf(
          "%.5f",
          mean_rate
        )
      ),
      colour = "red",
      fill = "white",
      fontface = "bold",
      size = 3.0,
      linewidth  = 0.20,
      hjust = -0.05
    ) +
    
    coord_cartesian(
      xlim = density_xlims,
      ylim = c(
        0,
        density_ymax
      ),
      clip = "off"
    ) +
    
    labs(
      title = method_name,
      subtitle = panel_subtitle,
      x = NULL,
      y =
        if(weighted)
          "Weighted Density"
      else
        "Density",
      tag = panel_letter
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      
      plot.title =
        element_text(
          face = "bold",
          size = 11,
          hjust = 0.5
        ),
      
      plot.subtitle =
        element_text(
          size = 8.5,
          hjust = 0.5
        ),
      
      plot.tag =
        element_text(
          face = "bold",
          size = 16
        ),
      
      axis.title =
        element_text(
          face = "bold"
        )
      
    )
  
}

############################################################
# COLUMN A
# UNWEIGHTED DISTRIBUTIONS
############################################################

A_linear <-
  
  make_density_panel(
    density_unweighted,
    "Linear",
    weighted = FALSE,
    panel_letter = "A",
    panel_subtitle = "Weights: none",
    fill_colour = "#C6DBEF"
  )

A_harmonic <-
  
  make_density_panel(
    density_unweighted,
    "Harmonic",
    weighted = FALSE,
    panel_letter = NULL,
    panel_subtitle = "Weights: none",
    fill_colour = "#C6DBEF"
  )

A_wgamm <-
  
  make_density_panel(
    density_unweighted,
    "wGAMM",
    weighted = FALSE,
    panel_letter = NULL,
    panel_subtitle = "Weights: none",
    fill_colour = "#C6DBEF"
  ) +
  
  labs(
    x =
      expression(
        paste(
          "Thermophilization Rate (",
          degree,
          "C yr"^-1,
          ")"
        )
      )
  )

############################################################
# COLUMN B
# PRECISION + EFFORT WEIGHTED
############################################################

B_linear <-
  
  make_density_panel(
    density_dat,
    "Linear",
    weighted = TRUE,
    panel_letter = "B",
    panel_subtitle =
      "Weights: survey effort (√total hauls)",
    fill_colour = "steelblue"
  )

B_harmonic <-
  
  make_density_panel(
    density_dat,
    "Harmonic",
    weighted = TRUE,
    panel_letter = NULL,
    panel_subtitle =
      "Weights: survey effort (√total hauls)",
    fill_colour = "steelblue"
  )

B_wgamm <-
  
  make_density_panel(
    density_dat,
    "wGAMM",
    weighted = TRUE,
    panel_letter = NULL,
    panel_subtitle =
      "Weights: survey effort (√total hauls)",
    fill_colour = "steelblue"
  ) +
  
  labs(
    x =
      expression(
        paste(
          "Thermophilization Rate (",
          degree,
          "C yr"^-1,
          ")"
        )
      )
  )

############################################################
# COLUMN C
# REGRESSION PANELS
############################################################

p_linear_single <-
  
  linear_prediction_data %>%
  
  filter(
    survey_unit == target_survey_unit
  ) %>%
  
  ggplot(
    aes(
      x = year,
      y = observed
    )
  ) +
  
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    fill = "lightblue",
    alpha = 0.35
  ) +
  
  geom_line(
    aes(
      y = fitted
    ),
    colour = "red",
    linewidth = 0.9
  ) +
  
  geom_point(
    size = 2
  ) +
  
  labs(
    title = linear_title,
    x = NULL,
    y = "Mean CTI (°C)"
  ) +
  
  theme_bw(base_size = 11)



p_harmonic_single <-
  
  harmonic_prediction_data %>%
  
  filter(
    survey_unit == target_survey_unit
  ) %>%
  
  ggplot(
    aes(
      x = year,
      y = observed
    )
  ) +
  
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    fill = "lightblue",
    alpha = 0.35
  ) +
  
  geom_line(
    aes(
      y = fitted
    ),
    colour = "red",
    linewidth = 0.9
  ) +
  
  geom_point(
    size = 2
  ) +
  
  labs(
    title = harmonic_title,
    x = NULL,
    y = "Mean CTI (°C)"
  ) +
  
  theme_bw(base_size = 11)



p_wgamm_single <-
  
  wgamm_prediction_data %>%
  
  filter(
    survey_unit == target_survey_unit
  ) %>%
  
  ggplot(
    aes(
      x = year,
      y = observed
    )
  ) +
  
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    fill = "lightblue",
    alpha = 0.35
  ) +
  
  geom_line(
    aes(
      y = fitted
    ),
    colour = "red",
    linewidth = 0.9
  ) +
  
  geom_point(
    size = 2
  ) +
  
  labs(
    title = wgamm_title,
    x = "Year",
    y = "Mean CTI (°C)"
  ) +
  
  theme_bw(base_size = 11)

p_linear_single <-
  
  p_linear_single +
  
  labs(
    tag = "C"
  ) +
  
  theme(
    plot.tag =
      element_text(
        face = "bold",
        size = 16
      )
  )

p_harmonic_single <-
  
  p_harmonic_single +
  
  labs(
    tag = NULL
  )

p_wgamm_single <-
  
  p_wgamm_single +
  
  labs(
    tag = NULL
  )

############################################################
# FINAL 3×3 FIGURE
############################################################

survey_model_diagnostic <-
  
  (
    A_linear |
      B_linear |
      p_linear_single
  ) /
  
  (
    A_harmonic |
      B_harmonic |
      p_harmonic_single
  ) /
  
  (
    A_wgamm |
      B_wgamm |
      p_wgamm_single
  ) +
  
  plot_layout(
    widths = c(
      1,
      1,
      1.20
    )
  )

############################################################
# DISPLAY
############################################################

survey_model_diagnostic

############################################################
# EXPORT
############################################################

ggsave(
  paste0(
    gsub(
      " ",
      "_",
      target_survey_unit
    ),
    "_Thermophilization_Diagnostic_Figure_3x3.png"
  ),
  survey_model_diagnostic,
  width = 16,
  height = 12,
  dpi = 600,
  bg = "white"
)


##########################################################
##### Q2: Regional drivers of thermophilizaiton ##########
##########################################################

colnames(global_spatially_valid_hex_stats)

#####################################################
#### return to building haul-hex-survey unit stats ##
#####################################################

calculate_yearly_hex_troph <- function(df) {
  
  yearly_hex_stats <- df %>%
    # 1. Filter out rows where Trophic data or SE is missing
    filter(!is.na(avg_troph), !is.na(total_se_troph), total_se_troph > 0) %>%
    
    group_by(survey_unit, year, hex_id) %>% # Updated to hexagon_id based on your colnames
    summarise(
      # Calculate weights (Inverse Variance: 1 / SE^2)
      # Summing the weights for the denominator
      total_weight = sum(1 / (total_se_troph^2), na.rm = TRUE),
      
      # 2. Calculate Weighted Average Trophic Level
      # Sum of (Value * Weight) / Total Weight
      avg_troph_hex = sum(avg_troph * (1 / (total_se_troph^2)), na.rm = TRUE) / total_weight,
      
      # 3. Calculate the Variance of the weighted mean
      # Propagated Variance = 1 / sum of weights
      var_troph_hex = 1 / total_weight,
      
      .groups = "drop"
    ) %>%
    # Convert Variance back to SE for easier interpretation later
    mutate(se_troph_hex = sqrt(var_troph_hex)) %>%
    select(-total_weight)
  
  return(yearly_hex_stats)
}

# Apply the function
global_troph_by_hex <- calculate_yearly_hex_troph(global_spatially_valid_haul_stats)

global_spatially_valid_hex_stats <- global_spatially_valid_hex_stats %>%
  left_join(
    global_troph_by_hex %>% select(-se_troph_hex), 
    by = c("survey_unit", "year", "hex_id")
  )

calculate_survey_yearly_troph <- function(hex_df) {
  
  survey_yearly_stats <- hex_df %>%
    # 1. Filter out invalid variances
    filter(!is.na(avg_troph_hex), !is.na(var_troph_hex), var_troph_hex > 0) %>%
    
    group_by(survey_unit, year) %>%
    summarise(
      # Weights derived from the hexagon-level variance
      total_survey_weight = sum(1 / var_troph_hex, na.rm = TRUE),
      
      # 2. Survey Unit Weighted Average Trophic Level
      avg_troph_survey = sum(avg_troph_hex * (1 / var_troph_hex), na.rm = TRUE) / total_survey_weight,
      
      # 3. Survey Unit Variance
      var_troph_survey = 1 / total_survey_weight,
      
      .groups = "drop"
    ) %>%
    # Convert Variance back to SE
    mutate(se_troph_survey = sqrt(var_troph_survey)) %>%
    select(-total_survey_weight)
  
  return(survey_yearly_stats)
}

# Now run the survey-level aggregation
global_troph_by_year <- calculate_survey_yearly_troph(global_troph_by_hex)

################################################################################
### survey unit stats dataframe                 ###
################################################################################
global_spatially_valid_surveyunit_stats <-global_spatially_valid_hex_stats %>%
  group_by(survey_unit, year) %>%
  summarise(
    lme_number = list(unique(lme_number)),
    lme_name = list(unique(lme_name)),
    lat = mean(lat, na.rm = TRUE),
    lon = mean(lon, na.rm = TRUE),
    depth = mean(depth, na.rm = TRUE),
    CV_sst = sd(sst, na.rm = TRUE) / mean(sst, na.rm = TRUE),
    sst = mean(sst, na.rm = TRUE),
    CV_sbt = sd(sbt, na.rm = TRUE) / mean(sbt, na.rm = TRUE),
    sbt = mean(sbt, na.rm = TRUE),
    glorys_bottomT = mean(glorys_bottomT, na.rm = TRUE),
    
    evenness_combined = mean(evenness_combined, na.rm = TRUE),
    CV_density_combined = sd(density_combined, na.rm = TRUE) / mean(density_combined, na.rm = TRUE),
    density_combined = mean(density_combined, na.rm = TRUE), 
    S_combined = mean(S_combined, na.rm = TRUE),
    body_weight_combined = mean(body_weight_combined, na.rm = TRUE),
    
    tl = mean(tl, na.rm = TRUE), 
    AR = mean(AR, na.rm = TRUE),
    offspring.size = mean(offspring.size, na.rm = TRUE),
    age.maturity = mean(age.maturity, na.rm = TRUE), 
    fecundity = mean(age.maturity, na.rm = TRUE), 
    length.infinity = mean(length.infinity, na.rm = TRUE), 
    growth.coefficient = mean(growth.coefficient, na.rm = TRUE),
    length_max = mean(length_max, na.rm = TRUE),
    age.max = mean(age.max, na.rm = TRUE),
    
    # Categorical traits (Dominant value across hauls)
    dom_habitat    = get_dominant(dom_habitat),
    dom_feeding    = get_dominant(dom_feeding),
    dom_spawning   = get_dominant(dom_spawning),
    dom_body.shape = get_dominant(dom_body.shape), 
    dom_fin.shape = get_dominant(dom_fin.shape),
    gear = list(unique(gear)),
    .groups = "drop"
  ) 

global_spatially_valid_surveyunit_stats <- global_spatially_valid_surveyunit_stats %>%
  mutate(across(where(is.numeric), ~ ifelse(is.finite(.), ., NA_real_)))

global_spatially_valid_surveyunit_stats <- global_spatially_valid_surveyunit_stats %>%
  left_join(
    global_MTC_by_year %>% select(c(survey_unit, year, MTC_survey_avg, MTC_survey_variance)),
    by = c("survey_unit", "year")
  )

global_spatially_valid_surveyunit_stats <- global_spatially_valid_surveyunit_stats %>%
  left_join(
    global_troph_by_year %>% select(-c(se_troph_survey)),
    by = c("survey_unit", "year")
  )

#####################################################
##### LME based fishing pressure estimates ##########
#####################################################
library(readr)
library(dplyr)
library(purrr)

# 1. Define the repair logic
repair_sau_data <- function(df) {
  # List of valid LME names found in your data
  valid_lmes <- c(
    "West Bering Sea", "Aleutian Islands", "East Bering Sea", 
    "California Current", "Southeast U.S. Continental Shelf", 
    "Northeast U.S. Continental Shelf", "Canary Current", 
    "Celtic-Biscay Shelf", "Greenland Sea", "Scotian Shelf", 
    "Hudson Bay Complex", "Northern Bering - Chukchi Seas", 
    "North Sea", "Newfoundland-Labrador Shelf", "Gulf of Alaska", 
    "Canadian Eastern Arctic -West Greenland", "Baltic Sea", 
    "Iceland Shelf and Sea", "Gulf of Mexico", "Norwegian Sea", 
    "Iberian Coastal", "Barents Sea", "Black Sea", "Mediterranean Sea", "Faroe Plateau"  
  )
  
  df %>%
    # Remove rows where area_name is a scientific name, a year, or starts with a comma
    #filter(area_name %in% valid_lmes) %>%
    # Standardize data types
    mutate(
      across(c(area_name, area_type, scientific_name, common_name, 
               functional_group, commercial_group, fishing_entity, 
               fishing_sector, catch_type, reporting_status, 
               gear_type, end_use_type), as.character),
      across(c(year, tonnes, landed_value), as.numeric)
    )
}


SAU_LME_21_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 21 v50-1.csv")
SAU_LME_21_v50_1      <- repair_sau_data(SAU_LME_21_v50_1)

SAU_LME_59_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 59 v50-1.csv")
SAU_LME_59_v50_1     <- repair_sau_data(SAU_LME_59_v50_1)

SAU_LME_25_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 25 v50-1.csv")
SAU_LME_25_v50_1 <- repair_sau_data(SAU_LME_25_v50_1)

SAU_LME_5_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 5 v50-1.csv")
SAU_LME_5_v50_1 <- repair_sau_data(SAU_LME_5_v50_1)

SAU_LME_20_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 20 v50-1.csv")
SAU_LME_20_v50_1 <- repair_sau_data(SAU_LME_20_v50_1)

SAU_LME_8_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 8 v50-1.csv")
SAU_LME_8_v50_1 <- repair_sau_data(SAU_LME_8_v50_1)

SAU_LME_54_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 54 v50-1.csv")
SAU_LME_54_v50_1 <- repair_sau_data(SAU_LME_54_v50_1)

SAU_LME_63_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 63 v50-1.csv")
SAU_LME_63_v50_1 <- repair_sau_data(SAU_LME_63_v50_1)

SAU_LME_19_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 19 v50-1.csv")
SAU_LME_19_v50_1 <- repair_sau_data(SAU_LME_19_v50_1)

SAU_LME_24_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 24 v50-1.csv")
SAU_LME_24_v50_1 <- repair_sau_data(SAU_LME_24_v50_1)

SAU_LME_22_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 22 v50-1.csv")
fix_northsea <- function(df){
  
  bad <- df$area_name == ",lme"
  
  fixed <- df
  
  fixed$landed_value[bad] <- fixed$tonnes[bad]
  fixed$tonnes[bad] <- NA
  
  fixed$tonnes[bad] <- suppressWarnings(
    as.numeric(fixed$end_use_type[bad])
  )
  
  fixed$end_use_type[bad] <- fixed$gear_type[bad]
  fixed$gear_type[bad] <- fixed$reporting_status[bad]
  fixed$reporting_status[bad] <- fixed$catch_type[bad]
  fixed$catch_type[bad] <- fixed$fishing_sector[bad]
  fixed$fishing_sector[bad] <- fixed$fishing_entity[bad]
  fixed$fishing_entity[bad] <- fixed$commercial_group[bad]
  fixed$commercial_group[bad] <- fixed$functional_group[bad]
  fixed$functional_group[bad] <- fixed$common_name[bad]
  fixed$common_name[bad] <- fixed$scientific_name[bad]
  fixed$scientific_name[bad] <- fixed$year[bad]
  
  fixed$year[bad] <- fixed$area_type[bad]
  fixed$area_type[bad] <- "lme"
  fixed$area_name[bad] <- "North Sea"
  
  fixed
  
}
SAU_LME_22_v50_1 <- fix_northsea(SAU_LME_22_v50_1)
SAU_LME_22_v50_1 <- repair_sau_data(SAU_LME_22_v50_1)

SAU_LME_9_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 9 v50-1.csv")
SAU_LME_9_v50_1 <- repair_sau_data(SAU_LME_9_v50_1)

SAU_LME_2_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 2 v50-1.csv")
SAU_LME_2_v50_1 <- repair_sau_data(SAU_LME_2_v50_1)

SAU_LME_18_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 18 v50-1.csv")
SAU_LME_18_v50_1 <- repair_sau_data(SAU_LME_18_v50_1)

SAU_LME_23_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 23 v50-1.csv")
SAU_LME_23_v50_1 <- repair_sau_data(SAU_LME_23_v50_1)

fix_shifted_block <- function(df, lme_name){
  
  bad <- df$area_name == ",lme"
  
  fixed <- df
  
  fixed$landed_value[bad] <- fixed$tonnes[bad]
  
  fixed$tonnes[bad] <- suppressWarnings(
    as.numeric(fixed$end_use_type[bad])
  )
  
  fixed$end_use_type[bad] <- fixed$gear_type[bad]
  fixed$gear_type[bad] <- fixed$reporting_status[bad]
  fixed$reporting_status[bad] <- fixed$catch_type[bad]
  fixed$catch_type[bad] <- fixed$fishing_sector[bad]
  fixed$fishing_sector[bad] <- fixed$fishing_entity[bad]
  fixed$fishing_entity[bad] <- fixed$commercial_group[bad]
  fixed$commercial_group[bad] <- fixed$functional_group[bad]
  fixed$functional_group[bad] <- fixed$common_name[bad]
  fixed$common_name[bad] <- fixed$scientific_name[bad]
  
  fixed$scientific_name[bad] <- NA
  
  fixed$year[bad] <- as.numeric(fixed$area_type[bad])
  fixed$area_type[bad] <- "lme"
  fixed$area_name[bad] <- lme_name
  
  fixed
}
SAU_LME_62_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 62 v50-1.csv")
repair_black <- fix_shifted_block(
  SAU_LME_62_v50_1,
  "Black Sea"
)
SAU_LME_62_v50_1 <- repair_sau_data(repair_black)

SAU_LME_26_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 26 v50-1.csv")

SAU_LME_26_v50_1 <- fix_shifted_block(
  SAU_LME_26_v50_1,
  "Mediterranean Sea"
)

SAU_LME_26_v50_1 <- repair_sau_data(SAU_LME_26_v50_1)

SAU_LME_6_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 6 v50-1.csv")
SAU_LME_6_v50_1 <- repair_sau_data(SAU_LME_6_v50_1)

SAU_LME_53_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 53 v50-1.csv")
SAU_LME_53_v50_1 <- repair_sau_data(SAU_LME_53_v50_1)

SAU_LME_7_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 7 v50-1.csv")
SAU_LME_7_v50_1 <- repair_sau_data(SAU_LME_7_v50_1)

SAU_LME_65_3_27_1_60_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 65,3,27,1,60 v50-1.csv")
SAU_LME_65_3_27_1_60_v50_1 <- repair_sau_data(SAU_LME_65_3_27_1_60_v50_1)

# 2. Standardize classes across the list
LME_fishing_pressure_SeaAroundUs <- bind_rows(
  SAU_LME_21_v50_1, SAU_LME_59_v50_1, SAU_LME_25_v50_1, SAU_LME_5_v50_1,
  SAU_LME_20_v50_1, SAU_LME_8_v50_1, SAU_LME_54_v50_1, SAU_LME_63_v50_1,
  SAU_LME_19_v50_1, SAU_LME_24_v50_1, SAU_LME_22_v50_1, SAU_LME_9_v50_1,
  SAU_LME_2_v50_1, SAU_LME_18_v50_1, SAU_LME_23_v50_1, SAU_LME_62_v50_1,
  SAU_LME_26_v50_1, SAU_LME_6_v50_1, SAU_LME_53_v50_1, SAU_LME_7_v50_1, SAU_LME_65_3_27_1_60_v50_1
)

"Record shifted columns of SeaAroundUs data for posterity"
SAU_LME_26_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 26 v50-1.csv")
SAU_LME_62_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 62 v50-1.csv")
SAU_LME_22_v50_1 <- read_csv("C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU LME 22 v50-1.csv")

shift_rows <- tibble(
  dataframe = c(
    "SAU_LME_22_v50_1",
    "SAU_LME_62_v50_1",
    "SAU_LME_26_v50_1"
  ),
  first_shifted_row = c(
    which(grepl("^[0-9]{4}$", SAU_LME_22_v50_1$scientific_name))[1],
    which(SAU_LME_62_v50_1$area_name == ",lme")[1],
    which(SAU_LME_26_v50_1$area_name == ",lme")[1]
  ), 
  total_rows = c(
    nrow(SAU_LME_22_v50_1),
    nrow(SAU_LME_62_v50_1),
    nrow(SAU_LME_26_v50_1)
  )
)
shift_rows

#########################################################
######### add high seas data for Rockall ################
#########################################################
library(readr)
library(dplyr)

# Import once to identify where parsing first breaks
SAU_HighSeas_27_v50_1 <- read_csv(
  "C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU HighSeas 27 v50-1.csv",
  show_col_types = FALSE
)

# Find first problematic row
first_bad_row <- min(problems(SAU_HighSeas_27_v50_1)$row)

# Read raw file
file_path <- "C:/Users/besch/Desktop/flagged FISHGLOB/LME fishing pressure data/SAU HighSeas 27 v50-1.csv"
lines <- readLines(file_path)

# Print diagnostics
cat("First problematic row:", first_bad_row, "\n")
cat("Total rows in file (excluding header):", length(lines) - 1, "\n")

# Header determines expected number of columns
expected_cols <- length(strsplit(lines[1], ",")[[1]])

# Fix all rows from the first bad row onward
for(i in first_bad_row:length(lines)) {
  
  fields <- strsplit(
    lines[i],
    ',(?=(?:[^"]*"[^"]*")*[^"]*$)',
    perl = TRUE
  )[[1]]
  
  if(length(fields) == expected_cols - 1) {
    fields <- c(fields, "")
  }
  
  if(length(fields) < expected_cols) {
    fields <- c(fields, rep("", expected_cols - length(fields)))
  }
  
  lines[i] <- paste(fields, collapse = ",")
}

# Write temporary corrected file
tmp_file <- tempfile(fileext = ".csv")
writeLines(lines, tmp_file)

# Re-import and overwrite original object
SAU_HighSeas_27_v50_1 <- read_csv(
  tmp_file,
  show_col_types = FALSE
)

# Verify fix
cat("Remaining parsing issues:", nrow(problems(SAU_HighSeas_27_v50_1)), "\n")

#########################################################
#########  ################
#########################################################

# Create the aggregated dataframe
LME_annual_tonnes <- LME_fishing_pressure_SeaAroundUs %>%
  group_by(area_name, year) %>%
  summarize(LME_total_tonnes_fished = sum(tonnes, na.rm = TRUE), .groups = "drop")
print(unique(LME_annual_tonnes$area_name))
print(unique(global_spatially_valid_hex_stats$lme_name))

High_Seas_annual_tonnes <- SAU_HighSeas_27_v50_1 %>%
  group_by(area_name, year) %>%
  summarize(High_Seas_total_tonnes_fished = sum(tonnes, na.rm = TRUE), .groups = "drop")

LME_annual_tonnes <- LME_annual_tonnes %>%
  mutate(
    area_name = recode(
      area_name,
      
      "Canadian Eastern Arctic -West Greenland" =
        "Canadian Eastern Arctic - West Greenland",
      
      "Newfoundland-Labrador Shelf" =
        "Labrador - Newfoundland"
    )
  )
print(unique(LME_annual_tonnes$area_name))
print(unique(global_spatially_valid_hex_stats$lme_name))
'
##############################################################################
# DIAGNOSTICS: CHECK COVERAGE OF LME-YEAR COMBINATIONS
##############################################################################

library(dplyr)

# ---------------------------------------------------------------------------
# 1. UNIQUE LME-YEAR COMBINATIONS AVAILABLE FROM FISHING PRESSURE DATA
# ---------------------------------------------------------------------------

source_combinations <- bind_rows(
  
  LME_annual_tonnes %>%
    transmute(
      lme_name = area_name,
      year
    ),
  
  High_Seas_annual_tonnes %>%
    transmute(
      lme_name = area_name,
      year
    )
  
) %>%
  distinct()

# ---------------------------------------------------------------------------
# 2. UNIQUE LME-YEAR COMBINATIONS PRESENT IN HEX DATA
# ---------------------------------------------------------------------------

hex_combinations <- global_spatially_valid_hex_stats %>%
  select(lme_name, year) %>%
  distinct()

# ---------------------------------------------------------------------------
# 3. WHICH HEX LME-YEAR COMBINATIONS ARE NOT COVERED?
# ---------------------------------------------------------------------------

missing_combinations <- anti_join(
  hex_combinations,
  source_combinations,
  by = c("lme_name", "year")
)

cat("\n====================================================\n")
cat("UNMATCHED LME-YEAR COMBINATIONS\n")
cat("====================================================\n")

print(
  missing_combinations %>%
    arrange(lme_name, year),
  n = Inf
)

# ---------------------------------------------------------------------------
# 4. COVERAGE SUMMARY
# ---------------------------------------------------------------------------

coverage_summary <- tibble(
  total_hex_combinations =
    nrow(hex_combinations),
  
  matched_combinations =
    nrow(
      semi_join(
        hex_combinations,
        source_combinations,
        by = c("lme_name", "year")
      )
    ),
  
  unmatched_combinations =
    nrow(missing_combinations)
) %>%
  mutate(
    pct_matched =
      round(
        100 * matched_combinations / total_hex_combinations,
        2
      )
  )

cat("\n====================================================\n")
cat("OVERALL COVERAGE SUMMARY\n")
cat("====================================================\n")

print(coverage_summary)

# ---------------------------------------------------------------------------
# 5. COUNT UNMATCHED COMBINATIONS BY LME
# ---------------------------------------------------------------------------

unmatched_by_lme <- missing_combinations %>%
  count(lme_name, sort = TRUE)

cat("\n====================================================\n")
cat("UNMATCHED COMBINATIONS BY LME\n")
cat("====================================================\n")

print(unmatched_by_lme, n = Inf)

# ---------------------------------------------------------------------------
# 6. WHICH SURVEYS ARE AFFECTED?
# ---------------------------------------------------------------------------

affected_surveys <- global_spatially_valid_hex_stats %>%
  distinct(survey_unit, lme_name, year) %>%
  semi_join(
    missing_combinations,
    by = c("lme_name", "year")
  ) %>%
  arrange(survey_unit, year)

cat("\n====================================================\n")
cat("SURVEYS CONTAINING UNMATCHED LME-YEAR COMBINATIONS\n")
cat("====================================================\n")

print(affected_surveys, n = Inf)

# ---------------------------------------------------------------------------
# 7. OPTIONAL: LIST THE YEAR RANGE AVAILABLE IN EACH SOURCE LME
# ---------------------------------------------------------------------------

source_year_ranges <- bind_rows(
  
  LME_annual_tonnes %>%
    transmute(
      lme_name = area_name,
      year,
      source = "LME"
    ),
  
  High_Seas_annual_tonnes %>%
    transmute(
      lme_name = area_name,
      year,
      source = "High Seas"
    )
  
) %>%
  group_by(lme_name, source) %>%
  summarise(
    min_year = min(year, na.rm = TRUE),
    max_year = max(year, na.rm = TRUE),
    n_years   = n_distinct(year),
    .groups = "drop"
  )

cat("\n====================================================\n")
cat("AVAILABLE YEAR RANGES IN SOURCE DATASETS\n")
cat("====================================================\n")

print(source_year_ranges, n = Inf)
'
##############################################################################
# JOIN LME FISHING PRESSURE
##############################################################################

global_spatially_valid_hex_stats <- global_spatially_valid_hex_stats %>%
  left_join(
    LME_annual_tonnes %>%
      rename(lme_name = area_name),
    by = c("lme_name", "year")
  )

##############################################################################
# JOIN HIGH SEAS FISHING PRESSURE (ROCKALL ONLY)
##############################################################################

global_spatially_valid_hex_stats <- global_spatially_valid_hex_stats %>%
  left_join(
    High_Seas_annual_tonnes %>%
      select(year, High_Seas_total_tonnes_fished),
    by = "year"
  ) %>%
  mutate(
    LME_total_tonnes_fished = case_when(
      survey_unit == "Rockall" &
        is.na(LME_total_tonnes_fished) &
        !is.na(High_Seas_total_tonnes_fished)
      ~ High_Seas_total_tonnes_fished,
      
      TRUE ~ LME_total_tonnes_fished
    )
  ) %>%
  select(-High_Seas_total_tonnes_fished)

'
##############################################################################
# QA CHECKS FOR FISHING-PRESSURE JOIN
##############################################################################

# ---------------------------------------------------------------------------
# 1. Overall missingness
# ---------------------------------------------------------------------------

global_spatially_valid_hex_stats %>%
  summarise(
    total_rows = n(),
    missing_fishing_pressure =
      sum(is.na(LME_total_tonnes_fished)),
    pct_missing =
      100 * mean(is.na(LME_total_tonnes_fished))
  )

# ---------------------------------------------------------------------------
# 2. Missingness by survey
# ---------------------------------------------------------------------------

global_spatially_valid_hex_stats %>%
  group_by(survey_unit) %>%
  summarise(
    pct_missing =
      100 * mean(is.na(LME_total_tonnes_fished)),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_missing)) %>%
  print(n = Inf)

# ---------------------------------------------------------------------------
# 3. Verify Rockall received High Seas fishing pressure
# ---------------------------------------------------------------------------

rockall_check <- global_spatially_valid_hex_stats %>%
  filter(survey_unit == "Rockall") %>%
  summarise(
    n_rows = n(),
    n_missing = sum(is.na(LME_total_tonnes_fished)),
    pct_missing = 100 * mean(is.na(LME_total_tonnes_fished))
  )

print(rockall_check)

# ---------------------------------------------------------------------------
# 4. Verify non-Rockall surveys are not missing due to LME assignment issues
# ---------------------------------------------------------------------------

global_spatially_valid_hex_stats %>%
  filter(
    survey_unit != "Rockall",
    is.na(LME_total_tonnes_fished)
  ) %>%
  distinct(survey_unit, lme_name, year) %>%
  arrange(survey_unit, year) %>%
  print(n = Inf)

# ---------------------------------------------------------------------------
# 5. Check whether all remaining missing values are simply years > 2019
# ---------------------------------------------------------------------------

global_spatially_valid_hex_stats %>%
  filter(is.na(LME_total_tonnes_fished)) %>%
  summarise(
    n_missing = n(),
    n_after_2019 = sum(year > 2019),
    pct_after_2019 =
      100 * n_after_2019 / n_missing
  )

# ---------------------------------------------------------------------------
# 6. Confirm all pre-2020 non-Rockall records have fishing pressure
# ---------------------------------------------------------------------------

pre2020_failures <- global_spatially_valid_hex_stats %>%
  filter(
    survey_unit != "Rockall",
    year <= 2019,
    is.na(LME_total_tonnes_fished)
  )

cat(
  "\nPre-2020 non-Rockall rows still missing fishing pressure:",
  nrow(pre2020_failures),
  "\n"
)

# ---------------------------------------------------------------------------
# 7. Confirm all pre-2020 Rockall records received High Seas values
# ---------------------------------------------------------------------------

rockall_failures <- global_spatially_valid_hex_stats %>%
  filter(
    survey_unit == "Rockall",
    year <= 2019,
    is.na(LME_total_tonnes_fished)
  )

cat(
  "Pre-2020 Rockall rows still missing fishing pressure:",
  nrow(rockall_failures),
  "\n"
)

# ---------------------------------------------------------------------------
# 8. PASS / FAIL SUMMARY
# ---------------------------------------------------------------------------

if(nrow(pre2020_failures) == 0 &
   nrow(rockall_failures) == 0) {
  
  cat(
    "\nSUCCESS: All non-Rockall pre-2020 records matched an LME fishing-pressure value,\n",
    "and all pre-2020 Rockall records matched a High Seas fishing-pressure value.\n"
  )
  
} else {
  
  cat(
    "\nWARNING: Some records expected to match still have missing fishing-pressure values.\n"
  )
  
}
'
##############################################################################
# Q2 prep determine how to treat richness in model of CTI
##############################################################################
##############################################################################
# RICHNESS AS A POTENTIAL GLOBAL AND REGIONAL DRIVER OF MTC
#
# Goal:
# Determine whether species richness should be treated as:
#
# 1. A nuisance covariate
# 2. A global driver
# 3. A regional driver
#
# Richness metric:
# S_combined
#
##############################################################################

library(tidyverse)
library(glmmTMB)
library(broom.mixed)
library(performance)

##############################################################################
# DATA PREPARATION
##############################################################################

rich_dat <-
  
  global_spatially_valid_hex_stats %>%
  
  mutate(
    
    ##########################################################################
    # RESPONSE
    ##########################################################################
    
    combined_T =
      coalesce(
        sbt,
        glorys_bottomT
      ),
    
    log_MTC =
      log(
        MTC_hex_avg
      )
    
  ) %>%
  
  ##########################################################################
# WITHIN / BETWEEN SURVEY RICHNESS DECOMPOSITION
##########################################################################

group_by(
  survey_unit
) %>%
  
  mutate(
    
    RichMean =
      mean(
        S_combined,
        na.rm = TRUE
      ),
    
    RichWithin =
      S_combined -
      RichMean
    
  ) %>%
  
  ungroup() %>%
  
  ##########################################################################
# STANDARDIZATION
##########################################################################

mutate(
  
  zT =
    as.numeric(
      scale(
        combined_T
      )
    ),
  
  zDepth =
    as.numeric(
      scale(
        log1p(depth)
      )
    ),
  
  zEven =
    as.numeric(
      scale(
        evenness_combined
      )
    ),
  
  zTroph =
    as.numeric(
      scale(
        avg_troph_hex
      )
    ),
  
  zFish =
    as.numeric(
      scale(
        log1p(
          LME_total_tonnes_fished
        )
      )
    ),
  
  zRichMean =
    as.numeric(
      scale(
        RichMean
      )
    ),
  
  zRichWithin =
    as.numeric(
      scale(
        RichWithin
      )
    ),
  
  wt =
    1 /
    MTC_hex_variance
  
) %>%
  
  filter(
    
    is.finite(log_MTC),
    
    is.finite(zT),
    is.finite(zDepth),
    is.finite(zEven),
    is.finite(zTroph),
    is.finite(zFish),
    
    is.finite(zRichMean),
    is.finite(zRichWithin),
    
    is.finite(wt),
    
    wt > 0,
    
    !is.na(survey_unit)
    
  )

##############################################################################
# CAP EXTREME WEIGHTS
##############################################################################

rich_dat$wt <-
  
  pmin(
    
    rich_dat$wt,
    
    quantile(
      rich_dat$wt,
      0.99,
      na.rm = TRUE
    )
    
  )

##############################################################################
# MODEL 1
#
# NO RICHNESS
#
# BASELINE DRIVER MODEL
##############################################################################

m0 <-
  
  glmmTMB(
    
    log_MTC ~
      
      zT +
      zDepth +
      zEven +
      zTroph +
      zFish +
      
      (
        1 +
          zT +
          zDepth +
          zEven +
          zTroph +
          zFish |
          
          survey_unit
      ),
    
    data = rich_dat,
    
    weights = wt,
    
    family = gaussian()
    
  )

##############################################################################
# MODEL 2
#
# GLOBAL RICHNESS EFFECTS
##############################################################################

m1 <-
  
  glmmTMB(
    
    log_MTC ~
      
      zT +
      zDepth +
      zEven +
      zTroph +
      zFish +
      
      zRichWithin +
      zRichMean +
      
      (
        1 +
          zT +
          zDepth +
          zEven +
          zTroph +
          zFish |
          
          survey_unit
      ),
    
    data = rich_dat,
    
    weights = wt,
    
    family = gaussian()
    
  )

##############################################################################
# MODEL 3
#
# ALLOW LOCAL RICHNESS EFFECT TO VARY AMONG SURVEYS
##############################################################################

m2 <-
  
  glmmTMB(
    
    log_MTC ~
      
      zT +
      zDepth +
      zEven +
      zTroph +
      zFish +
      
      zRichWithin +
      zRichMean +
      
      (
        1 +
          zT +
          zDepth +
          zEven +
          zTroph +
          zFish +
          
          zRichWithin |
          
          survey_unit
      ),
    
    data = rich_dat,
    
    weights = wt,
    
    family = gaussian()
    
  )

##############################################################################
# MODEL COMPARISON
##############################################################################

AIC(
  
  m0,
  m1,
  m2
  
)

anova(
  
  m0,
  m1,
  m2
  
)

##############################################################################
# GLOBAL COEFFICIENTS
##############################################################################

tidy(
  m1,
  effects = "fixed",
  conf.int = TRUE
)

##############################################################################
# RANDOM EFFECT VARIANCES
##############################################################################

VarCorr(
  m2
)

##############################################################################
# MODEL FIT
##############################################################################

performance::r2(
  m0
)

performance::r2(
  m1
)

performance::r2(
  m2
)

##############################################################################
# INTERPRETATION FRAMEWORK
##############################################################################

# Scenario A
# m1 ≈ m0
#
# Richness contributes little explanatory power.
#
# Treat richness as a control variable or remove it.

# Scenario B
# m1 substantially better than m0
# m2 not better than m1
#
# Richness is a global driver.
#
# Retain RichMean and RichWithin as fixed effects.
# Do not map richness.

# Scenario C
# m2 substantially better than m1
# Random slope variance for zRichWithin is non-trivial
#
# Richness is a regional driver.
#
# Extract survey-specific richness BLUPs
# and map richness exactly as temperature,
# depth, evenness, trophic level and fishing pressure.

##############################################################################
# OPTIONAL DIAGNOSTIC:
# EXTRACT SURVEY-SPECIFIC RICHNESS EFFECTS
##############################################################################

richness_blups <-
  
  ranef(
    m2
  )$cond$survey_unit %>%
  
  as.data.frame() %>%
  
  rownames_to_column(
    "survey_unit"
  )

head(
  richness_blups
)

##############################################################################
# Q2: GLOBAL AND REGIONAL DRIVERS OF MTC
#
# Response:
#   log(MTC_hex_avg)
#
# Regional drivers:
#   Temperature
#   Depth
#   Evenness
#   Trophic level
#   Fishing pressure
#
# Biodiversity context covariates:
#   RichnessWithin
#   RichnessMean
#
# Random effects:
#   Survey intercept
#   Survey-specific slopes:
#       Temperature
#       Depth
#       Evenness
#       Trophic level
#       Fishing pressure
#
##############################################################################
#
# Survey-specific uncertainty:
#   Conditional BLUP uncertainty from ranef(condVar = TRUE)
#
##############################################################################

library(tidyverse)
library(sf)
library(glmmTMB)
library(broom.mixed)
library(patchwork)
library(rnaturalearth)

##############################################################################
# 1. DATA PREPARATION
##############################################################################

hex_dat <- global_spatially_valid_hex_stats %>%
  
  mutate(
    
    ##########################################################################
    # RESPONSE
    ##########################################################################
    
    combined_T =
      coalesce(
        sbt,
        glorys_bottomT
      ),
    
    log_MTC =
      log(MTC_hex_avg)
    
  ) %>%
  
  ##########################################################################
# WITHIN- / AMONG-SURVEY RICHNESS DECOMPOSITION
##########################################################################

group_by(survey_unit) %>%
  
  mutate(
    
    RichMean =
      mean(
        S_combined,
        na.rm = TRUE
      ),
    
    RichWithin =
      S_combined -
      RichMean
    
  ) %>%
  
  ungroup() %>%
  
  ##########################################################################
# STANDARDIZE VARIABLES
##########################################################################

mutate(
  
  zT =
    as.numeric(
      scale(combined_T)
    ),
  
  zDepth =
    as.numeric(
      scale(
        log1p(depth)
      )
    ),
  
  zEven =
    as.numeric(
      scale(
        evenness_combined
      )
    ),
  
  zTroph =
    as.numeric(
      scale(
        avg_troph_hex
      )
    ),
  
  zFish =
    as.numeric(
      scale(
        log1p(
          LME_total_tonnes_fished
        )
      )
    ),
  
  zRichWithin =
    as.numeric(
      scale(
        RichWithin
      )
    ),
  
  zRichMean =
    as.numeric(
      scale(
        RichMean
      )
    ),
  
  wt =
    1 / MTC_hex_variance
  
) %>%
  
  filter(
    
    is.finite(log_MTC),
    
    is.finite(zT),
    is.finite(zDepth),
    is.finite(zEven),
    is.finite(zTroph),
    is.finite(zFish),
    
    is.finite(zRichWithin),
    is.finite(zRichMean),
    
    is.finite(wt),
    wt > 0,
    
    !is.na(survey_unit)
    
  )

##############################################################################
# CAP EXTREME WEIGHTS
##############################################################################

hex_dat$wt <- pmin(
  hex_dat$wt,
  quantile(
    hex_dat$wt,
    0.99,
    na.rm = TRUE
  )
)

##############################################################################
# RICHNESS DIAGNOSTIC
##############################################################################

richness_variation_table <- hex_dat %>%
  
  group_by(survey_unit) %>%
  
  summarise(
    
    mean_rich =
      mean(
        S_combined,
        na.rm = TRUE
      ),
    
    sd_rich =
      sd(
        S_combined,
        na.rm = TRUE
      ),
    
    n = n(),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(sd_rich)

print(
  richness_variation_table,
  n = 60
)

##############################################################################
# 2. HIERARCHICAL RANDOM-SLOPE MODEL
##############################################################################

rs_model <- glmmTMB(
  
  log_MTC ~
    
    zT +
    zDepth +
    zEven +
    zTroph +
    zFish +
    
    zRichWithin +
    zRichMean +
    
    (
      1 +
        zT +
        zDepth +
        zEven +
        zTroph +
        zFish |
        
        survey_unit
    ),
  
  data = hex_dat,
  
  weights = wt,
  
  family = gaussian()
  
)

##############################################################################
# MODEL OUTPUT
##############################################################################

summary(rs_model)

VarCorr(rs_model)

##############################################################################
# 3. GLOBAL FIXED EFFECTS
##############################################################################

global_coefficients <- tidy(
  rs_model,
  effects = "fixed",
  conf.int = TRUE
) %>%
  
  filter(
    term != "(Intercept)"
  ) %>%
  
  select(
    term,
    estimate,
    std.error,
    conf.low,
    conf.high,
    statistic,
    p.value
  )

print(
  global_coefficients,
  n = Inf
)

richness_summary <-
  
  global_coefficients %>%
  
  filter(
    term %in% c(
      "zRichWithin",
      "zRichMean"
    )
  ) %>%
  
  mutate(
    
    interpretation = case_when(
      
      term == "zRichWithin" ~
        "Local richness effect",
      
      term == "zRichMean" ~
        "Regional richness effect"
      
    )
    
  )
print(
  richness_summary,
  n = Inf
)

##############################################################################
# 4. SURVEY-SPECIFIC BLUPS + CONDITIONAL UNCERTAINTY
##############################################################################

re_obj <- ranef(
  rs_model,
  condVar = TRUE
)

re_mat <- re_obj$cond$survey_unit

##############################################################################
# SURVEY-SPECIFIC COEFFICIENTS
##############################################################################

survey_slopes <- tibble(
  
  survey_unit =
    rownames(re_mat),
  
  zT =
    fixef(rs_model)$cond["zT"] +
    re_mat[, "zT"],
  
  zDepth =
    fixef(rs_model)$cond["zDepth"] +
    re_mat[, "zDepth"],
  
  zEven =
    fixef(rs_model)$cond["zEven"] +
    re_mat[, "zEven"],
  
  zTroph =
    fixef(rs_model)$cond["zTroph"] +
    re_mat[, "zTroph"],
  
  zFish =
    fixef(rs_model)$cond["zFish"] +
    re_mat[, "zFish"]
  
)

##############################################################################
# FIXED-EFFECT UNCERTAINTY
##############################################################################

fixed_effects <- tidy(
  rs_model,
  effects = "fixed"
) %>%
  
  filter(
    term %in% c(
      "zT",
      "zDepth",
      "zEven",
      "zTroph",
      "zFish"
    )
  ) %>%
  
  select(
    predictor = term,
    fixed_estimate = estimate,
    fixed_se = std.error
  )

##############################################################################
# CONDITIONAL BLUP UNCERTAINTY
##############################################################################

cond_var <- attr(
  re_mat,
  "condVar"
)

survey_random_se <- purrr::map_dfr(
  
  seq_len(dim(cond_var)[3]),
  
  function(i){
    
    vmat <- cond_var[, , i]
    
    tibble(
      
      survey_unit =
        rownames(re_mat)[i],
      
      predictor = c(
        "zT",
        "zDepth",
        "zEven",
        "zTroph",
        "zFish"
      ),
      
      random_se = c(
        
        sqrt(vmat[2,2]),  # zT
        
        sqrt(vmat[3,3]),  # zDepth
        
        sqrt(vmat[4,4]),  # zEven
        
        sqrt(vmat[5,5]),  # zTroph
        
        sqrt(vmat[6,6])   # zFish
        
      )
      
    )
    
  }
  
)
##############################################################################
# SURVEY COEFFICIENT TABLE WITH CONFIDENCE INTERVALS
##############################################################################

survey_coefficients_long <- survey_slopes %>%
  
  pivot_longer(
    
    cols = c(
      zT,
      zDepth,
      zEven,
      zTroph,
      zFish
    ),
    
    names_to =
      "predictor",
    
    values_to =
      "survey_coef"
    
  ) %>%
  
  left_join(
    fixed_effects,
    by = "predictor"
  ) %>%
  
  left_join(
    survey_random_se,
    by = c(
      "survey_unit",
      "predictor"
    )
  ) %>%
  
  mutate(
    
    survey_se =
      
      sqrt(
        fixed_se^2 +
          random_se^2
      ),
    
    conf_low =
      survey_coef -
      1.96 * survey_se,
    
    conf_high =
      survey_coef +
      1.96 * survey_se,
    
    significant =
      conf_low > 0 |
      conf_high < 0
    
  )

##############################################################################
# EXPORT TABLE
##############################################################################
survey_coefficients_gt <- survey_coefficients_long %>%
  rename(
    `Survey Unit` = survey_unit,
    `Predictor` = predictor,
    `Survey Coefficient` = survey_coef,
    `Global Coefficient` = fixed_estimate,
    `Global SE` = fixed_se,
    `Survey Deviation SE` = random_se,
    `Combined SE` = survey_se,
    `95% CI Lower` = conf_low,
    `95% CI Upper` = conf_high,
    `Significant` = significant
  ) %>%
  gt() %>%
  fmt_number(
    columns = c(
      `Survey Coefficient`,
      `Global Coefficient`,
      `Global SE`,
      `Survey Deviation SE`,
      `Combined SE`,
      `95% CI Lower`,
      `95% CI Upper`
    ),
    decimals = 3
  ) %>%
  fmt_missing(
    columns = everything(),
    missing_text = "—"
  )

gtsave(
  survey_coefficients_gt,
  "Survey_Driver_Effects_With_CI.html"
)

##############################################################################
# SURVEY SUMMARY TABLE
##############################################################################

survey_unit_model_table <- hex_dat %>%
  
  count(
    survey_unit,
    name = "n_hexes"
  ) %>%
  
  left_join(
    survey_slopes,
    by = "survey_unit"
  )

##############################################################################
# BUILD SURVEY POLYGONS
##############################################################################

survey_polygons <- all_hex_geometries %>%
  
  mutate(
    survey_unit =
      as.character(
        survey_unit
      ),
    hex_id =
      as.character(
        hex_id
      )
  ) %>%
  
  semi_join(
    
    global_spatially_valid_hex_stats %>%
      
      distinct(
        survey_unit,
        hex_id
      ) %>%
      
      mutate(
        hex_id =
          as.character(
            hex_id
          )
      ),
    
    by = c(
      "survey_unit",
      "hex_id"
    )
    
  ) %>%
  
  group_by(
    survey_unit
  ) %>%
  
  summarise(
    geometry =
      st_union(
        geometry
      ),
    .groups = "drop"
  ) %>%
  
  st_make_valid()

##############################################################################
# WORLD MAP
##############################################################################

world_map_plot <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
) %>%
  st_make_valid()

##############################################################################
# JOIN COEFFICIENTS TO POLYGONS
##############################################################################

survey_polygon_map <- survey_polygons %>%
  
  left_join(
    survey_coefficients_long,
    by = "survey_unit"
  )

##############################################################################
# MAP LIMITS
##############################################################################

lon_lims <- c(
  -190,
  55
)

lat_lims <- c(
  st_bbox(survey_polygon_map)$ymin - 2,
  st_bbox(survey_polygon_map)$ymax + 2
)
##############################################################################
# SHARED MAP SETTINGS
##############################################################################

predictor_labels <- c(
  zT     = "Temperature",
  zDepth = "Depth",
  zEven  = "Evenness",
  zTroph = "Average Trophic Level",
  zFish  = "Total fish landings"
)

common_map_theme <-
  
  theme_bw(base_size = 12) +
  
  theme(
    
    plot.title =
      element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      ),
    
    legend.title =
      element_text(
        face = "bold"
      )
    
  )

common_annotation_theme <-
  
  theme(
    
    plot.title =
      element_text(
        size = 32,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = -15
        )
      ),
    
    plot.tag =
      element_text(
        size = 18,
        face = "bold"
      )
    
  )

##############################################################################
# DRIVER EFFECT MAPS
##############################################################################

make_predictor_map <- function(pred_name){
  
  pred_df <-
    
    survey_polygon_map %>%
    
    filter(
      predictor == pred_name
    ) %>%
    
    arrange(
      significant
    )
  
  coef_max <-
    
    max(
      abs(pred_df$survey_coef),
      na.rm = TRUE
    )
  
  ggplot() +
    
    ##########################################################################
  # BASEMAP
  ##########################################################################
  
  geom_sf(
    data = world_map_plot,
    fill = "grey92",
    colour = "grey50",
    linewidth = 0.15
  ) +
    
    ##########################################################################
  # NON-SIGNIFICANT SURVEYS
  ##########################################################################
  
  geom_sf(
    data =
      pred_df %>%
      filter(!significant),
    fill = "grey70",
    colour = "black",
    linewidth = 0.20
  ) +
    
    ##########################################################################
  # SIGNIFICANT SURVEYS
  ##########################################################################
  
  geom_sf(
    data =
      pred_df %>%
      filter(significant),
    aes(
      fill = survey_coef
    ),
    colour = "black",
    linewidth = 0.20
  ) +
    
    scale_fill_gradient2(
      
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      
      midpoint = 0,
      
      limits = c(
        -coef_max,
        coef_max
      ),
      
      oob = scales::squish,
      
      name = expression(beta)
      
    ) +
    
    coord_sf(
      xlim = lon_lims,
      ylim = lat_lims,
      expand = FALSE
    ) +
    
    labs(
      title =
        predictor_labels[pred_name]
    ) +
    
    common_map_theme
  
}

##############################################################################
# GENERATE DRIVER MAPS
##############################################################################

map_T     <- make_predictor_map("zT")
map_Depth <- make_predictor_map("zDepth")
map_Even  <- make_predictor_map("zEven")
map_Troph <- make_predictor_map("zTroph")
map_Fish  <- make_predictor_map("zFish")

##############################################################################
# COMBINED DRIVER FIGURE
##############################################################################

final_driver_maps <-
  
  (
    map_T |
      map_Depth |
      map_Even
  ) /
  
  (
    map_Troph |
      map_Fish |
      patchwork::plot_spacer()
  ) +
  
  plot_layout(
    heights = c(
      1,
      0.95
    )
  ) +
  
  plot_annotation(
    
    title =
      "Regional Drivers of Community Thermophilization",
    
    tag_levels = "A",
    
    theme =
      common_annotation_theme
    
  )

##############################################################################
# DRIVER FIGURE SAVE
##############################################################################

ggsave(
  "Regional_Driver_Effects_ConditionalBLUP_Significance.png",
  final_driver_maps,
  width = 18,
  height = 10,
  dpi = 600,
  bg = "white"
)
##############################################################################
# SHARED LEGEND SCALE
##############################################################################

deviation_fill_scale <-
  
  scale_fill_manual(
    
    values = c(
      
      "Stronger same sign" =
        "#B2182B",
      
      "Weaker same sign" =
        "#FDBE85",
      
      "Opposite sign" =
        "#92C5DE",
      
      "Strongly opposite sign" =
        "#2166AC",
      
      "Not different" =
        "grey70"
      
    ),
    
    breaks = c(
      "Stronger same sign",
      "Weaker same sign",
      "Opposite sign",
      "Strongly opposite sign",
      "Not different"
    ),
    
    limits = c(
      "Stronger same sign",
      "Weaker same sign",
      "Opposite sign",
      "Strongly opposite sign",
      "Not different"
    ),
    
    drop = FALSE,
    
    name = "Deviation From\nGlobal Effect"
  )

##############################################################################
# DEVIATION MAPS
##############################################################################

make_deviation_map <- function(pred_name){
  
  pred_df <-
    survey_deviation_map %>%
    filter(predictor == pred_name)
  
  ggplot() +
    
    geom_sf(
      data = world_map_plot,
      fill = "grey92",
      colour = "grey50",
      linewidth = 0.15
    ) +
    
    geom_sf(
      data = pred_df,
      aes(fill = deviation_class),
      colour = "black",
      linewidth = 0.20
    ) +
    
    deviation_fill_scale +
    
    coord_sf(
      xlim = lon_lims,
      ylim = lat_lims,
      expand = FALSE
    ) +
    
    labs(
      title = predictor_labels[pred_name],
      tag = NULL
    ) +
    
    common_map_theme
}

##############################################################################
# GENERATE DEVIATION MAPS
##############################################################################

survey_deviation_long <-
  
  survey_coefficients_long %>%
  
  mutate(
    
    deviation_from_global =
      survey_coef -
      fixed_estimate,
    
    deviation_low =
      conf_low -
      fixed_estimate,
    
    deviation_high =
      conf_high -
      fixed_estimate,
    
    deviation_significant =
      !(
        deviation_low <= 0 &
          deviation_high >= 0
      ),
    
    global_sign =
      sign(fixed_estimate),
    
    survey_sign =
      sign(survey_coef),
    
    same_sign =
      survey_sign == global_sign,
    
    deviation_class =
      case_when(
        
        ####################################################################
        # SURVEY EFFECT NOT DIFFERENT FROM GLOBAL EFFECT
        ####################################################################
        
        !deviation_significant ~
          
          "Not different",
        
        ####################################################################
        # SAME SIGN AS GLOBAL EFFECT
        ####################################################################
        
        same_sign &
          abs(survey_coef) >
          abs(fixed_estimate) ~
          
          "Stronger same sign",
        
        same_sign &
          abs(survey_coef) <=
          abs(fixed_estimate) ~
          
          "Weaker same sign",
        
        ####################################################################
        # OPPOSITE SIGN AS GLOBAL EFFECT
        ####################################################################
        
        !same_sign &
          abs(survey_coef) <
          abs(fixed_estimate) ~
          
          "Opposite sign",
        
        !same_sign &
          abs(survey_coef) >=
          abs(fixed_estimate) ~
          
          "Strongly opposite sign",
        
        TRUE ~
          
          "Not different"
        
      )
    
  ) %>%
  
  mutate(
    deviation_class =
      factor(
        deviation_class,
        levels = c(
          "Stronger same sign",
          "Weaker same sign",
          "Opposite sign",
          "Strongly opposite sign",
          "Not different"
        )
      )
  )

survey_deviation_map <-
  
  survey_polygons %>%
  
  left_join(
    survey_deviation_long,
    by = "survey_unit"
  )

##############################################################################
# INDIVIDUAL MAPS
##############################################################################

dev_T_raw     <- make_deviation_map("zT")
dev_Depth_raw <- make_deviation_map("zDepth")
dev_Even_raw  <- make_deviation_map("zEven")
dev_Troph_raw <- make_deviation_map("zTroph")
dev_Fish_raw  <- make_deviation_map("zFish")

##############################################################################
# EXTRACT SINGLE LEGEND
##############################################################################

library(cowplot)

shared_legend <-
  
  cowplot::get_legend(
    
    dev_T_raw +
      
      theme(
        legend.position = "right",
        legend.box = "vertical"
      )
    
  )

##############################################################################
# REMOVE LEGENDS FROM ALL MAPS
##############################################################################

dev_T <-
  dev_T_raw +
  labs(tag = "A") +
  theme(legend.position = "none")

dev_Depth <-
  dev_Depth_raw +
  labs(tag = "B") +
  theme(legend.position = "none")

dev_Even <-
  dev_Even_raw +
  labs(tag = "C") +
  theme(legend.position = "none")

dev_Troph <-
  dev_Troph_raw +
  labs(tag = "D") +
  theme(legend.position = "none")

dev_Fish <-
  dev_Fish_raw +
  labs(tag = "E") +
  theme(legend.position = "none")

##############################################################################
# CONVERT LEGEND TO PATCHWORK ELEMENT
##############################################################################

legend_plot <-
  patchwork::wrap_elements(shared_legend)

##############################################################################
# FINAL FIGURE
##############################################################################

driver_deviation_maps <-
  
  (
    dev_T |
      dev_Depth |
      dev_Even
  ) /
  
  (
    dev_Troph |
      dev_Fish |
      legend_plot
  ) +
  
  plot_layout(
    heights = c(1, 0.95),
    widths  = c(1, 1, 0.60)
  ) +
  
  plot_annotation(
    title = "Regional Deviations From Global Thermophilization Drivers",
    theme = common_annotation_theme
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 14
    ),
    plot.tag.position = c(0.01, 0.99)
  )
##############################################################################
# SAVE
##############################################################################

ggsave(
  "Regional_Driver_Deviation_From_Global_Effect.png",
  driver_deviation_maps,
  width = 18,
  height = 10,
  dpi = 600,
  bg = "white"
)


##############################################################################
# SUMMARY TABLE
##############################################################################

deviation_summary <- survey_deviation_long %>%
  
  count(
    predictor,
    deviation_class
  ) %>%
  
  arrange(
    predictor,
    desc(n)
  )

print(
  deviation_summary,
  n = Inf
)

##############################################################################
# OUTPUT OBJECTS
##############################################################################

driver_deviation_maps
survey_deviation_long
survey_deviation_map
deviation_summary

deviation_table <- survey_deviation_long %>%
  arrange(
    predictor,
    survey_unit
  ) %>%
  gt() %>%
  tab_header(
    title = md(
      "**Survey-Specific Deviations from Global Effects**"
    )
  ) %>%
  fmt_number(
    columns = c(
      survey_coef,
      fixed_estimate,
      survey_se,
      conf_low,
      conf_high,
      deviation_from_global,
      deviation_low,
      deviation_high
    ),
    decimals = 3
  ) %>%
  cols_label(
    survey_unit = "Survey",
    predictor = "Predictor",
    survey_coef = "Survey Effect",
    fixed_estimate = "Global Effect",
    significant = "Survey Sig.",
    deviation_from_global = "Deviation",
    deviation_significant = "Deviation Sig.",
    deviation_class = "Interpretation"
  )

gtsave(
  deviation_table,
  "Survey_Deviation_Table.html"
)

############################################################
# GLOBAL VS SURVEY-LEVEL DRIVERS OF MTC CHANGE
# Includes ALL predictors from the mixed model
############################################################

library(dplyr)
library(ggplot2)

predictor_order <- c(
  "zDepth",
  "zEven",
  "zFish",
  "zT",
  "zTroph"
)

predictor_labels <- c(
  zDepth      = "Depth",
  zEven       = "Evenness",
  zFish       = "Fishing Pressure",
  zT          = "Temperature",
  zTroph      = "Trophic Level"
)

plot_dat <-
  survey_coefficients_long %>%
  filter(predictor %in% predictor_order) %>%
  mutate(
    predictor = factor(
      predictor,
      levels = predictor_order
    )
  )

##############################################################################
# GLOBAL SIGNIFICANCE LABELS
##############################################################################

global_dat <-
  global_coefficients %>%
  rename(
    predictor   = term,
    global_coef = estimate,
    global_se   = std.error
  ) %>%
  filter(predictor %in% predictor_order) %>%
  mutate(
    predictor = factor(
      predictor,
      levels = predictor_order
    ),
    
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

##############################################################################
# LEGEND OBJECTS
##############################################################################

survey_legend_dat <- tibble(
  predictor = factor("zDepth", levels = predictor_order),
  x = -10,
  sig_class = c(
    "Not significant",
    "Significant"
  )
)

global_legend_dat <- tibble(
  predictor = factor("zDepth", levels = predictor_order),
  x = -10,
  global_class = c(
    "Global coefficient (blue point ±95% CI)",
    "Global significance (*, **, ***)"
  )
)

##############################################################################
# PLOT
##############################################################################

global_vs_survey_drivers <-
  
  ggplot(
    plot_dat,
    aes(
      x = survey_coef,
      y = predictor
    )
  ) +
  
  geom_vline(
    xintercept = 0,
    colour = "red",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  ##########################################################################
# BOXPLOTS
##########################################################################

geom_boxplot(
  aes(group = predictor),
  fill = "grey90",
  colour = "grey40",
  width = 0.45,
  outlier.shape = NA
) +
  
  ##########################################################################
# SURVEY COEFFICIENTS
##########################################################################

geom_point(
  aes(fill = significant),
  shape = 21,
  colour = "grey45",
  stroke = 0.4,
  size = 2.2,
  position = position_jitter(
    width = 0,
    height = 0.10
  )
) +
  
  ##########################################################################
# GLOBAL 95% CI
##########################################################################

geom_errorbar(
  data = global_dat,
  inherit.aes = FALSE,
  aes(
    y = predictor,
    xmin = global_coef - 1.96 * global_se,
    xmax = global_coef + 1.96 * global_se
  ),
  orientation = "y",
  height = 0.12,
  linewidth = 1.2,
  colour = "navy"
) +
  
  ##########################################################################
# GLOBAL POINT
##########################################################################

geom_point(
  data = global_dat,
  inherit.aes = FALSE,
  aes(
    x = global_coef,
    y = predictor
  ),
  shape = 21,
  size = 5,
  fill = "navy",
  colour = "white",
  stroke = 1.4
) +
  
  ##########################################################################
# GLOBAL SIGNIFICANCE STARS
##########################################################################

geom_text(
  data = global_dat,
  inherit.aes = FALSE,
  aes(
    x = global_coef,
    y = predictor,
    label = stars
  ),
  colour = "navy",
  fontface = "bold",
  size = 5,
  nudge_x = 0.06
) +
  
  ##########################################################################
# SURVEY LEGEND
##########################################################################

geom_point(
  data = survey_legend_dat,
  inherit.aes = FALSE,
  aes(
    x = x,
    y = predictor,
    fill = sig_class
  ),
  shape = 21,
  colour = "grey45",
  size = 3
) +
  
  ##########################################################################
# GLOBAL LEGEND
##########################################################################

geom_point(
  data = global_legend_dat %>%
    filter(global_class ==
             "Global coefficient (blue point ±95% CI)"),
  inherit.aes = FALSE,
  aes(
    x = x,
    y = predictor,
    shape = global_class
  ),
  colour = "navy",
  fill = "navy",
  size = 4
) +
  
  geom_point(
    data = global_legend_dat %>%
      filter(global_class ==
               "Global significance (*, **, ***)"),
    inherit.aes = FALSE,
    aes(
      x = x,
      y = predictor,
      shape = global_class
    ),
    colour = "navy",
    size = 4
  ) +
  
  ##########################################################################
# AXES
##########################################################################

scale_y_discrete(
  labels = predictor_labels
) +
  
  scale_x_continuous(
    limits = c(-0.5, 0.7),
    breaks = seq(-0.5, 0.7, by = 0.1)
  ) +
  
  ##########################################################################
# SURVEY LEGEND
##########################################################################

scale_fill_manual(
  name = "Survey effect",
  values = c(
    "TRUE" = "grey55",
    "FALSE" = NA,
    "Significant" = "grey55",
    "Not significant" = NA
  ),
  breaks = c(
    "Not significant",
    "Significant"
  )
) +
  
  ##########################################################################
# GLOBAL LEGEND
##########################################################################

scale_shape_manual(
  name = "Global effect",
  values = c(
    "Global coefficient (blue point ±95% CI)" = 16,
    "Global significance (*, **, ***)" = 8
  )
) +
  
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        size = 4,
        colour = "grey45"
      )
    ),
    shape = guide_legend(
      order = 2
    )
  ) +
  
  labs(
    x = "Standardized coefficient estimate",
    y = "CTI(hex) predictor"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    )
  )

print(global_vs_survey_drivers)

ggsave(
  "Global_vs_Survey_Driver_Effects.png",
  global_vs_survey_drivers,
  width = 11,
  height = 6.5,
  dpi = 600,
  bg = "white"
)


####################################################################################################################################################################################
######## Q3 Attempt, what species drive thermophilization of survey units, what families are dominant in these thermal classes, and what are their traits? ##
####################################################################################################################################################################################
global_spatially_valid_hex_stats$hex_id <- as.character(global_spatially_valid_hex_stats$hex_id)

temp_drivers_dataset <- Combined_FISH_optimized_spatially_valid %>%
  filter(!is.na(estim_thermal_mean) & wgt_combined > 0) %>%
  mutate(hexagon_id = as.character(hexagon_id)) %>%
  left_join(
    global_spatially_valid_hex_stats %>% 
      select(hex_id, year, survey_unit, MTC_hex_avg),
    by = c("hexagon_id" = "hex_id", "year", "survey_unit")
  )
'
Step 1: Calculate Species-Level Annual ContributionsFirst, we calculate how much each species pulls or pushes the $MTC$ away from the hexagons overall mean in any given year.Thermophilic Drivers (Warm Class): Species where $STI > \text{Hex Average } MTC$ whose biomass or contribution is increasing over time.Cold-Water Drivers (Cool Class): Species where $STI < \text{Hex Average } MTC$ whose biomass or contribution is increasing over time (or warm species that are collapsing).Here is the code to compute the annual baseline metrics:
'
library(dplyr)
library(tidyr)
library(purrr)

# 1. Calculate annual relative weights and absolute thermal contributions
species_contributions <- temp_drivers_dataset %>%
  group_by(hexagon_id, survey_unit, year) %>%
  mutate(
    # Total biomass in this hex/survey/year
    total_wgt = sum(wgt_combined, na.rm = TRUE),
    # Relative biomass weight of this species
    rel_wgt = wgt_combined / total_wgt,
    # This species' raw contribution to that year's MTC
    species_mtc_contrib = rel_wgt * estim_thermal_mean
  ) %>%
  ungroup() %>%
  # Bring in the baseline long-term average MTC for the Hex/Survey to classify species
  group_by(hexagon_id, survey_unit) %>%
  mutate(long_term_hex_MTC = mean(MTC_hex_avg, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    # Is the species thermal affinity above or below the long-term hex average?
    thermal_class = ifelse(estim_thermal_mean > long_term_hex_MTC, "Warm-Affinity", "Cold-Affinity")
  )


#Step 2: Estimate the Thermophilization Impact Over TimeTo find out which species are driving changes over time, we fit a linear regression ($y \sim \text{year}$) to each species annual thermal contribution (species_mtc_contrib) within each unique spatial unit (hexagon_id + survey_unit).The slope ($\beta$) of this line represents the annual thermophilization impact of that species.A positive slope means the species is actively driving the $MTC$ upward over time.A negative slope means the species is actively dragging the $MTC$ downward over time.

# 2. Estimate temporal trends per species per spatial unit
species_trends <- species_contributions %>%
  # Filter out species that only appear once or twice in a hex (need data to trend)
  group_by(hexagon_id, survey_unit, aphia_id, accepted_name, thermal_class) %>%
  filter(n_distinct(year) >= 3) %>% 
  summarise(
    # Slopes of relative weight and actual MTC impact over time
    wgt_trend = list(lm(rel_wgt ~ year)),
    mtc_impact_trend = list(lm(species_mtc_contrib ~ year)),
    .groups = "drop"
  ) %>%
  # Extract the slope coefficients using purrr
  mutate(
    biomass_shift_rate = map_dbl(wgt_trend, ~ coef(.x)["year"]),
    thermophilization_impact = map_dbl(mtc_impact_trend, ~ coef(.x)["year"])
  ) %>%
  select(-wgt_trend, -mtc_impact_trend)

# 3. Classify the drivers of directional MTC shifts
classified_drivers <- species_trends %>%
  mutate(
    driver_role = case_when(
      # Warm-water species expanding OR Cold-water species collapsing
      thermophilization_impact > 0 ~ "Warming Driver (Pushes MTC Up)",
      
      # Cold-water species expanding OR Warm-water species collapsing
      thermophilization_impact < 0 ~ "Cooling Driver (Pulls MTC Down)",
      
      TRUE ~ "Neutral / No Change"
    )
  ) %>%
  arrange(hexagon_id, survey_unit, desc(abs(thermophilization_impact)))

# 1. Isolate unique hex-to-LME mappings by taking the most frequent match
hex_lme_mapping <- Combined_FISH_optimized_spatially_valid %>%
  mutate(
    hexagon_id = as.character(hexagon_id)
  ) %>%
  count(hexagon_id, survey_unit, LME_NUMBER, LME_NAME) %>%
  group_by(hexagon_id, survey_unit) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>% # Keeps the dominant LME if a hex spans two
  ungroup() %>%
  select(-n)

# 2. Re-run your join safely
classified_drivers <- classified_drivers %>%
  left_join(hex_lme_mapping, by = c("hexagon_id", "survey_unit"))


'
To roll these hyper-local (hexagon-level) results up into regional, survey-wide impacts, you need to aggregate the species trends across all hexes within each survey_unit.

When averaging the thermophilization_impact across hexes, it is highly recommended to compute both the mean impact (to see the overall direction and magnitude) and the consistency/prevalence (in how many hexagons that species is actively driving change). This prevents a species that has a massive biomass spike in just one single hexagon from skewing your entire regional interpretation.

Here is the code to collate your results by survey_unit:
'

regional_species_impact <- classified_drivers %>%
  group_by(survey_unit, aphia_id, accepted_name, thermal_class) %>%
  summarise(
    # 1. Total number of hexes where this species was analyzed
    n_hexes_present = n(),
    
    # 2. Average thermophilization impact across all those hexes
    mean_regional_impact = mean(thermophilization_impact, na.rm = TRUE),
    
    # 3. Average rate of relative biomass change across the region
    mean_biomass_shift_rate = mean(biomass_shift_rate, na.rm = TRUE),
    
    # 4. Count how often they act as a warming vs cooling driver
    hexes_driving_warming = sum(driver_role == "Warming Driver (Pushes MTC Up)"),
    hexes_driving_cooling = sum(driver_role == "Cooling Driver (Pulls MTC Down)"),
    
    .groups = "drop"
  ) %>%
  # Calculate prevalence percentage to see how widespread their impact is
  mutate(
    regional_consistency_pct = (hexes_driving_warming + hexes_driving_cooling) / n_hexes_present * 100
  ) %>%
  # Rank by absolute regional impact to find the top heavy-hitters per survey unit
  group_by(survey_unit) %>%
  arrange(survey_unit, desc(abs(mean_regional_impact)))

library(dplyr)
library(tidyr)

# 1. Apply the 4-quadrant ecological labels to the regional dataset
labeled_regional_impact <- regional_species_impact %>%
  mutate(
    ecological_label = case_when(
      thermal_class == "Warm-Affinity" & mean_biomass_shift_rate > 0  ~ "Warm Expansion",
      thermal_class == "Cold-Affinity" & mean_biomass_shift_rate <= 0 ~ "Cold Contraction",
      thermal_class == "Cold-Affinity" & mean_biomass_shift_rate > 0  ~ "Cold Expansion",
      thermal_class == "Warm-Affinity" & mean_biomass_shift_rate <= 0 ~ "Warm Contraction",
      TRUE ~ "Stable / Neutral"
    ),
    # Define the high-level impact direction
    impact_direction = ifelse(mean_regional_impact > 0, "Warming Contribution", "Cooling Contribution")
  )

# 2. Extract Top 10 species pushing MTC UP (Warming) per survey unit
top_10_warming <- labeled_regional_impact %>%
  filter(impact_direction == "Warming Contribution") %>%
  group_by(survey_unit) %>%
  slice_max(order_by = mean_regional_impact, n = 10, with_ties = FALSE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# 3. Extract Top 10 species pulling MTC DOWN (Cooling) per survey unit
top_10_cooling <- labeled_regional_impact %>%
  filter(impact_direction == "Cooling Contribution") %>%
  group_by(survey_unit) %>%
  # Use abs() because we want the largest negative values
  slice_max(order_by = abs(mean_regional_impact), n = 10, with_ties = FALSE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# 4. Bind them into a final master table
final_top_drivers_table <- bind_rows(top_10_warming, top_10_cooling) %>%
  select(
    survey_unit, 
    impact_direction, 
    rank, 
    accepted_name, 
    aphia_id, 
    ecological_label, 
    mean_regional_impact, 
    mean_biomass_shift_rate,
    n_hexes_present
  ) %>%
  arrange(survey_unit, impact_direction, rank)

cold_hot_expansion_contraction <- final_top_drivers_table # TABLE THAT ANSWERS PART 1 OF H3, survey unit level identification of top temp infuencing taxa
print(cold_hot_expansion_contraction, n = 100)

#####################################################################
# PART 2: Which fish families dominate each thermal class?
#####################################################################
# Total sampled hexes per survey unit
survey_hex_totals <- classified_drivers %>%
  distinct(survey_unit, hexagon_id) %>%
  count(survey_unit, name = "total_hexes")

# Species-family lookup
species_family_lookup <- Combined_FISH_optimized_spatially_valid %>%
  select(
    aphia_id,
    accepted_name,
    family
  ) %>%
  distinct()

# Add family and occupancy information
family_species_impact <- labeled_regional_impact %>%
  left_join(
    species_family_lookup,
    by = c("aphia_id", "accepted_name")
  ) %>%
  left_join(
    survey_hex_totals,
    by = "survey_unit"
  ) %>%
  mutate(
    pct_hexes_present =
      100 * n_hexes_present / total_hexes
  )

dominant_families <- family_species_impact %>%
  group_by(
    survey_unit,
    ecological_label,
    family
  ) %>%
  summarise(
    total_family_impact =
      sum(abs(mean_regional_impact), na.rm = TRUE),
    
    n_species = n(),
    
    .groups = "drop"
  ) %>%
  group_by(
    survey_unit,
    ecological_label
  ) %>%
  slice_max(
    order_by = total_family_impact,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

dominant_family_table <- family_species_impact %>%
  inner_join(
    dominant_families %>%
      select(
        survey_unit,
        ecological_label,
        family
      ),
    by = c(
      "survey_unit",
      "ecological_label",
      "family"
    )
  ) %>%
  group_by(
    survey_unit,
    ecological_label
  ) %>%
  slice_max(
    order_by = abs(mean_regional_impact),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    survey_unit,
    ecological_label,
    family,
    aphia_id,
    accepted_name,
    mean_regional_impact,
    mean_biomass_shift_rate,
    pct_hexes_present
  )

gtsave(
  dominant_family_table %>%
    gt() %>%
    tab_header(
      title = md("**Dominant Families Driving Thermophilization by Survey Unit**")
    ) %>%
    fmt_number(
      columns = c(mean_regional_impact, mean_biomass_shift_rate, pct_hexes_present),
      decimals = 3
    ) %>%
    cols_label(
      survey_unit = "Survey Unit",
      ecological_label = "Thermal Class",
      family = "Dominant Family",
      aphia_id = "Aphia ID",
      accepted_name = "Most influential species member",
      mean_regional_impact = "Mean Regional CTI Impact",
      mean_biomass_shift_rate = "Mean Biomass Shift Rate",
      pct_hexes_present = "% Hexes Present"
    ),
  "Dominant_Families_Thermal_Classes.html"
)

family_taxonomy_lookup <- Combined_FISH_optimized_spatially_valid %>%
  select(
    family,
    order,
    class
  ) %>%
  distinct()
dominant_families_taxonomy <- dominant_families %>%
  left_join(
    family_taxonomy_lookup,
    by = "family"
  )
dominant_taxonomic_diversity <- dominant_families_taxonomy %>%
  
  summarise(
    
    n_unique_families =
      n_distinct(family),
    
    n_unique_orders =
      n_distinct(order),
    
    n_unique_classes =
      n_distinct(class)
    
  )

dominant_taxonomic_diversity
print(unique(dominant_families_taxonomy$class))

############################################################
# DOMINANT FAMILY COMPOSITION BY THERMAL CLASS
# TOP 6 FAMILIES + ALL OTHERS
# Ordered so:
#   - All Others centered at bottom
#   - Largest family upper-left
#   - Remaining families proceed clockwise by dominance
# Labels shown on every slice
############################################################

library(dplyr)
library(ggplot2)

thermal_levels <- c(
  "Warm Expansion",
  "Warm Contraction",
  "Cold Expansion",
  "Cold Contraction"
)############################################################
# TOP 4 FAMILIES + ALL OTHERS
############################################################

pie_dat <-
  
  dominant_family_table %>%
  
  filter(
    ecological_label %in% thermal_levels,
    !is.na(family)
  ) %>%
  
  count(
    ecological_label,
    family,
    name = "n_surveys"
  ) %>%
  
  group_by(ecological_label) %>%
  
  mutate(
    pct = 100 * n_surveys / sum(n_surveys)
  ) %>%
  
  arrange(
    ecological_label,
    desc(pct)
  ) %>%
  
  mutate(
    family_plot = if_else(
      row_number() <= 4,
      family,
      "All Others"
    )
  ) %>%
  
  group_by(
    ecological_label,
    family_plot
  ) %>%
  
  summarise(
    pct = sum(pct),
    .groups = "drop"
  ) %>%
  
  group_by(ecological_label) %>%
  
  mutate(
    pct = 100 * pct / sum(pct),
    
    label = paste0(
      family_plot,
      "\n",
      round(pct,1),
      "%"
    )
  ) %>%
  
  ungroup()

############################################################
# FACET TITLES
############################################################

facet_titles <-
  
  pie_dat %>%
  
  filter(
    family_plot != "All Others"
  ) %>%
  
  group_by(ecological_label) %>%
  
  slice_max(
    pct,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  mutate(
    facet_label =
      paste0(
        ecological_label,
        "\nDominant: ",
        family_plot,
        " (",
        round(pct,1),
        "% of survey units)"
      )
  ) %>%
  
  select(
    ecological_label,
    facet_label
  )

pie_dat <-
  
  left_join(
    pie_dat,
    facet_titles,
    by = "ecological_label"
  )

############################################################
# PIE FIGURE
############################################################

family_class_pies <-
  
  ggplot(
    pie_dat,
    aes(
      x = "",
      y = pct,
      fill = family_plot
    )
  ) +
  
  geom_col(
    width = 1,
    colour = "white",
    linewidth = 0.6
  ) +
  
  coord_polar(
    theta = "y"
  ) +
  
  facet_wrap(
    ~ facet_label,
    ncol = 2
  ) +
  
  geom_text(
    aes(
      label = label
    ),
    position = position_stack(
      vjust = 0.5
    ),
    size = 4,
    fontface = "bold"
  ) +
  
  labs(
    title =
      "Dominant Fish Families Across Thermophilization Classes",
    
    subtitle =
      "Top four dominant families per thermal class; remaining families pooled as 'All Others'"
  ) +
  
  theme_void(base_size = 14) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      ),
    
    plot.subtitle =
      element_text(
        hjust = 0.5
      ),
    
    strip.text =
      element_text(
        face = "bold",
        size = 12
      ),
    
    legend.position = "none"
    
  )

print(family_class_pies)

ggsave(
  "Thermal_Class_Dominant_Family_Pies.png",
  family_class_pies,
  width = 14,
  height = 16,
  dpi = 300
)

##############################################################################################
### H3: Retry trait analysis comparing thermal impacters to rest of taxa in each survey unit #####
##############################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)

# Helper function to find the most common categorical trait value (mode)
get_mode <- function(x) {
  x <- na.omit(x)
  if(length(x) == 0) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

################################################################################
### STEP 1: Process Survey Unit Level Driver Identifications
################################################################################

# Establish baseline survey-unit specific driver assignments 
survey_unit_species_classes <- labeled_regional_impact %>%
  mutate(
    archetype_class = case_when(
      thermal_class == "Warm-Affinity" & mean_biomass_shift_rate > 0  ~ "Warm Expansion",
      thermal_class == "Cold-Affinity" & mean_biomass_shift_rate <= 0 ~ "Cold Contraction",
      thermal_class == "Cold-Affinity" & mean_biomass_shift_rate > 0  ~ "Cold Expansion",
      thermal_class == "Warm-Affinity" & mean_biomass_shift_rate <= 0 ~ "Warm Contraction",
      TRUE ~ "Stable / Neutral"
    )
  ) %>%
  select(survey_unit, aphia_id, accepted_name, archetype_class)

################################################################################
### STEP 2: Compile Local Communities & Standardize Traits Within Survey Units
################################################################################

# 1. Pull unique species present in each survey unit from raw data
survey_communities <- Combined_FISH_optimized_spatially_valid %>%
  filter(!is.na(estim_thermal_mean) & wgt_combined > 0) %>%
  distinct(survey_unit, aphia_id, accepted_name, 
           age.maturity, length.max, growth.coefficient, Troph, AR, fecundity)

# 2. Join the classification labels. Species not classified become the "Remaining Community"
survey_community_classified <- survey_communities %>%
  left_join(survey_unit_species_classes, by = c("survey_unit", "aphia_id", "accepted_name")) %>%
  mutate(archetype_class = replace_na(archetype_class, "Remaining Community"))

# 3. Standardize continuous traits *WITHIN* each individual survey unit (Z-scoring)
# This removes biogeographic bias and makes across-survey comparisons valid.
survey_community_standardized <- survey_community_classified %>%
  group_by(survey_unit) %>%
  mutate(across(
    c(age.maturity, length.max, growth.coefficient, Troph, AR, fecundity),
    ~ as.vector(scale(.)),
    .names = "z_{.col}"
  )) %>%
  ungroup()

################################################################################
### STEP 3: Compute Marginal Deviations Against the Entire Community Pool
################################################################################

library(dplyr)
library(tidyr)
library(rstatix)
library(multcompView)

# 1. Calculate the baseline profile of the ENTIRE community (all species) per survey unit
entire_community_profiles <- survey_community_standardized %>%
  group_by(survey_unit) %>%
  summarise(across(starts_with("z_"), ~ median(., na.rm = TRUE), .names = "total_{.col}"), .groups = "drop")

# 2. Calculate the median standardized profile for the 4 target thermal classes per survey unit
thermal_class_profiles <- survey_community_standardized %>%
  filter(!archetype_class %in% c("Stable / Neutral", "Remaining Community")) %>%
  group_by(survey_unit, archetype_class) %>%
  summarise(across(starts_with("z_"), ~ median(., na.rm = TRUE)), .groups = "drop")

# 3. Calculate True Marginal Deviation: (Thermal Group Median - Total Community Median)
marginal_community_deviations <- thermal_class_profiles %>%
  inner_join(entire_community_profiles, by = "survey_unit") %>%
  mutate(
    `Age Maturity`     = z_age.maturity - total_z_age.maturity,
    `Length Max`       = z_length.max - total_z_length.max,
    `Growth Coef`      = z_growth.coefficient - total_z_growth.coefficient,
    `Trophic Level`    = z_Troph - total_z_Troph,
    `Aspect Ratio`     = z_AR - total_z_AR,
    `Fecundity`        = z_fecundity - total_z_fecundity
  ) %>%
  select(survey_unit, archetype_class, `Age Maturity`:`Fecundity`)

# 4. Melt data into long format for unified plotting and statistical groupings
plot_marginal_community <- marginal_community_deviations %>%
  pivot_longer(
    cols = c(`Age Maturity`:`Fecundity`),
    names_to = "Trait",
    values_to = "Marginal_Deviation"
  )

# 5. Run Pairwise Wilcoxon tests between thermal classes to group them by lack of difference
pairwise_results <- plot_marginal_community %>%
  group_by(Trait) %>%
  pairwise_wilcox_test(
    Marginal_Deviation ~ archetype_class, 
    p.adjust.method = "fdr"
  )

# 6. Build the Compact Letter Display (CLD)
cld_groupings <- pairwise_results %>%
  group_by(Trait) %>%
  do({
    p_vector <- .$p.adj
    names(p_vector) <- paste0(.$group1, "-", .$group2)
    
    levs <- multcompView::multcompLetters(p_vector)
    
    data.frame(
      archetype_class = names(levs$Letters),
      group_letter = levs$Letters
    )
  }) %>%
  ungroup()

# 🌟 NEW: Test if each group differs significantly from the Whole Community baseline (0)
# If p >= 0.05, we append ", 0" to show it's statistically tied with the baseline.
baseline_significance <- plot_marginal_community %>%
  group_by(Trait, archetype_class) %>%
  summarise(
    p_vs_baseline = suppressWarnings(wilcox.test(Marginal_Deviation, mu = 0)$p.value),
    .groups = "drop"
  )

# Merge baseline tests with pairwise letter groupings
cld_groupings <- cld_groupings %>%
  left_join(baseline_significance, by = c("Trait", "archetype_class")) %>%
  mutate(
    group_letter = ifelse(p_vs_baseline >= 0.05, paste0(group_letter, ", 0"), group_letter)
  )

################################################################################
### STEP 4: Plot Marginal Multi-Boxplots with Statistical Groupings
################################################################################

library(ggplot2)

# 1. Dynamically position labels above the upper whisker of each box to prevent overlap/clipping
label_positions <- plot_marginal_community %>%
  group_by(Trait, archetype_class) %>%
  summarise(
    # Calculates the upper whisker position (Q3 + 1.5 * IQR) bounded by max data point
    q3 = quantile(Marginal_Deviation, 0.75, na.rm = TRUE),
    iqr = IQR(Marginal_Deviation, na.rm = TRUE),
    max_val = max(Marginal_Deviation, na.rm = TRUE),
    upper_whisker = min(max_val, q3 + 1.5 * iqr),
    y_pos = upper_whisker + 0.06, # Lift just slightly above the whisker line
    .groups = "drop"
  )

cld_plot_data <- cld_groupings %>%
  left_join(label_positions, by = c("Trait", "archetype_class"))

# 2. Build the final visualization
marginal_boxplot_figure <- ggplot(plot_marginal_community, 
                                  aes(x = Trait, y = Marginal_Deviation, fill = archetype_class)) +
  
  # The essential reference anchor line: 0 means no deviation from the overall ecosystem pool
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  
  geom_boxplot(
    outlier.size = 1,
    alpha = 0.85,
    position = position_dodge(width = 0.8)
  ) +
  
  # ✅ Add sharing-letters layer. Shared letter = No difference. Includes ", 0" if tied with baseline.
  geom_text(
    data = cld_plot_data,
    aes(x = Trait, y = y_pos, label = group_letter, group = archetype_class),
    position = position_dodge(width = 0.8),
    size = 3.8,
    vjust = 0,
    fontface = "bold",
    color = "grey15"
  ) +
  
  scale_fill_brewer(palette = "Set1") +
  
  # Expanded upper limit slightly to 0.85 to allow breathing room for labels sitting over high whiskers
  coord_cartesian(ylim = c(-0.5, 0.85)) +
  scale_y_continuous(breaks = seq(-0.5, 0.75, by = 0.25)) +
  
  labs(
    title = "Marginal Trait Deviations Relative to Entire Survey Unit Community",
    subtitle = "Baseline (0) represents the full community median. Shared letters = groups match. Group includes '0' = group matches the full community baseline (FDR p >= 0.05).",
    x = "Functional Trait",
    y = "Marginal Deviation (Group Median - Total Community Median)",
    fill = "Thermal Class"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey30", size = 11, face = "italic")
  )

print(marginal_boxplot_figure)


###########################################################################
###### H3: Survey Unit Level Dynamic Animation (Full Boxplots) ############
###########################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(multcompView)
library(gganimate)
library(gifski)

################################################################################
### STEP 1: Compute Species-Level Marginal Deviations per Survey Unit
################################################################################

# 1. Melt data into long format at the individual species level
species_long_traits <- survey_community_standardized %>%
  select(survey_unit, archetype_class, aphia_id, accepted_name, starts_with("z_")) %>%
  pivot_longer(
    cols = starts_with("z_"),
    names_to = "Trait",
    values_to = "Z_Score"
  ) %>%
  mutate(
    Trait = gsub("z_", "", Trait),
    Trait = gsub("_|\\.", " ", Trait)
  )

# 2. Calculate the baseline median profile of the ENTIRE community pool per survey unit
su_community_medians <- species_long_traits %>%
  group_by(survey_unit, Trait) %>%
  summarise(su_total_median = median(Z_Score, na.rm = TRUE), .groups = "drop")

# 3. Compute species-level marginal deviations (Species Z-Score - SU Community Median)
# This retains all species data points so gganimate draws actual boxplots per frame!
species_marginal_deviations <- species_long_traits %>%
  filter(!archetype_class %in% c("Stable / Neutral", "Remaining Community")) %>%
  inner_join(su_community_medians, by = c("survey_unit", "Trait")) %>%
  mutate(Marginal_Deviation = Z_Score - su_total_median)

################################################################################
### STEP 2: Compute Survey Unit Pairwise & Baseline Groupings for Animation
################################################################################

# Get list of survey units that have a valid comparative framework
valid_su_list <- species_marginal_deviations %>%
  group_by(survey_unit) %>%
  summarise(n_classes = n_distinct(archetype_class), .groups = "drop") %>%
  filter(n_classes >= 2) %>%
  pull(survey_unit)

# 1. Generate Compact Letter Displays dynamically across all valid Survey Units and Traits
animated_cld_groupings <- species_marginal_deviations %>%
  filter(survey_unit %in% valid_su_list) %>%
  group_split(survey_unit, Trait) %>%
  purrr::map_df(function(sub_df) {
    current_su <- sub_df$survey_unit[1]
    current_trait <- sub_df$Trait[1]
    classes_present <- unique(sub_df$archetype_class)
    
    # Ensure groups have enough distinct data to test pairwise variations
    valid_groups <- sub_df %>%
      group_by(archetype_class) %>%
      summarise(n = sum(!is.na(Marginal_Deviation)), .groups = "drop") %>%
      filter(n >= 2)
    
    if(nrow(valid_groups) < 2) {
      return(data.frame(survey_unit = current_su, Trait = current_trait, archetype_class = classes_present, group_letter = "a"))
    }
    
    tryCatch({
      pairwise <- pairwise_wilcox_test(sub_df, Marginal_Deviation ~ archetype_class, p.adjust.method = "fdr")
      p_vector <- pairwise$p.adj
      names(p_vector) <- paste0(pairwise$group1, "-", pairwise$group2)
      
      levs <- multcompView::multcompLetters(p_vector)
      data.frame(survey_unit = current_su, Trait = current_trait, archetype_class = names(levs$Letters), group_letter = as.character(levs$Letters))
    }, error = function(e) {
      data.frame(survey_unit = current_su, Trait = current_trait, archetype_class = classes_present, group_letter = "a")
    })
  })

# 2. Test if each archetype significantly deviates from 0 (the total local baseline)
animated_baseline_sig <- species_marginal_deviations %>%
  filter(survey_unit %in% valid_su_list) %>%
  group_by(survey_unit, Trait, archetype_class) %>%
  summarise(
    p_vs_baseline = suppressWarnings(tryCatch({ wilcox.test(Marginal_Deviation, mu = 0)$p.value }, error = function(e) { 1 })),
    .groups = "drop"
  )

# 3. Combine baseline validations and append the ", 0" marker flag to shared categories
animated_cld_final <- animated_cld_groupings %>%
  left_join(animated_baseline_sig, by = c("survey_unit", "Trait", "archetype_class")) %>%
  mutate(
    group_letter = ifelse(p_vs_baseline >= 0.05, paste0(group_letter, ", 0"), group_letter)
  )

# 4. Dynamically anchor text label vertical positions based on local upper boxplot whiskers
animated_label_positions <- species_marginal_deviations %>%
  filter(survey_unit %in% valid_su_list) %>%
  group_by(survey_unit, Trait, archetype_class) %>%
  summarise(
    q3 = quantile(Marginal_Deviation, 0.75, na.rm = TRUE),
    iqr = IQR(Marginal_Deviation, na.rm = TRUE),
    max_val = max(Marginal_Deviation, na.rm = TRUE),
    upper_whisker = min(max_val, q3 + 1.5 * iqr),
    y_pos = upper_whisker + 0.06,
    .groups = "drop"
  )

animated_annotations <- animated_cld_final %>%
  left_join(animated_label_positions, by = c("survey_unit", "Trait", "archetype_class"))

# 5. Isolate final clean subset data to feed directly into the mapping engine
final_plot_data <- species_marginal_deviations %>%
  filter(survey_unit %in% valid_su_list)

################################################################################
### STEP 3: Build and Wrangle the Animated ggplot Architecture
################################################################################

su_animated_plot <- ggplot(final_plot_data, aes(x = Trait, y = Marginal_Deviation, fill = archetype_class)) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  
  geom_boxplot(
    outlier.size = 1, 
    alpha = 0.85, 
    position = position_dodge(width = 0.8)
  ) +
  
  # ✅ Dynamic Letter Annotation Layer — tracks the transition states flawlessly
  geom_text(
    data = animated_annotations,
    aes(x = Trait, y = y_pos, label = group_letter, group = archetype_class),
    position = position_dodge(width = 0.8),
    size = 3.8,
    vjust = 0,
    fontface = "bold",
    color = "grey15",
    na.rm = TRUE
  ) +
  
  scale_fill_brewer(palette = "Set1") +
  
  # Coordinate rules optimized across the combined dataset distribution range
  coord_cartesian(ylim = c(-0.6, 0.85)) +
  scale_y_continuous(breaks = seq(-0.5, 0.75, by = 0.25)) +
  
  labs(
    title = "Marginal Trait Deviations Across Thermal Classes",
    subtitle = "Survey Unit Location Context: {next_state}\n0 line = total local community median. Shared letters = groups match. Includes '0' = matches background baseline (FDR p >= 0.05).",
    x = "Functional Trait",
    y = "Species Marginal Deviation (Species Z-Score - SU Community Median)",
    fill = "Thermal Class Archetype"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "darkblue", size = 10.5, face = "italic")
  ) +
  
  # Instruct gganimate to loop through survey units sequentially
  transition_states(
    states = survey_unit,
    transition_length = 2,
    state_length = 5
  ) +
  enter_fade() +
  exit_fade()

################################################################################
### STEP 4: Render and Save the GIF
################################################################################

num_sus <- length(valid_su_list)

# Render block configuration
rendered_su_gif <- animate(
  su_animated_plot, 
  nframes = num_sus * 35, # Generates adequate frame counts per state change
  fps = 10,              # Display frame rendering velocity
  width = 900, 
  height = 750, 
  renderer = gifski_renderer()
)

# Render view in IDE panel
print(rendered_su_gif)

# Save tracking animation binary directly onto local working space
anim_save("Survey_Units_Trait_Anharmonic_Deviations.gif", animation = rendered_su_gif)

###############################################################################
# H3: THERMAL CLASS TRAIT PROFILES USING SPECIES-LEVEL MARGINAL DEVIATIONS
###############################################################################

'
Method: Question
Multinomial model: Which traits increase or decrease odds of class membership?
Random Forest: Which traits are most useful for predicting class membership or impact?
Indicator Trait Analysis: Which trait states are uniquely characteristic of each class? **** NOTE THIS ONLY HELPS EXPLAIN THE OTHER METHODS, NOT VERY USEFUL BY ITSELF ****
'

library(dplyr)
library(tidyr)
library(rstatix)
library(multcompView)
library(ggplot2)


# 1. Pull unique species present in each survey unit from raw data
full_trait_survey_communities <- Combined_FISH_optimized_spatially_valid %>%
  filter(!is.na(estim_thermal_mean) & wgt_combined > 0) %>%
  distinct(survey_unit, aphia_id, accepted_name, family, Troph, tl, habitat, feeding.mode, body.shape, fin.shape, AR,
           offspring.size, spawning.type, age.maturity, length.max, growth.coefficient, fecundity, trait_resolution)

# 2. Join the classification labels. Species not classified become the "Remaining Community"
full_trait_survey_communities_classified <- full_trait_survey_communities %>%
  left_join(survey_unit_species_classes, by = c("survey_unit", "aphia_id", "accepted_name")) %>%
  mutate(archetype_class = replace_na(archetype_class, "Remaining Community"))

# 3. Standardize continuous traits *WITHIN* each individual survey unit (Z-scoring)
# This removes biogeographic bias and makes across-survey comparisons valid.
full_trait_survey_communities_classified_standardized <- full_trait_survey_communities_classified %>%
  group_by(survey_unit) %>%
  mutate(across(
    c(Troph, tl, AR, offspring.size, age.maturity, fecundity, length.max, growth.coefficient),
    ~ as.vector(scale(.)),
    .names = "z_{.col}"
  )) %>%
  ungroup()

#_____________#
# Method 1: Multinomial Mixed Model (Inference)
# Answers: Which traits significantly increase/decrease the probability of belonging to a thermal class?
#_____________#
library(brms)

trait_membership_data <- full_trait_survey_communities_classified_standardized %>%
  filter(
    archetype_class %in%
      c(
        "Warm Expansion",
        "Warm Contraction",
        "Cold Expansion",
        "Cold Contraction"
      )
  )

multinom_model <- brm(
  
  archetype_class ~
    
    z_age.maturity +
    z_length.max +
    z_growth.coefficient +
    z_Troph +
    z_tl +
    z_AR +
    z_fecundity +
    z_offspring.size +
    
    habitat +
    feeding.mode +
    spawning.type +
    body.shape +
    fin.shape +
    
    (1|survey_unit) +
    (1|family),
  
  family = categorical(),
  
  data = trait_membership_data,
  
  chains = 4,
  iter = 4000
)
posterior_summary(multinom_model)

fx <- fixef(multinom_model)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)

fx_tbl <- as.data.frame(fx) |>
  rownames_to_column("parameter") |>
  mutate(
    thermal_class = str_extract(parameter,
                                "muColdExpansion|muWarmContraction|muWarmExpansion"),
    trait = str_remove(parameter,
                       "muColdExpansion_|muWarmContraction_|muWarmExpansion_"),
    effect = case_when(
      Q2.5 > 0  ~ "Positive",
      Q97.5 < 0 ~ "Negative",
      TRUE      ~ "NS"
    )
  ) |>
  select(thermal_class, trait, effect)

summary_table <- fx_tbl |>
  pivot_wider(
    names_from = thermal_class,
    values_from = effect
  ) |>
  rename(
    `Cold Expansion`   = muColdExpansion,
    `Warm Contraction` = muWarmContraction,
    `Warm Expansion`   = muWarmExpansion
  )

summary_table

supported_effects <- fx_tbl |>
  mutate(
    effect = case_when(
      effect == "Positive" ~ "+",
      effect == "Negative" ~ "-",
      TRUE ~ ""
    )
  ) |>
  pivot_wider(
    names_from = thermal_class,
    values_from = effect
  ) |>
  rename(
    `Cold Expansion` = muColdExpansion,
    `Warm Contraction` = muWarmContraction,
    `Warm Expansion` = muWarmExpansion
  ) |>
  filter(
    `Cold Expansion` != "" |
      `Warm Contraction` != "" |
      `Warm Expansion` != ""
  )

print(supported_effects, n = Inf)
library(dplyr)
library(gt)

supported_effects %>%
  filter(trait != "Intercept") %>%
  mutate(
    trait = recode(
      trait,
      "z_age.maturity" = "Age at maturity",
      "z_length.max" = "Maximum length",
      "z_growth.coefficient" = "Growth coefficient",
      "z_AR" = "Aspect ratio",
      "z_fecundity" = "Fecundity",
      "habitatbathypelagic" = "Bathypelagic habitat",
      "habitatpelagic" = "Pelagic habitat",
      "habitatreefMassociated" = "Reef-associated habitat",
      "feeding.modeherbivorous" = "Herbivorous feeding",
      "feeding.modeplanktivorous" = "Planktivorous feeding",
      "spawning.typebearer" = "Bearer spawning",
      "spawning.typeguarder" = "Guarder spawning",
      "spawning.typenonMguarder" = "Non-male guarder spawning",
      "body.shapeelongated" = "Elongated body",
      "body.shapeflat" = "Flat body",
      "body.shapefusiform" = "Fusiform body",
      "body.shapeshortandDordeep" = "Short/deep body",
      "fin.shapeheterocercal" = "Heterocercal fin",
      "fin.shapelunate" = "Lunate fin",
      "fin.shapepointed" = "Pointed fin",
      "fin.shaperounded" = "Rounded fin",
      "fin.shapetruncate" = "Truncate fin"
    )
  ) %>%
  gt() %>%
  cols_label(
    trait = "Trait"
  ) %>%
  tab_header(
    title = "Significant odds of thermal class membership given marginal trait variations from cold contracting class"
  ) %>%
  gtsave(
    filename = "Thermal_Class_Trait_Effects_Multinomial_Model.html"
  )

###############################################################################
# H3 RANDOM FOREST FRAMEWORK
#
# OPTION A
# Classification Forest
# Response = archetype_class
#
# OPTION B
# Regression Forest
# Response = mean_regional_impact
#
# Uses:
#   full_trait_survey_communities_classified_standardized
#   labeled_regional_impact
#
# Produces:
#   rf_class_model
#   rf_importance_df
#   rf_class_predictions
#   class_probability_df
#   rf_importance_plot
#   class_probability_plot
#
#   rf_impact_model
#   rf_impact_importance_df
#   rf_impact_predictions
#   rf_observed_vs_predicted
#   impact_importance_plot
#   observed_vs_predicted_plot
###############################################################################
'
set.seed(123)

###############################################################################
# OPTION A
# RANDOM FOREST CLASSIFICATION
#
# Which traits predict thermal class membership?
###############################################################################

rf_class_data <- full_trait_survey_communities_classified_standardized %>%
  filter(
    archetype_class %in%
      c(
        "Warm Expansion",
        "Warm Contraction",
        "Cold Expansion",
        "Cold Contraction"
      )
  ) %>%
  mutate(
    archetype_class = factor(archetype_class)
  ) %>%
  distinct()

###############################################################################
# FIT CLASSIFICATION FOREST
###############################################################################

rf_class_model <- ranger(
  
  archetype_class ~
    
    z_Troph +
    z_tl +
    z_AR +
    z_offspring.size +
    z_age.maturity +
    z_fecundity +
    z_length.max +
    z_growth.coefficient +
    
    habitat +
    feeding.mode +
    spawning.type +
    body.shape +
    fin.shape,
  
  data = rf_class_data,
  
  classification = TRUE,
  probability = TRUE,
  
  importance = "permutation",
  
  num.trees = 1000,
  
  seed = 123
  
)

###############################################################################
# GLOBAL TRAIT IMPORTANCE
###############################################################################

rf_importance_df <- data.frame(
  
  Trait =
    names(
      rf_class_model$variable.importance
    ),
  
  Importance =
    unname(
      rf_class_model$variable.importance
    )
  
) %>%
  arrange(
    desc(Importance)
  )

###############################################################################
# CLASSIFICATION PERFORMANCE
###############################################################################

rf_class_predictions <- predict(
  rf_class_model,
  data = rf_class_data
)

class_probability_df <- bind_cols(
  
  rf_class_data %>%
    select(
      survey_unit,
      aphia_id,
      accepted_name,
      archetype_class
    ),
  
  as.data.frame(
    rf_class_predictions$predictions
  )
  
)

###############################################################################
# MOST CHARACTERISTIC SPECIES OF EACH CLASS
###############################################################################

top_class_members <- class_probability_df %>%
  pivot_longer(
    cols = c(
      "Cold Expansion",
      "Cold Contraction",
      "Warm Expansion",
      "Warm Contraction"
    ),
    names_to = "Predicted_Class",
    values_to = "Probability"
  ) %>%
  group_by(Predicted_Class) %>%
  slice_max(
    Probability,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

###############################################################################
# FIGURE 1
# GLOBAL IMPORTANCE
###############################################################################

rf_importance_plot <- ggplot(
  
  rf_importance_df %>%
    slice_head(n = 15),
  
  aes(
    x = reorder(
      Trait,
      Importance
    ),
    y = Importance
  )
  
) +
  
  geom_col(
    fill = "steelblue"
  ) +
  
  coord_flip() +
  
  theme_minimal(base_size = 14) +
  
  labs(
    title =
      "Traits Predicting Thermophilization Class Membership",
    
    x = NULL,
    y = "Permutation Importance"
  )

print(rf_importance_plot)

ggsave(
  filename = "Thermophilization_Class_Importance_RF_model.png",
  plot = rf_importance_plot,
  width = 12,
  height = 6,
  dpi = 300
)
###############################################################################
# FIGURE 2
# CLASS PROBABILITY DISTRIBUTIONS
###############################################################################

class_probability_plot <- class_probability_df %>%
  
  pivot_longer(
    cols = c(
      "Cold Expansion",
      "Cold Contraction",
      "Warm Expansion",
      "Warm Contraction"
    ),
    names_to = "Class",
    values_to = "Probability"
  ) %>%
  
  ggplot(
    aes(
      Class,
      Probability,
      fill = Class
    )
  ) +
  
  geom_boxplot(
    alpha = 0.8
  ) +
  
  theme_minimal(base_size = 14) +
  
  labs(
    title =
      "Posterior Class Membership Probabilities",
    
    subtitle =
      "Random Forest Classification",
    
    x = NULL,
    y = "Predicted Probability"
  ) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )

print(class_probability_plot)

###############################################################################
# OPTION B
# RANDOM FOREST REGRESSION
#
# Which traits explain thermophilization magnitude?
###############################################################################

rf_impact_data <- labeled_regional_impact %>%
  
  ungroup() %>%
  
  left_join(
    
    full_trait_survey_communities_classified_standardized %>%
      select(
        survey_unit,
        aphia_id,
        z_Troph,
        z_tl,
        z_AR,
        z_offspring.size,
        z_age.maturity,
        z_fecundity,
        z_length.max,
        z_growth.coefficient,
        habitat,
        feeding.mode,
        spawning.type,
        body.shape,
        fin.shape
      ) %>%
      distinct(),
    
    by = c(
      "survey_unit",
      "aphia_id"
    )
    
  ) %>%
  
  filter(!is.na(mean_regional_impact)) %>%
  
  ungroup()

###############################################################################
# FIT REGRESSION FOREST
###############################################################################

rf_impact_model <- ranger(
  
  mean_regional_impact ~
    
    z_Troph +
    z_tl +
    z_AR +
    z_offspring.size +
    z_age.maturity +
    z_fecundity +
    z_length.max +
    z_growth.coefficient +
    
    habitat +
    feeding.mode +
    spawning.type +
    body.shape +
    fin.shape,
  
  data = rf_impact_data,
  
  importance = "permutation",
  
  num.trees = 1000,
  
  seed = 123
  
)

###############################################################################
# PREDICTIONS
###############################################################################

rf_impact_predictions <- predict(
  rf_impact_model,
  data = rf_impact_data
)

rf_observed_vs_predicted <- rf_impact_data %>%
  mutate(
    Predicted_Impact =
      rf_impact_predictions$predictions
  )

###############################################################################
# R²
###############################################################################

impact_r2 <- cor(
  
  rf_observed_vs_predicted$mean_regional_impact,
  rf_observed_vs_predicted$Predicted_Impact,
  
  use = "complete.obs"
  
)^2

print(impact_r2)

###############################################################################
# REGRESSION IMPORTANCE
###############################################################################

rf_impact_importance_df <- data.frame(
  
  Trait =
    names(
      rf_impact_model$variable.importance
    ),
  
  Importance =
    rf_impact_model$variable.importance
  
) %>%
  arrange(
    desc(Importance)
  )

###############################################################################
# FIGURE 3
# THERMOPHILIZATION IMPACT IMPORTANCE
###############################################################################

impact_importance_plot <- ggplot(
  
  rf_impact_importance_df %>%
    slice_head(n = 15),
  
  aes(
    x = reorder(
      Trait,
      Importance
    ),
    y = Importance
  )
  
) +
  
  geom_col(
    fill = "firebrick"
  ) +
  
  coord_flip() +
  
  theme_minimal(base_size = 14) +
  
  labs(
    title =
      "Traits Predicting Thermophilization Impact",
    
    subtitle =
      paste0(
        "Random Forest Regression   (R² = ",
        round(impact_r2, 3),
        ")"
      ),
    
    x = NULL,
    y = "Permutation Importance"
  )

print(impact_importance_plot)
ggsave(
  filename = "Thermophilization_Impact_Importance_RF_model.png",
  plot = impact_importance_plot,
  width = 12,
  height = 6,
  dpi = 300
)

###############################################################################
# FIGURE 4
# OBSERVED VS PREDICTED
###############################################################################

observed_vs_predicted_plot <- ggplot(
  
  rf_observed_vs_predicted,
  
  aes(
    x = mean_regional_impact,
    y = Predicted_Impact
  )
  
) +
  
  geom_point(
    alpha = 0.4,
    color = "steelblue"
  ) +
  
  geom_smooth(
    method = "lm",
    color = "red",
    se = FALSE
  ) +
  
  theme_minimal(base_size = 14) +
  
  labs(
    title =
      "Predicting Species Thermophilization Impact",
    
    subtitle =
      paste0(
        "Random Forest Regression (R² = ",
        round(impact_r2, 3),
        ")"
      ),
    
    x = "Observed Impact",
    
    y = "Predicted Impact"
  )

print(observed_vs_predicted_plot)
ggsave(
  filename = "Observed_vs_Predicted_Thermophilization_Impact_RF_model.png",
  plot = observed_vs_predicted_plot,
  width = 6,
  height = 6,
  dpi = 300
)

###############################################################################
# OPTIONAL:
# TOP TRAITS DRIVING THERMOPHILIZATION
###############################################################################

top_classification_traits <- rf_importance_df %>%
  slice_head(n = 10)

top_impact_traits <- rf_impact_importance_df %>%
  slice_head(n = 10)

###############################################################################
# KEY OUTPUTS
###############################################################################

# rf_class_model
# rf_importance_df
# class_probability_df
# top_class_members
# rf_importance_plot
# class_probability_plot
#
# rf_impact_model
# rf_impact_importance_df
# rf_observed_vs_predicted
# impact_r2
# impact_importance_plot
# observed_vs_predicted_plot
#
# top_classification_traits
# top_impact_traits
###############################################################################
'
'

When class membership RF and impact magnitude RF identify nearly identical predictors, it implies:

The traits that define thermal classes are also the traits that determine how strongly species influence thermophilization.

That is exactly the type of convergence reviewers tend to like.
'
###############################################################################
# RANDOM FOREST THERMOPHILIZATION ANALYSIS
# REPEATED TRAIN-TEST VALIDATION (100 REPEATS)
#
# Produces:
#   Figure 1:
#     Panel A = Relative predictor importance (95% CI)
#     Panel B = Class-specific prediction accuracy (95% CI)
#
#   Figure 2:
#     Panel A = Relative predictor importance (95% CI)
#     Panel B = Cross-validated predictive skill (R² and RMSE ±95% CI)
#
# Validation:
#   100 repeated 80/20 train-test splits
#
# Interpretation:
#   Classification model:
#     Which traits predict thermophilization archetype?
#
#   Regression model:
#     Which traits predict thermophilization impact magnitude?
###############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ranger)
library(purrr)
library(patchwork)
library(forcats)

set.seed(123)

n_repeats <- 100

###############################################################################
# PUBLICATION LABELS
###############################################################################

trait_labels <- c(
  z_Troph               = "Trophic Level",
  z_tl                  = "Total Length",
  z_AR                  = "Aspect Ratio",
  z_offspring.size      = "Offspring Size",
  z_age.maturity        = "Age at Maturity",
  z_fecundity           = "Fecundity",
  z_length.max          = "Maximum Length",
  z_growth.coefficient  = "Growth Coefficient",
  habitat               = "Habitat",
  feeding.mode          = "Feeding Mode",
  spawning.type         = "Spawning Strategy",
  body.shape            = "Body Shape",
  fin.shape             = "Fin Shape"
)

###############################################################################
# CLASSIFICATION DATA
###############################################################################

rf_class_data <- full_trait_survey_communities_classified_standardized %>%
  filter(
    archetype_class %in%
      c(
        "Warm Expansion",
        "Warm Contraction",
        "Cold Expansion",
        "Cold Contraction"
      )
  ) %>%
  mutate(
    archetype_class = factor(archetype_class)
  ) %>%
  distinct()


###############################################################################
# REPEATED CLASSIFICATION RF
###############################################################################

class_results <- map(
  seq_len(n_repeats),
  function(i){
    
    set.seed(1000 + i)
    
    train_index <- sample(
      seq_len(nrow(rf_class_data)),
      size = floor(0.80 * nrow(rf_class_data))
    )
    
    train_data <- rf_class_data[train_index, ]
    test_data  <- rf_class_data[-train_index, ]
    
    rf_model <- ranger(
      
      archetype_class ~
        
        z_Troph +
        z_tl +
        z_AR +
        z_offspring.size +
        z_age.maturity +
        z_fecundity +
        z_length.max +
        z_growth.coefficient +
        
        habitat +
        feeding.mode +
        spawning.type +
        body.shape +
        fin.shape,
      
      data = train_data,
      
      classification = TRUE,
      probability = TRUE,
      
      importance = "permutation",
      
      num.trees = 1000,
      
      seed = 123
      
    )
    
    pred <- predict(
      rf_model,
      data = test_data
    )
    
    pred_class <- colnames(
      pred$predictions
    )[max.col(pred$predictions)]
    
    overall_accuracy <- mean(
      pred_class ==
        test_data$archetype_class
    )
    
    importance_df <- tibble(
      Trait =
        names(
          rf_model$variable.importance
        ),
      Importance =
        rf_model$variable.importance
    ) %>%
      mutate(
        Importance =
          Importance /
          max(Importance),
        Repeat = i
      )
    
    prediction_df <- tibble(
      Observed =
        test_data$archetype_class,
      Predicted =
        factor(
          pred_class,
          levels =
            levels(
              rf_class_data$archetype_class
            )
        ),
      Repeat = i
    )
    
    list(
      importance = importance_df,
      predictions = prediction_df,
      overall_accuracy = overall_accuracy
    )
    
  }
)

###############################################################################
# CLASSIFICATION IMPORTANCE SUMMARY
###############################################################################

class_importance_summary <- bind_rows(
  map(class_results, "importance")
) %>%
  group_by(Trait) %>%
  summarise(
    
    Mean =
      mean(Importance),
    
    Lower =
      quantile(
        Importance,
        0.025
      ),
    
    Upper =
      quantile(
        Importance,
        0.975
      ),
    
    .groups = "drop"
    
  ) %>%
  mutate(
    Trait =
      recode(
        Trait,
        !!!trait_labels
      )
  ) %>%
  arrange(
    desc(Mean)
  )

###############################################################################
# ALL HELD-OUT CLASSIFICATION PREDICTIONS
###############################################################################

all_class_predictions <- bind_rows(
  map(
    class_results,
    "predictions"
  )
)

###############################################################################
# POOLED CROSS-VALIDATED CONFUSION MATRIX
###############################################################################

confusion_summary <- all_class_predictions %>%
  count(
    Observed,
    Predicted
  ) %>%
  group_by(
    Observed
  ) %>%
  mutate(
    Proportion =
      n / sum(n)
  ) %>%
  ungroup()

###############################################################################
# OVERALL CLASSIFICATION ACCURACY
###############################################################################

overall_class_accuracy <- map_dbl(
  class_results,
  "overall_accuracy"
)

overall_class_mean <- mean(
  overall_class_accuracy
)

overall_class_ci <- quantile(
  overall_class_accuracy,
  c(
    0.025,
    0.975
  )
)

###############################################################################
# PANEL A
# RF CLASSIFICATION IMPORTANCE
###############################################################################

classification_importance_panel <- ggplot(
  
  class_importance_summary %>%
    slice_head(n = 15),
  
  aes(
    x = reorder(
      Trait,
      Mean
    ),
    y = Mean
  )
  
) +
  
  geom_point(
    size = 3,
    color = "steelblue"
  ) +
  
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0
  ) +
  
  coord_flip() +
  
  theme_minimal(
    base_size = 13
  ) +
  
  labs(
    title =
      "A. RF Class Model Predictor Importance",
    x = NULL,
    y =
      "Relative Importance (0–1)"
  )

###############################################################################
# PANEL B
# RF CLASSIFICATION CONFUSION MATRIX
###############################################################################

classification_confusion_panel <- ggplot(
  
  confusion_summary,
  
  aes(
    x = Predicted,
    y = Observed,
    fill = Proportion
  )
  
) +
  
  geom_tile(
    colour = "white"
  ) +
  
  geom_text(
    aes(
      label =
        paste0(
          round(
            Proportion * 100,
            0
          ),
          "%"
        )
    ),
    size = 4
  ) +
  
  scale_fill_viridis_c(
    option = "C",
    limits = c(0,1),
    name = "Proportion"
  ) +
  
  coord_equal() +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
  ) +
  
  labs(
    title =
      paste0(
        "B. RF Class Model Confusion Matrix\n",
        "Mean Accuracy = ",
        round(
          overall_class_mean * 100,
          1
        ),
        "% (95% CI ",
        round(
          overall_class_ci[1] * 100,
          1
        ),
        "–",
        round(
          overall_class_ci[2] * 100,
          1
        ),
        "%)"
      ),
    x = "Predicted Class",
    y = "Observed Class"
  )

###############################################################################
# REGRESSION DATA
###############################################################################

rf_impact_data <- labeled_regional_impact %>%
  
  ungroup() %>%
  
  left_join(
    
    full_trait_survey_communities_classified_standardized %>%
      
      select(
        survey_unit,
        aphia_id,
        
        z_Troph,
        z_tl,
        z_AR,
        
        z_offspring.size,
        z_age.maturity,
        z_fecundity,
        z_length.max,
        z_growth.coefficient,
        
        habitat,
        feeding.mode,
        spawning.type,
        body.shape,
        fin.shape
      ) %>%
      
      distinct(),
    
    by = c(
      "survey_unit",
      "aphia_id"
    )
    
  ) %>%
  
  filter(
    !is.na(mean_regional_impact)
  )

###############################################################################
# REPEATED REGRESSION RF
###############################################################################

impact_results <- map(
  seq_len(n_repeats),
  function(i){
    
    set.seed(5000 + i)
    
    train_index <- sample(
      seq_len(nrow(rf_impact_data)),
      size = floor(0.80 * nrow(rf_impact_data))
    )
    
    train_data <- rf_impact_data[train_index, ]
    test_data  <- rf_impact_data[-train_index, ]
    
    rf_model <- ranger(
      
      mean_regional_impact ~
        
        z_Troph +
        z_tl +
        z_AR +
        
        z_offspring.size +
        z_age.maturity +
        z_fecundity +
        z_length.max +
        z_growth.coefficient +
        
        habitat +
        feeding.mode +
        spawning.type +
        body.shape +
        fin.shape,
      
      data = train_data,
      
      importance = "permutation",
      
      num.trees = 1000,
      
      seed = 123
      
    )
    
    pred <- predict(
      rf_model,
      data = test_data
    )$predictions
    
    r2 <- cor(
      pred,
      test_data$mean_regional_impact,
      use = "complete.obs"
    )^2
    
    rmse <- sqrt(
      mean(
        (
          pred -
            test_data$mean_regional_impact
        )^2
      )
    )
    
    spearman <- cor(
      pred,
      test_data$mean_regional_impact,
      method = "spearman",
      use = "complete.obs"
    )
    
    importance_df <- tibble(
      Trait =
        names(
          rf_model$variable.importance
        ),
      Importance =
        rf_model$variable.importance
    ) %>%
      mutate(
        Importance =
          Importance /
          max(Importance)
      ) %>%
      mutate(
        Repeat = i
      )
    
    metrics <- tibble(
      Metric = c(
        "R²",
        "RMSE",
        "Spearman"
      ),
      Value = c(
        r2,
        rmse,
        spearman
      ),
      Repeat = i
    )
    
    prediction_df <- tibble(
      Observed =
        test_data$mean_regional_impact,
      Predicted =
        pred,
      Repeat = i
    )
    
    list(
      importance  = importance_df,
      metrics     = metrics,
      predictions = prediction_df
    )
    
  }
)

###############################################################################
# REGRESSION IMPORTANCE SUMMARY
###############################################################################

impact_importance_summary <- bind_rows(
  map(
    impact_results,
    "importance"
  )
) %>%
  group_by(Trait) %>%
  summarise(
    
    Mean =
      mean(Importance),
    
    Lower =
      quantile(
        Importance,
        0.025
      ),
    
    Upper =
      quantile(
        Importance,
        0.975
      ),
    
    .groups = "drop"
    
  ) %>%
  mutate(
    Trait =
      recode(
        Trait,
        !!!trait_labels
      )
  ) %>%
  arrange(
    desc(Mean)
  )

###############################################################################
# ALL HELD-OUT PREDICTIONS
###############################################################################

all_rf_predictions <- bind_rows(
  map(
    impact_results,
    "predictions"
  )
)

###############################################################################
# SPEARMAN SUMMARY
###############################################################################

bootstrap_spearman <- all_rf_predictions %>%
  group_by(Repeat) %>%
  summarise(
    Spearman = cor(
      Observed,
      Predicted,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  )

global_spearman <- cor(
  all_rf_predictions$Observed,
  all_rf_predictions$Predicted,
  method = "spearman",
  use = "complete.obs"
)

mean_spearman <- mean(
  bootstrap_spearman$Spearman
)

spearman_ci <- quantile(
  bootstrap_spearman$Spearman,
  c(
    0.025,
    0.975
  )
)

###############################################################################
# BOOTSTRAP PREDICTION ENVELOPE
###############################################################################

n_bins <- 50

prediction_envelope <- all_rf_predictions %>%
  
  mutate(
    
    Observed_bin = ntile(
      Observed,
      n_bins
    )
    
  ) %>%
  
  group_by(
    Observed_bin
  ) %>%
  
  summarise(
    
    Observed =
      median(
        Observed
      ),
    
    Mean =
      mean(
        Predicted
      ),
    
    Lower =
      quantile(
        Predicted,
        0.025
      ),
    
    Upper =
      quantile(
        Predicted,
        0.975
      ),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(
    Observed
  )

###############################################################################
# PANEL A
# IMPORTANCE
###############################################################################

impact_importance_panel <- ggplot(
  
  impact_importance_summary %>%
    slice_head(n = 15),
  
  aes(
    x = reorder(
      Trait,
      Mean
    ),
    y = Mean
  )
  
) +
  
  geom_point(
    size = 3,
    color = "firebrick"
  ) +
  
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0
  ) +
  
  coord_flip() +
  
  theme_minimal(
    base_size = 13
  ) +
  
  labs(
    title = "A. Relative Predictor Importance",
    x = NULL,
    y = "Relative Importance (0–1)"
  )

###############################################################################
# PANEL B
# OBSERVED VS PREDICTED
###############################################################################

impact_skill_panel <- ggplot() +
  
  geom_hex(
    data = all_rf_predictions,
    aes(
      x = Observed,
      y = Predicted
    ),
    bins = 40
  ) +
  
  scale_fill_viridis_c(
    option = "C",
    trans = "sqrt",
    name = "Density"
  ) +
  
  geom_ribbon(
    data = prediction_envelope,
    aes(
      x = Observed,
      ymin = Lower,
      ymax = Upper
    ),
    fill = "dodgerblue3",
    alpha = 0.25
  ) +
  
  geom_line(
    data = prediction_envelope,
    aes(
      x = Observed,
      y = Mean
    ),
    colour = "dodgerblue4",
    linewidth = 1.2
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    colour = "black"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  labs(
    title =
      "B. Cross-Validated Prediction",
    
    subtitle =
      paste0(
        "Pooled Spearman \u03c1 = ",
        round(
          global_spearman,
          2
        ),
        "; Mean bootstrap \u03c1 = ",
        round(
          mean_spearman,
          2
        ),
        " (95% CI ",
        round(
          spearman_ci[1],
          2
        ),
        "–",
        round(
          spearman_ci[2],
          2
        ),
        ")"
      ),
    
    x =
      "Observed Thermophilization Impact",
    
    y =
      "Random Forest Predicted Thermophilization Impact"
  )

###############################################################################
# COMBINED FIGURE
###############################################################################

impact_rf_figure <-
  
  impact_importance_panel +
  impact_skill_panel +
  
  plot_layout(
    ncol = 2,
    widths = c(
      1,
      1.25
    )
  )

print(
  impact_rf_figure
)

###############################################################################
# SAVE
###############################################################################

ggsave(
  "RF_Impact_RepeatedValidation.png",
  impact_rf_figure,
  width = 16,
  height = 7,
  dpi = 300
)

###############################################################################
# PANEL C
# RF IMPACT IMPORTANCE
###############################################################################

impact_importance_panel <- impact_importance_panel +
  
  labs(
    title =
      "C. RF Impact Model Predictor Importance"
  )

###############################################################################
# PANEL D
# RF IMPACT PREDICTION
###############################################################################

impact_skill_panel <- impact_skill_panel +
  
  labs(
    title =
      "D. RF Impact Model Prediction"
  )

###############################################################################
# FINAL 2x2 MULTI-PANEL FIGURE
###############################################################################

rf_combined_figure <-
  
  (
    classification_importance_panel +
      classification_confusion_panel
  ) /
  
  (
    impact_importance_panel +
      impact_skill_panel
  ) +
  
  plot_layout(
    widths = c(1, 1),
    heights = c(1, 1)
  )

print(
  rf_combined_figure
)

###############################################################################
# SAVE
###############################################################################

ggsave(
  "RF_Thermophilization_Combined.png",
  rf_combined_figure,
  width = 18,
  height = 18,
  dpi = 300
)

###############################################################################
# OUTPUT OBJECTS
###############################################################################

rf_combined_figure

confusion_summary

all_class_predictions

overall_class_mean
overall_class_ci



'

###############################################################################
# REPEATED CLASSIFICATION RF
###############################################################################

class_results <- map(
  seq_len(n_repeats),
  function(i){
    
    set.seed(1000 + i)
    
    train_index <- sample(
      seq_len(nrow(rf_class_data)),
      size = floor(0.80 * nrow(rf_class_data))
    )
    
    train_data <- rf_class_data[train_index, ]
    test_data  <- rf_class_data[-train_index, ]
    
    rf_model <- ranger(
      
      archetype_class ~
        
        z_Troph +
        z_tl +
        z_AR +
        z_offspring.size +
        z_age.maturity +
        z_fecundity +
        z_length.max +
        z_growth.coefficient +
        
        habitat +
        feeding.mode +
        spawning.type +
        body.shape +
        fin.shape,
      
      data = train_data,
      
      classification = TRUE,
      probability = TRUE,
      
      importance = "permutation",
      
      num.trees = 1000,
      
      seed = 123
      
    )
    
    pred <- predict(
      rf_model,
      data = test_data
    )
    
    pred_class <- colnames(
      pred$predictions
    )[max.col(pred$predictions)]
    
    class_accuracy <- tibble(
      Observed = test_data$archetype_class,
      Predicted = pred_class
    ) %>%
      group_by(Observed) %>%
      summarise(
        Accuracy =
          mean(
            Observed == Predicted
          ),
        .groups = "drop"
      ) %>%
      mutate(
        Repeat = i
      )
    
    overall_accuracy <- mean(
      pred_class ==
        test_data$archetype_class
    )
    
    importance_df <- tibble(
      Trait =
        names(
          rf_model$variable.importance
        ),
      Importance =
        rf_model$variable.importance
    ) %>%
      mutate(
        Importance =
          Importance /
          max(Importance)
      ) %>%
      mutate(
        Repeat = i
      )
    
    list(
      importance = importance_df,
      class_accuracy = class_accuracy,
      overall_accuracy = overall_accuracy
    )
    
  }
)

###############################################################################
# CLASSIFICATION IMPORTANCE SUMMARY
###############################################################################

class_importance_summary <- bind_rows(
  map(class_results, "importance")
) %>%
  group_by(Trait) %>%
  summarise(
    
    Mean =
      mean(Importance),
    
    Lower =
      quantile(
        Importance,
        0.025
      ),
    
    Upper =
      quantile(
        Importance,
        0.975
      ),
    
    .groups = "drop"
    
  ) %>%
  mutate(
    Trait =
      recode(
        Trait,
        !!!trait_labels
      )
  ) %>%
  arrange(desc(Mean))

###############################################################################
# CLASSIFICATION ACCURACY SUMMARY
###############################################################################

class_accuracy_summary <- bind_rows(
  map(class_results, "class_accuracy")
) %>%
  group_by(Observed) %>%
  summarise(
    
    Mean =
      mean(Accuracy),
    
    Lower =
      quantile(
        Accuracy,
        0.025
      ),
    
    Upper =
      quantile(
        Accuracy,
        0.975
      ),
    
    .groups = "drop"
    
  )

overall_class_accuracy <- map_dbl(
  class_results,
  "overall_accuracy"
)

overall_class_mean <- mean(
  overall_class_accuracy
)

overall_class_ci <- quantile(
  overall_class_accuracy,
  c(0.025, 0.975)
)

###############################################################################
# CLASSIFICATION FIGURE
###############################################################################

classification_importance_panel <- ggplot(
  
  class_importance_summary %>%
    slice_head(n = 15),
  
  aes(
    x = reorder(Trait, Mean),
    y = Mean
  )
  
) +
  
  geom_point(
    size = 3,
    color = "steelblue"
  ) +
  
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0
  ) +
  
  coord_flip() +
  
  theme_minimal(base_size = 13) +
  
  labs(
    title = "A. Relative Predictor Importance",
    x = NULL,
    y = "Relative Importance (0–1)"
  )

classification_accuracy_panel <- ggplot(
  
  class_accuracy_summary,
  
  aes(
    x = fct_reorder(
      Observed,
      Mean
    ),
    y = Mean
  )
  
) +
  
  geom_point(
    size = 4,
    color = "darkorange"
  ) +
  
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0.15
  ) +
  
  coord_flip() +
  
  ylim(0,1) +
  
  theme_minimal(base_size = 13) +
  
  labs(
    title =
      paste0(
        "B. Prediction Accuracy\nMean Overall Accuracy = ",
        round(
          overall_class_mean * 100,
          1
        ),
        "% (95% CI ",
        round(
          overall_class_ci[1] * 100,
          1
        ),
        "-",
        round(
          overall_class_ci[2] * 100,
          1
        ),
        "%)"
      ),
    x = NULL,
    y = "Prediction Accuracy"
  )

classification_rf_figure <-
  
  classification_importance_panel +
  classification_accuracy_panel +
  plot_layout(ncol = 2)

print(classification_rf_figure)

ggsave(
  "RF_Classification_RepeatedValidation.png",
  classification_rf_figure,
  width = 15,
  height = 7,
  dpi = 300
)

'

'


###############################################################################
# OUTPUT OBJECTS
###############################################################################

classification_rf_figure
impact_rf_figure

class_importance_summary
impact_importance_summary

class_accuracy_summary
impact_metric_summary

overall_class_mean
overall_class_ci

###############################################################################
'



###############################################################################
# H3: INDICATOR TRAIT ANALYSIS FOR THERMAL CLASS INCLUSION
#
# Goal:
# Identify trait states that are disproportionately associated with:
#
#   • Cold Contraction
#   • Cold Expansion
#   • Warm Contraction
#   • Warm Expansion
#
# Uses:
#   full_trait_survey_communities_classified_standardized
#
# Method:
#   1. Convert continuous traits into ordered trait states (tertiles)
#   2. Combine with categorical traits
#   3. Create trait-state presence/absence matrix
#   4. Run IndVal indicator analysis
#   5. Extract significant indicator traits
#   6. Plot strongest thermal-class signatures
###############################################################################

library(dplyr)
library(tidyr)
library(indicspecies)
library(ggplot2)
library(tibble)
library(permute)

###############################################################################
# STEP 1
# RETAIN THERMAL CLASSES ONLY
###############################################################################

indicator_data <- full_trait_survey_communities_classified_standardized %>%
  filter(
    archetype_class %in%
      c(
        "Cold Contraction",
        "Cold Expansion",
        "Warm Contraction",
        "Warm Expansion"
      )
  )

###############################################################################
# STEP 2
# CONVERT CONTINUOUS TRAITS TO ECOLOGICALLY INTERPRETABLE STATES
###############################################################################

indicator_data <- indicator_data %>%
  mutate(
    
    Age_Maturity_State =
      paste0(
        "AgeMaturity_",
        ntile(z_age.maturity, 3)
      ),
    
    Length_State =
      paste0(
        "LengthMax_",
        ntile(z_length.max, 3)
      ),
    
    Growth_State =
      paste0(
        "Growth_",
        ntile(z_growth.coefficient, 3)
      ),
    
    Trophic_State =
      paste0(
        "Trophic_",
        ntile(z_Troph, 3)
      ),
    
    Thermal_State =
      paste0(
        "Thermal_",
        ntile(z_tl, 3)
      ),
    
    AR_State =
      paste0(
        "AspectRatio_",
        ntile(z_AR, 3)
      ),
    
    Fecundity_State =
      paste0(
        "Fecundity_",
        ntile(z_fecundity, 3)
      ),
    
    Offspring_State =
      paste0(
        "OffspringSize_",
        ntile(z_offspring.size, 3)
      )
    
  )
indicator_data_complete <- indicator_data %>%
  select(
    archetype_class,
    
    Age_Maturity_State,
    Length_State,
    Growth_State,
    Trophic_State,
    Thermal_State,
    AR_State,
    Fecundity_State,
    Offspring_State,
    
    habitat,
    feeding.mode,
    spawning.type,
    body.shape,
    fin.shape
  ) %>%
  drop_na()

###############################################################################
# STEP 3
# BUILD TRAIT STATE MATRIX
###############################################################################

trait_state_matrix <- indicator_data_complete %>%
  select(
    
    Age_Maturity_State,
    Length_State,
    Growth_State,
    Trophic_State,
    Thermal_State,
    AR_State,
    Fecundity_State,
    Offspring_State,
    
    habitat,
    feeding.mode,
    spawning.type,
    body.shape,
    fin.shape
    
  )

trait_state_matrix <- model.matrix(
  ~ . - 1,
  data = trait_state_matrix
)

###############################################################################
# STEP 4
# RUN INDICATOR TRAIT ANALYSIS
###############################################################################

indicator_model <- multipatt(
  
  x = trait_state_matrix,
  
  cluster =
    indicator_data_complete$archetype_class,
  
  func = "IndVal.g",
  
  duleg = FALSE,
  
  control =
    how(
      nperm = 999
    )
  
)

###############################################################################
# STEP 5
# EXTRACT RESULTS
###############################################################################

indicator_traits_df <- indicator_model$sign %>%
  as.data.frame() %>%
  rownames_to_column("Trait") %>%
  arrange(p.value)

###############################################################################
# STEP 6
# DETERMINE ASSOCIATED THERMAL CLASS
###############################################################################
group_cols <- setdiff(
  colnames(indicator_traits_df),
  c("index", "stat", "p.value")
)

indicator_traits_df$thermal_class <-
  apply(
    indicator_traits_df[, group_cols, drop = FALSE],
    1,
    function(x){
      
      classes <- group_cols[x == 1]
      
      if(length(classes) == 0){
        return(NA_character_)
      }
      
      paste(classes, collapse = " + ")
      
    }
  )


###############################################################################
# STEP 7
# RETAIN SIGNIFICANT INDICATORS
###############################################################################

significant_indicator_traits_df <- indicator_traits_df %>%
  filter(
    p.value <= 0.05
  ) %>%
  arrange(
    p.value
  )

###############################################################################
# STEP 8
# TOP INDICATOR TRAITS WITHIN EACH THERMAL CLASS
###############################################################################

top_indicator_traits_df <- significant_indicator_traits_df %>%
  group_by(
    thermal_class
  ) %>%
  slice_max(
    stat,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

###############################################################################
# STEP 9
# CLEAN LABELS FOR FIGURES
###############################################################################

plot_indicator_traits <- top_indicator_traits_df %>%
  mutate(
    
    Trait = gsub(
      "Age_Maturity_State",
      "",
      Trait
    ),
    
    Trait = gsub(
      "Length_State",
      "",
      Trait
    ),
    
    Trait = gsub(
      "Growth_State",
      "",
      Trait
    )
    
  )

###############################################################################
# STEP 10
# VISUALIZATION
###############################################################################

indicator_trait_plot <- ggplot(
  
  plot_indicator_traits,
  
  aes(
    
    x =
      reorder(
        Trait,
        stat
      ),
    
    y = stat,
    
    fill = thermal_class
    
  )
  
) +
  
  geom_col() +
  
  coord_flip() +
  
  facet_wrap(
    ~ thermal_class,
    scales = "free_y"
  ) +
  
  theme_minimal(
    base_size = 14
  ) +
  
  labs(
    
    title =
      "Indicative Traits (IndVal) of Thermophilization Classes",
    
    x = NULL,
    
    y =
      "Indicator Value"
    
  ) +
  
  theme(
    legend.position = "none"
  )

print(indicator_trait_plot)
ggsave("Indicator_Trait_Plot.png", indicator_trait_plot, width = 12, height = 3)

###############################################################################
# STEP 11
# HUMAN-READABLE SUMMARY TABLE
###############################################################################

thermal_class_signatures_df <- significant_indicator_traits_df %>%
  group_by(
    thermal_class
  ) %>%
  summarise(
    
    top_indicator_traits =
      paste(
        head(
          Trait[order(-stat)],
          10
        ),
        collapse = "; "
      ),
    
    .groups = "drop"
    
  )

print(thermal_class_signatures_df)
gtsave(
  gt(
    thermal_class_signatures_df
  ),
  filename = "Thermal_Class_Indicator_Traits_Summary.png"
)

###############################################################################
# OUTPUTS
###############################################################################

# indicator_model
#
# indicator_traits_df
#
# significant_indicator_traits_df
#
# top_indicator_traits_df
#
# thermal_class_signatures_df
#
# indicator_trait_plot
###############################################################################
significant_indicator_traits_df

###############################################################################
######## last minute checks for thesis writing ################################
###############################################################################
number_of_reportable_surveys <- length(unique(Combined_FISH_optimized_spatially_valid$survey_unit))
number_of_reportable_hauls <- length(unique(Combined_FISH_optimized_spatially_valid$haul_id))
number_of_unique_taxa <- length(unique(Combined_FISH_optimized_spatially_valid$aphia_id))
number_of_animals <- sum(Combined_FISH_optimized_spatially_valid$num,na.rm = TRUE)

'
total_area_trawled <- Combined_FISH_optimized_spatially_valid %>%
  distinct(haul_id, .keep_all = TRUE) %>%
  summarise(total_area_swept = sum(area_swept, na.rm = TRUE))
'

troph_coverage <- Combined_FISH_optimized_spatially_valid %>%
  distinct(aphia_id, Troph) %>%
  summarise(
    total_unique_aphia = n_distinct(aphia_id),
    aphia_with_troph = n_distinct(aphia_id[!is.na(Troph)]),
    percent_taxa_with_troph = 100 * aphia_with_troph / total_unique_aphia
  )

troph_coverage_species <- Combined_FISH_optimized_spatially_valid %>%
  filter(rank == "Species") %>%
  summarise(
    total_unique_aphia = n_distinct(aphia_id),
    aphia_with_troph = n_distinct(aphia_id[!is.na(Troph)]),
    percent_species_with_troph = 100 * aphia_with_troph / total_unique_aphia
  )

print(troph_coverage_species$percent_species_with_troph)

missing_depth_haul_percentage <- length(unique(missing_depth_hauls$haul_id)) / length(unique(Combined_FISH_optimized_spatially_valid$haul_id))

sbt_validation <- global_spatially_valid_haul_stats %>%
  filter(
    !is.na(sbt),
    !is.na(glorys_bottomT)
  ) %>%
  distinct(haul_id, sbt, glorys_bottomT)

glorys_sbt_agreement_estimate <- sbt_validation %>%
  summarise(
    n_hauls = n(),
    spearman_rho = cor(sbt, glorys_bottomT, method = "spearman"),
  )
  
glorys_reliance_set <- global_spatially_valid_haul_stats %>%
  filter(
    is.na(sbt),
    !is.na(glorys_bottomT)
  ) 

glorys_reliance <- length(unique(glorys_reliance_set$haul_id)) / number_of_reportable_hauls

thermal_coverage <- Combined_FISH_optimized_spatially_valid %>%
  filter(!is.na(aphia_id)) %>%
  summarise(
    total_unique_taxa = n_distinct(aphia_id),
    taxa_with_thermal_mean = n_distinct(aphia_id[!is.na(estim_thermal_mean)]),
    percent_with_thermal_mean = 100 * taxa_with_thermal_mean / total_unique_taxa
  )

print(thermal_coverage)

thermal_coverage_species <- Combined_FISH_optimized_spatially_valid %>%
  filter(
    rank == "Species",
    !is.na(aphia_id)
  ) %>%
  summarise(
    total_unique_species = n_distinct(aphia_id),
    species_with_thermal_mean = n_distinct(aphia_id[!is.na(estim_thermal_mean)]),
    percent_with_thermal_mean = 100 * species_with_thermal_mean / total_unique_species
  )

print(thermal_coverage_species)

hex_summary_gt # for num hexes increasing CTI trend

survey_trend_summary <-
  
  survey_results_all %>%
  
  group_by(method, trend_cat) %>%
  
  summarise(
    n_surveys = n(),
    .groups = "drop"
  ) %>%
  
  group_by(method) %>%
  
  mutate(
    total_surveys = sum(n_surveys),
    percent_surveys = 100 * n_surveys / total_surveys
  ) %>%
  
  ungroup() %>%
  
  arrange(
    method,
    desc(percent_surveys)
  )

survey_trend_summary

richness_effects <- broom.mixed::tidy(
  rs_model,
  effects = "fixed",
  conf.int = TRUE
) %>%
  filter(
    term %in% c(
      "zRichWithin",
      "zRichMean"
    )
  ) %>%
  select(
    term,
    estimate,
    std.error,
    conf.low,
    conf.high,
    statistic,
    p.value
  )

richness_effects

fixed_effect_summary <- broom.mixed::tidy(
  rs_model,
  effects = "fixed",
  conf.int = TRUE
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Predictor = recode(
      term,
      zT = "Temperature",
      zDepth = "Depth",
      zEven = "Evenness",
      zTroph = "Trophic level",
      zFish = "Fishing pressure",
      zRichWithin = "Local richness",
      zRichMean = "Regional richness"
    ),
    
    abs_effect = abs(estimate),
    
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  select(
    Predictor,
    estimate,
    std.error,
    conf.low,
    conf.high,
    p.value,
    significance,
    abs_effect
  ) %>%
  arrange(desc(abs_effect))

fixed_effect_summary

dominant_taxonomic_diversity
print(unique(dominant_families_taxonomy$class))

depth_profile <- hist(global_spatially_valid_haul_stats$depth)
