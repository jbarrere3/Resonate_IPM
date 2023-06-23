#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien Barrere
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to get the path of a file, and create directories if they don't exist
#' @param file.in character: path of the file, filename included (ex: "plot/plot.png")
create_dir_if_needed <- function(file.in){
  
  path.in <- strsplit(file.in, "/")[[1]]
  if(length(path.in) > 1){
    for(i in 1:(length(path.in)-1)){
      if(i == 1) path.in_i <- path.in[i]
      else path.in_i <- paste(path.in_i, path.in[i], sep = "/")
      if(!dir.exists(path.in_i)) dir.create(path.in_i)
    }
  }
}


#' Format the raw dataset to get the proportion of each IPM species per case study
#' @param species_share_raw raw data of species share per CS
format_species_share_raw = function(species_share_raw){
  
  # Species included in the IPM
  species.IPM = gsub("\\_", "\\ ", fit_species)
  
  # Change format of the table, and get species name and numeric share
  data = species_share_raw %>%
    # countries from columns to lines
    gather(key = "CS", value = "sp_raw", colnames(.)) %>%
    # Remove nas
    drop_na() %>%
    # Extract species and relative share
    mutate(species = gsub("\\ \\(.+", "", sp_raw), 
           share = gsub("\\%.+", "", gsub(".+\\(", "", sp_raw)), 
           share = as.numeric(gsub("\\,", "\\.", share))) %>%
    # Standardize other broadleaf and conifer categories
    mutate(species = case_when(
      substr(species, 1, 8) %in% c("Other br", "Other an", "other br") ~ "Other broadleaf", 
      substr(species, 1, 8) %in% c("Other co", "Other gy") ~ "Other conifer", 
      TRUE ~ species))
  
  # Table linking family to order
  data(vascular.families)
  
  # Extract genus and family of all species present
  data.species <- cbind(data.frame(species.original = unique(data$species)), 
                        TPL(splist = unique(data$species)))
  data.species = data.species %>%
    dplyr::select(species = species.original, 
                  genus = New.Genus) %>% 
    left_join((data.species %>%
                 dplyr::select(genus = New.Genus, family = Family) %>%
                 filter(!is.na(family)) %>%
                 filter(family != "") %>%
                 distinct()), 
              by = "genus") %>%
    left_join((vascular.families %>%
                 dplyr::select(family = Family, group = Group) %>%
                 rbind(data.frame(family = "Leguminosae", group = "angiosperms"))), 
              by = "family") %>%
    filter(!is.na(genus)) %>%
    mutate(group = ifelse(is.na(group), group, paste(
      toupper(substr(group, 1, 1)), substr(group, 2, nchar(group)), sep="")))
  
  # Final dataset
  data.out = data %>%
    left_join(data.species, by = "species") %>%
    mutate(inIPM = ifelse(species %in% species.IPM, TRUE, FALSE)) %>%
    mutate(sp_final = case_when(
      (!inIPM & group == "Angiosperms") ~ "Other broadleaf",
      (!inIPM & group == "Gymnosperms") ~ "Other conifer", 
      species %in% c("Other broadleaf", "Other conifer", species.IPM) ~ species, 
      TRUE ~ "Other unknown"
    )) %>%
    group_by(CS, sp_final) %>%
    summarize(share = sum(share)) %>%
    rename(species = sp_final) %>%
    arrange(CS, desc(share)) %>%
    # Adjust so that the sum of species share equals 100 in each CS
    ungroup() %>% group_by(CS) %>%
    mutate(share = round(share*100/sum(share), digits = 2))
  
  # return the dataset
  return(data.out)
  
}