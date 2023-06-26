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
    mutate(share = round(share*100/sum(share), digits = 2)) %>%
    mutate(species = gsub("\\ ", "\\_", species))
  
  # return the dataset
  return(data.out)
  
}



#' Function to generate a list with climate with species combinations
#' @param FUNDIV_climate_species data with climate and sp presence per plot
#' @param n_per_richness number of sp combinations to select per sp richness
#' @param richness_max Maximum level of species richness to explore in each CS
#' @param disturbance.in which disturbance do we plan to apply ?
#' @param exclude.in vector of species to exclude if bad estimation or IPM fit
make_CS <- function(species_share, n_per_richness = 10, richness_max = 4, 
                    disturbance.in = "storm",
                    exclude.in = c("Carpinus betulus", "Quercus ilex")){
  
  # Initialize output list
  list.out = list()
  
  # vector containing the names of all case studies
  vec.cs = unique(species_share$CS)
  
  # Load data from matreex
  data("climate_species")
  data("disturb_coef")
  
  # Species that can in the end be included
  species.in = (disturb_coef %>%
                  # We need the parameter of the disturbance
                  filter(disturbance %in% disturbance.in) %>%
                  # Included in the IPM
                  filter(species %in% fit_species) %>%
                  # Not in the species we choose to exclude
                  filter(!(species %in% exclude.in)))$species
  
  
  # Loop on all case studies
  for(i in 1:length(vec.cs)){
    
    # Initialize the output list for case study i
    list.i = list()
    
    # Name of the case study 
    CS.i = vec.cs[i]
    
    # Identify the species of interest and their relative share
    species.i = species_share %>%
      ungroup() %>%
      filter(CS == CS.i) %>%
      filter(species %in% species.in) %>%
      dplyr::select(-CS)
    
    # Loop on all level of richness
    for(j in 1:min(dim(species.i)[1], richness_max)){
      # All possible combinations
      combi.ij = combn(species.i$species, m = j)
      # Initialize data for csi richness j
      data.ij = data.frame(CS = CS.i, Richness = j, 
                           Assemblage = array(NA_character_, dim = dim(combi.ij)[2]), 
                           share = NA_real_)
      # Loop on all combinations
      for(k in 1:dim(data.ij)[1]){
        # Sum of the share of each species in the assemblage
        data.ij$share[k] = sum(filter(species.i, species %in% combi.ij[, k])$share)
        # ID in species i of the species included in the assemblage k
        data.ij$Assemblage[k] = paste(combi.ij[, k], collapse = ".")
      }
      
      # Arrange by descending share
      data.ij = data.ij %>% arrange(desc(share))
      # Remove less likely assemblages if more than n_per_richness
      if(dim(data.ij)[1] > n_per_richness) data.ij = data.ij[c(1:n_per_richness), ]
      # Add to final CS dataset 
      if(j == 1) data.i = data.ij
      else data.i = rbind(data.i, data.ij)
    }
    
    # Add the selected assemblages to the output list
    list.i$assemblages = data.i$Assemblage
    
    # Climate for the case study: mean of sgdd and wai weighted 
    list.i$climate = setNames(
      object = as.numeric(
        species.i %>%
          left_join(climate_species %>%
                      filter(N == 2) %>%
                      dplyr::select(species = sp, sgdd, wai, PC1, PC2), 
                    by = "species") %>%
          summarize(sgdd = weighted.mean(sgdd, w = share), 
                    wai = weighted.mean(wai, w = share), 
                    PC1 = weighted.mean(PC1, w = share), 
                    PC2 = weighted.mean(PC2, w = share)) %>%
          mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
                 waib = 1/(1 + wai), N = 2, SDM = 0) %>%
          dplyr::select(sgdd, wai, sgddb, waib, wai2, sgdd2, PC1, PC2, N, SDM)
      ), c("sgdd", "wai", "sgddb", "waib", "wai2", "sgdd2", "PC1", "PC2", "N", "SDM")
    )
    
    # Add list of species included in the output list
    list.i$species = unique(unlist(strsplit(data.i$Assemblage, split = "\\."))) 
    
    # Save the information on species combinations simulated
    list.i$data = data.i
    
    # Add to the final list
    eval(parse(text = paste0("list.out$", CS.i, " = list.i")))
    
  }
  
  # Return output
  return(list.out)
}



#' Function to create a list of IPM to run
#' @param CS list generate by make_CS function
make_species_list = function(CS){
  
  # Loop on all climates
  for(i in 1:length(names(CS))){
    # Dataframe for climate i
    out.i = data.frame(
      ID.CS = i, 
      CS = names(CS)[i],
      species = CS[[i]]$species
    )
    
    # Add to final dataset
    if(i == 1) out = out.i
    else out = rbind(out, out.i)
    
  }
  
  # Final formatting
  out = out %>%
    mutate(ID.species = c(1:dim(.)[1])) %>%
    dplyr::select(ID.species, ID.CS, CS, species)
  
  # Return output
  return(out)
  
}



#' Function to make a species object, save it as rds and return filename
#' @param CS list generate by make_CS function
#' @param species_list table containing all species to create
#' @param ID.species.in ID of the species to make in species_list
make_species_rds = function(CS, species_list, ID.species.in){
  
  # Load demographic parameter of the species 
  eval(parse(text=paste0("fit.in <- fit_", species_list$species[ID.species.in])))
  
  # Make IPM
  IPM.in = make_IPM(
    species = species_list$species[ID.species.in], 
    climate = CS[[species_list$ID.CS[ID.species.in]]]$climate, 
    fit =  fit.in, 
    clim_lab = paste0("clim_", species_list$CS[ID.species.in]),
    mesh = c(m = 700, L = 100, U = as.numeric(fit.in$info[["max_dbh"]]) * 1.1),
    BA = 0:200, verbose = TRUE, correction = "none"
  )
  
  # Create species object 
  species.in = species(IPM.in, init_pop = def_initBA(20), harvest_fun = def_harv, 
                       disturb_fun = def_disturb)
  
  # Name of the file to save
  file.in = paste0("rds/", species_list$CS[ID.species.in], "/species/", 
                   species_list$species[ID.species.in], ".rds")
  
  # Save species object in a rdata
  create_dir_if_needed(file.in)
  saveRDS(species.in, file.in)
  
  # Return output list
  return(file.in)
  
}



#' Function to create a list of forest for the simulations
#' @param CS list generate by make_CS function
make_forest_list = function(CS){
  
  # Loop on all climates
  for(i in 1:length(names(CS))){
    # Dataframe for climate i
    out.i = CS[[i]]$data %>%
      mutate(ID.CS = i)
    # Add to final dataset
    if(i == 1) out = out.i
    else out = rbind(out, out.i)
    
  }
  
  # Final formatting
  out = out %>%
    mutate(ID.forest = c(1:dim(.)[1])) %>%
    dplyr::select("ID.forest", "ID.CS", "CS", "Richness", "Assemblage")
  
  # Return output
  return(out)
  
}




#' Function to make a list of simulations till equilibrium
#' @param CS list generate by make_CS function
#' @param species_list df with information on all species object
#' @param forest_list df with information on all forest to simulate 
#' @param species vector containing all species rds files created
#' @param ID.forest.in ID of the forest to simulate in forest_list
make_simulations_equilibrium = function(
  CS, species_list, forest_list, species, ID.forest.in){
  
  # Identify ID of the climate
  ID.CS.in = forest_list$ID.CS[ID.forest.in]
  
  # Identify species combination i
  assemblage.in = forest_list$Assemblage[ID.forest.in]
  
  # vector of species in forest i
  species.in = unlist(strsplit(assemblage.in, "\\."))
  
  # Initialize list of species to create
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  # Loop on all species
  for(i in 1:length(species.in)){
    
    # Identify the file in species containing species i
    species.file.i = species[(species_list %>%
                                filter(species == species.in[i]) %>%
                                filter(ID.CS == ID.CS.in))$ID.species]
    
    # Store the file in the list
    list.species[[i]] = readRDS(species.file.i)
    
  }
  
  
  # Make forest
  forest.in = new_forest(species = list.species, harv_rules = c(
    Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1))
  
  # Run simulation till equilibrium
  sim.in = sim_deter_forest(
    forest.in, tlim = 4000, equil_time = 50000, equil_dist = 2000, 
    equil_diff = 0.5, harvest = "default", SurfEch = 0.03, verbose = TRUE)
  
  # Name of the file to save
  file.in = paste0("rds/", forest_list$CS[ID.forest.in], "/sim_equilibrium/", 
                   gsub("\\.", "\\-", forest_list$Assemblage[ID.forest.in]), 
                   ".rds")
  
  # Save species object in a rdata
  create_dir_if_needed(file.in)
  saveRDS(sim.in, file.in)
  
  # Return output list
  return(file.in)
}

