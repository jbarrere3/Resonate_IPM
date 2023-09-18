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


#' Get coordinate in first pca traits axis per species
#' @param traits dataframe containing trait value per species
get_pc1_per_species <- function(traits){
  
  # - Make PCA 
  pca <- prcomp((traits %>% dplyr::select(-species)), 
                center = T, scale = T)
  
  # - Extract the coordinates of the ndividuals on pca axis
  out <- data.frame(species = traits$species, 
                    pca1 = get_pca_ind(pca)[[1]][, 1]) 
  
  # return the output
  return(out)
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


#' Format the raw climate data to extract the initial climate per case study
#' @param climate_raw Raw climate data
#' @param ssp.in ssp scenario to choose the data
#' @param year.in year to choose for the reference climate
get_climate_ref = function(climate_raw, ssp.in = "ssp126", year.in = 2015){
  
  # Data frame to join case studies name
  df.cs = data.frame(
    cs = c("Bauges", "Catalonia",  "Galicia / Northern Portugal", 
           "Ireland", "Istria", "Kostelek", "New Forest", 
           "South Western Finland", "Upper Rhine Valley and Foothills"), 
    CS = c("FRANCE", "CATALONIA", "GALICIA", "REPUBLIC_OF_IRELAND", 
           "CROATIA", "CZECH_REPUBLIC", "UK", "FINLAND", "GERMANY")
  )
  
  # Format the output dataset
  out = climate_raw %>%
    # Only keep the right year and ssp scenario
    filter(year == year.in & ssp == ssp.in) %>%
    # Add the case study names used for the simulations
    left_join(df.cs, by = "cs") %>%
    # Remove errors in the calculation
    filter(pet > 0) %>%
    # Calculate average climate per case study 
    group_by(CS) %>%
    summarize(wai = mean(wai), sgdd = mean(sgdd))
  
  # Return formatted data set
  return(out)
}


#' Function to generate a list with climate with species combinations
#' @param species_share data with share of each species per cs
#' @param climate_ref data frame containing climate per case study
#' @param n_per_richness number of sp combinations max to select per sp richness
#' @param richness_max Maximum level of species richness to explore in each CS
#' @param share_min_coef Vector of length two giving intercept and slope of 
#'                      the relation btw richness and minimum share of all sp
#' @param disturbance.in which disturbance do we plan to apply ?
#' @param exclude.in vector of species to exclude if bad estimation or IPM fit
make_CS <- function(species_share, climate_ref, n_per_richness = 20, richness_max = 7, 
                    share_min_coef = c(-10, 15), disturbance.in = "storm",
                    exclude.in = c("Carpinus_betulus", "Quercus_ilex")){
  
  # Initialize output list
  list.out = list()
  
  # vector containing the names of all case studies
  vec.cs = unique(species_share$CS)
  
  # Load data from matreex
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
      
      # Ensure that we select combinations that represent a sufficeintly high share
      data.ij = data.ij %>% 
        filter(share > (share_min_coef[1] + share_min_coef[2]*j))
      
      # Keep the combinations selected below n_per_richness
      data.ij = data.ij %>% arrange(desc(share))
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
        climate_ref %>% 
          filter(CS == CS.i) %>%
          mutate(PC1 = 0, PC2 = 0, sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
                 waib = 1/(1 + wai), N = 2, SDM = 0) %>%
          dplyr::select(sgdd, wai, sgddb, waib, wai2, sgdd2, PC1, PC2, N, SDM)
      ), c("sgdd", "wai", "sgddb", "waib", "wai2", "sgdd2", "PC1", "PC2", "N", "SDM")
    )
    
    # Add list of species included in the output list
    list.i$species = unique(unlist(strsplit(data.i$Assemblage, split = "\\."))) 
    
    # Save the information on species combinations simulated
    list.i$data = data.i
    
    # Add to the final list if enough data
    if(dim(data.i)[1] >= 10){
      eval(parse(text = paste0("list.out$", CS.i, " = list.i")))}
    
    
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


#' Disturbance function
#'
#' @param x population state distribution at time t
#' @param species The species class object of interest to get mesh and RDIcoef
#' values from. RDIcoef is a one line dataframe with RDI coefficient for one
#' species.
#' @param disturb Disturbance parameters. Highly depend on the disturbance
#' impact parameters given to the species.
#' @param ... Not used in this case.
#' \describe{
#' \item{qmd}{Forest Quadratic Mean Diameter}
#' }
#' @author Maxime Jeaunatre
#'
disturb_fun <- function(x, species, disturb = NULL, ...){
  
  dots <- list(...)
  qmd <- dots$qmd 
  size <- species$IPM$mesh
  coef <- species$disturb_coef
  if(any(disturb$type %in% coef$disturbance)){
    coef <- subset(coef, disturbance == disturb$type)
  } else {
    stop(sprintf("The species %s miss this disturbance type (%s) parameters",
                 sp_name(species), disturb$type))
  }
  
  # edits for delay
  size[size == 0] <- min(size[size !=0])
  
  logratio <-  log(size / qmd)
  dbh.scaled = coef$dbh.intercept + size * coef$dbh.slope
  logratio.scaled = coef$logratio.intercept + logratio * coef$logratio.slope
  Pkill <- plogis(coef$a0 + coef$a1 * logratio.scaled + 
                    coef$b * disturb$intensity ^(coef$c * dbh.scaled))
  
  return(x* Pkill) # always return the mortality distribution
}





#' Function to make a list of simulations with disturbance
#' @param CS list generate by make_CS function
#' @param species_list df with information on all species object
#' @param forest_list df with information on all forest to simulate 
#' @param species vector containing all species rds files created
#' @param sim_equilibrium Vector containing file names of simulations till equil
#' @param ID.forest.in ID of the forest to simulate in forest_list
#' @param disturbance.df disturbance dataframe
make_simulations_disturbance = function(
  CS, species_list, forest_list, species, sim_equilibrium, 
  ID.forest.in, disturbance.df){
  
  # Identify ID of the climate
  ID.CS.in = forest_list$ID.CS[ID.forest.in]
  
  # Identify species combination i
  assemblage.in = forest_list$Assemblage[ID.forest.in]
  
  # vector of species in forest i
  species.in = unlist(strsplit(assemblage.in, "\\."))
  
  # Initialize list of species to create
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  # Read the simulation at equilibrium
  sim_equilibrium.in = readRDS(sim_equilibrium[ID.forest.in])
  
  # Checked that the population reached equilibrium
  reached_equil = ifelse(
    is.na(sum((sim_equilibrium.in %>%
                 filter(var == "BAsp") %>%
                 filter(time == max(.$time) - 1))$value)), 
    FALSE, TRUE
  )
  
  # Basal area at equilibrium
  BA_eq = ifelse(!reached_equil, 200, 
                 sum((sim_equilibrium.in %>% 
                        filter(var == "BAsp" & equil))$value))
  
  # Only make the simulation id ==f population reached an equilibrium
  if(reached_equil & BA_eq < 120){
    # Loop on all species
    for(i in 1:length(species.in)){
      
      # Identify the file in species containing species i
      species.file.i = species[(species_list %>%
                                  filter(species == species.in[i]) %>%
                                  filter(ID.CS == ID.CS.in))$ID.species]
      
      # Store the file in the list
      list.species[[i]] = readRDS(species.file.i)
      
      # Extract the equilibrium for species i
      equil.i = sim_equilibrium.in %>%
        filter(var == "n", equil, species == species.in[i]) %>% 
        pull(value)
      
      # Initiate the population at equilibrium
      list.species[[i]]$init_pop <- def_init_k(equil.i*0.03)
      
      # Update disturbance function
      list.species[[i]]$disturb_fun <- disturb_fun
      
      # Add disturbance coefficients
      list.species[[i]]$disturb_coef <- filter(matreex::disturb_coef, 
                                               species == species.in[i])
    }
    
    
    # Make forest
    forest.in = new_forest(species = list.species, harv_rules = c(
      Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1))
    
    # Run simulation till equilibrium
    sim.in = sim_deter_forest(
      forest.in, tlim = 4000, equil_time = 4000, disturbance = disturbance.df, 
      SurfEch = 0.03, verbose = TRUE)
    
  } else {
    sim.in = matrix()
  }
  
  
  # Name of the file to save
  file.in = paste0("rds/", forest_list$CS[ID.forest.in], "/sim_disturbance/", 
                   gsub("\\.", "\\-", forest_list$Assemblage[ID.forest.in]), 
                   ".rds")
  
  # Save species object in a rdata
  create_dir_if_needed(file.in)
  saveRDS(sim.in, file.in)
  
  # Return output list
  return(file.in)
}




#' Get resilience, resistance and recovery from simulations with disturbance
#' @param sim_disturbance vector containing the file names of sim with disturbance
#' @param disturbance.df disturbance dataset used to generate the disturbance
#' @param forest_list Table giving the information on each forest generated
get_resilience_metrics <- function(sim_disturbance, disturbance.df, forest_list){
  
  # Initialize the output
  out <- forest_list %>%
    dplyr::select(ID.forest, ID.CS, Assemblage) %>%
    mutate(resistance = NA_real_, recovery = NA_real_, resilience = NA_real_, 
           t0 = NA_real_, thalf = NA_real_, SD = NA_real_, BA_diff = NA_real_, 
           BA_eq = NA_real_, dbh_mean = NA_real_, dbh_q10 = NA_real_, 
           dbh_q90 = NA_real_, dbh_mean_postdist = NA_real_, 
           dbh_q10_postdist = NA_real_, dbh_q90_postdist = NA_real_)
  
  # Identify disturbance time
  tdist = min(disturbance.df$t)
  
  # Loop on all species combination
  for(i in 1:length(sim_disturbance)){
    
    # Printer
    print(paste0("Reading simulation ", i, "/", length(sim_disturbance)))
    
    # Read simulation i
    sim.i = readRDS(sim_disturbance[i])
    
    # First, verify that equilibrium was reached
    if(!is.na(sim.i[1, 1])){
      
      # Also verify that there was not a problem after the disturbance
      if(!any(is.na(filter(sim.i, var == "BAsp" & 
                           time == (max(disturbance.df$t)+1))$value))){
        
        # mean dbh at equilibrium and after disturbance
        dbh_i = sim.i %>%
          filter(var == "n") %>%
          filter(time %in% c(1, (max(disturbance.df$t)+1))) %>%
          group_by(size, time) %>%
          summarize(ntot = sum(value)) %>%
          ungroup() %>% group_by(time) %>%
          filter(size > 0) %>%
          mutate(ntot_size = ntot*size) %>%
          summarize(mean_dbh = weighted.mean(size, w = ntot), 
                    q10_dbh = weighted.quantile(size, w = ntot, prob = 0.1), 
                    q90_dbh = weighted.quantile(size, w = ntot, prob = 0.9))
        out$dbh_mean[i] <- subset(dbh_i, time == 1)$mean_dbh
        out$dbh_q10[i] <- subset(dbh_i, time == 1)$q10_dbh
        out$dbh_q90[i] <- subset(dbh_i, time == 1)$q90_dbh
        out$dbh_mean_postdist[i] <- subset(dbh_i, time != 1)$mean_dbh
        out$dbh_q10_postdist[i] <- subset(dbh_i, time != 1)$q10_dbh
        out$dbh_q90_postdist[i] <- subset(dbh_i, time != 1)$q90_dbh
        
        # Format the output
        data.i <- sim.i %>%
          filter(var == "BAsp") %>%
          filter(!equil) %>%
          group_by(time) %>%
          summarize(BA = sum(value))
        
        ## Calculate stability before disturbance (to check equilibrium)
        out$SD[i] = sd(subset(data.i, time %in% c(1:(tdist-1)))$BA)
        out$BA_diff[i] = diff(range(subset(data.i, time %in% c(1:(tdist-1)))$BA))
        
        ## Calculate resistance
        #  - Basal area at equilibrium
        Beq.i = mean((data.i %>% filter(time < min(disturbance.df$t)))$BA)
        out$BA_eq[i] = Beq.i
        # - Basal area after disturbance
        Bdist.i = (data.i %>% filter(time == max(disturbance.df$t)+1))$BA
        # - Resistance : logit of the percentage of basal area that survived 
        #out$resistance[i] = Beq.i/(Beq.i - Bdist.i)
        out$resistance[i] = log((Bdist.i/Beq.i)/(1 - (Bdist.i/Beq.i)))
        
        ## Calculate recovery
        #  - Time at which population recovered fully
        Rec.time.i = min((data.i %>% 
                            filter(time > max(disturbance.df$t)) %>%
                            filter(BA > Beq.i))$time)
        # - Basal area 20 years after disturbance
        Bdist20.i = (data.i %>% filter(time == max(disturbance.df$t)+21))$BA
        # - Recovery = slope of BA increase in teh 20 years after disturbance
        out$recovery[i] = abs(Bdist20.i - Bdist.i)/20
        
        ## Calculate resilience
        out$resilience[i] <- 1/sum((data.i %>%
                                      mutate(BA0 = .[which(.$time == 1), "BA"]) %>%
                                      mutate(diff = abs(BA - BA0)))$diff)
        
        ## Calculate t0
        #  - Time at which population recovered to 5% of the basal area lost
        Rec.0.time.i = min((data.i %>% 
                              filter(time > max(disturbance.df$t)) %>%
                              filter(BA > (Beq.i + 19*Bdist.i)/20))$time)
        # - Recovery = time to recover minus time of disturbance
        out$t0[i] = Rec.0.time.i - max(disturbance.df$t)
        
        ## Calculate thalf
        #  - Time at which population recovered to 50% of the basal area lost
        Rec.half.time.i = min((data.i %>% 
                                 filter(time > max(disturbance.df$t)) %>%
                                 filter(BA > (Beq.i + Bdist.i)/2))$time)
        # - Recovery = time to recover minus time of disturbance
        out$thalf[i] = Rec.half.time.i - max(disturbance.df$t)
        
      }
      
    }
    
  }
  
  # Return output
  return(out)
}



#' Get functional diversity from simulations with FD package
#' @param forest_list df containing info on each forest simulated
#' @param sim_disturbance vector containing the file names of sim with disturbance
#' @param pc1_per_species Position of each species along the growth-mortality trade-off
get_FD <- function(forest_list, sim_disturbance, pc1_per_species){
  
  # Initialize vector of all successful simulations
  vec.sim = c()
  
  # Initialize the original fd data
  data.fd.original <- forest_list %>%
    dplyr::select(ID.forest, ID.CS, Assemblage, nsp = Richness) %>%
    mutate(FD = NA_real_, H = NA_real_, D = NA_real_, Nha = NA_real_)
  
  # Loop on all species combination
  for(i in 1:length(sim_disturbance)){
    
    # Printer
    print(paste0(i, "/", length(sim_disturbance)))
    
    # Read simulation i
    sim.i = readRDS(sim_disturbance[i])
    
    # First, verify that equilibrium was reached
    if(!is.na(sim.i[1, 1])){
      
      ## - Calculate the number of trees per ha at equilibrium
      data.fd.original$Nha[i] = sum((sim.i %>%
                                       filter(var == "N") %>%
                                       filter(time == 1))$value)
      
      ## - FD with the original approach
      data.fd.original.i <- sim.i %>%
        filter(var == "BAsp") %>%
        filter(time == 1) %>%
        left_join(pc1_per_species, by = "species") %>%
        mutate(p = value/sum(.$value), 
               plnp = p*log(p), 
               p2 = p^2) %>%
        summarise(FD = weighted.var(pca1, w = value), 
                  H = -sum(plnp), 
                  D = 1/sum(p2)) 
      data.fd.original$FD[i] <- data.fd.original.i$FD
      data.fd.original$H[i] <- data.fd.original.i$H
      data.fd.original$D[i] <- data.fd.original.i$D
      
      ## - FD with the FD package (dimension approach)
      data.fd.dimension.i <- sim.i %>%
        mutate(ID.CS = forest_list$ID.CS[i], 
               ID.forest = forest_list$ID.forest[i], 
               ID.community = paste(ID.CS, ID.forest, sep = ".")) %>%
        filter(var == "BAsp") %>%
        filter(time == 1) %>%
        dplyr::select(ID.CS, ID.forest, ID.community, species, value)
      
      # Also check we have trait value for all species in the community
      if(all(data.fd.dimension.i$species %in% pc1_per_species$species)){
        # Add to the final dataframe
        if(length(vec.sim) == 0) data.fd.dimension = data.fd.dimension.i
        else data.fd.dimension = rbind(data.fd.dimension, data.fd.dimension.i)
        
        # Increment the counter
        vec.sim = c(vec.sim, i)
      }
    }
    
  }
  
  # Replace NA by 0 in original approach
  data.fd.original <- data.fd.original %>% mutate(FD = ifelse(is.na(FD), 0, FD))
  
  # Abundance dataframe
  abun.df = data.fd.dimension %>%
    spread(key = "species", value = "value") %>% 
    replace(is.na(.), 0)
  
  # Abundance matrix
  abun.matrix = as.matrix(
    abun.df %>% dplyr::select(-ID.CS, -ID.forest, -ID.community))
  rownames(abun.matrix) = abun.df$ID.community
  
  # Trait df
  trait.df = data.frame(species = colnames(abun.matrix)) %>%
    left_join(pc1_per_species, by = "species") %>%
    dplyr::select(-species)
  rownames(trait.df) = colnames(abun.matrix)
  
  # Calculate FD per community
  fd.raw <- dbFD(trait.df, abun.matrix, w.abun = TRUE)
  
  # Final dataframe
  data.fd.dimension.final <- data.frame(
    ID.forest = forest_list$ID.forest[vec.sim], 
    ID.CS = forest_list$ID.CS[vec.sim], 
    Assemblage = forest_list$Assemblage[vec.sim], 
    FRic = fd.raw$FRic, 
    FDis = fd.raw$FDis, 
    CWM = fd.raw$CWM$pca1
  )
  
  # Join the two datasets
  out = left_join(data.fd.original, data.fd.dimension.final, 
                  by = c("ID.forest", "ID.CS", "Assemblage")) %>%
    filter(!is.na(FRic))
  
  # Return output
  return(out)}


#' Format data for the models
#' @param climate list of climate objects used for the simulations
#' @param resilience df with resilience per forest ID
#' @param FD df with FD and CWM per forest ID
get_data_model = function(climate, resilience, FD){
  
  # Build a dataset to associate ID climate with sgdd, wai and pca1
  data.climate = data.frame(
    ID.climate = c(1:length(names(climate))), 
    pca1 = NA_real_, 
    sgdd = NA_real_, 
    wai = NA_real_
  )
  for(i in 1:dim(data.climate)[1]){
    data.climate$pca1[i] = climate[[i]]$climate[7]
    data.climate$sgdd[i] = climate[[i]]$climate[1]
    data.climate$wai[i] = climate[[i]]$climate[2]
  }
  
  # Format final dataset
  out = resilience %>%
    left_join(FD, by = c("ID.forest", "ID.climate", "combination")) %>%
    left_join(data.climate, by = "ID.climate") %>%
    rename(forest.composition = combination)
  
  # Return output
  return(out)
}



#' Format future climatic data to keep only information of interest for simul
#' @param future_climate_all Raw data on future climatic oconditions
format_future_climate = function(future_climate_all){
  
  # -- data to join case studies name
  df.cs = data.frame(
    cs = c("Bauges", "Catalonia",  "Galicia / Northern Portugal", 
           "Ireland", "Istria", "Kostelek", "New Forest", 
           "South Western Finland", "Upper Rhine Valley and Foothills"), 
    CS = c("FRANCE", "CATALONIA", "GALICIA", "REPUBLIC_OF_IRELAND", 
           "CROATIA", "CZECH_REPUBLIC", "UK", "FINLAND", "GERMANY")
  )
  
  # -- homogenize case study in climate file 
  out = future_climate_all %>%
    left_join(df.cs, by = "cs") %>%
    dplyr::select(-cs)
  
  # -- average climate over case studies
  out = out %>%
    group_by(CS, year, ssp) %>%
    filter(pet > 0) %>%
    summarize(sgdd = mean(sgdd, na.rm = TRUE), 
              wai = mean(wai, na.rm = TRUE))
  
  # -- Return output
  return(out)
}


#' Function to make a species object, save it as rds and return filename
#' @param future_climate future climate data, formatted
#' @param species_list table containing all species to create
#' @param ID.species.in ID of the species to make in species_list
make_species_mu_rds = function(future_climate, species_list, ID.species.in){
  
  # Identify the species and case study for iteration i
  species.in = species_list$species[ID.species.in]
  CS.in = species_list$CS[ID.species.in]
  
  # Load demographic parameter of the species 
  eval(parse(text=paste0("fit.in <- fit_", species.in)))
  
  # Restrict climate file to the case study of interest
  clim.in = subset(future_climate, CS == CS.in) 
  
  # From clim.in, build the margins file for iteration i
  margins.in = data.frame(
    sgdd = c(max(clim.in$sgdd), mean(clim.in$sgdd), min(clim.in$sgdd)), 
    wai = c(min(clim.in$wai), mean(clim.in$wai), max(clim.in$wai)), 
    N = c(1, 2, 3)) %>%
    mutate(PC1 = 0, PC2 = 0, sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
           waib = 1/(1 + wai), SDM = 0) %>%
    dplyr::select(sgdd, wai, sgddb, waib, wai2, sgdd2, PC1, PC2, N, SDM)
  
  # Make the mu matrix
  mu.in <- make_mu_gr(
    species = species.in, fit = fit.in, climate = margins.in, 
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit.in) * 1.1),
    verbose = TRUE, stepMu = 0.001)
  
  # Create species object from random distribution
  sp.in = species(IPM = mu.in, init_pop = def_initBA(20),
                  harvest_fun = def_harv)
  # Update disturbance function
  sp.in$disturb_fun = disturb_fun
  # Add disturbance coefficients
  sp.in$disturb_coef  <- filter(
    matreex::disturb_coef, species == species.in)
  
  # Name of the file to save
  file.in = paste0("rds/", CS.in, "/species_mu/", species.in, ".rds")
  
  # Save species object in a rdata
  create_dir_if_needed(file.in)
  saveRDS(sp.in, file.in)
  
  # Return output list
  return(file.in)
  
}



#' Function to make a list of simulations till equilibrium via mu integration
#' @param future_climate future climate data, formatted
#' @param species_list df with information on all species object
#' @param species_mu vector containing all species mu rds files created
#' @param ID.species.in ID of the species to simulate in forest_list
make_simulations_equilibrium_mu = function(
  future_climate, species_list, species_mu, ID.species.in){
  
  
  # Identify the species and case study for iteration i
  species.in = species_list$species[ID.species.in]
  CS.in = species_list$CS[ID.species.in]
  
  # Climate for the simulation : climate of the first year
  climate.in = setNames(
    object = as.numeric(
      future_climate %>% 
        ungroup() %>%
        filter(CS == CS.in & year == min(future_climate$year)) %>%
        summarize(sgdd = mean(sgdd), wai = mean(wai)) %>%
        mutate(PC1 = 0, PC2 = 0, sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
               waib = 1/(1 + wai), N = 2, SDM = 0) %>%
        dplyr::select(sgdd, wai, sgddb, waib, wai2, sgdd2, PC1, PC2, N, SDM)
    ), c("sgdd", "wai", "sgddb", "waib", "wai2", "sgdd2", "PC1", "PC2", "N", "SDM")
  )
  
  
  # Load species object
  sp.in = readRDS(species_mu[ID.species.in])
  
  # Make forest
  forest.in = forest(species = list(mu = sp.in), harv_rules = c(
    Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1))
  
  # Run simulation till equilibrium
  sim.in = sim_deter_forest(
    forest.in, tlim = 4000, equil_time = 50000, equil_dist = 2000, climate = climate.in,
    equil_diff = 0.5, harvest = "default", SurfEch = 0.03, verbose = TRUE)
  
  # Name of the file to save
  file.in = paste0("rds/", CS.in, "/sim_equilibrium_mu/", species.in, ".rds")
  
  # Save species object in a rdata
  create_dir_if_needed(file.in)
  saveRDS(sim.in, file.in)
  
  # Return output list
  return(file.in)
}



#' Make mu simulations with dist and changing climate
#' @param future_climate future climate data, formatted
#' @param species_list df with information on all species object
#' @param species_mu vector containing all species mu rds files created
#' @param sim_equilibrium_mu Vector containing file names of simulations till equil
#' @param ID.scenarios ID of the climatic - dist scenario
#' @param df_scenarios disturbance and climatic scenarios dataframe
make_simulations_disturbance_mu = function(
  species_list, future_climate, species_mu, sim_equilibrium_mu, ID.scenarios, 
  df_scenarios){
  
  print(ID.scenarios)
  # Identify the species, case studies, climate scenarios, etc.
  ID.species.in = df_scenarios$ID.species[ID.scenarios]
  species.in = species_list$species[ID.species.in]
  CS.in = species_list$CS[ID.species.in]
  ssp.in = df_scenarios$ssp[ID.scenarios]
  dist.n.in = df_scenarios$dist.n[ID.scenarios]
  dist.type.in = df_scenarios$dist.type[ID.scenarios]
  dist.Iref.in = df_scenarios$dist.Iref[ID.scenarios]
  dist.Islope.in = df_scenarios$dist.Islope[ID.scenarios]
  
  # Build climate file for the simulations 
  #   /!\ CHANGE TIME WHEN CLIMATE PROBLEM SOLVED
  climate.in = future_climate %>% 
    ungroup() %>%
    filter(CS == CS.in & ssp == ssp.in) %>%
    mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
           waib = 1/(1 + wai), t = c(1:dim(.)[1])) %>% 
    dplyr::select(sgdd, wai, sgddb, waib, wai2, sgdd2, t)
  
  # Build disturbance df
  disturbance.df = data.frame(type = dist.type.in, IsSurv = FALSE, t = as.numeric(
    trunc(quantile(climate.in$t, probs = c(1:(dist.n.in))/(dist.n.in+1))))) %>%
    mutate(intensity = dist.Iref.in + dist.Islope.in*t)
  
  # Get equilibrium simulation
  sim.equil.in = readRDS(sim_equilibrium_mu[ID.species.in])
  
  # Checked that the population did reach equilibrium
  reached_equil = !is.na(sum((sim.equil.in %>%
                                filter(var == "BAsp") %>%
                                filter(time == max(.$time) - 1))$value))
  
  # Continue only if reached equilibrium
  if(reached_equil){
    
    # Get the distribution at equilibrium
    distrib.in <- filter(sim.equil.in, equil, var == "n") %>% pull(value) * 0.03
    
    # Load species mu object
    sp.in = readRDS(species_mu[ID.species.in])
    
    # Change initial distribution
    sp.in$init_pop <- def_init_k(distrib.in*0.03)
    
    # Make forest
    forest.in = forest(species = list(mu = sp.in), harv_rules = c(
      Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1))
    
    # Run simulation with a changing climate and disturbance
    sim.dist.in = sim_deter_forest(
      forest.in, tlim = max(climate.in$t), climate = climate.in, 
      equil_dist = max(climate.in$t),  equil_time = max(climate.in$t), 
      verbose = TRUE, correction = "cut", disturbance = disturbance.df)
    
  } else {
    sim.dist.in = matrix()
  }
  
  # Name of the file to save
  file.in = paste0("rds/", CS.in, "/sim_disturbance_mu/", species.in, "_", ssp.in, 
                   "_", dist.n.in, dist.type.in, "dist_", dist.Iref.in, "Iref_", 
                   dist.Islope.in, "Islope.rds")
  
  
  # Save species object in a rdata
  create_dir_if_needed(file.in)
  saveRDS(sim.dist.in, file.in)
  
  # Return output list
  return(file.in)
  
}




