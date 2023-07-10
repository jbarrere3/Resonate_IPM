#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name _targets.R  
#' @description R script to launch the target pipeline
#' @author Julien BARRERE
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options and packages ----------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load targets
library(targets)
# Load functions
lapply(grep("R$", list.files("R"), value = TRUE), function(x) source(file.path("R", x)))
# install if needed and load packages
packages.in <- c("dplyr", "ggplot2", "matreex", "tidyr", "readxl", "modi", "cowplot",
                 "Taxonstand", "WorldFlora", "data.table", "factoextra", "FD", 
                 "xtable")
for(i in 1:length(packages.in)){
  if(!(packages.in[i] %in% rownames(installed.packages()))){
    install.packages(packages.in[i])
  }
}  
# Targets options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in,
               memory = "transient")
future::plan(future::multisession, workers = 6)
set.seed(2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  # Load data file
  tar_target(species_share_file, "data/species_share_CS.xlsx", format = "file"), 
  tar_target(species_share_raw, read_xlsx(species_share_file)), 
  
  # Load functional traits file
  tar_target(traits_file, "data/traits.csv", format = "file"),
  tar_target(traits, fread(traits_file)),
  tar_target(pc1_per_species, get_pc1_per_species(traits)),
  
  # Format species data
  tar_target(species_share, format_species_share_raw(species_share_raw)), 
  
  # Create list of case studies (objects with information about clim and sp)
  tar_target(CS, make_CS(species_share)), 
  
  # Make species objects
  # -- List of species objects (i.e., IPM) to make
  tar_target(species_list, make_species_list(CS)), 
  # -- Vector from 1 to number of species to make (useful for parallel computing)
  tar_target(ID.species, species_list$ID.species), 
  # -- Make species via branching over ID.species
  tar_target(species, make_species_rds(CS, species_list, ID.species), 
    pattern = map(ID.species), iteration = "vector", format = "file"),
  
  # Make simulations till equilibrium
  # -- Start with a list of forest to simulate
  tar_target(forest_list, make_forest_list(CS)), 
  # -- Make a vector of ID for each forest to simulate
  tar_target(ID.forest, forest_list$ID.forest),
  # -- Make simulations till equilibrium
  tar_target(sim_equilibrium, make_simulations_equilibrium(
    CS, species_list, forest_list, species, ID.forest), 
    pattern = map(ID.forest), iteration = "vector", format = "file"), 
  
  # Make simulations with disturbance
  # -- Data.frame containing storm disturbance to apply
  tar_target(disturbance.df, data.frame(
    type = rep("storm", 3), intensity = rep(0.5, 3), 
    IsSurv = rep(FALSE, 3), t = c(500:502))), 
  # -- Make simulations with disturbance
  tar_target(sim_disturbance, make_simulations_disturbance(
    CS, species_list, forest_list, species, sim_equilibrium, 
    ID.forest.in = ID.forest, disturbance.df), 
    pattern = map(ID.forest), iteration = "vector", format = "file"),
  
  # Extract results
  # -- Get functional diversity
  tar_target(FD, get_FD(forest_list, sim_disturbance, pc1_per_species)),
  # -- Get resilience metrics
  tar_target(resilience, get_resilience_metrics(
    sim_disturbance, disturbance.df, forest_list)), 
  
  # Plot results
  # -- Vector of color for plotting
  tar_target(CS_color.vec, create_CS_color.vec(CS)),
  # -- Plot results of the simulations
  tar_target(fig_CS_metrics, plot_div_struct_resil_CS(
    FD, resilience, CS_color.vec, "output/fig/metrics_per_CS.jpg"), format = "file"), 
  # -- Plot the effect of species composition on resilience
  tar_target(fig_spcompo_vs_resilience, plot_FD_effect_resilience(
    FD, resilience, CS_color.vec, "output/fig/spcompo_vs_resilience"), format = "file"), 
  # -- Plot the climate in each case study
  tar_target(fig_CS_climate, plot_CS_climate(
    CS, CS_color.vec, "output/fig/CS_climate.jpg"), format = "file")
  
)

