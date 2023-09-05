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
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- CURRENT CLIMATE ---- 
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Load and format species share data file
  tar_target(species_share_file, "data/species_share_CS.xlsx", format = "file"), 
  tar_target(species_share_raw, read_xlsx(species_share_file)), 
  tar_target(species_share, format_species_share_raw(species_share_raw)), 
  
  # Load and format functional traits file
  tar_target(traits_file, "data/traits.csv", format = "file"),
  tar_target(traits, fread(traits_file)),
  tar_target(pc1_per_species, get_pc1_per_species(traits)),
  
  # Load and format climate data file
  tar_target(climate_file, "data/RESONATE_climate.csv", format = "file"),
  tar_target(climate_raw, fread(climate_file)), 
  tar_target(climate_ref, get_climate_ref(climate_raw)),

  # Create list of case studies (objects with information about clim and sp)
  tar_target(CS, make_CS(species_share, climate_ref)),
  
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
  
  # Plot info on case studies
  tar_target(fig_nperrichness, plot_n_per_richness(
    CS, "output/fig/fig_nperrichness.jpg"), format = "file"),

  # Plot results
  # -- Vector of color for plotting
  tar_target(CS_color.vec, create_CS_color.vec(CS)),
  # -- Plot results of the simulations
  tar_target(fig_CS_metrics, plot_div_struct_resil_CS(
    FD, resilience, CS_color.vec, "output/fig/metrics_per_CS.jpg"), format = "file"),
  # -- Plot the effect of species composition on resilience
  tar_target(fig_spcompo_vs_resilience, plot_FD_effect_resilience(
    FD, resilience, CS_color.vec, "output/fig/spcompo_vs_resilience"), format = "file"), 
  
  
  
  
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- FUTURE CLIMATE ---- 
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Load and format future climate data
  tar_target(future_climate_all_file, "data/future_clim_all.csv", format = "file"), 
  tar_target(future_climate_all, fread(future_climate_all_file)), 
  tar_target(future_climate, format_future_climate(future_climate_all)),
  
  
  # -- Make species mu objects via branching over ID.species
  tar_target(species_mu, make_species_mu_rds(future_climate, species_list, ID.species),
             pattern = map(ID.species), iteration = "vector", format = "file"),
  
  # -- Make mu simulation till equil via branching over ID.species
  tar_target(sim_equilibrium_mu, make_simulations_equilibrium_mu(
    future_climate, species_list, species_mu, ID.species),
             pattern = map(ID.species), iteration = "vector", format = "file"),
  
  # -- Build crossed disturbance and climate scenarios
  tar_target(df_scenarios, expand.grid(
    ID.species = species_list$ID.species, ssp = c("ssp126", "ssp370", "ssp585"), 
    dist.n = c(1:3), dist.type = "storm", dist.Iref = 0.3, dist.Islope = 0.005*c(0:2))),
  tar_target(ID.scenarios, c(1:dim(df_scenarios)[1])),
  
  # Make simulations with disturbance and changing climate
  tar_target(sim_disturbance_mu, make_simulations_disturbance_mu(
    species_list, future_climate, species_mu, sim_equilibrium_mu, ID.scenarios, 
    df_scenarios), pattern = map(ID.scenarios), iteration = "vector", format = "file"),
  
  # Plot species margins vs climate in case studies
  tar_target(fig_margin_vs_climate_2020, plot_margin_vs_climate(
    future_climate_all, species_list, 2020, "output/fig/fig_margin_vs_clim2020.jpg"), 
    format = "file"), 
  tar_target(fig_margin_vs_climate_2100, plot_margin_vs_climate(
    future_climate_all, species_list, 2100, "output/fig/fig_margin_vs_clim2100.jpg"), 
    format = "file")
  
  
)

