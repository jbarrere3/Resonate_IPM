#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_plot.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien Barrere
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Create vector of color, one color per CS
#' @param CS
create_CS_color.vec = function(CS){
  
  # Vector of color
  color.vec = c("#FFB703", "#DC2F02", "#81B29A", "#48CAE4", "#005F73", 
                "#6930C3", "#588157", "#55A930", "#57CC99")
  
  # Named vector to use it in ggplot functions scale_color/fill
  names(color.vec) = names(CS)
  
  # Return vector
  return(color.vec)
}


#' Function to plot between-CS differences in composition, structure and resilience
#' @param FD dataset with species_composition data
#' @param resilience dataset with data on resilience and structure at euilibrium
#' @param CS_color.vec Vector of colors
#' @param file.in name of the file to save, including path
plot_div_struct_resil_CS = function(FD, resilience, CS_color.vec, file.in){
  
  # Create directory of the file to save if needed
  create_dir_if_needed(file.in)
  
  # Format data
  data = FD %>%
    # Add name of the case studies
    left_join(data.frame(ID.CS = c(1:length(CS_color.vec)), 
                         CS = names(CS_color.vec)), 
              by = "ID.CS") %>%
    # Add structure and resilience metrics
    left_join((resilience %>% dplyr::select(
      ID.forest, ID.CS, resistance, recovery, resilience, BA_eq, dbh_mean)), 
      by = c("ID.forest", "ID.CS")) %>%
    # Remove outliers
    filter(BA_eq < 100)
  
  
  
  
  #' Internal function to plot boxplot of three variables per case study
  #' @param data df with a column CS (case study) and columns to plot
  #' @param var names of the columns to plot
  #' @param var.names name of the variables to include in the plot
  #' @param title.in title of the plot
  #' @param color.vec names vector containing color codes for each CS
  plot_var = function(data, var, var.names, title.in, color.vec){
    # Rank countries by ascending order of the first var
    CS.in = (data %>%
               dplyr::select("CS", "var" = var[1]) %>%
               group_by(CS) %>%
               summarize(mean = mean(var, na.rm = TRUE)) %>%
               arrange(desc(mean)))$CS
    # Select variables of interest and change CS into a factor
    data.plot = data %>% 
      dplyr::select("CS", var) %>%
      mutate(CS = factor(CS, levels = CS.in))
    # Change column names
    colnames(data.plot)[1+c(1:length(var))] = var.names
    # Make plot
    data.plot %>%
      gather(key = "var", value = "value", var.names) %>%
      mutate(var = factor(var, levels = var.names)) %>%
      ggplot(aes(x = CS, y = value)) + 
      geom_boxplot(color = "black", aes(fill = CS), 
                   outlier.shape = 21, outlier.color = "black") + 
      facet_wrap(~ var, scales = "free_x") + 
      scale_fill_manual(values = color.vec) +
      coord_flip() + 
      ggtitle(title.in) + 
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            strip.background = element_blank(), 
            plot.title = element_text(hjust = 0.5), 
            axis.title = element_blank(), 
            legend.position = "none")
  }
  
  # Plot species composition at equilibrium 
  plot.speciescompo = plot_var(
    data, var = c("H", "FD", "CWM"),
    var.names = c("Shannon\nindex (H)", "Functional\ndiversity (FD)", 
                  "Community\nweighted\nmean (CWM)"),
    title.in = "Species composition at equilibrium", CS_color.vec
  )
  
  # Plot forest structure at equilibrium
  plot.structure = plot_var(
    data, var = c("Nha", "BA_eq", "dbh_mean"),
    var.names = c("Density\n(# trees per ha)", "Basal area\nat equilibrium\n(m2/ha)", 
                  "Mean diameter\nat breast\nheight (mm)"),
    title.in = "Forest structure at equilibrium", CS_color.vec
  )
  
  # Plot forest structure at equilibrium
  plot.resilience = plot_var(
    data, var = c("resistance", "recovery", "resilience"),
    var.names = c("Resistance\n(logratio of BA\nafter/before)", 
                  "Recovery\n(slope of post-dist\nincrease in BA)", 
                  "Resilience\n(integral of difference\nBA(t) - BAeq) "),
    title.in = "Demographic response to storm disturbance", CS_color.vec
  )
  
  # Final plot
  plot.out = plot_grid(plot.speciescompo, plot.structure, plot.resilience, 
                       ncol = 1, align = "v", labels = c("(a)", "(b)", "(c)"), 
                       scale = 0.9)
  
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 17, height = 17, units = "cm", 
         dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
  
  
}



#' Function to estimate of FD metrics on resilience
#' @param FD Data on species composition
#' @param resilience Data on demographic response to disturbance
#' @param CS_color.vec vector of colors for each case study
#' @param dir.in Name of the directory where to save files
plot_FD_effect_resilience = function(FD, resilience, CS_color.vec, dir.in){
  
  # Name of the files
  fig.file.in = paste0(dir.in, "/fig_composition_vs_resilience.jpg")
  table.file.stats = paste0(dir.in, "/stats_composition_vs_resilience.tex")
  
  # create output directory if it doesn't exist
  create_dir_if_needed(fig.file.in)
  
  # Vector of response variables for which to run models
  response.vec = c("resilience", "resistance", "recovery")
  
  # Merge the two dataset
  data_model = FD %>% 
    left_join(resilience, by = c("ID.forest", "ID.CS", "Assemblage")) %>%
    left_join(data.frame(CS = names(CS_color.vec), 
                         ID.CS = c(1:length(CS_color.vec))), 
              by = "ID.CS")
  
  # Data to fit the models
  data.in = cbind(data_model[, c("ID.CS", "CS", response.vec, "H", "FD", "CWM")], 
                  scale((data_model %>% dplyr::select(
                    "H.scaled" = "H", "FD.scaled" = "FD", "CWM.scaled" = "CWM")), 
                    center = TRUE, scale = TRUE)) %>%
    mutate(resilience = log(resilience), 
           recovery = log(recovery))
  
  # Models to unscale explanatory variables
  unscale.H = lm(H ~ H.scaled, data.in)
  unscale.FD = lm(FD ~ FD.scaled, data.in)
  unscale.CWM = lm(CWM ~ CWM.scaled, data.in)
  
  
  # Initialize the table with statistics
  table.stats = data.frame(col1 = "", col2 = c("Resistance", "Est (se)"), col3 = c("", "F"), 
                           col4 = c("", "p"), col5 = "", col6 = c("Recovery", "Est (se)"), 
                           col7 = c("", "F"), col8 = c("", "p"), col9 = "", 
                           col10 = c("Resilience", "Est (se)"), col11 = c("", "F"), 
                           col12 = c("", "p"))
  
  # Loop on all case studies
  for(i in 1:length(names(CS_color.vec))){
    
    # Initialize statistics table for this case study
    table.stats.i = data.frame(col1 = "", col2 = "", col3 = "", col4 = "", col5 = "", 
                               col6 = c("", names(CS_color.vec)[i], "", "", ""), 
                               col7 = "", col8 = "", col9 = "", col10 = "", 
                               col11 = "", col12 = "")
    
    # Loop on all response variables
    for(j in 1:length(response.vec)){
      
      # Fit model
      eval(parse(text = paste0(
        "model.ij = lm(", response.vec[j], 
        " ~ H.scaled + FD.scaled + CWM.scaled, data = subset(data.in, ID.CS == ", 
        i, "))"
      )))
      
      # Output data set for model i j 
      data.out.ij = data.frame(
        ID.CS = i, 
        CS = names(CS_color.vec)[i],
        var.resp = response.vec[j], 
        var.exp = c("H", "FD", "CWM"), 
        var.pos = c(1:3),
        est = as.numeric(coef(model.ij)[-1]), 
        est.low = as.numeric(confint(model.ij)[-1, 1]), 
        est.high = as.numeric(confint(model.ij)[-1, 2])
      )
      
      # Statistics of model ij
      stats.ij = data.frame(
        est = paste0(
          round(summary(model.ij)$coefficients[-1, 1], digits = 2), " (", 
          round(summary(model.ij)$coefficients[-1, 2], digits = 2), ")"
        ), 
        Fval = round(anova(model.ij)[-dim(anova(model.ij))[1], 4], digits = 2), 
        pval = scales::pvalue(anova(model.ij)[-dim(anova(model.ij))[1], 5])
      )
      
      # Complete the statistics table
      if(response.vec[j] == "resistance"){
        table.stats.i[3:5, paste0("col", c(2:4))] = stats.ij}
      if(response.vec[j] == "recovery"){
        table.stats.i[3:5, paste0("col", c(6:8))] = stats.ij}
      if(response.vec[j] == "resilience"){
        table.stats.i[3:5, paste0("col", c(10:12))] = stats.ij}
      
      # Add to the final output dataset
      if(j == 1 & i == 1) data.out = data.out.ij
      else data.out = rbind(data.out, data.out.ij)
      
    }
    
    # Complete statistics table
    table.stats = rbind(table.stats, table.stats.i)
  }
  
  
  # Plot the estimates
  plot.out = data.out %>%
    mutate(significance = ifelse(est.low > 0 | est.high < 0, "yes", "no")) %>%
    mutate(var.resp = factor(var.resp, levels = c("resistance", "recovery", "resilience")), 
           var.exp = factor(var.exp, levels = c("H", "FD", "CWM"))) %>%
    ggplot(aes(x = CS, y = est, color = significance, fill = CS)) + 
    geom_errorbar(aes(ymin = est.low, ymax = est.high),
                  width = 0) + 
    geom_point(shape = 21, size = 3) +
    xlab("") + ylab("Effect on response to\ndisturbance metric") +
    scale_color_manual(values = c(`no` = "gray", `yes` = "black")) +
    scale_fill_manual(values = CS_color.vec) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(var.exp ~ var.resp, scales = "free") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.position = "none") + 
    coord_flip()
  
  # Save plot 
  ggsave(fig.file.in, plot.out, width = 21, height = 11, units = "cm", 
         dpi = 600, bg = "white")
  
  # Save tex file
  print(xtable(table.stats, type = "latex", 
               caption = "Statistics of the models used to test the effect of species composition on forest response to disturbances", 
               label = "table_stat_spcompo_resilience"), 
        include.rownames=FALSE, hline.after = c(0, 1, dim(table.stats)[1]), 
        include.colnames = FALSE, caption.placement = "top", 
        file = table.file.stats)
  
  
  # Return name of the file saved
  return(c(fig.file.in, table.file.stats))
}

#' Plot the climate in each case study
#' @param CS case study
#' @param CS_color.vec vector of color for each case study
#' @param file.in Name of the file to save, including path
plot_CS_climate = function(CS, CS_color.vec, file.in){
  
  # Create directory if necessary
  create_dir_if_needed(file.in)
  
  # Format data for plotting
  data = data.frame(CS = names(CS), sgdd = NA_real_, wai = NA_real_)
  for(i in 1:length(names(CS))) data[i, c("sgdd", "wai")] = CS[[i]]$climate[c("sgdd", "wai")]
  
  # Make the plot
  plot.out = data %>% 
    ggplot(aes(x = wai, y = sgdd)) + 
    geom_point(aes(fill = CS), shape = 21, color = "black", size = 2) +
    scale_fill_manual(values = CS_color.vec) +
    scale_color_manual(values = CS_color.vec) +
    geom_text(aes(label = CS, color = CS), nudge_y = 50, nudge_x = 0.01) + 
    xlim(range(data$wai)*1.4) +
    theme(legend.position = "none", 
          panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank())
  
  # Save plot 
  ggsave(file.in, plot.out, width = 12, height = 9, units = "cm", 
         dpi = 600, bg = "white")
  
  # return name of the file to save
  return(file.in)
  
}