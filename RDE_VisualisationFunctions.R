
#################################################
### VISUALISATION FUNCTIONS FOR SEIR EPIDEMICS ##
#################################################

library(gridExtra)
library(MASS)

#~~~~~~~~~~~~~#
# TRACE PLOTS #
#~~~~~~~~~~~~~#


gg_trace_plot_inf <- function(results, other_results, params_true, burn_in, b_c_y_limits, g_y_limits, F_y_limits, e_y_limits){
  
  # create data frame of results
  
  df = data.frame(samples = other_results[, 11], b_c = results[, 1], g = results[, 2], F_ = results[, 3], e = results[, 4])
  
  # Calc number of samples
  
  num_samples = nrow(df)
  
  # Highlight burn in
  
  df_burn_in <- df %>% 
    mutate(burnin = c(rep("yes", burn_in), rep("no", (num_samples - burn_in))))
  
  
  #Create beta_c plot
  
  beta_c_plot = df_burn_in %>% 
    ggplot(aes(x = samples, y = b_c)) +
    geom_line(size  = 0.2, aes(colour = burnin)) + 
    scale_colour_manual(values = c("dodgerblue", "#fd8f24")) +
    geom_hline(aes(yintercept = params_true[1]), size = 0.7, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(b_c_y_limits[1], b_c_y_limits[2])) +
    labs(
      x = "Sample index", 
      y = "Value", 
      title = expression("Traceplot of parameter samples" ~ beta[c])) + 
    guides(colour = "none") 
  
  
  # Create delta plot
  
  g_plot = df_burn_in %>% 
    ggplot(aes(x = samples, y = g)) +
    geom_line(size  = 0.2, aes(colour = burnin)) + 
    scale_colour_manual(values = c("dodgerblue", "#fd8f24")) +
    geom_hline(aes(yintercept = params_true[3]), size = 0.7, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(g_y_limits[1], g_y_limits[2])) +
    labs(
      x = "Sample index", 
      y = "Value", 
      title = expression("Traceplot of parameter samples" ~ delta)) + 
    guides(colour = "none") 
  
  
  #Create F plot
  
  F_plot = df_burn_in %>% 
    ggplot(aes(x = samples, y = F_)) +
    geom_line(size  = 0.2, aes(colour = burnin)) + 
    scale_colour_manual(values = c("dodgerblue", "#fd8f24")) +
    geom_hline(aes(yintercept = params_true[4]), size = 0.7, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(F_y_limits[1], F_y_limits[2])) +
    labs(
      x = "Sample index", 
      y = "Value", 
      title = expression("Traceplot of parameter samples" ~ F)) + 
    guides(colour = "none") 
  
  
  #Create e plot
  
  e_plot = df_burn_in %>% 
    ggplot(aes(x = samples, y = e)) +
    geom_line(size  = 0.2, aes(colour = burnin)) + 
    scale_colour_manual(values = c("dodgerblue", "#fd8f24")) +
    geom_hline(aes(yintercept = params_true[5]), size = 0.7, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(e_y_limits[1], e_y_limits[2])) +
    labs(
      x = "Sample index", 
      y = "Value", 
      title = expression("Traceplot of parameter samples" ~ epsilon)) + 
    guides(colour = "none") 
  
  
  
  
  # Arrange plots into a grid
  
  trace_plots = grid.arrange(beta_c_plot, g_plot, F_plot, e_plot, nrow = 4)
  
  # use ggsave() to save it, e.g.
  # ggsave(filename = "F_augs_trace.png", width = 10, height = 8) ##### HERE IS CHANGE NAME
  
  return(trace_plots)
}





gg_trace_plot_det <- function(results, other_results, params_true, burn_in, r_y_limits, re_y_limits){
  
  # create data frame of results
  
  df = data.frame(samples = other_results[, 11], r = results[, 5], re = results[, 6])
  
  # Calc number of samples
  
  num_samples = nrow(df)
  
  # Highlight burn in
  
  df_burn_in <- df %>% 
    mutate(burnin = c(rep("yes", burn_in), rep("no", (num_samples - burn_in))))
  
  
  #Create r plot
  
  r_plot = df_burn_in %>% 
    ggplot(aes(x = samples, y = r)) +
    geom_line(size  = 0.2, aes(colour = burnin)) + 
    scale_colour_manual(values = c("dodgerblue", "#fd8f24")) +
    geom_hline(aes(yintercept = params_true[6]), size = 0.7, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(r_y_limits[1], r_y_limits[2])) +
    labs(
      x = "Sample index", 
      y = "Value", 
      title = expression("Traceplot of parameter samples" ~ rho)) + 
    guides(colour = "none") 
  
  
  
  # Create re plot
  
  re_plot = df_burn_in %>% 
    ggplot(aes(x = samples, y = re)) +
    geom_line(size  = 0.2, aes(colour = burnin)) + 
    scale_colour_manual(values = c("dodgerblue", "#fd8f24")) +
    geom_hline(aes(yintercept = params_true[7]), size = 0.7, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(re_y_limits[1], re_y_limits[2])) +
    labs(
      x = "Sample index", 
      y = "Value", 
      title = expression("Traceplot of parameter samples" ~ rho[E])) + 
    guides(colour = "none") 
  
  
  
  
  # Arrange plots into a grid
  
  trace_plots = grid.arrange(r_plot, re_plot, nrow = 2)
  
  # use ggsave() to save it, e.g.
  # ggsave(filename = "F_augs_trace.png", width = 10, height = 8) ##### HERE IS CHANGE NAME
  
  return(trace_plots)
}





#~~~~~~~~~~~~#
# HISTOGRAMS #
#~~~~~~~~~~~~#

gg_hist_plot_inf <- function(results, other_results, params_true, burn_in, b_c_x_limits, g_x_limits, F_x_limits, e_x_limits){
  
  # create data frame of results
  
  df = data.frame(samples = other_results[, 11], b_c = results[, 1], g = results[, 2], F_ = results[, 3], e = results[, 4])
  
  # Calc number of samples
  
  num_samples = nrow(df)
  
  # Highlight burn in
  
  df_burn_in <- df[-(1:burn_in), ]
  
  
  #Create beta_c plot
  
  beta_c_plot = df_burn_in %>% 
    ggplot(aes(x = b_c, y = ..density..)) +
    geom_histogram(fill = "#f5c04a", colour = "grey15", alpha = 0.85) +
    geom_vline(aes(xintercept = params_true[1]), size = 1.1, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), 
                       limits = c(b_c_x_limits[1], b_c_x_limits[2])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)
    ) +
    labs(
      x = expression("Samples of parameter" ~ beta[c]), ##### HERE IS CHANGE SYMBOL
      y = "Density") 
  
  
  #Create delta plot
  
  g_plot = df_burn_in %>% 
    ggplot(aes(x = g, y = ..density..)) +
    geom_histogram(fill = "#f5c04a", colour = "grey15", alpha = 0.85) +
    geom_vline(aes(xintercept = params_true[3]), size = 1.1, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), 
                       limits = c(g_x_limits[1], g_x_limits[2])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)
    ) +
    labs(
      x = expression("Samples of parameter" ~ delta), ##### HERE IS CHANGE SYMBOL
      y = "Density") 
  
  
  #Create F plot
  
  F_plot = df_burn_in %>% 
    ggplot(aes(x = F_, y = ..density..)) +
    geom_histogram(fill = "#f5c04a", colour = "grey15", alpha = 0.85) +
    geom_vline(aes(xintercept = params_true[4]), size = 1.1, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), 
                       limits = c(F_x_limits[1], F_x_limits[2])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)
    ) +
    labs(
      x = expression("Samples of parameter" ~ F), ##### HERE IS CHANGE SYMBOL
      y = "Density") 
  
  
  #Create e plot
  
  e_plot = df_burn_in %>% 
    ggplot(aes(x = e, y = ..density..)) +
    geom_histogram(fill = "#f5c04a", colour = "grey15", alpha = 0.85) +
    geom_vline(aes(xintercept = params_true[5]), size = 1.1, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), 
                       limits = c(e_x_limits[1], e_x_limits[2])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)
    ) +
    labs(
      x = expression("Samples of parameter" ~ epsilon), ##### HERE IS CHANGE SYMBOL
      y = "Density") 
  
  
  hist_plots = grid.arrange(beta_c_plot, g_plot, F_plot, e_plot, nrow = 2)
  
  return(hist_plots)
}



gg_hist_plot_det <- function(results, other_results, params_true, burn_in, r_x_limits, re_x_limits){
  
  # create data frame of results
  
  df = data.frame(samples = other_results[, 11], r = results[, 5], re = results[, 6])
  
  # Calc number of samples
  
  num_samples = nrow(df)
  
  # Highlight burn in
  
  df_burn_in <- df[-(1:burn_in), ]
  
  
  #Create beta_c plot
  
  r_plot = df_burn_in %>% 
    ggplot(aes(x = r, y = ..density..)) +
    geom_histogram(fill = "#f5c04a", colour = "grey15", alpha = 0.85) +
    geom_vline(aes(xintercept = params_true[6]), size = 1.1, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), 
                       limits = c(r_x_limits[1], r_x_limits[2])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)
    ) +
    labs(
      x = expression("Samples of parameter" ~ rho), ##### HERE IS CHANGE SYMBOL
      y = "Density") 
  
  
  #Create beta_b plot
  
  re_plot = df_burn_in %>% 
    ggplot(aes(x = re, y = ..density..)) +
    geom_histogram(fill = "#f5c04a", colour = "grey15", alpha = 0.85) +
    geom_vline(aes(xintercept = params_true[7]), size = 1.1, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), 
                       limits = c(re_x_limits[1], re_x_limits[2])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)
    ) +
    labs(
      x = expression("Samples of parameter" ~ rho[E]), ##### HERE IS CHANGE SYMBOL
      y = "Density") 
  
  
  
  hist_plots = grid.arrange(r_plot, re_plot, nrow = 1)
  
  return(hist_plots)
}






#~~~~~~~~~~#
# CONTOURS #
#~~~~~~~~~~#


gg_contour_plot_inf <- function(results, other_results, params_true, burn_in, b_c_limits, g_limits, F_limits, e_limits){
  
  # create data frame of results
  
  df = data.frame(samples = other_results[, 11], b_c = results[, 1], g = results[, 2], F_ = results[, 3], e = results[, 4])
  
  # Calc number of samples
  
  num_samples = nrow(df)
  
  # Highlight burn in
  
  df_burn_in <- df[-(1:burn_in), ]
  
  
  
  
  # Beta_c vs delta posterior max
  
  
  hexplot_init = ggplot(df_burn_in, aes(x=b_c, y=g) ) +
    geom_point(color = NA, fill = NA) +
    geom_hex(bins = 50) +
    scale_fill_continuous(type = "viridis")
  
  meta_data = ggplot_build(hexplot_init)$data
  
  mean_x2 = meta_data[[2]][which.max(meta_data[[2]]$density),]$x
  mean_y2 = meta_data[[2]][which.max(meta_data[[2]]$density),]$y
  
  
  
  # Beta_c vs delta plot
  
  bc_g_plot = df_burn_in %>% 
    ggplot(aes(x = b_c, y = g)) +
    stat_density2d(aes(fill = ..level..), geom = "polygon", h = 0.0005) +
    scale_fill_continuous(type = "viridis") +
    geom_vline(aes(xintercept = mean_x2), size = 0.7, linetype = 2, 
               colour = 'yellow')+
    geom_hline(aes(yintercept = mean_y2), size = 0.7, linetype = 2, 
               colour = 'yellow') +  
    geom_vline(aes(xintercept = params_true[1]), size = 0.7, linetype = 2, 
               colour = 'red')+
    geom_hline(aes(yintercept = params_true[3]), size = 0.7, linetype = 2, 
               colour = 'red') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(min(df_burn_in$b_c), max(df_burn_in$b_c))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(min(df_burn_in$g), max(df_burn_in$g))) +
    coord_cartesian(xlim = c(b_c_limits[1], b_c_limits[2]),  ylim = c(g_limits[1], g_limits[2])) +
    labs(
      x = as.expression(bquote(beta[c])), 
      y = as.expression(bquote(delta)), 
      title = expression(beta[c] ~ "vs." ~ delta)) + 
    guides(colour = "none") 
  
  
  
  
  
  # Beta_c vs F posterior max
  
  
  hexplot_init = ggplot(df_burn_in, aes(x=b_c, y=F_) ) +
    geom_point(color = NA, fill = NA) +
    geom_hex(bins = 50) +
    scale_fill_continuous(type = "viridis")
  
  meta_data = ggplot_build(hexplot_init)$data
  
  mean_x3 = meta_data[[2]][which.max(meta_data[[2]]$density),]$x
  mean_y3 = meta_data[[2]][which.max(meta_data[[2]]$density),]$y
  
  
  
  # Beta_c vs F plot
  
  bc_F_plot = df_burn_in %>% 
    ggplot(aes(x = b_c, y = F_)) +
    stat_density2d(aes(fill = ..level..), geom = "polygon", h = 0.0004) +
    scale_fill_continuous(type = "viridis") +
    geom_vline(aes(xintercept = mean_x3), size = 0.7, linetype = 2, 
               colour = 'yellow')+
    geom_hline(aes(yintercept = mean_y3), size = 0.7, linetype = 2, 
               colour = 'yellow') +  
    geom_vline(aes(xintercept = params_true[1]), size = 0.7, linetype = 2, 
               colour = 'red')+
    geom_hline(aes(yintercept = params_true[4]), size = 0.7, linetype = 2, 
               colour = 'red') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(min(df_burn_in$b_c), max(df_burn_in$b_c))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(min(df_burn_in$F_), max(df_burn_in$F_))) +
    coord_cartesian(xlim = c(b_c_limits[1], b_c_limits[2]),  ylim = c(F_limits[1], F_limits[2])) +
    labs(
      x = as.expression(bquote(beta[c])), 
      y = as.expression(bquote(F)), 
      title = expression(beta[c] ~ "vs." ~ F)) + 
    guides(colour = "none") 
  
  
  
  
  
  # Beta_c vs epsilon posterior max
  
  
  hexplot_init = ggplot(df_burn_in, aes(x=b_c, y=e) ) +
    geom_point(color = NA, fill = NA) +
    geom_hex(bins = 50) +
    scale_fill_continuous(type = "viridis")
  
  meta_data = ggplot_build(hexplot_init)$data
  
  mean_x4 = meta_data[[2]][which.max(meta_data[[2]]$density),]$x
  mean_y4 = meta_data[[2]][which.max(meta_data[[2]]$density),]$y
  
  
  
  # Beta_c vs epsilon plot
  
  bc_e_plot = df_burn_in %>% 
    ggplot(aes(x = b_c, y = e)) +
    stat_density2d(aes(fill = ..level..), geom = "polygon", h = 0.0007) +
    scale_fill_continuous(type = "viridis") +
    geom_vline(aes(xintercept = mean_x4), size = 0.7, linetype = 2, 
               colour = 'yellow')+
    geom_hline(aes(yintercept = mean_y4), size = 0.7, linetype = 2, 
               colour = 'yellow') +  
    geom_vline(aes(xintercept = params_true[1]), size = 0.7, linetype = 2, 
               colour = 'red')+
    geom_hline(aes(yintercept = params_true[5]), size = 0.7, linetype = 2, 
               colour = 'red') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(min(df_burn_in$b_c), max(df_burn_in$b_c))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(min(df_burn_in$e), max(df_burn_in$e))) +
    coord_cartesian(xlim = c(b_c_limits[1], b_c_limits[2]),  ylim = c(e_limits[1], e_limits[2])) +
    labs(
      x = as.expression(bquote(beta[c])), 
      y = as.expression(bquote(epsilon)), 
      title = expression(beta[c] ~ "vs." ~ epsilon)) + 
    guides(colour = "none") 
  
  
  
  

  # delta vs F posterior max
  
  
  hexplot_init = ggplot(df_burn_in, aes(x=g, y=F_) ) +
    geom_point(color = NA, fill = NA) +
    geom_hex(bins = 50) +
    scale_fill_continuous(type = "viridis")
  
  meta_data = ggplot_build(hexplot_init)$data
  
  mean_x2 = meta_data[[2]][which.max(meta_data[[2]]$density),]$x
  mean_y2 = meta_data[[2]][which.max(meta_data[[2]]$density),]$y
  
  
  
  # delta vs F plot
  
  g_F_plot = df_burn_in %>% 
    ggplot(aes(x = g, y = F_)) +
    stat_density2d(aes(fill = ..level..), geom = "polygon", h = 0.0003) +
    scale_fill_continuous(type = "viridis") +
    geom_vline(aes(xintercept = mean_x2), size = 0.7, linetype = 2, 
               colour = 'yellow')+
    geom_hline(aes(yintercept = mean_y2), size = 0.7, linetype = 2, 
               colour = 'yellow') +  
    geom_vline(aes(xintercept = params_true[3]), size = 0.7, linetype = 2, 
               colour = 'red')+
    geom_hline(aes(yintercept = params_true[4]), size = 0.7, linetype = 2, 
               colour = 'red') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(min(df_burn_in$g), max(df_burn_in$g))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(min(df_burn_in$F_), max(df_burn_in$F_))) +
    coord_cartesian(xlim = c(g_limits[1], g_limits[2]),  ylim = c(F_limits[1], F_limits[2])) +
    labs(
      x = as.expression(bquote(delta)), 
      y = as.expression(bquote(F)), 
      title = expression(delta ~ "vs." ~ F)) + 
    guides(colour = "none") 
  
  
  
  
  
  # delta vs E posterior max
  
  
  hexplot_init = ggplot(df_burn_in, aes(x=g, y=e) ) +
    geom_point(color = NA, fill = NA) +
    geom_hex(bins = 50) +
    scale_fill_continuous(type = "viridis")
  
  meta_data = ggplot_build(hexplot_init)$data
  
  mean_x3 = meta_data[[2]][which.max(meta_data[[2]]$density),]$x
  mean_y3 = meta_data[[2]][which.max(meta_data[[2]]$density),]$y
  
  
  
  # delta vs E plot
  
  g_e_plot = df_burn_in %>% 
    ggplot(aes(x = g, y = e)) +
    stat_density2d(aes(fill = ..level..), geom = "polygon", h = 0.0006) +
    scale_fill_continuous(type = "viridis") +
    geom_vline(aes(xintercept = mean_x3), size = 0.7, linetype = 2, 
               colour = 'yellow')+
    geom_hline(aes(yintercept = mean_y3), size = 0.7, linetype = 2, 
               colour = 'yellow') +  
    geom_vline(aes(xintercept = params_true[3]), size = 0.7, linetype = 2, 
               colour = 'red')+
    geom_hline(aes(yintercept = params_true[5]), size = 0.7, linetype = 2, 
               colour = 'red') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(min(df_burn_in$g), max(df_burn_in$g))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(min(df_burn_in$e), max(df_burn_in$e))) +
    coord_cartesian(xlim = c(g_limits[1], g_limits[2]),  ylim = c(e_limits[1], e_limits[2])) +
    labs(
      x = as.expression(bquote(delta)), 
      y = as.expression(bquote(epsilon)), 
      title = expression(delta ~ "vs." ~ epsilon)) + 
    guides(colour = "none") 
  
  
  
  
  
  # F vs E posterior max
  
  
  hexplot_init = ggplot(df_burn_in, aes(x=F_, y=e) ) +
    geom_point(color = NA, fill = NA) +
    geom_hex(bins = 50) +
    scale_fill_continuous(type = "viridis")
  
  meta_data = ggplot_build(hexplot_init)$data
  
  mean_x4 = meta_data[[2]][which.max(meta_data[[2]]$density),]$x
  mean_y4 = meta_data[[2]][which.max(meta_data[[2]]$density),]$y
  
  
  
  # F vs epsilon plot
  
  F_e_plot = df_burn_in %>% 
    ggplot(aes(x = F_, y = e)) +
    stat_density2d(aes(fill = ..level..), geom = "polygon", h = 0.0003) +
    scale_fill_continuous(type = "viridis") +
    geom_vline(aes(xintercept = mean_x4), size = 0.7, linetype = 2, 
               colour = 'yellow')+
    geom_hline(aes(yintercept = mean_y4), size = 0.7, linetype = 2, 
               colour = 'yellow') +  
    geom_vline(aes(xintercept = params_true[4]), size = 0.7, linetype = 2, 
               colour = 'red')+
    geom_hline(aes(yintercept = params_true[5]), size = 0.7, linetype = 2, 
               colour = 'red') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(min(df_burn_in$F_), max(df_burn_in$F_))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(min(df_burn_in$e), max(df_burn_in$e))) +
    coord_cartesian(xlim = c(F_limits[1], F_limits[2]),  ylim = c(e_limits[1], e_limits[2])) +
    labs(
      x = as.expression(bquote(F)), 
      y = as.expression(bquote(epsilon)), 
      title = expression(F ~ "vs." ~ epsilon)) + 
    guides(colour = "none") 
  
  
  
  trace_plots = grid.arrange(bc_g_plot, bc_F_plot, bc_e_plot, g_F_plot, g_e_plot, F_e_plot, nrow = 3)
  
  # use ggsave() to save it, e.g.
  # ggsave(filename = "F_augs_trace.png", width = 10, height = 8) ##### HERE IS CHANGE NAME
  
  
  return(trace_plots)
}









gg_contour_plot_det <- function(results, other_results, params_true, burn_in, r_limits, re_limits){
  
  # create data frame of results
  
  df = data.frame(samples = other_results[, 11], r = results[, 5], re = results[, 6])
  
  # Calc number of samples
  
  num_samples = nrow(df)
  
  # Highlight burn in
  
  df_burn_in <- df[-(1:burn_in), ]
  
  
  
  
  
  
  # r vs re posterior max
  
  
  hexplot_init = ggplot(df_burn_in, aes(x=r, y=re) ) +
    geom_point(color = NA, fill = NA) +
    geom_hex(bins = 50) +
    scale_fill_continuous(type = "viridis")
  
  meta_data = ggplot_build(hexplot_init)$data
  
  mean_x5 = meta_data[[2]][which.max(meta_data[[2]]$density),]$x
  mean_y5 = meta_data[[2]][which.max(meta_data[[2]]$density),]$y
  
  
  
  # r vs re plot
  
  r_re_plot = df_burn_in %>% 
    ggplot(aes(x = r, y = re)) +
    stat_density2d(aes(fill = ..level..), geom = "polygon", h = 0.005) +
    scale_fill_continuous(type = "viridis") +
    geom_vline(aes(xintercept = mean_x5), size = 0.7, linetype = 2, 
               colour = 'yellow')+
    geom_hline(aes(yintercept = mean_y5), size = 0.7, linetype = 2, 
               colour = 'yellow') +  
    geom_vline(aes(xintercept = params_true[6]), size = 0.7, linetype = 2, 
               colour = 'red')+
    geom_hline(aes(yintercept = params_true[7]), size = 0.7, linetype = 2, 
               colour = 'red') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(min(df_burn_in$r), max(df_burn_in$r))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(min(df_burn_in$re), max(df_burn_in$re))) +
    coord_cartesian(xlim = c(r_limits[1], r_limits[2]),  ylim = c(re_limits[1], re_limits[2])) +
    labs(
      x = as.expression(bquote(rho)), 
      y = as.expression(bquote(rho[E])), 
      title = expression(rho ~ "vs." ~ rho[E])) + 
    guides(colour = "none") 

  

  
  trace_plots = grid.arrange(r_re_plot, nrow = 1)
  
  # use ggsave() to save it, e.g.
  # ggsave(filename = "F_augs_trace.png", width = 10, height = 8) ##### HERE IS CHANGE NAME
  
  
  return(trace_plots)
}











