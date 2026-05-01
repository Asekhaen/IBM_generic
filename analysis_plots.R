#Library

# library(lme4)
# library(lmerTest)
# library(dplyr)
# library(ggplot2)
# library(patchwork)
# library(purrr)
# library(data.table)


# customised colour scheme
my_colors <- c("No effect" ="#482173", "Lethal effect" = "#bdb726", "Sterile effect" = "#29af7f")


#Source files to call some simulation parameters for analysis
source("R/dependencies.R")



#===============================================
#  FOUNDER EVENT AND ALLEE EFFECT 
#===============================================

## Load simulation output and clean data 
allee_df <- readRDS("~/Documents/Curtin-PhD/R_and_IBM/phd_codes/0.allee_effect/R/output/allee_effect.rds")

# rename scenarios, remove redundant columns and source patch, rename patch to founder pop size
allee_data <- allee_df |>
  mutate(
    status = case_when(
      lethal_effect == "FALSE" & complete_sterile == "TRUE" ~ "Sterile effect",
      lethal_effect == "TRUE" & complete_sterile == "FALSE" ~ "Lethal effect",
      lethal_effect == "FALSE" & complete_sterile == "FALSE" ~ "No effect",
      TRUE          ~ NA_character_
    )
  ) |> 
  select(-complete_sterile, -lethal_effect, -init_frequency, -time_half_K, -scenario, -patch_occupied) |>
  filter(patch != 1) |>
  mutate(patch = case_when(
    patch == 2 ~ 5,
    patch == 3 ~ 10,
    patch == 4 ~ 15,
    patch == 5 ~ 25,
    patch == 6 ~ 50,
    patch == 7 ~ 75,
    patch == 8 ~ 100,
    TRUE ~ patch   # keep other values unchanged
  ))


#re-order status level as required
allee_data$status <- factor(allee_data$status,
                            levels = c("No effect", "Lethal effect", "Sterile effect"))




## ALLEE EFFECT GROWTH CURVE

#cal mean values; pop, g_rate and their 5% & 95% quantiles
allee_summary <- allee_data |>
  group_by(year, status, patch, n_load, n_loci) |>
  summarise(
    mean_pop = mean(pop_size),
    pop_q5 = quantile(pop_size, 0.05),
    pop_q95 = quantile(pop_size, 0.95),
    pcg = mean(g_rate),
    q5 = quantile(g_rate, 0.05),
    q95 = quantile( g_rate, 0.95),
    .groups = "drop"
  )


#extract some timesteps/generations to observe Allee effect dynamics per generations)
allee_plot_gen <- allee_summary |>
  filter(year %in% c(1, 2, 3, 4, 5, 7, 10, 15, 25, 50),
         n_load == 0.25,
         n_loci == 1000)


  #Figure 3:
  allee_plot_load <- allee_summary |>
  filter(year == 4)
  
  ggplot(allee_plot_load, aes(x = patch, y = pcg, colour = factor(status))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = q5, ymax = q95, fill = factor(status)),
                alpha = 0.2, color = NA, show.legend = FALSE) +
    facet_grid(n_load~n_loci) +
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    scale_x_continuous(breaks = seq(0,100, by = 20))+
    labs(x = "Founder size",
         y = "Per capita growth rate",
         color = "Load effect") +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank())
  
  
  
#Allee effect at generation 4
allee_plot_25 <- allee_summary |>
  filter(n_load == "0.25")



#Figure 4 A and B

#Allee effect curve: growth rate ~ population size

#plot for population dynamics 
A <- ggplot(allee_plot_25, aes(x = year, y = mean_pop, colour = factor(patch), group = factor(patch))) +
  geom_line(size = 1) +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = pop_q5, ymax = pop_q95, fill = factor(patch)),
              alpha = 0.2, color = NA) +
  facet_grid(n_loci ~ status) +
  # facet_wrap (~ status, ncol = 10) +
  scale_colour_viridis_d(option = "viridis") +
  scale_fill_viridis_d(option = "viridis") +
  scale_y_continuous(breaks = seq(100,1000, by = 150))+
  labs(x = "Generations",
       y = "Population size",
       colour = "Founder size",
       fill = "Founder size"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        legend.position = "none")


#plot for growth dynamics 

B <- ggplot(allee_plot_25, aes(x = year, y = pcg, colour = factor(patch), group = factor(patch))) +
  geom_line(size = 1) +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = q5, ymax = q95, fill = factor(patch)),
              alpha = 0.2, color = NA) +
  facet_grid(n_loci ~ status) +
  # facet_wrap (~ status, ncol = 10) +
  scale_colour_viridis_d(option = "viridis") +
  scale_fill_viridis_d(option = "viridis") +
  labs(x = "Generations",
       y = "Growth rate",
       colour = "Founder size",
       fill = "Founder size"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank())




## EXTINCTION PROBABILITY

allee_ext <- allee_data %>%
  mutate(extinct = case_when(
    pop_size == 0 ~ 0,
       TRUE ~ 1   # keep other values unchanged
  ))


#Extinction probability was measured at the end of the simulation (i.e. generation = 50 or 100)
# as proportion of the simulation where population is = 0 by the end of the simulation


allee_ext_summary <- allee_ext |>
  group_by(year, status, patch, n_load, n_loci) |>
  summarise(
    extinct_prob = mean(extinct == 0),
    pop = mean(pop_size),
    g_rate = mean(g_rate),
    .groups = "drop"
  )


#filter the last generation
allee_ext_df <- allee_ext_summary |>
  filter(year == 50)



#Extinct probability Figure S1
ggplot(allee_ext_df, aes(x = patch, y = extinct_prob, colour = factor(status))) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  facet_grid(n_load~n_loci) +
  # facet_wrap(~ n_loci, ncol = 4) +
  scale_color_manual(values = my_colors) +
  #scale_fill_manual(values = my_colors) +
  labs(x = "Founder population",
       y = "Extinction Probability",
       color = "Load effect") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())





#___________________________________
# Genetics: Allee effect
#___________________________________

allee_gen <- readRDS("~/Documents/Curtin-PhD/R_and_IBM/phd_codes/0.allee_effect/R/output/allee_genetic.rds")


allee_gen_scenarios <- allee_gen |>
  mutate(
    status = case_when(
      lethal_effect == "FALSE" & complete_sterile == "TRUE" ~ "Sterile effect",
      lethal_effect == "TRUE" & complete_sterile == "FALSE" ~ "Lethal effect",
      lethal_effect == "FALSE" & complete_sterile == "FALSE" ~ "No effect",
      TRUE          ~ NA_character_
    )
  )|> 
  select(-complete_sterile, -lethal_effect, -init_frequency, -scenario)



#allee_gen_df <- allee_gen_scenarios |> select(-complete_sterile, -lethal_effect, -init_frequency, -scenario)



allee_gen_df$status <- factor(allee_gen_df$status,
                            levels = c("No effect", "Lethal effect", "Sterile effect"))



#estimate homozygosisty proportion and mean frequency of the deleterious allele across population 
allee_div <- allee_gen_df |>
  group_by(patch, year, replicate, n_loci, n_load, status) |>
  summarise(
    load = 1 - prod(1 - freq_a^2),
    .groups = "drop"
  )

#extract load size 0.25 (25%), exclude source population, and rename patches to 
# founder population sizes for analysis

load_gen <- allee_div |>
  filter(
    # status != "No effect",
    n_load ==  0.25,
    patch %in% 2:8
  ) |>
  mutate(patch = case_when(
    patch == 2 ~ 5,
    patch == 3 ~ 10,
    patch == 4 ~ 15,
    patch == 5 ~ 25,
    patch == 6 ~ 50,
    patch == 7 ~ 75,
    patch == 8 ~ 100,
    TRUE ~ patch   # keep other values unchanged
  ))


#mean homozygosisty and 95% quantile
plot_allee_df <- load_gen |>
  group_by(patch, year, n_loci, n_load, status) |>
  summarise(
    mean_load = mean(load),
    load_q5 = quantile(load, 0.05),
    load_q95 = quantile(load, 0.95),
    .groups = "drop"
  )


# Figure 4 C
#plot for homozygosity   (Re run this simulation and plot 1, 3, 5, 10, 25, 50 generations)
C <- ggplot(plot_allee_df, aes(x = year, y = mean_load, colour = factor(patch), group = factor(patch))) +
  geom_line(size = 1) +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = load_q5, ymax = load_q95, fill = factor(patch)),
              alpha = 0.2, color = NA) +
  facet_grid(n_loci ~ status) +
  # facet_wrap (~ status, ncol = 10) +
  scale_colour_viridis_d(option = "viridis") +
  scale_fill_viridis_d(option = "viridis") +
  labs(x = "Generations",
       y = "Homozygous probability",
       colour = "Founder size",
       fill = "Founder size"
       ) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
       legend.position = "none")




#Plot
(A/B/C)



#===============================================
#  TWO-PATCH SOURCE-SINK ANALYSIS AND PLOTS
#===============================================


# patch population statistics

# import data from saved folder
two_patch_data <- readRDS("~/Documents/Curtin-PhD/R_and_IBM/phd_codes/ibm_generic/R/output/two_patch.rds")

# #select only relevant scenarios
# patch_data <- two_patch_data |> 
#   filter(scenario %in% 2:4)

patch_data_df <- two_patch_data |>
  mutate(
      status = case_when(
        lethal_effect == "FALSE" & complete_sterile == "TRUE" ~ "Sterile effect",
        lethal_effect == "TRUE" & complete_sterile == "FALSE" ~ "Lethal effect",
        lethal_effect == "FALSE" & complete_sterile == "FALSE" ~ "No effect",
        TRUE          ~ NA_character_
      )
    )

#remove redundant column
pop_data <- patch_data_df |> select(-scenario, -complete_sterile, -lethal_effect, -patch_occupied)


#Extract sink patch for downstream analysis
patch_df <- pop_data |>
  filter(patch == 2)

patch_df$status <- factor(patch_df$status,
                               levels = c("No effect", "Lethal effect", "Sterile effect"))


#estimate time to establishment 
time_to_K_df <- patch_df |>
  group_by(status, replicate, n_load, n_loci) |>
  mutate(
    time_est = {
      if (any(pop_size >= establishment_threshold, na.rm = TRUE)) {
        min(year[pop_size >= establishment_threshold], na.rm = TRUE)
      } else {
        NA_integer_
      }
    }
  ) |>
  ungroup()

#aggregation of replicates
patch_summary <- time_to_K_df |>
  group_by(year, status, init_frequency, n_load, n_loci) |>
  summarise(mean_pop = mean(pop_size),
            pop5 = quantile(pop_size, 0.05),
            pop95= quantile(pop_size, 0.95),
            mean_g = mean(g_rate),
            g5 = quantile(g_rate, 0.05),
            g95 = quantile(g_rate, 0.95),
            .groups = "drop")



patch_25 <- patch_summary |>
  filter(n_load == 0.25)


#------------------------------------
#Plots
#------------------------------------

#growth 

A <- ggplot(patch_25, aes(x = year, y = mean_g, colour = factor(status))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = g5, ymax = g95, fill = factor(status)), alpha = 0.25, colour = NA,  show.legend = FALSE) +
  scale_colour_manual(values = c(my_colors)) +
  scale_fill_manual(values = c(my_colors)) +
  facet_wrap(~ n_loci, ncol = 4) +
  #facet_grid(n_load ~ n_loci) +
  labs(x = "Generations)", y = "Growth rate", color = "Load effect") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), legend.position = "none")



#population size

B <- ggplot(patch_25, aes(x = year, y = mean_pop, colour = factor(status))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = pop5, ymax = pop95, fill = factor(status)), alpha = 0.25, colour = NA,  show.legend = FALSE) +
  geom_hline(yintercept = 1000, linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = 750,  linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_text(x = max(patch_25$year) * 0.98, y = 1075, label = "100%", 
            size = 3, inherit.aes = FALSE) +
  geom_text(x = max(patch_25$year) * 0.98, y = 800, label = "75%", 
            size = 3, inherit.aes = FALSE) +
  scale_colour_manual(values = c(my_colors)) +
  scale_fill_manual(values = c(my_colors)) +
  facet_wrap(~ n_loci, ncol = 4) +
  #facet_grid(n_load ~ n_loci) +
  scale_y_continuous(breaks = seq(100,1000, by = 100))+
  labs(x = "Generations", y = "Population size", color = "Load effect") +
  coord_cartesian(ylim = c(0, 1150)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


#Time to carrying capacity 


#This removes  NAs and selects only the unique values i.e.
#year the population reach/crossed establishment threshold

time_to_half <- time_to_K_df |>
  group_by(status, replicate, init_frequency, n_load, n_loci) |>
  reframe(
    time_est = first(na.omit(time_half_K)),
    mean_pop = mean(pop_size),
    .groups = "drop"
  )

patch_establish <- time_to_half |>
  filter(n_load == 0.25)

C <- ggplot(patch_establish, aes(x = time_est, y = status, colour = status)) +
  scale_x_continuous(limits = c(0, 100)) +
  geom_boxplot(aes(fill = status), outlier.shape = NA, width = 0.6, alpha = 0.6, show.legend = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  #facet_grid(n_load ~ n_loci) +
  facet_wrap(~ n_loci, ncol = 4) +
  labs(x = "Generations", y = "", color = "Load effect") +
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), legend.position = "none" )



(A/B/C)



# Supplementary figures (all load scenarios)

  ggplot(time_to_half, aes(x = time_est, y = status, colour = status)) +
  scale_x_continuous(limits = c(0, 100)) +
  geom_boxplot(aes(fill = status), outlier.shape = NA, width = 0.6, alpha = 0.6, show.legend = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_grid(n_load ~ n_loci) +
  #facet_wrap(~ n_loci, ncol = 4) +
  labs(x = "Generations", y = "Genetic load", color = "Genetic load") +
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), legend.position = "none" )


  ggplot(patch_summary, aes(x = year, y = mean_pop, colour = factor(status))) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = pop5, ymax = pop95, fill = factor(status)), alpha = 0.25, colour = NA,  show.legend = FALSE) +
    geom_hline(yintercept = 1000, linetype = "dashed", colour = "black", linewidth = 0.5) +
    geom_hline(yintercept = 750,  linetype = "dashed", colour = "black", linewidth = 0.5) +
    geom_text(x = max(patch_25$year) * 0.98, y = 1075, label = "100%", 
              size = 3, inherit.aes = FALSE) +
    geom_text(x = max(patch_25$year) * 0.98, y = 800, label = "75%", 
              size = 3, inherit.aes = FALSE) +
    scale_colour_manual(values = c(my_colors)) +
    scale_fill_manual(values = c(my_colors)) +
    # facet_wrap(~ n_loci, ncol = 4) +
    facet_grid(n_load ~ n_loci) +
    labs(x = "Generations", y = "Population size", color = "Load effect") +
    coord_cartesian(ylim = c(0, 1150)) +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank())
  
  
#------------------------------------
# 2. Genetic structure
#------------------------------------


# # import data from folder
# 
# two_genetic <- readRDS("~/Documents/Curtin-PhD/R_and_IBM/Generic_IBM_Proj/IBM_generic/R/output/two_patch_genetic.rds")
# 
# 
# # #select only relevant scenarios if neccessary
# # two_genetic_df <- two_genetic |> 
# #   filter(scenario %in% 2:4)
# 
# two_genetic_data <- two_genetic |>
#   mutate(
#     status = case_when(
#       scenario == 2 ~ "Sterile",
#       scenario == 3 ~ "Lethal",
#       scenario == 4 ~ "No_load",
#       TRUE          ~ NA_character_
#     )
#   )
# 
# two_genetic_data <- two_genetic_data |> select(-lethal_effect, -complete_sterile)
# 
# gen_df <- two_genetic_data |>
#   filter(patch == 2)
# 
# # # frequency of homozygous deleterious allele per locus 
# 
# # HWE p^2 + 2pq^ + q^2 = 1
# # where p is the frequency of the wildtype allele, q is the frequency of the 
# # deleterious allele, p^2 is the frequency of homozygous wildtype individuals, 
# # 2pq is the frequency of heterozygous recessive individuals, and q^2 is the frequency of 
# # homozygous deleterious individuals (p + q = 1)
#  
# 
# 
# #summary per locus
# per_locus_summary <- gen_df |>
#   group_by(year, status, locus) |>
#   summarise(AA = mean(AA),
#             Aa = mean(Aa),
#             aa = mean(aa),
#             A_freq = mean(freq_A),
#             a_freq = mean(freq_a),
#             homo = mean(a_freq^2),
#             not_homo = ifelse(homo > 0, 1-homo, 0),
#             .groups = "drop")
# 
# per_locus_summary_df <- per_locus_summary |>
#   filter(status != "No_load")
# 
# 
# 
# # #probability that any locus is homozygous deleterious
# # 
# # all_loci <- per_locus_summary_df |>
# #   group_by(year,status) |>
# #   summarise(any_loci_homo = (1 - prod(not_homo)),
# #   .groups = "drop")
# 
# 
# # A. per locus all series visualisation
# ggplot(per_locus_summary_df, aes(x = year, y = a_freq, color = factor(status), group = interaction(status, locus))) +
#   geom_line(size = 0.1) +
#   scale_fill_manual(
#     values = c(my_colors)
#   ) +
#   labs(
#     x = "Generation",
#     y = "Mean deleterious allele frequency", color = "Genetic load") +
#   #facet_wrap(~patch) +
#   theme_minimal()
# 
# 
# # averaged acros loci 
# average_del_freq <- per_locus_summary_df |>
#   group_by(year, status) |>
#   summarize(
#     mean_a = mean(a_freq),
#     min_a = min(a_freq),
#     max_a = max(a_freq),
#     .groups = "drop")
# 
# # combine deleterious alleles with mean population for visualisation 
# 
# avg_pop <- patch_summary |> 
#   filter(status != "No_load")
# 
#  avg_pop_df <- avg_pop|>
#    select(mean_pop)
# 
# combine_df <- bind_cols(average_del_freq, avg_pop_df)
#   
# #B. mean visualisation with minimum and maximum values
# F3 <- ggplot(combine_df, aes(x = mean_pop, y = mean_a, color = factor(status), group = interaction(status))) +
#   geom_line(size = 0.1) +
#   scale_color_manual(
#     values = c(my_colors)) +
#   geom_ribbon(aes(ymin = min_a, ymax = max_a, fill = factor(status)), alpha = 0.2) +
#   labs(
#     x = "Generation",
#     y = "Deleterious allele frequency", color = "Genetic load", fill = "Genetic load") +
#   theme_classic(base_size = 12)+
#   theme(legend.position = "none")
# 
# (F1|F2)/(F3|F4)
# 
# 


#===============================================
#  STEPPING-STONE ANALYSIS AND PLOTS
#===============================================

step_data <- readRDS("~/Documents/Curtin-PhD/R_and_IBM/phd_codes/ibm_generic/R/output/step_data.rds")  
  # (need to remove source patch from all analysis and plot)

  
  # remove source patch and renumber patches
  step_data_df <- step_data |>
    filter(patch != 1) |>
    mutate(patch = patch -1,
           status = case_when(
             lethal_effect == "FALSE" & complete_sterile == "TRUE" ~ "Sterile effect",
             lethal_effect == "TRUE" & complete_sterile == "FALSE" ~ "Lethal effect",
             lethal_effect == "FALSE" & complete_sterile == "FALSE" ~ "No effect",
             TRUE          ~ NA_character_
           )) |>
    select(-scenario,-lethal_effect, -time_half_K, -g_rate, -complete_sterile)
  
  
  #reorder status
  step_data_df$status <- factor(step_data_df$status,
                                levels = c("No effect", "Lethal effect", "Sterile effect"))
  


#Visualise invasion process using mean population size across replicates
step_df <- step_data_df |>
  group_by(year, status, patch, dispersal_frac) |>
  summarise(
    mean_pop = mean(pop_size),
    pop_05 = quantile(pop_size, 0.05),
    pop_95 = quantile(pop_size, 0.95),
    col_pop = ifelse(mean_pop <= colonisation_threshold, NA, mean_pop),
    log_pop = log10(col_pop),
    .groups = "drop"
  )


# Varying dispersal rates
ggplot(step_df, aes( year, patch, fill = col_pop)) +
  geom_tile() +
  facet_grid(status~dispersal_frac) +
  scale_fill_viridis_c(option = "viridis", 
                       limits = c(0, 1050),
                       na.value = "white") +
  scale_y_continuous(breaks = 1:10) +
  scale_x_continuous(breaks = seq(0,100, by = 20)) +
  theme_bw(base_size = 14) +
  labs(
    x = "Generations",
    y = "Patch",
    fill = "Population size"
  )



#Estimate speed using arrival time

   arrival_df <- step_data_df |>
     group_by(status, dispersal_frac, replicate, patch) |>
     summarise(
       arrival = if (any(pop_size >= colonisation_threshold)) {
         min(year[pop_size >= colonisation_threshold])
       } else {
         NA_real_
       },
       .groups = "drop"
     )
   
   

   speed_df <- arrival_df |>
     group_by(status, dispersal_frac, replicate) |>
     summarise(
       speed = coef(lm(patch ~ arrival))[2],          # Patches per time
       r_squared = summary(lm(patch ~ arrival))$r.squared,
       .groups = "drop"
     )
   
   speed_summary <- speed_df |>
     group_by(status, dispersal_frac) |>
     summarise(
       mean_speed = mean(speed, na.rm = TRUE),
       q5 = quantile(speed, 0.05, na.rm = TRUE),
       q95 = quantile(speed, 0.95, na.rm = TRUE),
       .groups = "drop"
     )
   
   
   ggplot(speed_summary,
          aes(dispersal_frac, mean_speed, colour = status)) +
     geom_line(linewidth = 1) +
     geom_point(size = 1) +
     geom_ribbon(
       aes(ymin = q5, ymax = q95, fill = status),
       alpha = 0.2,
       colour = NA
     ) +
     scale_fill_manual(values = my_colors) +
     scale_color_manual(values = my_colors) +
     theme_bw(base_size = 14) +
     labs(
       x = "Dispersal probability",
       y = "Invasion speed",
       fill = "Load effect",
       color = "Load effect") 
   
   
   
   # Compare the invasion speed between load status

     library(lme4)
     library(lmerTest)


     test <- lm(
       speed ~ dispersal_frac * status,
       data = speed_df
     )

     anova(test)
     summary(test)

     plot(test)
     shapiro.test(residuals(test))
   
   
   
#===================================
# Genetics 
#===================================

step_genetic <-  readRDS("~/Documents/Curtin-PhD/R_and_IBM/phd_codes/ibm_generic/R/output/step_gene_df.rds")
   
   
step_gen <- step_genetic |>
  mutate(
    status = case_when(
      lethal_effect == "FALSE" & complete_sterile == "TRUE" ~ "Sterile effect",
      lethal_effect == "TRUE" & complete_sterile == "FALSE" ~ "Lethal effect",
      lethal_effect == "FALSE" & complete_sterile == "FALSE" ~ "No effect",
      TRUE          ~ NA_character_
    )
  )

#remove redundant columns
step_gen <- step_gen |> select(-scenario, -lethal_effect, -complete_sterile)


step_gen_df <- step_gen |>
  mutate(pop = AA + Aa +aa)


#remove source patch and renumber the patches
step_gen_data <- step_gen_df |>
  filter(patch != 1) |>
  mutate(
    patch = patch -1
  )



step_gen_data$status <- factor(step_gen_data$status,
                            levels = c("No effect", "Lethal effect", "Sterile effect"))


# # #save this edited file for reuse to avoid cleaning the heavy data all over
# saveRDS(step_gen_data, file = file.path("R/output", "step_gen_data.rds"))
# step_gen_data <-  readRDS("~/Documents/Curtin-PhD/R_and_IBM/phd_codes/ibm_generic/R/output/step_gen_data.rds")



#calculate time to colonisation for each replicates
step_gen_occ <- step_gen_data |>
  group_by(patch, status, locus, replicate, dispersal_frac) |>
  mutate(
    arrival = {
      if (any(pop >= colonisation_threshold)) {
        min(year[pop >= colonisation_threshold], na.rm = TRUE)
      } else {
        NA_integer_
      }
    }
  ) |>
  ungroup()


# #memory issues, so I store the output for reuse
# saveRDS(step_gen_occ, file = file.path("R/output", "step_gen_occ.rds"))
# step_gen_occ <-  readRDS("~/Documents/Curtin-PhD/R_and_IBM/phd_codes/ibm_generic/R/output/step_gen_occ.rds")


# Retain  95% of the simulations above the 5% cut off 
  step_95_data <- step_gen_occ |>
    filter(arrival >= quantile(arrival, 0.05, na.rm = TRUE))


 # #memory issues, so I store the output for reuse
 # saveRDS(step_95_data, file = file.path("R/output", "step_95_data.rds"))
 # step_95_data <-  readRDS("~/Documents/Curtin-PhD/R_and_IBM/phd_codes/ibm_generic/R/output/step_95_data.rds")
  

 # ## earliest arrival for 95% simulations
 #  step_95_minimum <- step_95_data |>
 #    group_by(patch, year, status, dispersal_frac) |>
 #    summarise(min_arrival = min(arival, na.rm = TRUE), 
 #              .groups = "drop") |>
 #    filter(min_arrival >= quantile(min_arrival, 0.05, na.rm = TRUE))

    
    
## calculate the mean homozygosity for the entire 95% data and mean_L

 homo_df <- step_95_data |>
   group_by(replicate, year, patch, dispersal_frac, status) |>
   summarise(
     pop = first(pop),
     arrival = first(arrival),
     L = 1 - prod(1 - freq_a^2),
     .groups = "drop"
   ) #|>
   # mutate(
   #   col_L = if_else(pop <= colonisation_threshold, NA_real_, L),
   #   col_pop = if_else(pop <= colonisation_threshold, NA_real_, pop)
   # )
 
 
 
 # ## calculate the mean homozygosity for the entire 95% data and mean_L
 # 
 # homo_df <- step_gen_occ |>
 #   group_by(replicate, year, patch, dispersal_frac, status) |>
 #   summarise(
 #     pop = first(pop),
 #     L = 1 - prod(1 - freq_a^2),
 #     .groups = "drop"
 #   )
 
 
 mean_homo <- homo_df |>
  group_by(year, patch, dispersal_frac, status) |>
  summarise(
    mean_L = mean(L),
    mean_pop = mean(pop),
    .groups = "drop"
  )|>
   mutate(
     col_pop = ifelse(mean_pop <= colonisation_threshold, NA, mean_pop),
     col_L = ifelse(mean_pop >= colonisation_threshold, mean_L, NA)
   )
 
 #re-order the status
 mean_homo$status <- factor(mean_homo$status,
                                levels = c("No effect", "Lethal effect", "Sterile effect"))
 
 
#homozygosity
ggplot(mean_homo, aes( year, patch, fill = col_L)) +
  geom_tile() +
  facet_grid(status~dispersal_frac) +
  scale_fill_viridis_c(option = "cividis", na.value = "white") +
  scale_y_continuous(breaks = 1:10) +
  theme_bw(base_size = 14) +
  labs(
    x = "Generations",
    y = "Patch",
    fill = "Homozygosity"
  )


## Cal CoV with the 95% data

cov_L <- homo_df |>
  group_by(year, patch, dispersal_frac, status) |>
  summarise(
    mean_L = mean(L, na.rm = TRUE),
    sd_L   = sd(L, na.rm = TRUE),
    cov_A  = if_else(mean_L == 0, NA_real_, sd_L / mean_L),
    .groups = "drop"
  )

cov_L$status <- factor(cov_L$status,
                           levels = c("No effect", "Lethal effect", "Sterile effect"))

# Plot coefficient of Variance 
ggplot(cov_L, aes(x = year, y = patch, fill = cov_A)) +
  geom_tile() +
  scale_fill_viridis_c(
    option = "cividis",
    na.value = "white"
    # name = "Homozygosity"
  ) +
  facet_grid(status~dispersal_frac) +
  labs(x = "Time (in generations)", y = "Patch", fill = "CoV") +
  theme_bw(base_size = 14)



cov_L |>
  group_by(year, status, dispersal_frac) |>
  summarise(mean_cv = mean(cov_A, na.rm = TRUE)) |>
  ggplot(aes(year, mean_cv)) +
  geom_line() +
  facet_grid(status ~ dispersal_frac)
