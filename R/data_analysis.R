# Data Analysis 




#load <- readRDS("C:\\Users\\22181916\\Documents\\Curtin-PhD\\R_and_IBM\\Generic_IBM_Proj\\IBM_generic\\output\\patch_load.rds")



#The codes below where used to analyse the data from the Latin hypercube sampling
# data aggregation

# Aggregated by simulation scenarios
load_df2 <- patch_stats %>%
  group_by(scenario, init_frequency, n_loci, fecundity, dispersal_prob) %>%
  summarise(
    mean_pop = mean(pop_size, na.rm = TRUE),
    sd_pop = sd(pop_size, na.rm = TRUE),
    mean_growth = mean(growth_rate, na.rm = TRUE),
    sd_growth = sd(growth_rate, na.rm = TRUE),
    sd_pop = sd(pop_size, na.rm = TRUE),
    mean_load = mean(load, na.rm = TRUE),
    sd_load = sd(load, na.rm = TRUE),
    mean_speed = mean(speed, na.rm = TRUE),
    sd_speed = sd(speed, na.rm = TRUE),
    .groups = "drop"
  )


# aggregated by patch and timestep
load_df <- patch_stats %>%
  group_by(patch, year, scenario, init_frequency, n_loci, dispersal_prob, arrival, time_k) %>%
  summarise(
    mean_pop = mean(pop_size, na.rm = TRUE),
    sd_pop = sd(pop_size, na.rm = TRUE),
    mean_growth = mean(growth_rate, na.rm = TRUE),
    sd_growth = sd(growth_rate, na.rm = TRUE),
    sd_pop = sd(pop_size, na.rm = TRUE),
    mean_load = mean(load, na.rm = TRUE),
    sd_load = sd(load, na.rm = TRUE),
    mean_speed = mean(speed, na.rm = TRUE),
    sd_speed = sd(speed, na.rm = TRUE),
    .groups = "drop"
  )


# Boosted regression tree (BRT) using gbm package 

# BRT Speed
brt_speed <- gbm(
  formula = mean_speed ~ init_frequency + fecundity + n_loci + dispersal_prob,
  data = load_df,
  distribution = "gaussian",  # use "bernoulli" for binary outcomes
  n.trees = 50000,
  interaction.depth = 3,
  shrinkage = 0.01,
  bag.fraction = 0.5,
  train.fraction = 0.8,
  cv.folds = 5,
  verbose = FALSE
)

# shows relative influence of each predictor
summary(brt_speed) 


# use emulator to get conditional effect plot
pred_df <- tibble(
  init_frequency = 0.1,
  fecundity = 5,
  n_loci = 100,
  dispersal_prob = seq(0.0001, 0.1, length.out = 100)
)

pred_df$pred_speed <- predict(brt_speed, newdata = pred_df)
plot(pred_speed ~ dispersal_prob, data = pred_df, type = "l")



# use LHS samples (without emulator) to get marginal effect plot   (From Nick G.)

load %>%
  # filter(
  #   init_frequency > 0.05 & init_frequency < 0.15,
  #   fecundity < 10 & fecundity > 2
  # ) %>%
  plot(speed ~ jitter(dispersal_prob),
       data = .,
       bg = patch,
       pch = 16,
       type = "p")


summary(load_df)



#multi-scenario
ggplot(all_patch_stats, aes(x = year, y = pop_size, color = factor(patch))) +
  geom_line(aes(group = patch), alpha = 0.6) +
  facet_wrap(~ dispersal_prob, ncol = 2, scales = "free_y") +
  labs(
    title = "population dynamics",
    x = "year",
    y = "pop",
    color = "Patch"
  ) +
  theme_minimal()







no_load_df <- no_load %>%
  group_by(patch, year, init_frequency, fecundity, n_loci, dispersal_prob) %>%
  summarise(
    mean_pop = mean(pop_size, na.rm = TRUE),
    sd_pop = sd(pop_size, na.rm = TRUE),
    mean_growth = mean(growth_rate, na.rm = TRUE),
    sd_growth = sd(growth_rate, na.rm = TRUE),
    mean_load = mean(load, na.rm = TRUE),
    sd_load = sd(load, na.rm = TRUE),
    mean_speed = mean(speed, na.rm = TRUE),
    sd_speed = sd(speed, na.rm = TRUE),
    .groups = "drop"
  )



load_growth <- load_df |> select(mean_growth)
no_load_growth <- no_load_df |> select(mean_growth)
load_growth <- unlist(load_growth$mean_growth)
no_load_growth <- unlist(no_load_growth$mean_growth)
growth_diff <- no_load_growth - load_growth


load_speed <- load_df |> select(mean_speed)
no_load_speed <- no_load_df |> select(mean_speed)
load_speed <- unlist(load_speed$mean_speed)
no_load_speed <- unlist(no_load_speed$mean_speed)
speed_diff <- no_load_speed - load_speed




clean_data <- load_df |> select(scenario, init_frequency, fecundity, n_loci, dispersal_prob, mean_load)

clean_data$speed <- speed_diff
clean_data$growth <- growth_diff





# Boosted regression tree (BRT) using gbm package 

# BRT Speed
brt_speed <- gbm(
  formula = mean_speed ~ init_frequency + fecundity + n_loci + dispersal_prob,
  data = load_df,
  distribution = "gaussian",  # use "bernoulli" for binary outcomes
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.01,
  bag.fraction = 0.5,
  train.fraction = 0.8,
  cv.folds = 5,
  verbose = FALSE
)

# shows relative influence of each predictor
summary(brt_speed) 


# BRT Speed
brt_speed <- gbm(
  formula = speed ~ init_frequency + fecundity + n_loci + dispersal_prob,
  data = clean_data,
  distribution = "gaussian",  # use "bernoulli" for binary outcomes
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.01,
  bag.fraction = 0.5,
  train.fraction = 0.8,
  cv.folds = 5,
  verbose = FALSE
)

# shows relative influence of each predictor
summary(brt_speed) 


plot(brt_speed, i.var = "init_frequency")
plot(brt_speed, i.var = "fecundity")


#Sobol analysis

# get data for emulator
emul_data <- patch_stats %>%
  ungroup() %>%
  select(speed, init_frequency, fecundity, n_loci, dispersal_prob, prob_survival) 

# emulate speed
ref_emul <- randomForest(
  speed ~ init_frequency + fecundity + n_loci + dispersal_prob + prob_survival,
  data = emul_data,
  n.tree = 1000,
  importance = TRUE)


# print(ref_emul)
# importance(ref_emul)

# -----------------------------
# 2. Generate parameters with LHS Parameter Sets ---
# -----------------------------

N <- 1000

param_ranges <- data.frame(
  init_frequency = c(0.1, 0.5),
  fecundity       = c(2, 30),
  n_loci          = c(5, 30),
  dispersal_prob  = c(0.0001, 0.1),
  prob_survival = c(0.6, 1)
)

param_names <- names(param_ranges)

lhs_sampleA <- randomLHS(N, ncol(param_ranges))
lhs_sampleB <- randomLHS(N, ncol(param_ranges))


A <- data.frame(
  init_frequency = qunif(lhs_sampleA[,1], param_ranges[1,1], param_ranges[2,1]),
  fecundity = qinteger(lhs_sampleA[,2], param_ranges[1,2], param_ranges[2,2]),
  n_loci = qinteger(lhs_sampleA[,3], param_ranges[1,3], param_ranges[2,3]),
  dispersal_prob = qunif(lhs_sampleA[,4], param_ranges[1,4], param_ranges[2,4]),
  prob_survival = qunif(lhs_sampleA[,5], param_ranges[1,5], param_ranges[2,5])
)


B <- data.frame(
  init_frequency = qunif(lhs_sampleB[,1], param_ranges[1,1], param_ranges[2,1]),
  fecundity = qinteger(lhs_sampleB[,2], param_ranges[1,2], param_ranges[2,2]),
  n_loci = qinteger(lhs_sampleB[,3], param_ranges[1,3], param_ranges[2,3]),
  dispersal_prob = qunif(lhs_sampleB[,4], param_ranges[1,4], param_ranges[2,4]),
  prob_survival = qunif(lhs_sampleB[,5], param_ranges[1,5], param_ranges[2,5])
)



# 3. define model function
emulator <- function(X) {
  newdata <- as.data.frame(X)
  colnames(newdata) <- param_names
  predict(ref_emul, newdata)
}


# 4. run sobol on the emulator
sob <- soboljansen(
  model = emulator,
  X1 = as.matrix(A),
  X2 = as.matrix(B),
  nboot = 100
)


print(sob)


# plot Sobol Analysis

si <- sob$S   # first-order indices
ti <- sob$T   # total-order indices



sob_df <- data.frame(
  parameter = rownames(sob$S),
  S  = sob$S$original,
  S_low  = sob$S$`min. c.i.`,
  S_high = sob$S$`max. c.i.`,
  T  = sob$T$original,
  T_low  = sob$T$`min. c.i.`,
  T_high = sob$T$`max. c.i.`
)



library(ggplot2)

# first-order sensitivity analysis 

ggplot(sob_df, aes(x = reorder(parameter, S), y = S)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = S_low, ymax = S_high), width = 0.2) +
  coord_flip() +
  labs(
    x = "Parameter",
    y = "First-order Sobol Index",
    title = "Sensitivity Analysis (First-order)"
  ) +
  theme_minimal(base_size = 14)


# first-order and total sensitivity analysis 

sob_long <- sob_df %>%
  select(parameter, S, T) %>%
  pivot_longer(cols = c(S, T), names_to = "index", values_to = "value")

ggplot(sob_long, aes(x = reorder(parameter, value), y = value, fill = index)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    x = "Parameter",
    y = "Sobol Index",
    title = "First-order vs Total Sobol Indices"
  ) +
  scale_fill_manual(values = c("S" = "skyblue", "T" = "orange")) +
  theme_minimal(base_size = 14)






# Example: marginal relationship init_frequency vs mean_pop
ggplot(clean_data, aes(x = growth, y = mean_load)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_minimal() +
  labs(x = "difference in growth rate", y = "genetic load",
       title = "Relationship growth rate and genetic load")



ggplot(clean_data, aes(x = speed, y = mean_load)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_minimal() +
  labs(x = "Invasion speed", y = "genetic load",
       title = "Relationship between fecundity and speed")


# Example: marginal relationship init_frequency vs mean_pop
ggplot(patch_df, aes(x = mean_speed, y = dispersal_prob)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_minimal() +
  labs(x = "Invasion Speed", y = "dispersal",
       title = "Relationship between dispersal and speed")
