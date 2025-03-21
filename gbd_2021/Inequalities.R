# Visualization of Health Inequalities 
# EXAMPLE


# Load required packages ------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(car) # Heteroscedasticity test
library(MASS) # Robust regression
# install.packages("devtools")
# devtools::install_github("nicolash2/ggbrace")
library(ggbrace) # Draw braces
library(mgcv) # Provides Lorenz curve fitting (spline functions, etc.)
library(splines) # Fit spline functions
library(broom)

# Data Preparation --------------------------------------------------------------------
# Read and merge three datasets
setwd("C:\\Users\\35111\\Desktop\\Visualization of Health Inequalities")
# Disease burden data
burden <- read.csv("YourDatFromGBDExchange.csv", header = T)
burden <- burden %>% 
  filter(age != "Age-standardized")
# SDI data
sdi <- read.csv("SDIDatFromGBDExchange.csv", header = T, check.names = F)
sdi <- sdi %>% # Convert wide data to long format
  pivot_longer(cols = `1990`:`2019`, names_to = "year") %>%
  rename(sdi = value) %>%
  dplyr::select(location, year, sdi)
sdi$year <- as.integer(sdi$year)
# Merge SDI and disease burden, create data
data <- left_join(burden, sdi, by = c("location", "year"))
# Note: one country is Georgia
# Read population data
dir <- "/PopDir/"
files <- list.files(dir,
                    pattern = ".CSV",
                    full.names = T)
pop <- map_dfr(files, read.csv)
pop <- pop %>%
  dplyr::select(location_name, sex_name, age_group_name, year_id, metric_name, val)
colnames(pop) <- c("location", "sex", "age", "year", "metric", "val")
pop$sex[pop$sex == "male"] <- "Male"
pop$sex[pop$sex == "female"] <- "Female"
pop$sex[pop$sex == "both"] <- "Both"
pop <- pop %>%
  filter(age == "All Ages") %>%
  dplyr::select("location", "sex", "year", "val") %>%
  rename(pop = val)
# Merge population data, create mydata
mydata <- left_join(data, pop,
                    by = c("location", "sex", "year"))

# Visualization of the Slope Index ----------------------------------------------------------------

## 1. Prepare data for plotting -----------------------------------------------------------------
# Calculate total population
a <- mydata %>%
  filter(metric == "Number") %>%
  group_by(year) %>%
  summarise(sum = sum(pop))
pop1990 <- a$sum[1]
pop2019 <- a$sum[2]
# Calculate weighted rank
rank <- mydata %>%
  mutate(pop_global = ifelse(year == 1990, pop1990, pop2019)) %>%
  group_by(year, metric) %>%
  arrange(sdi) %>%
  mutate(cummu = cumsum(pop)) %>%       # Calculate cumulative population
  mutate(half = pop / 2) %>%              # Calculate half of each country's population
  mutate(midpoint = cummu - half) %>%     # Cumulative population minus half of the country's population gives the population midpoint
  mutate(weighted_order = midpoint / pop_global)  # The relative position of the country is the ratio of its population midpoint to the total population
# Convert year to a factor
rank$year <- factor(rank$year)
# Select data
temp1 <- rank %>%
  filter(metric == "Rate") %>%
  filter(year == 1990)
temp2 <- rank %>%
  filter(metric == "Rate") %>%
  filter(year == 2019)
# Model to calculate the slope index
fit1 <- lm(data = temp1, val ~ weighted_order)
fit2 <- lm(data = temp2, val ~ weighted_order)
# Check for heteroscedasticity (heteroscedasticity exists)
ncvTest(fit1)
ncvTest(fit2)
# Use robust regression: iteratively reweighted least squares
r.huber1 <- rlm(data = temp1, val ~ weighted_order)
r.huber2 <- rlm(data = temp2, val ~ weighted_order)
# Get coefficients and intercepts
coef(r.huber1)
coef(r.huber2)
# Calculate the 95% confidence intervals for the robust regression
confint.default(r.huber1)
confint.default(r.huber2)

# 2. Plotting ----------------------------------------------------------------------
color <- c("#6699FF", "#990000")
colnames(rank)
p1 <- rank %>%
  filter(metric == "Rate") %>%
  ggplot(aes(x = weighted_order, y = val, fill = year, group = year, color = year)) +
  geom_point(aes(color = year, size = pop / 1e6), alpha = 0.8, shape = 21) +
  scale_size_area("Population\n(million)", breaks = c(200, 400, 600, 800, 1000, 1200)) +
  geom_smooth(method = "rlm", size = 0.6, alpha = 0.1) +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) +
  # Add braces: the positions of the braces need to be predetermined based on the intercept and slope
  geom_brace(aes(x = c(1.003, 1.103), y = c(87548, 87548 - 70136)),
             inherit.data = F, size = 0.6,
             rotate = 90, color = "#6699FF") +
  geom_brace(aes(c(1, 1.1), c(40734, 40734 - 13805)),
             inherit.data = F, size = 0.6,
             rotate = 90, color = "#990000") +
  # Add horizontal dashed lines
  geom_segment(x = 0.02, xend = 0.99,
               y = 87548, yend = 87548,   # Position of the intercept
               color = "#6699FF", linetype = 2, size = 0.4, alpha = 0.4) +
  geom_segment(x = 0.02, xend = 0.99,
               y = 40734, yend = 40734,   # Position of the intercept
               color = "#990000", linetype = 2, size = 0.4, alpha = 0.4) +
  # Add labels for certain countries: for example, China and India
  geom_text(aes(label = ifelse(location == "China" | location == "India", as.character(location), ""),
                color = year),
            hjust = 0, vjust = 1.7,  # Avoid overlapping points and text
            size = 3) +
  # Add slope index labels
  annotate("text", label = "Slope Index of Inequality", x = 1.22, y = 40670 + 6751 / 2, size = 4, angle = 90) +
  annotate("text", label = "-70136", x = 1.15, y = 87548 - 70136 / 2, size = 3.5) +   # Position
  annotate("text", label = "-13805", x = 1.15, y = 40734 - 13805 / 2, size = 3.5) +
  scale_x_continuous(limits = c(0, 1.22), labels = c("0", "0.25", "0.50", "0.75", "1.00", "")) +
  xlab("Relative rank by SDI") +
  ylab("Crude DALY rate (per 100,000)") +
  theme_bw()

p1

# Visualization of the Concentration Index ----------------------------------------------------------------

# 1. Prepare data for plotting ------------------------------------------------------------------
a <- mydata %>%
  filter(metric == "Number") %>%
  group_by(year) %>%
  summarise(sum = sum(val))
daly1990 <- a$sum[1]
daly2019 <- a$sum[2]

ci <- rank %>%
  filter(metric == "Number") %>%
  mutate(total_daly = ifelse(year == 1990, daly1990, daly2019)) %>%
  group_by(year) %>%
  arrange(sdi) %>%
  mutate(cummu_daly = cumsum(val)) %>%         # Calculate cumulative DALY
  mutate(frac_daly = cummu_daly / total_daly) %>% # Calculate the fraction of cumulative DALY relative to the total
  mutate(frac_population = cummu / pop_global)    # Calculate the fraction of cumulative population relative to the total population

# 2. Plotting --------------------------------------------------------------------
p2 <- ci %>%
  ggplot(aes(x = frac_population, y = frac_daly, fill = year, color = year, group = year)) +
  # Add two line segments at X = 0, y = 0
  geom_segment(x = 0, xend = 1,
               y = 0, yend = 0,
               linetype = 1, size = 1, color = "gray") +
  geom_segment(x = 1, xend = 1,
               y = 0, yend = 1,
               linetype = 1, size = 1, color = "gray") +
  # Diagonal line
  geom_segment(x = 0, xend = 1,
               y = 0, yend = 1,
               color = "#CD853F", linetype = 1, size = 0.7, alpha = 1) +
  # Scatter points
  geom_point(aes(fill = year, size = pop / 1e6), alpha = 0.75, shape = 21) +
  scale_fill_manual(values = color) +
  scale_size_area("Population\n(million)", breaks = c(200, 400, 600, 800, 1000, 1200)) +
  # Fit Lorenz curve using cubic spline (setting knots and boundary conditions)
  geom_smooth(method = "gam", # You could also simply use geom_line to connect the points
              formula = y ~ ns(x,
                               knots = c(0.0000000001, 0.25, 0.5, 0.75, 0.9999999),  # Set knots to
                               # df = 3, not needed here since already set above
                               Boundary.knots = c(0, 1)),
              linetype = 1, size = 0.1, alpha = 0.6, se = T) +
  scale_color_manual(values = color) +
  # Add concentration index for the two years
  annotate("text", label = "Concentration Index", x = 0.75, y = 0.35, size = 5) +
  annotate("text", label = "1990: 0.1?", x = 0.75, y = 0.3, size = 4, color = "#6699FF") +
  annotate("text", label = "2019: 0.11?", x = 0.75, y = 0.25, size = 4, color = "#990000") +
  # Add labels for certain countries in 1990
  geom_text(aes(label = ifelse(location == "China" & year == "1990" | location == "India" & year == "1990",
                               as.character(location), "")),
            hjust = -0.6, vjust = 0.8,
            size = 3) +
  # Add labels for certain countries in 2019
  geom_text(aes(label = ifelse(location == "China" & year == "2019" | location == "India" & year == "2019",
                               as.character(location), "")),
            hjust = 1.8, vjust = -0.0,
            size = 3) +
  # Add labels for major population countries
  geom_text(aes(label = ifelse(location %in% a & year == "1990",
                               as.character(location), "")),
            hjust = -0.6, vjust = 0.8,
            size = 3) +
  geom_text(aes(label = ifelse(location %in% a & year == "2019",
                               as.character(location), "")),
            hjust = 1.8, vjust = -0.0,
            size = 3) +
  # x and y labels
  xlab("Cumulative fraction of population ranked by SDI") +
  ylab("Cumulative fraction of DALY") +
  theme_bw()
p2
