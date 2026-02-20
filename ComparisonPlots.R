#####################################################################################
#####################################################################################
##**                        Lamprey Comparison Plots                             **##
##**                             February 2025                                   **##
##**                            Mikayla Stinson                                  **##
#####################################################################################
#####################################################################################

#This code is solely for interest in looking at the comparison between years of
#         different programs. This is not to go with the data request.


#####################################################################################
#####################################################################################
##**                           Packages to Load                                  **##
#####################################################################################
#####################################################################################

library(tidyverse)
library(readr)
library(ggplot2)
library(purrr)
library(stringr)
install.packages("janitor")
library(janitor)
install.packages("readxl")
library(readxl)
install.packages("viridis")
library(viridis)
set.seed(42)

#####################################################################################
#####################################################################################

#95% confidence interval (Bootstrap  per-year resampling of fish)
bootstrap_rate_ci <- function(x, n_boot = 5000, alpha = 0.05) {
  # x = vector of per-fish total_wounds for a given year
  n <- length(x)
  means <- replicate(n_boot, mean(sample(x, size = n, replace = TRUE)))
  rate_per100 <- 100 * means
  quantile(rate_per100, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
}

#95% confidence interval (Poisson exact on total counts)
poisson_count_ci <- function(k, alpha = 0.05) {
  lo <- if (k == 0) 0 else 0.5 * qchisq(alpha / 2, 2 * k)
  hi <- 0.5 * qchisq(1 - alpha / 2, 2 * (k + 1))
  c(lo = lo, hi = hi)
}

#####################################################################################
#####################################################################################
##**                                  GL1                                        **##
#####################################################################################
#####################################################################################

###############                GL1 - Read in Files              #####################

#2025
GL1 = read.csv(file = "GL1/LOA_IA25_GL1_bio.csv")

#2024
GL1_2024 = read.csv(
  file = "G:/Fisheries Management/Assessment/Lamprey/Lamprey Data Request/2024/GL1/LOA_IA24_GL1_bio.csv"
)

#2023
GL1_2023 = read_excel(
  "G:/Fisheries Management/Assessment/Lamprey/Lamprey Data Request/2023/LOA_IA23_TW1_GL1_lamprey.xlsx",
  sheet = 1
)

#2022
GL1_2022 = read_excel(
  "G:/Fisheries Management/Assessment/Lamprey/Lamprey Data Request/2022/GL1lamprey_biodata_export.xlsx",
  sheet = 1
)


##########    Reduce the tables down to only the lamprey wound columns      #########

GL1_2025_lamprey = GL1[18:25]
GL1_2024_lamprey = GL1_2024[20:25]
GL1_2023_lamprey = GL1_2023[4:11]
GL1_2022_lamprey = GL1_2022[18:25]

##########################     Add year column     ##################################


GL1_2025 = GL1_2025_lamprey %>%
  mutate(year = 2025)

GL1_2024 = GL1_2024_lamprey %>%
  mutate(year = 2024)

GL1_2023 = GL1_2023_lamprey %>%
  mutate(year = 2023)

GL1_2022 = GL1_2022_lamprey %>%
  mutate(year = 2022)


###############              Remainder of code           ############################


#combine all the rows from all lamprey tables
GL1_allyears = bind_rows(GL1_2025, GL1_2024, GL1_2023, GL1_2022)

#turn NA values into 0
GL1_allyears[is.na(GL1_allyears)] = 0

#add a column of total wounds = sum across all columns
Clean_GL1 = GL1_allyears %>%
  mutate(total_wounds = A1 + A2 + A3 + A4 + B1 + B2 + B3 + B4)


#yearly summary including rate per 100 fish
yearly = Clean_GL1 %>%
  group_by(year) %>%
  summarise(
    n_fish = n(),
    Wounds_total = sum(total_wounds, na.rm = TRUE),
    total_per100 = 100 * Wounds_total / n_fish,
    .groups = "drop"
  ) %>%
  arrange(year)


#long format for plotting all the years
by_year_long = yearly %>%
  select(year, Wounds_total) %>%
  pivot_longer(
    cols = Wounds_total,
    names_to = "category",
    values_to = "count"
  ) %>%
  mutate(category = toupper(category))

#Bootstrap
ci_boot <- Clean_GL1 |>
  group_by(year) |>
  reframe({
    qs <- bootstrap_rate_ci(total_wounds)
    tibble(total_per100_lo = qs[[1]], total_per100_hi = qs[[2]])
  })

yearly_boot <- yearly |>
  left_join(ci_boot, by = "year")

#Poisson
by_year_ci <- by_year_long |>
  rowwise() |>
  mutate(
    ci = list(poisson_count_ci(count)),
    total_lo = ci["lo"],
    total_hi = ci["hi"]
  ) |>
  select(-ci) |>
  ungroup()


#Ploting rate per 100 fish
vcol <- viridis(1, option = "magma", begin = 0.8) # default
# options: "A", "B", "C", "D", "E", "F", "G", "H", "magma", "plasma", "inferno", "cividis"

ggplot(yearly_boot, aes(x = year, y = total_per100)) +
  geom_ribbon(
    aes(ymin = total_per100_lo, ymax = total_per100_hi),
    fill = vcol,
    alpha = 0.20
  ) +
  geom_line(color = vcol, linewidth = 1) +
  geom_point(color = vcol, size = 2) +
  scale_x_continuous(
    breaks = seq(min(yearly_boot$year), max(yearly_boot$year), 1)
  ) +
  labs(
    title = "Total Lamprey Wounds per 100 Fish by Year",
    subtitle = "Community Index Gillnetting Program",
    x = "Year",
    y = "Wounds per 100 fish"
  ) +
  theme_minimal()

#bar plot - rate per 100 fish
ggplot(yearly_boot, aes(x = year, y = total_per100)) +
  geom_col(aes(fill = "bar"), position = "dodge", width = 0.7) +
  scale_fill_viridis_d(option = "turbo", direction = 1, begin = 1) +
  geom_errorbar(
    aes(ymin = total_per100_lo, ymax = total_per100_hi),
    width = 0.25,
    linewidth = 0.7
  ) +
  scale_x_continuous(
    breaks = seq(min(yearly_boot$year), max(yearly_boot$year), 1)
  ) +
  labs(
    title = "Total Lamprey Wounds per 100 Fish by Year",
    subtitle = "Community Index Gillnetting Program",
    x = "Year",
    y = "Wounds per 100 fish"
  ) +
  theme_classic(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")


#line plot of all counts per year
ggplot(by_year_long, aes(x = year, y = count)) +
  geom_line(color = vcol, linewidth = 1) +
  geom_point(color = vcol, size = 2) +
  scale_x_continuous(
    breaks = seq(
      min(by_year_long$year, na.rm = TRUE),
      max(by_year_long$year, na.rm = TRUE),
      by = 1
    )
  ) +
  labs(
    title = "Total Lamprey Wounds (raw counts) by Year",
    subtitle = "Community Index Gillnetting Program",
    x = "Year",
    y = "Total wound count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )


# Stacked bar showing all counts per year
ggplot(by_year_ci, aes(x = factor(year), y = count, fill = category)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_viridis_d(option = "inferno", direction = 1, begin = 0.5) + # options: "viridis","magma","plasma","inferno","cividis","turbo"
  coord_cartesian(ylim = c(0, NA)) + # Avoid clipping CI tops
  geom_errorbar(
    aes(ymin = pmax(0, total_lo), ymax = total_hi),
    width = 0.18,
    linewidth = 0.6
  ) +
  labs(
    title = "Lamprey Wound Counts by Year",
    subtitle = "Community Index Gillnetting Program",
    x = "Year",
    y = "Total wounds count",
    fill = "Category"
  ) +
  theme_classic(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")


#####################################################################################
#####################################################################################
##**                                   EC1                                       **##
#####################################################################################
#####################################################################################

################              EC1 - Read in Files            ########################

EC1_2025 = read.csv(
  file = "G:/Fisheries Management/Assessment/Lamprey/Lamprey Data Request/2025/EC1_EC2/LOA_SF25_EC1_bio.csv"
)
EC1_2024 = read.csv(
  file = "G:/Fisheries Management/Assessment/Lamprey/Lamprey Data Request/2024/EC1/LOA_SF24_EC1_bio.csv"
)
EC1_2023 = read_excel(
  "G:/Fisheries Management/Assessment/Lamprey/Lamprey Data Request/2023/LOA_SF23_EC1.xlsx",
  sheet = 2
)
EC1_2022 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2022"
)
EC1_2021 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2021"
)
EC1_2020 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2020"
)
EC1_2019 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2019"
)
EC1_2018 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2018"
)
EC1_2017 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2017"
)
EC1_2016 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2016"
)
EC1_2015 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2015"
)
EC1_2014 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2014"
)
EC1_2013 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2013"
)
EC1_2012 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2012"
)
EC1_2011 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2011"
)
EC1_2010 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2010"
)
EC1_2009 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2009"
)
EC1_2008 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2008"
)
EC1_2007 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2007"
)
EC1_2006 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2006"
)
EC1_2005 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2005"
)
EC1_2004 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2004"
)
EC1_2003 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2003"
)
EC1_2002 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2002"
)
EC1_2001 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2001"
)
EC1_2000 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "2000"
)
EC1_1999 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "1999"
)
EC1_1998 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "1998"
)
EC1_1997 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "1997"
)
EC1_1996 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "1996"
)
EC1_1995 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "1995"
)
EC1_1994 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "1994"
)
EC1_1993 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "1993"
)
EC1_1992 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA92-22_EC1-edited.xlsx",
  sheet = "1992"
)


############     Reduce the tables down to only the lamprey wound columns     #######

EC1_2025_lamprey = EC1_2025[14:21]
EC1_2024_lamprey = EC1_2024[14:21]
EC1_2023_lamprey = EC1_2023[19:26]
EC1_2022_lamprey = EC1_2022[19:26]
EC1_2021_lamprey = EC1_2021[19:26]
EC1_2020_lamprey = EC1_2020[19:26]
EC1_2019_lamprey = EC1_2019[19:26]
EC1_2018_lamprey = EC1_2018[19:26]
EC1_2017_lamprey = EC1_2017[19:26]
EC1_2016_lamprey = EC1_2016[19:26]
EC1_2015_lamprey = EC1_2015[19:26]
EC1_2014_lamprey = EC1_2014[19:26]
EC1_2013_lamprey = EC1_2013[19:26]
EC1_2012_lamprey = EC1_2012[19:26]
EC1_2011_lamprey = EC1_2011[19:26]
EC1_2010_lamprey = EC1_2010[19:26]
EC1_2009_lamprey = EC1_2009[19:26]
EC1_2008_lamprey = EC1_2008[19:26]
EC1_2007_lamprey = EC1_2007[19:26]
EC1_2006_lamprey = EC1_2006[19:26]
EC1_2005_lamprey = EC1_2005[19:26]
EC1_2004_lamprey = EC1_2004[19:26]
EC1_2003_lamprey = EC1_2003[19:26]
EC1_2002_lamprey = EC1_2002[19:26]
EC1_2001_lamprey = EC1_2001[19:26]
EC1_2000_lamprey = EC1_2000[19:26]
EC1_1999_lamprey = EC1_1999[19:26]
EC1_1998_lamprey = EC1_1998[19:26]
EC1_1997_lamprey = EC1_1997[19:26]
EC1_1996_lamprey = EC1_1996[19:26]
EC1_1995_lamprey = EC1_1995[19:26]
EC1_1994_lamprey = EC1_1994[19:26]
EC1_1993_lamprey = EC1_1993[19:26]
EC1_1992_lamprey = EC1_1992[19:26]


##################             Add year column                   ####################

EC1_2025 = EC1_2025_lamprey %>%
  mutate(year = 2025)

EC1_2024 = EC1_2024_lamprey %>%
  mutate(year = 2024)

EC1_2023 = EC1_2023_lamprey %>%
  mutate(year = 2023)

EC1_2022 = EC1_2022_lamprey %>%
  mutate(year = 2022)

EC1_2021 = EC1_2021_lamprey %>%
  mutate(year = 2021)

EC1_2020 = EC1_2020_lamprey %>%
  mutate(year = 2020)

EC1_2019 = EC1_2019_lamprey %>%
  mutate(year = 2019)

EC1_2018 = EC1_2018_lamprey %>%
  mutate(year = 2018)

EC1_2017 = EC1_2017_lamprey %>%
  mutate(year = 2017)

EC1_2016 = EC1_2016_lamprey %>%
  mutate(year = 2016)

EC1_2015 = EC1_2015_lamprey %>%
  mutate(year = 2015)

EC1_2014 = EC1_2014_lamprey %>%
  mutate(year = 2014)

EC1_2013 = EC1_2013_lamprey %>%
  mutate(year = 2013)

EC1_2012 = EC1_2012_lamprey %>%
  mutate(year = 2012)

EC1_2011 = EC1_2011_lamprey %>%
  mutate(year = 2011)

EC1_2010 = EC1_2010_lamprey %>%
  mutate(year = 2010)

EC1_2009 = EC1_2009_lamprey %>%
  mutate(year = 2009)

EC1_2008 = EC1_2008_lamprey %>%
  mutate(year = 2008)

EC1_2007 = EC1_2007_lamprey %>%
  mutate(year = 2007)

EC1_2006 = EC1_2006_lamprey %>%
  mutate(year = 2006)

EC1_2005 = EC1_2005_lamprey %>%
  mutate(year = 2005)

EC1_2004 = EC1_2004_lamprey %>%
  mutate(year = 2004)

EC1_2003 = EC1_2003_lamprey %>%
  mutate(year = 2003)

EC1_2002 = EC1_2002_lamprey %>%
  mutate(year = 2002)

EC1_2001 = EC1_2001_lamprey %>%
  mutate(year = 2001)

EC1_2000 = EC1_2000_lamprey %>%
  mutate(year = 2000)

EC1_1999 = EC1_1999_lamprey %>%
  mutate(year = 1999)

EC1_1998 = EC1_1998_lamprey %>%
  mutate(year = 1998)

EC1_1997 = EC1_1997_lamprey %>%
  mutate(year = 1997)

EC1_1996 = EC1_1996_lamprey %>%
  mutate(year = 1996)

EC1_1995 = EC1_1995_lamprey %>%
  mutate(year = 1995)

EC1_1994 = EC1_1994_lamprey %>%
  mutate(year = 1994)

EC1_1993 = EC1_1993_lamprey %>%
  mutate(year = 1993)

EC1_1992 = EC1_1992_lamprey %>%
  mutate(year = 1992)

###################           Remainder of code                ######################

#combine all the rows from all lamprey tables
EC1_allyears = bind_rows(
  EC1_2025,
  EC1_2024,
  EC1_2023,
  EC1_2022,
  EC1_2021,
  EC1_2020,
  EC1_2019,
  EC1_2018,
  EC1_2017,
  EC1_2016,
  EC1_2015,
  EC1_2014,
  EC1_2013,
  EC1_2012,
  EC1_2011,
  EC1_2010,
  EC1_2009,
  EC1_2008,
  EC1_2007,
  EC1_2006,
  EC1_2005,
  EC1_2004,
  EC1_2003,
  EC1_2002,
  EC1_2001,
  EC1_2000,
  EC1_1999,
  EC1_1998,
  EC1_1997,
  EC1_1996,
  EC1_1995,
  EC1_1994,
  EC1_1993,
  EC1_1992
)


#turn NA values into 0
EC1_allyears[is.na(EC1_allyears)] = 0

#add a column of total wounds = sum across all columns
Clean_EC1 = EC1_allyears %>%
  mutate(total_wounds = A1 + A2 + A3 + A4 + B1 + B2 + B3 + B4)


#yearly summary including rate per 100 fish
yearly2 = Clean_EC1 %>%
  group_by(year) %>%
  summarise(
    n_fish = n(),
    Wounds_total = sum(total_wounds, na.rm = TRUE),
    total_per100 = 100 * Wounds_total / n_fish,
    .groups = "drop"
  ) %>%
  arrange(year)


#long format for plotting all the years
by_year_long2 = yearly2 %>%
  select(year, Wounds_total) %>%
  pivot_longer(
    cols = Wounds_total,
    names_to = "category",
    values_to = "count"
  ) %>%
  mutate(category = toupper(category))

#Bootstrap
ci_boot2 <- Clean_EC1 |>
  group_by(year) |>
  reframe({
    qs <- bootstrap_rate_ci(total_wounds)
    tibble(total_per100_lo = qs[[1]], total_per100_hi = qs[[2]])
  })

yearly_boot2 <- yearly2 |>
  left_join(ci_boot2, by = "year")


#Poisson
by_year_ci2 <- by_year_long2 |>
  rowwise() |>
  mutate(
    ci = list(poisson_count_ci(count)),
    total_lo = ci["lo"],
    total_hi = ci["hi"]
  ) |>
  select(-ci) |>
  ungroup()


#Ploting rate per 100 fish
vcol <- viridis(1, option = "magma", begin = 0.8) # default

ggplot(yearly_boot2, aes(x = year, y = total_per100)) +
  geom_ribbon(
    aes(ymin = total_per100_lo, ymax = total_per100_hi),
    fill = vcol,
    alpha = 0.20
  ) +
  geom_line(color = vcol, linewidth = 1) +
  geom_point(color = vcol, size = 2) +
  scale_x_continuous(
    breaks = seq(1992, max(yearly_boot2$year, na.rm = TRUE), by = 5)
  ) +
  labs(
    title = "Total Lamprey Wounds per 100 Fish by Year",
    subtitle = "Credit River Chinook salmon Egg Collection",
    x = "Year",
    y = "Wounds per 100 fish"
  ) +
  theme_bw()

#bar plot - rate per 100 fish
ggplot(yearly_boot2, aes(x = year, y = total_per100)) +
  geom_col(aes(fill = "bar")) +
  scale_fill_viridis_d(option = "turbo", direction = 1, begin = 1) +
  geom_errorbar(
    aes(ymin = total_per100_lo, ymax = total_per100_hi),
    width = 0.25,
    linewidth = 0.7
  ) +
  scale_x_continuous(
    breaks = seq(1992, max(yearly_boot2$year, na.rm = TRUE), by = 5)
  ) +
  labs(
    title = "Total Lamprey Wounds per 100 Fish by Year",
    subtitle = "Credit River Chinook salmon Egg Collection",
    x = "Year",
    y = "Wounds per 100 fish"
  ) +
  theme_classic(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")


#line plot of all counts per year
ggplot(by_year_long2, aes(x = year, y = count)) +
  geom_line(color = vcol, linewidth = 1) +
  geom_point(color = vcol, size = 2) +
  scale_x_continuous(
    breaks = seq(1992, max(yearly_boot2$year, na.rm = TRUE), by = 5)
  ) +
  labs(
    title = "Total Lamprey Wounds (raw counts) by Year",
    subtitle = "Credit River Chinook salmon Egg Collection",
    x = "Year",
    y = "Total wound count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )


# Stacked bar showing all counts per year
ggplot(by_year_ci2, aes(x = year, y = count, fill = category)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_x_continuous(
    breaks = seq(1992, max(yearly_boot2$year, na.rm = TRUE), by = 5)
  ) +
  scale_fill_viridis_d(option = "inferno", direction = 1, begin = 0.5) + # options: "viridis","magma","plasma","inferno","cividis","turbo"
  coord_cartesian(ylim = c(0, NA)) + # Avoid clipping CI tops
  geom_errorbar(
    aes(ymin = pmax(0, total_lo), ymax = total_hi),
    width = 0.18,
    linewidth = 0.6
  ) +
  labs(
    title = "Lamprey Wound Counts by Year",
    subtitle = "Credit River Chinook salmon Egg Collection",
    x = "Year",
    y = "Total wounds count",
    fill = "Category"
  ) +
  theme_classic(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")



#####################################################################################
#####################################################################################
##**                                    WGS                                      **##
#####################################################################################
#####################################################################################

################              WGS - Read in Files              ######################

WGS_2025 = read.csv(
  file = "G:/Fisheries Management/Assessment/Lamprey/Lamprey Data Request/2025/WGS/LOA_IA25_WGS_bio.csv"
)
WGS_2024 = read.csv(
  file = "G:/Fisheries Management/Assessment/Lamprey/Lamprey Data Request/2024/WGS/LOA_IA24_WGS_bio.csv"
)

WGS_2022 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2022"
)
WGS_2021 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2021"
)

WGS_2019 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2019"
)
WGS_2018 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2018"
)
WGS_2017 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2017"
)
WGS_2016 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2016"
)
WGS_2015 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2015"
)

WGS_2013 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2013"
)

WGS_2011 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2011"
)
WGS_2010 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2010"
)
WGS_2009 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2009"
)
WGS_2008 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2008"
)
WGS_2007 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2007"
)
WGS_2006 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2006"
)
WGS_2005 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2005"
)
WGS_2004 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2004"
)
WGS_2003 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2003"
)

WGS_2001 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2001"
)
WGS_2000 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "2000"
)
WGS_1999 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "1999"
)
WGS_1998 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "1998"
)
WGS_1997 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "1997"
)
WGS_1996 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "1996"
)
WGS_1995 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "1995"
)
WGS_1994 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "1994"
)
WGS_1993 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "1993"
)
WGS_1992 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "1992"
)
WGS_1991 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "1991"
)
WGS_1990 = read_excel(
  "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lamprey/LOA_IA90-22_WGS.xlsx",
  sheet = "1990"
)


############     Reduce the tables down to only the lamprey wound columns    ########
WGS_2025_lamprey = WGS_2025[19:26]
WGS_2024_lamprey = WGS_2024[19:26]

WGS_2022_lamprey = WGS_2022[19:26]
WGS_2021_lamprey = WGS_2021[19:26]

WGS_2019_lamprey = WGS_2019[19:26]
WGS_2018_lamprey = WGS_2018[19:26]
WGS_2017_lamprey = WGS_2017[19:26]
WGS_2016_lamprey = WGS_2016[19:26]
WGS_2015_lamprey = WGS_2015[19:26]

WGS_2013_lamprey = WGS_2013[19:26]

WGS_2011_lamprey = WGS_2011[19:26]
WGS_2010_lamprey = WGS_2010[19:26]
WGS_2009_lamprey = WGS_2009[19:26]
WGS_2008_lamprey = WGS_2008[19:26]
WGS_2007_lamprey = WGS_2007[19:26]
WGS_2006_lamprey = WGS_2006[19:26]
WGS_2005_lamprey = WGS_2005[19:26]
WGS_2004_lamprey = WGS_2004[19:26]
WGS_2003_lamprey = WGS_2003[19:26]

WGS_2001_lamprey = WGS_2001[19:26]
WGS_2000_lamprey = WGS_2000[19:26]
WGS_1999_lamprey = WGS_1999[19:26]
WGS_1998_lamprey = WGS_1998[19:26]
WGS_1997_lamprey = WGS_1997[19:26]
WGS_1996_lamprey = WGS_1996[19:26]
WGS_1995_lamprey = WGS_1995[19:26]
WGS_1994_lamprey = WGS_1994[19:26]
WGS_1993_lamprey = WGS_1993[19:26]
WGS_1992_lamprey = WGS_1992[19:26]
WGS_1991_lamprey = WGS_1991[19:26]
WGS_1990_lamprey = WGS_1990[19:26]


#################              Add year column               ########################

WGS_2025 = WGS_2025_lamprey %>%
  mutate(year = 2025)

WGS_2024 = WGS_2024_lamprey %>%
  mutate(year = 2024)


WGS_2022 = WGS_2022_lamprey %>%
  mutate(year = 2022)

WGS_2021 = WGS_2021_lamprey %>%
  mutate(year = 2021)


WGS_2019 = WGS_2019_lamprey %>%
  mutate(year = 2019)

WGS_2018 = WGS_2018_lamprey %>%
  mutate(year = 2018)

WGS_2017 = WGS_2017_lamprey %>%
  mutate(year = 2017)

WGS_2016 = WGS_2016_lamprey %>%
  mutate(year = 2016)

WGS_2015 = WGS_2015_lamprey %>%
  mutate(year = 2015)


WGS_2013 = WGS_2013_lamprey %>%
  mutate(year = 2013)


WGS_2011 = WGS_2011_lamprey %>%
  mutate(year = 2011)

WGS_2010 = WGS_2010_lamprey %>%
  mutate(year = 2010)

WGS_2009 = WGS_2009_lamprey %>%
  mutate(year = 2009)

WGS_2008 = WGS_2008_lamprey %>%
  mutate(year = 2008)

WGS_2007 = WGS_2007_lamprey %>%
  mutate(year = 2007)

WGS_2006 = WGS_2006_lamprey %>%
  mutate(year = 2006)

WGS_2005 = WGS_2005_lamprey %>%
  mutate(year = 2005)

WGS_2004 = WGS_2004_lamprey %>%
  mutate(year = 2004)

WGS_2003 = WGS_2003_lamprey %>%
  mutate(year = 2003)


WGS_2001 = WGS_2001_lamprey %>%
  mutate(year = 2001)

WGS_2000 = WGS_2000_lamprey %>%
  mutate(year = 2000)

WGS_1999 = WGS_1999_lamprey %>%
  mutate(year = 1999)

WGS_1998 = WGS_1998_lamprey %>%
  mutate(year = 1998)

WGS_1997 = WGS_1997_lamprey %>%
  mutate(year = 1997)

WGS_1996 = WGS_1996_lamprey %>%
  mutate(year = 1996)

WGS_1995 = WGS_1995_lamprey %>%
  mutate(year = 1995)

WGS_1994 = WGS_1994_lamprey %>%
  mutate(year = 1994)

WGS_1993 = WGS_1993_lamprey %>%
  mutate(year = 1993)

WGS_1992 = WGS_1992_lamprey %>%
  mutate(year = 1992)

WGS_1991 = WGS_1991_lamprey %>%
  mutate(year = 1991)

WGS_1990 = WGS_1990_lamprey %>%
  mutate(year = 1990)


##################            Remainder of code              ########################

#combine all the rows from all lamprey tables
WGS_allyears = bind_rows(
  WGS_2025,
  WGS_2024,
  WGS_2022,
  WGS_2021,
  WGS_2019,
  WGS_2018,
  WGS_2017,
  WGS_2016,
  WGS_2015,
  WGS_2013,
  WGS_2011,
  WGS_2010,
  WGS_2009,
  WGS_2008,
  WGS_2007,
  WGS_2006,
  WGS_2005,
  WGS_2004,
  WGS_2003,
  WGS_2001,
  WGS_2000,
  WGS_1999,
  WGS_1998,
  WGS_1997,
  WGS_1996,
  WGS_1995,
  WGS_1994,
  WGS_1993,
  WGS_1992,
  WGS_1991,
  WGS_1990
)


#turn NA values into 0
WGS_allyears[is.na(WGS_allyears)] = 0

#add a column of total wounds = sum across all columns
Clean_WGS = WGS_allyears %>%
  mutate(total_wounds = A1 + A2 + A3 + A4 + B1 + B2 + B3 + B4)


#yearly summary including rate per 100 fish
yearly3 = Clean_WGS %>%
  group_by(year) %>%
  summarise(
    n_fish = n(),
    Wounds_total = sum(total_wounds, na.rm = TRUE),
    total_per100 = 100 * Wounds_total / n_fish,
    .groups = "drop"
  ) %>%
  arrange(year)


#long format for plotting all the years
by_year_long3 = yearly3 %>%
  select(year, Wounds_total) %>%
  pivot_longer(
    cols = Wounds_total,
    names_to = "category",
    values_to = "count"
  ) %>%
  mutate(category = toupper(category))


#Bootstrap
ci_boot3 <- Clean_WGS |>
  group_by(year) |>
  reframe({
    qs <- bootstrap_rate_ci(total_wounds)
    tibble(total_per100_lo = qs[[1]], total_per100_hi = qs[[2]])
  })

yearly_boot3 <- yearly3 |>
  left_join(ci_boot3, by = "year")

#Poisson
by_year_ci3 <- by_year_long3 |>
  rowwise() |>
  mutate(
    ci = list(poisson_count_ci(count)),
    total_lo = ci["lo"],
    total_hi = ci["hi"]
  ) |>
  select(-ci) |>
  ungroup()

#Ploting rate per 100 fish
vcol <- viridis(1, option = "magma", begin = 0.8) # default

ggplot(yearly_boot3, aes(x = year, y = total_per100)) +
  geom_ribbon(
    aes(ymin = total_per100_lo, ymax = total_per100_hi),
    fill = vcol,
    alpha = 0.20
  ) +
  geom_line(color = vcol, linewidth = 1) +
  geom_point(color = vcol, size = 2) +
  scale_x_continuous(
    breaks = seq(1990, max(yearly_boot3$year, na.rm = TRUE), by = 5)
  ) +
  labs(
    title = "Total Lamprey Wounds per 100 Fish by Year",
    subtitle = "Rainbow Trout Spring Assessment",
    x = "Year",
    y = "Wounds per 100 fish"
  ) +
  theme_bw()

#bar plot - rate per 100 fish
ggplot(yearly_boot3, aes(x = year, y = total_per100)) +
  geom_col(aes(fill = "bar")) +
  scale_fill_viridis_d(option = "turbo", direction = 1, begin = 1) +
  geom_errorbar(
    aes(ymin = total_per100_lo, ymax = total_per100_hi),
    width = 0.25,
    linewidth = 0.7
  ) +
  scale_x_continuous(
    breaks = seq(1990, max(yearly_boot3$year, na.rm = TRUE), by = 5)
  ) +
  labs(
    title = "Total Lamprey Wounds per 100 Fish by Year",
    subtitle = "Rainbow Trout Spring Assessment",
    x = "Year",
    y = "Wounds per 100 fish"
  ) +
  theme_classic(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")


#line plot of all counts per year
ggplot(by_year_long3, aes(x = year, y = count)) +
  geom_line(color = vcol, linewidth = 1) +
  geom_point(color = vcol, size = 2) +
  scale_x_continuous(
    breaks = seq(1990, max(yearly_boot3$year, na.rm = TRUE), by = 5)
  ) +
  labs(
    title = "Total Lamprey Wounds (raw counts) by Year",
    subtitle = "Rainbow Trout Spring Assessment",
    x = "Year",
    y = "Total wound count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )


# Stacked bar showing all counts per year
ggplot(by_year_ci3, aes(x = year, y = count, fill = category)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_x_continuous(
    breaks = seq(1990, max(yearly_boot3$year, na.rm = TRUE), by = 5)
  ) +
  scale_fill_viridis_d(option = "inferno", direction = 1, begin = 0.5) + # options: "viridis","magma","plasma","inferno","cividis","turbo"
  coord_cartesian(ylim = c(0, NA)) + # Avoid clipping CI tops
  geom_errorbar(
    aes(ymin = pmax(0, total_lo), ymax = total_hi),
    width = 0.18,
    linewidth = 0.6
  ) +
  labs(
    title = "Lamprey Wound Counts by Year",
    subtitle = "Rainbow Trout Spring Assessment",
    x = "Year",
    y = "Total wounds count",
    fill = "Category"
  ) +
  theme_classic(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")



#####################################################################################
#####################################################################################
##**                                  CF25_001                                   **##
#####################################################################################
#####################################################################################

################              CF1 - Read in Files             #######################

CF1 = read.csv(
  file = "G:/Fisheries Management/Assessment/Lamprey/Lamprey Data Request/2025/CF_001/LOA_CF25_001_bio.csv"
)

#2024
CF1_2024 = read.csv(
  file = "G:/Fisheries Management/Assessment/Lamprey/Lamprey Data Request/2024/CF24_001/LOA_CF24_001_bio.csv"
)

#2023
CF1_2023 = read_excel(
  "G:/Fisheries Management/Assessment/Lamprey/Lamprey Data Request/2023/LOA_CF23_001_lamprey.xlsx",
  sheet = 1
)

#2022
CF1_2022 = read_excel(
  "G:/Fisheries Management/Assessment/Lamprey/Lamprey Data Request/2022/CF1lamprey_biodata_export.xlsx",
  sheet = 1
)

#############    Reduce the tables down to only the lamprey wound columns     #######

CF1_2025_lamprey = CF1[18:25]
CF1_2024_lamprey = CF1_2024[20:27]
CF1_2023_lamprey = CF1_2023[19:26]
CF1_2022_lamprey = CF1_2022[19:26]

#################              Add year column                #######################

CF1_2025 = CF1_2025_lamprey %>%
  mutate(year = 2025)

CF1_2024 = CF1_2024_lamprey %>%
  mutate(year = 2024)

CF1_2023 = CF1_2023_lamprey %>%
  mutate(year = 2023)

CF1_2022 = CF1_2022_lamprey %>%
  mutate(year = 2022)


##################            Remainder of code               #######################

#combine all the rows from all lamprey tables
CF1_allyears = bind_rows(CF1_2025, CF1_2024, CF1_2023, CF1_2022)

#turn NA values into 0
CF1_allyears[is.na(CF1_allyears)] = 0

#add a column of total wounds = sum across all columns
Clean_CF1 = CF1_allyears %>%
  mutate(total_wounds = A1 + A2 + A3 + A4 + B1 + B2 + B3 + B4)


#yearly summary including rate per 100 fish
yearly4 = Clean_CF1 %>%
  group_by(year) %>%
  summarise(
    n_fish = n(),
    Wounds_total = sum(total_wounds, na.rm = TRUE),
    total_per100 = 100 * Wounds_total / n_fish,
    .groups = "drop"
  ) %>%
  arrange(year)


#long format for plotting all the years
by_year_long4 = yearly4 %>%
  select(year, Wounds_total) %>%
  pivot_longer(
    cols = Wounds_total,
    names_to = "category",
    values_to = "count"
  ) %>%
  mutate(category = toupper(category))



#Bootstrap
ci_boot4 <- Clean_CF1 |>
  group_by(year) |>
  reframe({
    qs <- bootstrap_rate_ci(total_wounds)
    tibble(total_per100_lo = qs[[1]], total_per100_hi = qs[[2]])
  })

yearly_boot4 <- yearly4 |>
  left_join(ci_boot4, by = "year")

#Poisson
by_year_ci4 <- by_year_long4 |>
  rowwise() |>
  mutate(
    ci = list(poisson_count_ci(count)),
    total_lo = ci["lo"],
    total_hi = ci["hi"]
  ) |>
  select(-ci) |>
  ungroup()


#Ploting rate per 100 fish
vcol <- viridis(1, option = "magma", begin = 0.8) # default
# options: "A", "B", "C", "D", "E", "F", "G", "H", "magma", "plasma", "inferno", "cividis"

ggplot(yearly_boot4, aes(x = year, y = total_per100)) +
  geom_ribbon(
    aes(ymin = total_per100_lo, ymax = total_per100_hi),
    fill = vcol,
    alpha = 0.20
  ) +
  geom_line(color = vcol, linewidth = 1) +
  geom_point(color = vcol, size = 2) +
  scale_x_continuous(
    breaks = seq(min(yearly_boot4$year), max(yearly_boot4$year), 1)
  ) +
  labs(
    title = "Total Lamprey Wounds per 100 Fish by Year",
    subtitle = "Commercial Catch Program",
    x = "Year",
    y = "Wounds per 100 fish"
  ) +
  theme_minimal()

#bar plot - rate per 100 fish
ggplot(yearly_boot4, aes(x = year, y = total_per100)) +
  geom_col(aes(fill = "bar"), position = "dodge", width = 0.7) +
  scale_fill_viridis_d(option = "turbo", direction = 1, begin = 1) +
  geom_errorbar(
    aes(ymin = total_per100_lo, ymax = total_per100_hi),
    width = 0.25,
    linewidth = 0.7
  ) +
  scale_x_continuous(
    breaks = seq(min(yearly_boot4$year), max(yearly_boot4$year), 1)
  ) +
  labs(
    title = "Total Lamprey Wounds per 100 Fish by Year",
    subtitle = "Commerical Catch Program",
    x = "Year",
    y = "Wounds per 100 fish"
  ) +
  theme_classic(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")


#line plot of all counts per year
ggplot(by_year_long4, aes(x = year, y = count)) +
  geom_line(color = vcol, linewidth = 1) +
  geom_point(color = vcol, size = 2) +
  scale_x_continuous(
    breaks = seq(
      min(by_year_long4$year, na.rm = TRUE),
      max(by_year_long4$year, na.rm = TRUE),
      by = 1
    )
  ) +
  labs(
    title = "Total Lamprey Wounds (raw counts) by Year",
    subtitle = "Commerical Catch Program",
    x = "Year",
    y = "Total wound count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )


# Stacked bar showing all counts per year
ggplot(by_year_ci4, aes(x = factor(year), y = count, fill = category)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_viridis_d(option = "inferno", direction = 1, begin = 0.5) + # options: "viridis","magma","plasma","inferno","cividis","turbo"
  coord_cartesian(ylim = c(0, NA)) + # Avoid clipping CI tops
  geom_errorbar(
    aes(ymin = pmax(0, total_lo), ymax = total_hi),
    width = 0.18,
    linewidth = 0.6
  ) +
  labs(
    title = "Wound Counts by Year",
    subtitle = "Commerical Catch Program",
    x = "Year",
    y = "Total wounds count",
    fill = "Category"
  ) +
  theme_classic(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

##########################END##############################################################################################
