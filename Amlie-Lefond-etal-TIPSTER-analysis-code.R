############################################################################## #
#
# Code to reproduce analysis in Amlie-Lefond et al. in review. 
# Risk of Symptomatic Intracranial Hemorrhage Following Intravenous Tissue 
# Plasminogen Activator for Acute Stroke is Not Increased in Children
# 
# Code by: Andy Cooper, Dwight Barry
# Date: 16 August 2019
# Version: 1.0
#
# URL: https://github.com/Rmadillo/Amlie-Lefondetal-TIPSTERS
# Contact: dwight [dot] barry {at} seattlechildrens [dot] org
# PRs, suggestions, improvements, and corrections welcome!
#
# License: This code is released under GNU GENERAL PUBLIC LICENSE Version 3
#
# This code is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# It is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this code. If not, see <https://www.gnu.org/licenses/>.
#
############################################################################## #



############################################################################## #
# A. SET UP                                                                 ####
############################################################################## #


#### Load packages ####
library(tidyverse)
library(viridis)
library(forcats)



############################################################################## #
# B. GET DATA                                                               ####
############################################################################## #


#### Set up Bayesian model data ####


## Bayesian beta-binomial model for probability of risk of SICH after tPA


# SICH incidents

obs = 0


# Number of tPA patients in TIPS

npats = 26


# Set up x axis points for posterior plot

x = seq(0, 1, 0.005)


# Primary prior assumption: same as young adult (18-40) risk, 1.7% (mode)
#   Dodds, J. A., Y. Xian, S. Sheng, G. C. Fonarow, D. L. Bhatt, R. Matsouaka, 
#       L. H. Schwamm, E. D. Peterson and E. E. Smith (2019). "Thrombolysis in  
#       young adults with stroke: Findings from Get With The Guidelines-Stroke." 
#       Neurology 92(24): e2784-e2792.

alpha_17 = 1.5
beta_17 = (alpha_17 - 1 - 0.017 * alpha_17 + 2 * 0.017) / 0.017
df_17 = data.frame(x, y = pbeta(x, alpha_17 + obs, beta_17 + npats - obs)) %>% 
          mutate(prior = "1.7% Risk Prior")


# No apriori assumption of risk, uniform prior

alpha_unif = 1
beta_unif = 1
df_unif = data.frame(x, y = pbeta(x, alpha_unif + obs, beta_unif + npats - obs)) %>% 
        mutate(prior = "Uniform Prior")


# Conservative prior assumption: risk is similar to old adult estimate, 6.4% (mode)
#   National Institute of Neurological Disorders and Stroke rt-PA Stroke Study
#       Group (1995). Tissue plasminogen activator for acute ischemic stroke. 
#       New England Journal of Medicine, 333(24), 1581-1588.

alpha_64 = 1.5
beta_64 = (alpha_64 - 1 - 0.064 * alpha_64 + 2 * 0.064) / 0.064
df_64 = data.frame(x, y = pbeta(x, alpha_64 + obs, beta_64 + npats - obs)) %>%
        mutate(prior = "6.4% Risk Prior")


# Put all three together for plotting, order by assumption

df_figure1 = bind_rows(df_17, df_unif, df_64) %>%
    mutate(prior = factor(prior, levels = c("1.7% Risk Prior",
        "Uniform Prior", "6.4% Risk Prior")))



#### Data for Figures ####

# Data to build Figure 2

figure2_data = read_csv("Amlie-Lefond-etal-TIPSTER-Figure2-data.csv") %>%
    mutate(NIHSS_decrease = NIHSS_624 - NIHSS_initial)


#### Tables from paper, for Supplemental Figures 2 and 3 ####

table1 = read_csv("Amlie-Lefond-etal-TIPSTER-Table1.csv")

table2 = read_csv("Amlie-Lefond-etal-TIPSTER-Table2.csv")





############################################################################## #
#### C. ANALYSIS                                                            ####
############################################################################## #


#### Risk values and Figure 1: Posterior cdf plot ####

## Risk values

# Posterior medians

qbeta(0.5, alpha_17 + obs, beta_17 + npats - obs)       # 1.7% risk prior
qbeta(0.5, alpha_unif + obs, beta_unif + npats - obs)   # uniform prior
qbeta(0.5, alpha_64 + obs, beta_64 + npats - obs)       # 6.4% risk prior


# Posterior modes

(alpha_17 + obs - 1)/(alpha_17 + obs + beta_17 + npats - obs - 2)         # 1.7% risk prior
(alpha_unif + obs - 1)/(alpha_unif + obs + beta_unif + npats - obs - 2)   # uniform prior
(alpha_64 + obs - 1)/(alpha_64 + obs + beta_64 + npats - obs - 2)         # 6.4% risk prior


# HDPIs from grid approximation for 1.7% risk prior

df_grid_17 = tibble(p_grid = seq(from = 0, to = 1, length.out = 10001),
                    prior = dbeta(p_grid, alpha_17, beta_17),
                    likelihood = dbinom(obs, size = npats, prob = p_grid), 
                    posterior = (likelihood * prior) / sum(likelihood * prior))


samples = df_grid_17 %>% 
    sample_n(size = 1e4, weight = posterior, replace = TRUE)


rethinking::HPDI(samples$p_grid, 0.95)



## Figure 1 

ggplot(df_figure1, aes(x = x, y = y, color = prior)) +
    
    geom_path(aes(linetype = prior)) +
    scale_linetype_manual(values=c("solid", "longdash", "dotdash")) +
    
    geom_point(data = df_17, show.legend = F) +
    
    scale_color_manual(values = c("black", "gray30", "gray60")) +
    
    scale_x_continuous(
        name = "Percent Risk of SICH", 
            labels = scales::percent, limits = c(0, 0.16), 
            breaks = seq(0, 1, 0.05), minor_breaks = seq(0, 1, 0.01),
            expand = c(0.001, 0)) +
    
    scale_y_continuous(
        name = expression(paste("Probability that risk of SICH is below ", 
               italic("x"))),
            labels = scales::percent, limits = c(0, 1),
            breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0.001, 0.001)) +
    
    labs(subtitle = paste0("Total IV tPA Patients = ", npats, ", Observed SICH = ", obs), 
         color = "Prior Assumption: ", 
         linetype = "Prior Assumption: ") +

    theme_bw() +
    theme(legend.position = "bottom")
      

# Save to local directory
ggsave("TIPSTER_Figure_1_Bayesian_Regression.png", width = 6.8, height = 4, units = "in", dpi = 320)



#### NIHSS/PSOM results and Figure 2: Outcomes ####


## NIHSS and PSOM results

NIHSS = figure2_data %>%
    select(RecordID, NIHSS_initial, NIHSS_624) %>%
    gather(SCORE, VALUE, 2:3) %>%
    mutate(SCORE = ordered(SCORE, levels = c("NIHSS_initial", "NIHSS_624")),
           RecordID = factor(RecordID)) %>%
    as_tibble()


NIHSS_counts = NIHSS %>%
    filter(complete.cases(.)) %>%
    group_by(SCORE) %>%
    summarise(COUNT = n(),
              MED = round(median(VALUE, na.rm = T), 0),
              PCT = round(COUNT/npats, 2),
              LABEL = paste0("n = ", COUNT))


# NIHSS Medians, n, pct
NIHSS_counts


# Compare the two groups, regardless of pre-post combination
psych::describeBy(NIHSS$VALUE, NIHSS$SCORE, quant = c(0.25, 0.5, 0.75), 
                  digits = 0, mat = T)


# Those that have both pre- and post NIHSS scores
NIHSS_both = filter(figure2_data, !is.na(NIHSS_initial) & !is.na(NIHSS_624))
psych::describe(NIHSS_both$NIHSS_624, quant = c(0.25, 0.5, 0.75))


# Look at pattern of decreases
NIHSS_decreases = filter(figure2_data, NIHSS_decrease < 0)
median(NIHSS_decreases$NIHSS_decrease, na.rm = T)




## Figure 2 

# Figure 2A

NIHSS_plot = ggplot(NIHSS, aes(SCORE, VALUE)) + 
    geom_line(aes(group = RecordID),  alpha = 0.5) + 
    geom_point(alpha = 0.5) +
    geom_line(data = NIHSS_counts, aes(SCORE, MED, group = 1), 
              color = "black", size = 3, alpha = 0.4) + 
    scale_y_continuous(limits = c(-0.5, 44.5), breaks = seq(0, 40, 5), 
                       minor_breaks = NULL, expand = c(0.01, 0)) +
    scale_x_discrete(labels = c("NIHSS_initial" = "Before tPA", 
                              "NIHSS_624" = "6-24 hours\nlater")) +
    labs(x = "", y = "NIHSS Score") +
    theme_minimal() +
    theme(legend.position="none") + 
    annotate(geom = "text", x = 1:2, y = 42.5, 
             label = paste0(NIHSS_counts$LABEL),
             size = 2.5, fontface = 3, hjust = 0.5, vjust = 0.5) 


# Figure 2B

PSOM = figure2_data %>%
    select(RecordID, PSOM_7dc, PSOM_3months, PSOM_1year) %>%
    gather(SCORE, VALUE, 2:4) %>%
    mutate(SCORE = ordered(SCORE, 
                levels = c("PSOM_7dc", "PSOM_3months", "PSOM_1year")),
           RecordID = factor(RecordID)) %>%
    as_tibble()


PSOM_counts = PSOM %>%
    filter(complete.cases(.)) %>%
    group_by(SCORE) %>%
    summarise(COUNT = n(),
              MED = round(median(VALUE, na.rm = T), 0),
              PCT = round(COUNT/26, 2),
              LABEL = paste0("n = ", COUNT))


# PSOM medians, n, pct
PSOM_counts


PSOM_plot = ggplot(PSOM, aes(SCORE, VALUE, group = RecordID)) + 
    geom_line(alpha = 0.5) + 
    geom_point(alpha = 0.5) +
    geom_line(data = PSOM_counts, aes(SCORE, MED, group = 1), 
              color = "black", size = 3, alpha = 0.4) + 
    scale_y_continuous(limits = c(-0.1, 11.1), breaks = c(0, 2, 4, 6, 8, 10),
                       minor_breaks = NULL, expand = c(0.01, 0)) +
    scale_x_discrete(labels=c("PSOM_7dc" = "Discharge or\n7 days later", 
                              "PSOM_3months" = "3 months\nlater",
                              "PSOM_1year" = "1 year\n later")) +
    expand_limits(y = 0) +
    labs(x = "", y = "PSOM Score") +
    theme_minimal() +
    theme(legend.position="none") + 
    annotate(geom = "text", x = 1:3, y = 10.6, 
        label = paste0(PSOM_counts$LABEL),
        size = 2.5, fontface = 3, hjust = 0.5, vjust = 0.5) 


cowplot::plot_grid(NIHSS_plot, PSOM_plot, labels = c("A", "B"),
                   align = "h", nrow = 1, rel_widths = c(3.025, 4.1))


ggsave("TIPSTER_Figure_2_NIHSS_PSOM.png", width = 6.8, height = 4, units = "in", dpi = 320)



#### Suppl Figure 1: distributions of priors, likelihood, and posteriors #### 


# Prior distributions

df_prior_17 = tibble(x, y = dbeta(x, alpha_17, beta_17)) %>% 
    mutate(prior = "1.7% Risk Prior")

df_prior_unif = tibble(x, y = dbeta(x, alpha_unif, beta_unif)) %>% 
    mutate(prior = "Uniform Prior")

df_prior_64 = tibble(x, y = dbeta(x, alpha_64, beta_64)) %>% 
    mutate(prior = "6.4% Risk Prior")

df_priors = bind_rows(df_prior_17, df_prior_unif, df_prior_64) %>%
    mutate(prior = factor(prior, levels = c("1.7% Risk Prior",
                                            "Uniform Prior", 
                                            "6.4% Risk Prior")))


# Likelihood function

likelihood = tibble(x = x, like_density = dbinom(x = obs, size = npats, 
                                    prob = x, log = FALSE))

likelihood_groups = tibble(group = rep(c(" 0 SICH in 26 patients", 
                                            " 0 SICH in 26 patients ", 
                                            "0 SICH in 26 patients "), each = 201),
                         x = rep(likelihood$x, 3),
                         like_density = rep(likelihood$like_density, 3))


# Posterior distributions

df_post_17 = data.frame(x, y = dbeta(x, alpha_17 + obs, beta_17 + npats - obs)) %>% 
    mutate(prior = "1.7% Risk Prior")

df_post_unif = data.frame(x, y = dbeta(x, alpha_unif + obs, beta_unif + npats - obs)) %>% 
    mutate(prior = "Uniform Prior")

df_post_64 = data.frame(x, y = dbeta(x, alpha_64 + obs, beta_64 + npats - obs)) %>% 
    mutate(prior = "6.4% Risk Prior")

df_posts = bind_rows(df_post_17, df_post_unif, df_post_64) %>%
    mutate(prior = factor(prior, levels = c("1.7% Risk Prior",
                                            "Uniform Prior", 
                                            "6.4% Risk Prior")))


## Figure SM.I

priors_plot = ggplot(df_priors, aes(x = x, y = y, color = prior)) +
  
    geom_path(size = 1, aes(linetype = prior)) +
    
    scale_linetype_manual(values=c("solid", "longdash", "dotdash")) +
    scale_color_manual(values = c("black", "#BCAF6FFF", "#61646FFF")) +
    
    scale_x_continuous(name = "Probability of SICH", expand = c(0.001, 0.0025)) +
    scale_y_continuous(name = "", labels = NULL, 
                     breaks = NULL, expand = c(0.01, 0.01)) +
    
    expand_limits(x = 1.005) +
    
    labs(title = "Prior Probability", color = "Prior Assumption: ", 
         linetype = "Prior Assumption: ") +
    
    facet_wrap(~prior, ncol = 1, scales = "free") +
    
    theme_bw(base_size = 10) + 
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))


like_plot = ggplot(likelihood_groups, aes(x, like_density)) + 
    
    geom_line(size = 1) +
    
    scale_x_continuous(expand = c(0.001, 0.0025)) +
    scale_y_continuous(labels = NULL, breaks = NULL, expand = c(0.01, 0.01)) +
    
    expand_limits(x = 1.005) +
    
    labs(x = "Likelihood", y = "", title = "Likelihood Function") +
    
    facet_wrap(~group, ncol = 1, scales = "free") +
    
    theme_bw(base_size = 10) + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          strip.text = element_text(hjust = 0.5))


post_plot = ggplot(df_posts, aes(x = x, y = y, color = prior)) +
    
    geom_path(size = 1, aes(linetype = prior)) +
    
    scale_linetype_manual(values=c("solid", "longdash", "dotdash")) +
    scale_color_manual(values = c("black", "#BCAF6FFF", "#61646FFF")) +
    
    scale_x_continuous(name = "Probability of SICH", expand = c(0.001, 0.0025)) +
    scale_y_continuous(name = "", labels = NULL, breaks = NULL, expand = c(0.01, 0.01)) +
    
    expand_limits(x=1.005) +
  
    labs(title = "Posterior Probability", color = "Prior Assumption: ") +
    
    facet_wrap(~prior, ncol = 1, scales = "free") +
    
    theme_bw(base_size = 10) + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))


# Legend and spacers

priors_plot_nolege = priors_plot + theme(legend.position = "none")


priors_plot_leggie = ggplot(df_priors, aes(x = x, y = y, color = prior)) +
    geom_path(size = 1) +
    scale_color_manual(values = c("black", "#BCAF6FFF", "#61646FFF")) +
    scale_linetype_manual(values=c("solid", "longdash", "dotdash")) +
    theme_bw(base_size = 10) + 
    labs( color = "Prior Assumption: ", linetype = "Prior Assumption: ") +
    theme(legend.position = "bottom")


leggie = cowplot::get_legend(priors_plot_leggie) 


times = ggplot() + 
    labs(title = " * ") + 
    theme_void() + 
    theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"))


proportional = ggplot() + 
    labs(title = expression(" " %prop% " ")) + 
    theme_void() + 
    theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"))


# Plot Supplemental Figure I

cowplot::plot_grid(priors_plot_nolege, times, like_plot, proportional, post_plot, 
                   nrow = 1, scale = 0.95, rel_widths = c(2, 0.3, 2, 0.3, 2),
                   hjust = c(-1, -1, -1))


ggsave("TIPSTER_Supplemental_Figure_I.png", width = 6.5, height = 5, dpi = 320)



#### Suppl Figure 2: No patterns in demographics, etc. ####


ggplot(table1, aes(Age, Time_to_tPA, 
                   fill = Sex, 
                   size = forcats::fct_rev(Location),
                   shape = forcats::fct_rev(Location))) +
    
    geom_point(alpha = 0.5, color = "black") +
    geom_point(data = filter(table1, Sex == "Female"), alpha = 0.5, color = "black") +
    
    scale_y_continuous(limits = c(0, 5), expand = c(0.01, 0)) +
    scale_x_continuous(limits = c(0, 20), minor_breaks = seq(1:19), 
                       expand = c(0.01, 0)) +
    
    scale_fill_viridis_d() +
    coord_flip() +
    scale_size_discrete(range = c(4, 2)) +
    scale_shape_manual(values = c(21, 23)) +
    labs(x = "Patient Age (years)", y = "Time to tPA (hours)", 
         fill = "Sex", size = "Location",
         shape = "Location") +
    
    theme_minimal() +
    guides(fill = guide_legend(nrow = 2, byrow = F, 
        override.aes = list(shape = 22, size = 6, alpha = 0.8)), 
        shape = guide_legend(nrow = 2, byrow = T, override.aes=list(size = c(4, 2))),
        size = guide_legend(nrow = 2, byrow = T)) 


ggsave("TIPSTER_Supplemental_Figure_II.png", dpi = 320, width = 6.5, height = 4, units = "in")



#### Suppl Figure 3: scanning time ####


table1a = select(table1, ID, Age, Sex, Location)

table2a = inner_join(table2, table1a) 


ggplot(table2a, aes(Time, ASPECTS_score, 
                   fill = Sex, 
                   size = Age, 
                   shape = forcats::fct_rev(Location))) +
    
    geom_point(alpha = 0.5) +
    geom_point(data = filter(table2a, Sex == "Female"), alpha = 0.25, 
               color = "black") +
    
    scale_y_continuous(limits = c(1, 10.5), breaks = seq(1:10), minor_breaks = NULL, 
                       expand = c(0.01, 0)) +
    scale_x_continuous(limits = c(0, 200), breaks = c(0, 60, 120, 180), 
                       expand = c(0.01, 0)) +
    
    scale_fill_viridis_d() +
    scale_shape_manual(values = c(21, 23)) +
    scale_size(range = c(0, 4)) +
    labs(y = "ASPECTS Score", x = "Time to CT (minutes)", fill = "Sex", size = "Age",
         shape = "Location") +
    
    theme_minimal() +
    guides(fill = guide_legend(nrow = 2, byrow = F, 
           override.aes = list(shape = 22, size = 6, alpha = 0.8)), 
           shape = guide_legend(nrow = 2, byrow = T, override.aes=list(size = 4)),
           size = guide_legend(nrow = 2, byrow = T)) 


ggsave("TIPSTER_Supplemental_Figure_III.png", dpi = 320, width = 6.5, height = 4, units = "in")



#### End of file ####
