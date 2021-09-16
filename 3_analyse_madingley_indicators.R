
## REPOSITORY: https://github.com/conservationscience/madingley_terrestrial_indicators

rm(list = ls())

# Directory path to git repo

# Simone's deakin laptop

# cd "C:\\Users\\ssteven\\OneDrive - Deakin University\\Deakin\\Chapter_3_indicator_testing\\madingley_terrestrial_indicators"

# Data structure

#' Input data is structured hierarchically as follows:
#' # Location (1)
#' ## Scenarios (3)
#' ### Replicates (25)
#' 
#' Output indicator data is calculated at the same level, so each indicator should
#' have values for:
#' # Location (1)
#' ## Scenarios (3)
#' ### Replicates (25) 


# TODO LIST ----

# Libraries ----

## Data wrangling
library(tidyverse)
library(tidylog)

## Analysis
library(mgcv)
library(nlme)
library(changepoint)
library(MASS)
library(strucchange)

## Viz
library(ggpubr)
library(cowplot)

# Functions ----

#' Description

#' @param 
#' @param 
#' @param 
#' @return 

function_name <- function(param){

}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

input <- data

smooth_gam <- function(input, smoothing_param) {
  
  input <- input %>% 
    mutate(abundance = abundance +
             (mean(abundance)*0.01))
  
  # mod <- gam(log10(abundance) ~ s(annual_time_step), sp = smoothing_param,  
  #            data = input, method = "REML")
  
  mod <- gam(log10(abundance) ~ s(annual_time_step, k = 4), sp = smoothing_param,  
             data = input, method = "REML")
  
  # plot(mod)
  # summary(mod)
  # gam.check(mod)
  
  modelled_abundance <- predict(mod, pred_annual_time_step = unique(input$annual_time_step))
  
  modelled_abundance <- as.data.frame(cbind(annual_time_step =unique(input$annual_time_step),
                                            transformed_abundance = exp(modelled_abundance)))
  
  newdata <- input %>% 
    merge(modelled_abundance, by = "annual_time_step")
  
  return(newdata)
  
  
}

## Model Checking function
# https://raw.githubusercontent.com/gavinsimpson/random_code/master/tsDiagGamm.R

tsDiagGamm <- function(x, timevar, observed, f = 0.3, type = "normalized") {
  resi <- resid(x$lme, type = type)
  fits <- fitted(x$lme)
  on.exit(layout(1))
  layout(matrix(1:6, ncol = 3, byrow = TRUE))
  plot(resi ~ fits, ylab = "Normalized Residuals",
       xlab = "Fitted Values", main = "Fitted vs. Residuals")
  lines(lowess(x = fits, y = resi, f = f), col = "blue",
        lwd = 2)
  plot(resi ~ timevar, ylab = "Normalized Residuals",
       xlab = "Time", main = "Time series of residuals")
  lines(lowess(x = timevar, y = resi, f = f), col = "blue", lwd = 2)
  plot(observed ~ fits, ylab = "Observed",
       xlab = "Fitted Values", main = "Fitted vs. Observed",
       type = "n")
  abline(a = 0, b = 1, col = "red")
  points(observed ~ fits)
  lines(lowess(x = fits, y = observed, f = f), col = "blue",
        lwd = 2)
  hist(resi, freq = FALSE, xlab = "Normalized Residuals")
  qqnorm(resi)
  qqline(resi)
  acf(resi, main = "ACF of Residuals")
}

# Derivative gam functions

# https://raw.githubusercontent.com/gavinsimpson/random_code/master/derivFun.R

Deriv <- function(mod, n = 300, eps = 1e-7, newdata) {
  
    #if(isTRUE(all.equal(class(mod), "list")))
    mod <- mod$gam
    m.terms <- attr(terms(mod), "term.labels")
  
    if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x), length = n))
    names(newD) <- m.terms
  
    } else {
    newD <- newdata
  }
  X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
  newD <- newD + eps
  X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  # number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  for(i in seq_len(nt)) {
    Xi <- Xp * 0
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(mod)
    df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
    ## Xi <- Xp * 0 ##matrix(0, nrow = Xp.r, ncol = Xp.c)
    ## J <- bs.dims[i]
    ## Xi[,(i-1) * J + 1:J + 1] <- Xp[,(i-1) * J + 1:J +1]
    ## df <- Xi %*% coef(mod)
    ## df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    ## lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  return(lD)
}

confint.Deriv <- function(object, term, alpha = 0.05, ...) {
  l <- length(object) - 3
  term.labs <- names(object[seq_len(l)])
  if(missing(term))
    term <- term.labs
  Term <- match(term, term.labs)
  ##term <- term[match(term, term.labs)]
  if(any(miss <- is.na(Term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  ## if(is.na(term))
  ##     stop("'term' not a valid model term.")
  res <- vector(mode = "list", length = length(term))
  names(res) <- term
  residual.df <- length(object$gamModel$y) - sum(object$gamModel$edf)
  tVal <- qt(1 - (alpha/2), residual.df)
  ## tVal <- qt(1 - (alpha/2), object$gamModel$df.residual)
  for(i in seq_along(term)) {
    upr <- object[[term[i]]]$deriv + tVal * object[[term[i]]]$se.deriv
    lwr <- object[[term[i]]]$deriv - tVal * object[[term[i]]]$se.deriv
    res[[term[i]]] <- list(upper = drop(upr), lower = drop(lwr))
  }
  res$alpha = alpha
  res
}

signifD <- function(x, d, upper, lower, eval = 0) {
  miss <- upper > eval & lower < eval
  incr <- decr <- x
  want <- d > eval
  incr[!want | miss] <- NA
  want <- d < eval
  decr[!want | miss] <- NA
  list(incr = incr, decr = decr)
}

plot.Deriv <- function(x, alpha = 0.05, polygon = TRUE,
                       sizer = FALSE, term, eval = 0, lwd = 3,
                       col = "lightgrey", border = col,
                       ylab, xlab, ...) {
  l <- length(x) - 3
  ## get terms and check specified (if any) are in model
  term.labs <- names(x[seq_len(l)])
  if(missing(term))
    term <- term.labs
  Term <- match(term, term.labs)
  if(any(miss <- is.na(Term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  if(all(is.na(Term)))
    stop("All terms in 'term' not found in model.")
  l <- sum(!miss)
  nplt <- n2mfrow(l)
  ## tVal <- qt(1 - (alpha/2), x$gamModel$df.residual)
  residual.df <- length(x$gamModel$y) - sum(x$gamModel$edf)
  tVal <- qt(1 - (alpha/2), residual.df)
  if(missing(ylab))
    ylab <- expression(italic(hat(f)*"'"*(x)))
  if(missing(xlab)) {
    xlab <- attr(terms(x$gamModel), "term.labels")[Term]
    names(xlab) <- xlab
  }
  layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
  CI <- confint(x, term = term, alpha = alpha)
  for(i in seq_along(term)) {
    ## for(i in seq_len(l)) {
    upr <- CI[[term[i]]]$upper
    lwr <- CI[[term[i]]]$lower
    ylim <- range(upr, lwr)
    plot(x$eval[,term[i]], x[[term[i]]]$deriv, type = "n",
         ylim = ylim, ylab = ylab, xlab = xlab[term[i]], ...)
    if(isTRUE(polygon)) {
      polygon(c(x$eval[,term[i]], rev(x$eval[,term[i]])),
              c(upr, rev(lwr)), col = col, border = border)
    } else {
      lines(x$eval[,term[i]], upr, lty = "dashed")
      lines(x$eval[,term[i]], lwr, lty = "dashed")
    }
    abline(h = 0, ...)
    if(isTRUE(sizer)) {
      lines(x$eval[,term[i]], x[[term[i]]]$deriv, lwd = 1)
      S <- signifD(x[[term[i]]]$deriv, x[[term[i]]]$deriv, upr, lwr,
                   eval = eval)
      lines(x$eval[,term[i]], S$incr, lwd = lwd, col = "blue")
      lines(x$eval[,term[i]], S$decr, lwd = lwd, col = "red")
    } else {
      lines(x$eval[,term[i]], x[[term[i]]]$deriv, lwd = 2)
    }
  }
  layout(1)
  invisible(x)
}


# Set up paths ----

indicators_project <- "N:/Quantitative-Ecology/Indicators-Project"
  
# Inputs ----

# Get date to label outputs

today <- Sys.Date()

location <- 'Serengeti'

scenarios <- list("baseline", "land use", 
                  "carnivore harvesting", 
                  "herbivore harvesting")

disturbance_string <- c("pre-disturbance", "disturbance", "post-disturbance")

disturbance_factor <- factor(dist, ordered = TRUE, 
                             levels = c("pre-disturbance", 
                                        "disturbance", 
                                        "post-disturbance"))

# Set up output folders

analysis_inputs_folder <- file.path(indicators_project,  
                           "/Serengeti/Outputs_from_analysis_code/Analysis_inputs",
                           today)

if( !dir.exists( file.path(analysis_inputs_folder) ) ) {
  dir.create( file.path(analysis_inputs_folder), recursive = TRUE )
  
}

analysis_outputs_folder <- file.path(indicators_project, 
                                     "/Serengeti/Outputs_from_analysis_code/Analysis_outputs",
                                     today)

if( !dir.exists( file.path(analysis_outputs_folder) ) ) {
  dir.create( file.path(analysis_outputs_folder), recursive = TRUE )
  
}

analysis_plots_folder <- file.path(indicators_project, 
                                     "/Serengeti/Outputs_from_analysis_code/Analysis_plots_folder",
                                   today)

if( !dir.exists( file.path(analysis_plots_folder) ) ) {
  dir.create( file.path(analysis_plots_folder), recursive = TRUE )
  
}

# Load indicator > scenario level data ----

indicators_all_df <- readRDS("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_indicator_code/Indicator_outputs/general/2021-09-13/2021-09-14_all_indicators_output_dataframe.rds")
head(indicators_all_df)
unique(indicators_all_df$indicator)

# Indicator time series ----

indicators_all_df_list <- split(indicators_all_df, 
                             indicators_all_df$indicator)

indicator_time_series_plots <- list()

for (i in seq_along(indicators_all_df_list)) {
  
  indicator_name <- indicators_all_df_list[[i]]$indicator[1]
  
  if(indicator_name == "abundance harvested groups") {

    indicator_data <- indicators_all_df_list[[i]] %>%
                      filter(scenario != "000_Baseline") %>%
                      mutate(scenario = ifelse(scenario == "100_Land_use",
                             "1 - Land use scenario",
                             ifelse(scenario == "200_Harvesting_carnivores",
                             "2 - Carnivore harvesting scenario",
                             ifelse(scenario == "300_Harvesting_herbivores",
                             "3 - Herbivore harvesting scenario", NA))))


    indicator_time_series_plots[[i]] <- ggplot(data = indicator_data,
                                               aes(x = annual_time_step,
                                                   y = log(indicator_score),
                                                   col = scenario)) +
      geom_line() +
      scale_color_manual(values = c("#440154FF",
                                    "#39568CFF",
                                    "#1F968BFF")) +
      facet_wrap( ~ scenario, ncol = 3) +
      theme_classic() +
      theme(panel.background = element_rect(fill = "snow2"),
            legend.position = "none",
            strip.text = element_text(face="bold", size=9)) +
      geom_vline(xintercept = 100, linetype = "dashed") +
      annotate(x=100,y=+Inf,label="Impact start",vjust=2,geom="label",
               size = 3) +
      geom_vline(xintercept = 200, linetype = "dashed") +
      annotate(x=200,y=+Inf,label="Impact end",vjust=2,geom="label",
               size = 3) +
      labs(x = "Annual time step",
           y = paste(indicator_name, "(log)", sep = " "))

  } else {
  
  indicator_data <- indicators_all_df_list[[i]] %>% 
    filter(scenario != "000_Baseline") %>% 
    mutate(scenario = ifelse(scenario == "100_Land_use",
                             "1 - Land use scenario",
                             ifelse(scenario == "200_Harvesting_carnivores",
                                    "2 - Carnivore harvesting scenario",
                                    ifelse(scenario == "300_Harvesting_herbivores",
                                           "3 - Herbivore harvesting scenario", NA))))
  
indicator_time_series_plots[[i]] <- ggplot(data = indicator_data, 
                                           aes(x = annual_time_step, 
                                    y = indicator_score,
                                    col = scenario)) +
    geom_line() +
    scale_color_manual(values = c("#440154FF",
                         "#39568CFF",
                         "#1F968BFF")) +
    facet_wrap( ~ scenario, ncol = 3) + 
    theme_classic() +
    theme(panel.background = element_rect(fill = "snow2"),
          legend.position = "none",
          strip.text = element_text(face="bold", size=9)) +
    geom_vline(xintercept = 100, linetype = "dashed") +
    annotate(x=100,y=+Inf,label="Impact start",vjust=2,geom="label",
             size = 3) +
    geom_vline(xintercept = 200, linetype = "dashed") +
    # annotate(x=200,y=+Inf,label="Impact end",vjust=2,geom="label",
    #          size = 3) +
    labs(x = "Annual time step",
         y = indicator_name)
  }
}


RLI <- indicator_time_series_plots[[1]]
RLI_large <- indicator_time_series_plots[[2]]
LPI <- indicator_time_series_plots[[3]]
harvested <- indicator_time_series_plots[[4]]
  


fig <- plot_grid(harvested, RLI, RLI_large, LPI, align = "v", 
                 nrow = 4, rel_heights = c(1/4, 1/4, 1/4, 1/4))

fig

ggsave(file.path(analysis_plots_folder, 
                 paste(today,
                       "indicator_time_series.png",
                       sep = "_")), fig,  device = "png")

# Distribution histograms ----

indicators_all_list <- readRDS("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_indicator_code/Indicator_outputs/general/2021-09-13/2021-09-14_all_indicators_output_list.rds")

# indicator_histograms <- list()
# 
# for (i in seq_along(indicators_all)) {
#   
#   single_indicator <- indicators_all[[i]]
#  
#   scenario_histograms <- list()
#   
#   for (j in seq_along(single_indicator)) {
#     
#     scenario_name <- str_to_title(scenarios[[j]])
#     
#     indicator_name <- single_indicator[[j]]$indicator[1]
#    
#     
# 
#     
#     ggsave(file.path(analysis_plots_folder,
#                      paste(indicator_name,"_", scenario_name, ".png", sep = "")),
#            scenario_scatterplots[[j]],  device = "png")
#     
#   }
#   
#   names(scenario_scatterplots) <- scenarios
#   
#   indicator_scatterplots[[i]] <- scenario_scatterplots
#   
# }

# Correlation analysis ----

indicator_scatterplots <- list()

for (i in seq_along(indicators_all_list)) {
  
  single_indicator <- indicators_all_list[[i]]
  harvest_indicator <- indicators_all_list[[length(indicators_all_list)]]
  
  scenario_scatterplots <- list()
  
  for (j in seq_along(single_indicator)) {
    
    scenario_name <- str_to_title(scenarios[[j]])
    
    indicator_name <- single_indicator[[j]]$indicator[1]
    harvest_name <- "Harvested group abundance"
    
    scatterplot_data <- single_indicator[[j]][c("annual_time_step",
                                      "indicator_score")] %>% 
                        rename(indicator = indicator_score) %>% 
                        merge(harvest_indicator[[j]][c("annual_time_step",
                                "indicator_score")],
                        by = "annual_time_step") %>% 
                        rename(harvested = indicator_score) %>% 
                        mutate(disturbance = ifelse(annual_time_step < 100, 
                                  "pre-disturbance",
                                  ifelse(annual_time_step >= 100 &
                                         annual_time_step < 200, 
                                         "disturbance", 
                                         "post-disturbance"))) 
    head(scatterplot_data)
                      
    
    scenario_scatterplots[[j]] <- ggplot(scatterplot_data, 
                                         aes(x = range01(harvested), 
                                             y = range01(indicator),
                           col = disturbance)) +
      geom_point() +
      labs(x = harvest_name,
           y = indicator_name,
           title = paste(indicator_name, scenario_name, sep = " ")) + 
      stat_cor(method = "spearman")
    
    ggsave(file.path(analysis_plots_folder,
                     paste(indicator_name,"_", scenario_name, ".png", sep = "")),
                     scenario_scatterplots[[j]],  device = "png")
    
  }
  
  names(scenario_scatterplots) <- scenarios
  
  indicator_scatterplots[[i]] <- scenario_scatterplots
  
}

names(indicator_scatterplots) <- names(indicators_all_list)

indicator_scatterplots[["LPI"]][[1]]
indicator_scatterplots[["LPI"]][[2]]
indicator_scatterplots[["LPI"]][[3]]
indicator_scatterplots[["LPI"]][[4]]

indicator_scatterplots[["RLI all"]][[1]]
indicator_scatterplots[["RLI all"]][[2]]
indicator_scatterplots[["RLI all"]][[3]]
indicator_scatterplots[["RLI all"]][[4]]

indicator_scatterplots[["RLI large spp"]][[1]]
indicator_scatterplots[["RLI large spp"]][[2]]
indicator_scatterplots[["RLI large spp"]][[3]]
indicator_scatterplots[["RLI large spp"]][[4]]

## * Calculate correlation coefficients ----

indicator_cor_scores <- list()

for (i in seq_along(indicators_all_list)) {
  
  # Get the indicator (all scenarios)
  
  indicator_scenarios <- indicators_all_list[[i]]
  
  print(indicator_scenarios[[1]]$indicator[1])
  
  harvested_scenarios <- indicators_all_list[["abundance harvested groups"]]
  
  # Make a list to hold scenario correlation coefficients
  
  scenario_cor_scores <- list()
  
  for (j in seq_along(indicator_scenarios)) {
    
  print(scenarios[i])
    
  indicator <- indicator_scenarios[[j]] %>% 
               dplyr::select(annual_time_step, indicator_score) %>% 
               rename(indicator = indicator_score)
  
  harvested<- harvested_scenarios[[j]] %>% 
              dplyr::select(annual_time_step, indicator_score) %>% 
              rename(harvested = indicator_score)
  
  comparison <- indicator %>% 
                merge(harvested, by = "annual_time_step")
    
  cor <- cor(comparison$indicator, 
             comparison$harvested, method = "spearman")
  
  scenario_cor_scores[[j]] <- data.frame(indicator = indicator_scenarios[[j]]$indicator[1],
                                         correlation = cor,
                                         scenario = scenarios[[j]])
  
  }
  
  scenario_cor_df <- do.call(rbind, scenario_cor_scores)
  
  indicator_cor_scores[[i]] <- scenario_cor_df

}

correlation_dataframe <- do.call(rbind, indicator_cor_scores) %>% 
                         arrange(cor)

# Breakpoint analyais ----


scenario_indicator_names <- names(indicators_all_list)

scenario_changepoint_summaries <- list()
scenario_variance_summaries <- list()

for (i in seq_along(indicators_all_list)) {

# Get a single indicator data
  
single_indicator <- indicators_all_list[[i]]

indicator_changepoint_summaries <- list()
indicator_variance_summaries <- list()

  for (j in seq_along(single_indicator)) {

  # Get a single scenario for the indicator
    
  data <- single_indicator[[j]]
    
  # Convert into a time series
  
  indicator_ts <- as.ts(data$indicator_score)
  
  # Calculate change points
  
  cpt <- cpt.mean(indicator_ts, method = "PELT", 
                   penalty = "AIC", 
                   pen.value = c(1,25))
  
  cptvar <- cpt.var(indicator_ts, method = "PELT", 
                  penalty = "AIC", 
                  pen.value = c(1,25))
  
  # Save break points
  
  indicator_changepoint_summaries[[j]] <- cpt
  
  indicator_variance_summaries[[j]] <- cptvar
  
  }

names(indicator_changepoint_summaries) <- scenarios

scenario_changepoint_summaries[[i]] <- indicator_changepoint_summaries
scenario_variance_summaries[[i]] <- indicator_variance_summaries

}

names(scenario_changepoint_summaries) <- names(all_indicators_list)
# Harvested groups

i <- i + 1
plot(scenario_changepoint_summaries[["RLI large spp"]][[i]])
summary(scenario_changepoint_summaries[["RLI all"]][[i]])
i <- i + 1
plot(scenario_changepoint_summaries[["LPI"]][[i]])
summary(scenario_variance_summaries[[2]][[i]])

# LPI
plot(scenario_changepoint_summaries[[15]][[i]])
summary(scenario_changepoint_summaries[[15]][[i]])
plot(scenario_variance_summaries[[15]][[i]])
summary(scenario_variance_summaries[[15]][[i]])

# RLI

i <- i + 1
plot(scenario_changepoint_summaries[[23]][[i]])
summary(scenario_changepoint_summaries[[23]][[i]])
plot(scenario_variance_summaries[[23]][[i]])
summary(scenario_variance_summaries[[23]][[i]])

# Test strucchange ----
library(strucchange) # Seems to only find first breakpoint?

lpi_landuse <- indicators_all_df %>% 
  filter(indicator == "RLI 5y all spp") %>%
  filter(scenario == "100_Land_use") 

dat <- tibble(ylag0 = lpi_landuse$indicator_score,
              ylag1 = lag(lpi_landuse$indicator_score)
) %>%
  drop_na()

qlr <- Fstats(ylag0 ~ ylag1, data = dat)

breakpoints(ylag0 ~ ylag1, data = dat, breaks = 2)

plot.ts(lpi_landuse$indicator_score)

lpi_landuse_ts <- as.ts(lpi_landuse$indicator_score)

m_binseg <- cpt.mean(lpi_landuse_ts, penalty = "BIC", method = "PELT", Q = 5)
plot(m_binseg, type = "l", xlab = "Index", cpt.width = 4)

# Generalise additive models ----

## Following this tutorial https://fromthebottomoftheheap.net/2011/06/12/additive-modelling-and-the-hadcrut3v-global-mean-temperature-series/
## And this: https://www.youtube.com/watch?v=sgw4cu8hrZM&ab_channel=BottomoftheHeap
## (this might also be useful https://fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/)

# * HARVESTED ABUNDANCE ----

# ** Land use scenario ----

gam_inputs <- indicators_all_df %>% 
              filter(indicator == "total abundance harvested") %>% 
              filter(scenario == "100_Land_use") 

ylabel <- gam_inputs$indicator[1]
timesteps <- max(gam_inputs$annual_time_step)

row.names(gam_inputs) <- gam_inputs$annual_time_step
head(gam_inputs)

plot(indicator_score ~ annual_time_step, 
     data = gam_inputs, type = "o", ylab = ylabel)

# *** Fit a few different models ----

term <- 50
## Fit a smoother for Year to the data
m1 <- gamm(log(indicator_score+ 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs)
summary(m1$gam)

## look at autocorrelation in residuals:
acf(resid(m1$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m1$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best

## ...so fit the AR1
m2 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 1))
## ...and fit the AR2
m3 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 2))

# *** Model selection ----

anova(m1$lme, m2$lme, m3$lme)

# m2 looks best for our data too, lets have a look ...
selected_mod <- m2

gam.check(selected_mod$gam)

plot(selected_mod$gam, residuals = TRUE, pch = 19, cex = 0.75)

summary(selected_mod$gam)


# *** Plot the fitted trend ----

with(gam_inputs, tsDiagGamm(selected_mod, timevar = annual_time_step, 
                       observed = indicator_score))

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, 
     type = "p", ylab = ylabel)

pdat <- with(gam_inputs,
             data.frame(annual_time_step = seq(min(annual_time_step), 
                                               max(annual_time_step),
                                   length = 300)))

p1 <- predict(m1$gam, newdata = pdat)
p2 <- predict(selected_mod$gam, newdata = pdat)
lines(p1 ~ annual_time_step, data = pdat, col = "red")
lines(p2 ~ annual_time_step, data = pdat, col = "blue")
legend("topleft",
       legend = c("Uncorrelated Errors","AR(1) Errors"),
       bty = "n", col = c("red","blue"), lty = 1)

# *** Plot the derivatives and periods of change ----

selected_mod.d <- Deriv(selected_mod, n = timesteps)
plot(selected_mod.d, sizer = TRUE, alpha = 0.01)

# Add periods of change to time series

plot(log(indicator_score) ~ annual_time_step, 
     data = gam_inputs, type = "p", ylab = ylabel) 

lines(p2 ~ annual_time_step, data = pdat) 
CI <- confint(selected_mod.d, alpha = 0.01)
S <- signifD(p2, selected_mod.d$annual_time_step$deriv, 
             CI$annual_time_step$upper, 
             CI$annual_time_step$lower,
             eval = 0)
lines(S$incr ~ annual_time_step, data = pdat, lwd = 3, col = "blue")
lines(S$decr ~ annual_time_step, data = pdat, lwd = 3, col = "red")

# Recreate with ggplot

deriv_plot <- ggplot() +
  geom_point(data = gam_inputs,
             aes(x = annual_time_step, y = log(indicator_score)),
             alpha = 0.3) +
  geom_line(aes(x = gam_inputs$annual_time_step,
                y = p2)) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$incr), col = "blue", size = 2) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$decr), col = "red", size = 2) +
  labs(title = paste(gam_inputs$scenario[1],
                     gam_inputs$indicator[1], sep = " "))

deriv_plot

ggsave(file.path(analysis_plots_folder,
                 paste(gam_inputs$indicator[1],"_", gam_inputs$scenario[1], 
                       "derivative_plot.png", sep = "")),
       deriv_plot,  device = "png")

# *** Plot uncertainty ----

## simulate from posterior distribution of beta
Rbeta <- mvrnorm(n = 1000, coef(m2$gam), vcov(m2$gam))
Xp <- predict(m2$gam, newdata = pdat, type = "lpmatrix")
sim1 <- Xp %*% t(Rbeta)

# plot the observation and 25 of the 1000 trends
set.seed(321)
want <- sample(1000, 25)
# ylim <- range(sim1[,want], gam_inputs$annual_time_step)
ylim <- c(-1.2,1.2)
plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, 
     ylim = ylim, ylab = ylabel)
matlines(pdat$annual_time_step, sim1[,want], col = "black", lty = 1, pch = NA)

# ** Carnivore harvesting scenario ----

gam_inputs <- indicators_all_df %>% 
  filter(indicator == "total abundance harvested") %>% 
  filter(scenario == "200_Harvesting_carnivores") 

ylabel <- gam_inputs$indicator[1]
timesteps <- max(gam_inputs$annual_time_step)

row.names(gam_inputs) <- gam_inputs$annual_time_step
head(gam_inputs)

plot(indicator_score ~ annual_time_step, 
     data = gam_inputs, type = "o", ylab = ylabel)

# *** Fit a few different models ----

term <- 100
## Fit a smoother for Year to the data
m1 <- gamm(log(indicator_score+ 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs)
summary(m1$gam)

## look at autocorrelation in residuals:
acf(resid(m1$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m1$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best

## ...so fit the AR1
m2 <- gamm(indicator_score ~ s(annual_time_step, k = term), data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 1))
## ...and fit the AR2
m3 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 2))

# *** Model selection ----

anova(m1$lme, m2$lme, m3$lme)

# m2 looks best for our data too, lets have a look ...

plot(m3$gam, residuals = TRUE, pch = 19, cex = 0.75)

summary(m3$gam)


# *** Plot the fitted trend ----

with(gam_inputs, tsDiagGamm(m2, timevar = annual_time_step, 
                            observed = indicator_score))

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, 
     type = "p", ylab = ylabel)

pdat <- with(gam_inputs,
             data.frame(annual_time_step = seq(min(annual_time_step), 
                                               max(annual_time_step),
                                               length = 300)))

p1 <- predict(m1$gam, newdata = pdat)
p2 <- predict(m2$gam, newdata = pdat)
lines(p1 ~ annual_time_step, data = pdat, col = "red")
lines(p2 ~ annual_time_step, data = pdat, col = "blue")
legend("topleft",
       legend = c("Uncorrelated Errors","AR(1) Errors"),
       bty = "n", col = c("red","blue"), lty = 1)

# *** Plot the derivatives and periods of change ----

m2.d <- Deriv(m3, n = timesteps)
plot(m2.d, sizer = TRUE, alpha = 0.01)

# Add periods of change to time series

plot(log(indicator_score) ~ annual_time_step, 
     data = gam_inputs, type = "p", ylab = ylabel)

lines(p2 ~ annual_time_step, data = pdat)
CI <- confint(m2.d, alpha = 0.01)
S <- signifD(p2, m2.d$annual_time_step$deriv, 
             CI$annual_time_step$upper, 
             CI$annual_time_step$lower,
             eval = 0)
lines(S$incr ~ annual_time_step, data = pdat, lwd = 3, col = "blue")
lines(S$decr ~ annual_time_step, data = pdat, lwd = 3, col = "red")

# Recreate with ggplot

deriv_plot <- ggplot() +
  geom_point(data = gam_inputs,
             aes(x = annual_time_step, y = log(indicator_score)),
             alpha = 0.3) +
  geom_line(aes(x = gam_inputs$annual_time_step,
                y = p2)) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$incr), col = "blue", size = 2) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$decr), col = "red", size = 2) +
  labs(title = paste(gam_inputs$scenario[1],
                     gam_inputs$indicator[1], sep = " "))

deriv_plot

ggsave(file.path(analysis_plots_folder,
                 paste(gam_inputs$indicator[1],"_", gam_inputs$scenario[1], 
                       "derivative_plot.png", sep = "")),
       deriv_plot,  device = "png")


# *** Plot uncertainty ----

## simulate from posterior distribution of beta
Rbeta <- mvrnorm(n = 1000, coef(m2$gam), vcov(m2$gam))
Xp <- predict(m2$gam, newdata = pdat, type = "lpmatrix")
sim1 <- Xp %*% t(Rbeta)

# plot the observation and 25 of the 1000 trends
set.seed(321)
want <- sample(1000, 25)
# ylim <- range(sim1[,want], gam_inputs$annual_time_step)
ylim <- c(-1.2,1.2)
plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, ylim = ylim, ylab = ylabel)
matlines(pdat$annual_time_step, sim1[,want], col = "black", lty = 1, pch = NA)

# ** Herbivore harvesting scenario ----

gam_inputs <- indicators_all %>% 
  filter(indicator == "abundance harvested groups") %>% 
  filter(scenario == "300_Harvesting_herbivores") 

ylabel <- gam_inputs$indicator[1]
timesteps <- max(gam_inputs$annual_time_step)

row.names(gam_inputs) <- gam_inputs$annual_time_step
head(gam_inputs)

plot(indicator_score ~ annual_time_step, 
     data = gam_inputs, type = "o", ylab = ylabel)

# *** Fit a few different models ----

term <- 20
## Fit a smoother for Year to the data
m1 <- gamm(log(indicator_score+ 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs)
summary(m1$gam)

## look at autocorrelation in residuals:
acf(resid(m1$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m1$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best

## ...so fit the AR1
m2 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 1))
## ...and fit the AR2
m3 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 2))

# *** Model selection ----

anova(m1$lme, m2$lme, m3$lme)

# m2 best?

plot(m2$gam, residuals = TRUE, pch = 19, cex = 0.75)

summary(m2$gam)


# *** Plot the fitted trend ----

with(gam_inputs, tsDiagGamm(m2, timevar = annual_time_step, 
                            observed = indicator_score))

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, 
     type = "p", ylab = ylabel)

pdat <- with(gam_inputs,
             data.frame(annual_time_step = seq(min(annual_time_step), 
                                               max(annual_time_step),
                                               length = 300)))

p1 <- predict(m1$gam, newdata = pdat)
p2 <- predict(m2$gam, newdata = pdat)
lines(p1 ~ annual_time_step, data = pdat, col = "red")
lines(p2 ~ annual_time_step, data = pdat, col = "blue")
legend("topleft",
       legend = c("Uncorrelated Errors","AR(1) Errors"),
       bty = "n", col = c("red","blue"), lty = 1)

# *** Plot the derivatives and periods of change ----

m2.d <- Deriv(m2, n = 300)
plot(m2.d, sizer = TRUE, alpha = 0.01)

# Add periods of change to time series

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, 
     type = "p", ylab = ylabel)
lines(p2 ~ annual_time_step, data = pdat)
CI <- confint(m2.d, alpha = 0.01)
S <- signifD(p2, m2.d$annual_time_step$deriv, 
             CI$annual_time_step$upper, 
             CI$annual_time_step$lower,
             eval = 0)
lines(S$incr ~ annual_time_step, data = pdat, lwd = 3, col = "blue")
lines(S$decr ~ annual_time_step, data = pdat, lwd = 3, col = "red")

# Recreate with ggplot

deriv_plot <- ggplot() +
  geom_point(data = gam_inputs,
             aes(x = annual_time_step, y = log(indicator_score)),
             alpha = 0.3) +
  geom_line(aes(x = gam_inputs$annual_time_step,
                y = p2)) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$incr), col = "blue", size = 2) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$decr), col = "red", size = 2) +
  labs(title = paste(gam_inputs$scenario[1],
                     gam_inputs$indicator[1], sep = " "))

deriv_plot

ggsave(file.path(analysis_plots_folder,
                 paste(gam_inputs$indicator[1],"_", gam_inputs$scenario[1], 
                       "derivative_plot.png", sep = "")),
       deriv_plot,  device = "png")

# *** Plot uncertainty ----

## simulate from posterior distribution of beta
Rbeta <- mvrnorm(n = 1000, coef(m2$gam), vcov(m2$gam))
Xp <- predict(m2$gam, newdata = pdat, type = "lpmatrix")
sim1 <- Xp %*% t(Rbeta)

# plot the observation and 25 of the 1000 trends
set.seed(321)
want <- sample(1000, 25)
# ylim <- range(sim1[,want], gam_inputs$annual_time_step)
ylim <- c(-1.2,1.2)
plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, ylim = ylim, ylab = ylabel)
matlines(pdat$annual_time_step, sim1[,want], col = "black", lty = 1, pch = NA)

# * RED LIST INDEX ----

# ** Land use scenario ----

unique(indicators_all_df$indicator)

response <- indicators_all_df %>% 
  filter(indicator == "LPI") %>% 
  filter(scenario == "100_Land_use") %>% 
  dplyr::select(annual_time_step, indicator_score) %>% 
  rename(LPI = indicator_score)

predictor <- indicators_all_df %>% 
  filter(indicator == "total abundance harvested") %>% 
  filter(scenario == "100_Land_use") %>% 
  dplyr::select(annual_time_step, indicator_score) %>% 
  rename(harvested = indicator_score)

gam_inputs <- response %>% 
              merge(predictor, by = "annual_time_step") %>% 
              mutate(LPI_scaled = scale(LPI),
                     harvested_scaled = scale(harvested))

head(gam_inputs)
summary(gam_inputs)

# Looking for trends ----
gam_inputs <- indicators_all_df %>% 
  filter(indicator == "RLI 5y all spp") %>% 
  filter(scenario == "100_Land_use")

ylabel <- gam_inputs$indicator[1]
timesteps <- max(gam_inputs$annual_time_step)

row.names(gam_inputs) <- gam_inputs$annual_time_step
head(gam_inputs)

plot(indicator_score ~ annual_time_step, 
     data = gam_inputs, type = "o", ylab = ylabel)

# *** Fit a few different models ----

## Have a look at the distribution (positively skewed)
hist(gam_inputs$indicator_score, breaks = 20)
hist(log(gam_inputs$indicator_score), breaks = 20)
range(gam_inputs$indicator_score)

term <- 30

## Test some different options for family and temporal autocorrelation

scatmod <- gamm(indicator_score ~ s(annual_time_step, k = term), 
           data = gam_inputs, family = "scat")

summary(scatmod$gam)
gam.check(scatmod$gam)

with(gam_inputs, tsDiagGamm(scatmod, timevar = annual_time_step, 
                            observed = indicator_score))

betamod <- gamm(indicator_score ~ s(annual_time_step, k = term), 
           data = gam_inputs, family = "betar")

summary(betamod$gam)
gam.check(betamod$gam)

with(gam_inputs, tsDiagGamm(betamod, timevar = annual_time_step, 
                            observed = indicator_score))

#####
logmod <- gamm(log(indicator_score) ~ s(annual_time_step, k = term), 
                data = gam_inputs)

summary(logmod$gam)
gam.check(logmod$gam)

with(gam_inputs, tsDiagGamm(logmod, timevar = annual_time_step, 
                            observed = indicator_score))

sqrtmod <- gamm(sqrt(indicator_score) ~ s(annual_time_step, k = term), 
               data = gam_inputs)

summary(sqrtmod$gam)
gam.check(sqrtmod$gam)

with(gam_inputs, tsDiagGamm(logmod, timevar = annual_time_step, 
                            observed = indicator_score))

binommod <- gamm(indicator_score ~ s(annual_time_step, k = term), 
                data = gam_inputs, family = "gaulss")

summary(binommod$gam)
gam.check(binommod$gam)

with(gam_inputs, tsDiagGamm(logmod, timevar = annual_time_step, 
                            observed = indicator_score))
#####

bm <- gam(indicator_score ~s(annual_time_step, k = term),
          family=betar(link="cauchit"),
          data= gam_inputs)

bm
plot(bm,pages=1)
summary(bm)
gam.check(bm)

bmm <- gamm(indicator_score ~s(annual_time_step, k = 40, bs = 'cs'),
          family="betar",
          data= gam_inputs)

bmm
plot(bmm$gam,pages=1)
summary(bmm$gam)
gam.check(bmm$gam)

## ...so fit the AR1
betamodar <- gamm(indicator_score ~ s(annual_time_step,k = 50, bs = 'cs'), 
                  data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 2),
           family= "betar")

betamodar
plot(betamodar$gam,pages=1)
summary(betamodar$gam)
gam.check(betamodar$gam)

tweediemod <- gamm(indicator_score ~ s(annual_time_step, k = 50, bs = 'ts'), 
                   data = gam_inputs,
                  correlation = corARMA(form = ~ annual_time_step, p = 2),
                  family= "tw")

tweediemod
plot(tweediemod$gam,pages=1)
gam.check(tweediemod$gam)
summary(tweediemod$gam)

nullmod <- gamm(indicator_score ~ s(1), family = betar, data = gam_inputs)

# *** Model selection ----

## Using ANOVA
anova(m1$lme, m2$lme, m3$lme)
anova(bm, bmm$lme, m1$lme)

## Using AIC
mod_sel_table <- as.data.frame(AIC(tweediemod$lme, bmm$lme, betamodar$lme))
mod_sel_table

# Define selected model

selected_mod <- tweediemod

plot(selected_mod$gam, residuals = TRUE, pch = 19, cex = 0.75)

summary(selected_mod$gam)
gam.check(selected_mod$gam)

# *** Plot the fitted trend ----

## Only works for GAMM
# with(gam_inputs, tsDiagGamm(m2, timevar = annual_time_step, 
#                             observed = indicator_score))

plot(indicator_score ~ annual_time_step, data = gam_inputs, 
     type = "p", ylab = ylabel)

pdat <- with(gam_inputs,
             data.frame(annual_time_step = seq(5, 295, 5)))

# p1 <- predict(m1$gam, newdata = pdat)
# p2 <- predict(m2$gam, newdata = pdat)
p <- predict(selected_mod$gam, newdata = pdat, type = "response")
p
lines(p ~ annual_time_step, data = pdat, col = "pink")
# lines(p1 ~ annual_time_step, data = pdat, col = "red")
# lines(p2 ~ annual_time_step, data = pdat, col = "blue")
legend("topleft",
       legend = c("Uncorrelated Errors","AR(1) Errors"),
       bty = "n", col = c("red","blue"), lty = 1)

# *** Plot the derivatives and periods of change ----

m2.d <- Deriv(selected_mod, n = 59)
plot(m2.d, sizer = TRUE, alpha = 0.01)

# Add periods of change to time series

plot(indicator_score ~ annual_time_step, data = gam_inputs, 
     type = "p", ylab = ylabel)
lines(p ~ annual_time_step, data = pdat)
CI <- confint(m2.d, alpha = 0.01)
S <- signifD(p, m2.d$annual_time_step$deriv, 
             CI$annual_time_step$upper, 
             CI$annual_time_step$lower,
             eval = 0)
lines(S$incr ~ annual_time_step, data = pdat, lwd = 3, col = "blue")
lines(S$decr ~ annual_time_step, data = pdat, lwd = 3, col = "red")
# Recreate with ggplot

deriv_plot <- ggplot() +
  geom_point(data = gam_inputs,
             aes(x = annual_time_step, y = indicator_score),
             alpha = 0.3) +
  geom_line(aes(x = gam_inputs$annual_time_step,
                y = p2)) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$incr), col = "blue", size = 2) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$decr), col = "red", size = 2) +
  labs(title = paste(gam_inputs$scenario[1],
                     gam_inputs$indicator[1], sep = " "))

deriv_plot

ggsave(file.path(analysis_plots_folder,
                 paste(gam_inputs$indicator[1],"_", gam_inputs$scenario[1], 
                       "derivative_plot.png", sep = "")),
       deriv_plot,  device = "png")

# *** Plot uncertainty ----

## simulate from posterior distribution of beta
Rbeta <- mvrnorm(n = 1000, coef(m2$gam), vcov(m2$gam))
Xp <- predict(m2$gam, newdata = pdat, type = "lpmatrix")
sim1 <- Xp %*% t(Rbeta)

# plot the observation and 25 of the 1000 trends
set.seed(321)
want <- sample(1000, 25)
# ylim <- range(sim1[,want], gam_inputs$annual_time_step)
ylim <- c(-1.2,1.2)
plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, ylim = ylim, ylab = ylabel)
matlines(pdat$annual_time_step, sim1[,want], col = "black", lty = 1, pch = NA)

# ** Carnivore harvesting scenario ----

gam_inputs <- indicators_all %>% 
  filter(indicator == "abundance harvested groups") %>% 
  filter(scenario == "200_Harvesting_carnivores") 

ylabel <- gam_inputs$indicator[1]
timesteps <- max(gam_inputs$annual_time_step)

row.names(gam_inputs) <- gam_inputs$annual_time_step
head(gam_inputs)

plot(indicator_score ~ annual_time_step, 
     data = gam_inputs, type = "o", ylab = ylabel)

# *** Fit a few different models ----

term <- 20
## Fit a smoother for Year to the data
m1 <- gamm(log(indicator_score+ 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs)
summary(m1$gam)

## look at autocorrelation in residuals:
acf(resid(m1$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m1$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best

## ...so fit the AR1
m2 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), 
           data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 1))
## ...and fit the AR2
m3 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term),
           data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 2))

# *** Model selection ----

anova(m1$lme, m2$lme, m3$lme)

# m2 looks best for our data too, lets have a look ...
m2 <- m3

plot(m2$gam, residuals = TRUE, pch = 19, cex = 0.75)

summary(m2$gam)

# *** Plot the fitted trend ----

with(gam_inputs, tsDiagGamm(m2, timevar = annual_time_step, 
                            observed = indicator_score))

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, 
     type = "p", ylab = ylabel)

pdat <- with(gam_inputs,
             data.frame(annual_time_step = seq(min(annual_time_step), 
                                               max(annual_time_step),
                                               length = 300)))

p1 <- predict(m1$gam, newdata = pdat)
p2 <- predict(m2$gam, newdata = pdat)
lines(p1 ~ annual_time_step, data = pdat, col = "red")
lines(p2 ~ annual_time_step, data = pdat, col = "blue")
legend("topleft",
       legend = c("Uncorrelated Errors","AR(1) Errors"),
       bty = "n", col = c("red","blue"), lty = 1)

# *** Plot the derivatives and periods of change ----

m2.d <- Deriv(m2, n = 300)
plot(m2.d, sizer = TRUE, alpha = 0.01)

# Add periods of change to time series

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, type = "p", ylab = ylabel)
lines(p2 ~ annual_time_step, data = pdat)
CI <- confint(m2.d, alpha = 0.01)
S <- signifD(p2, m2.d$annual_time_step$deriv, 
             CI$annual_time_step$upper, 
             CI$annual_time_step$lower,
             eval = 0)
lines(S$incr ~ annual_time_step, data = pdat, lwd = 3, col = "blue")
lines(S$decr ~ annual_time_step, data = pdat, lwd = 3, col = "red")

# Recreate with ggplot

deriv_plot <- ggplot() +
  geom_point(data = gam_inputs,
             aes(x = annual_time_step, y = log(indicator_score)),
             alpha = 0.3) +
  geom_line(aes(x = gam_inputs$annual_time_step,
                y = p2)) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$incr), col = "blue", size = 2) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$decr), col = "red", size = 2) +
  labs(title = paste(gam_inputs$scenario[1],
                     gam_inputs$indicator[1], sep = " "))

deriv_plot

ggsave(file.path(analysis_plots_folder,
                 paste(gam_inputs$indicator[1],"_", gam_inputs$scenario[1], 
                       "derivative_plot.png", sep = "")),
       deriv_plot,  device = "png")

# *** Plot uncertainty ----

## simulate from posterior distribution of beta
Rbeta <- mvrnorm(n = 1000, coef(m2$gam), vcov(m2$gam))
Xp <- predict(m2$gam, newdata = pdat, type = "lpmatrix")
sim1 <- Xp %*% t(Rbeta)

# plot the observation and 25 of the 1000 trends
set.seed(321)
want <- sample(1000, 25)
# ylim <- range(sim1[,want], gam_inputs$annual_time_step)
ylim <- c(-1.2,1.2)
plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, ylim = ylim, ylab = ylabel)
matlines(pdat$annual_time_step, sim1[,want], col = "black", lty = 1, pch = NA)

# ** Herbivore harvesting scenario ----

gam_inputs <- indicators_all %>% 
  filter(indicator == "abundance harvested groups") %>% 
  filter(scenario == "300_Harvesting_herbivores") 

ylabel <- gam_inputs$indicator[1]
timesteps <- max(gam_inputs$annual_time_step)

row.names(gam_inputs) <- gam_inputs$annual_time_step
head(gam_inputs)

plot(indicator_score ~ annual_time_step, 
     data = gam_inputs, type = "o", ylab = ylabel)

# *** Fit a few different models ----

term <- 20
## Fit a smoother for Year to the data
m1 <- gamm(log(indicator_score+ 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs)
summary(m1$gam)

## look at autocorrelation in residuals:
acf(resid(m1$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m1$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best

## ...so fit the AR1
m2 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 1))
## ...and fit the AR2
m3 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 2))

# *** Model selection ----

anova(m1$lme, m2$lme, m3$lme)

# m2 looks best for our data too, lets have a look ...
m2 <- m3

plot(m2$gam, residuals = TRUE, pch = 19, cex = 0.75)

summary(m2$gam)


# *** Plot the fitted trend ----

with(gam_inputs, tsDiagGamm(m2, timevar = annual_time_step, 
                            observed = indicator_score))

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, 
     type = "p", ylab = ylabel)

pdat <- with(gam_inputs,
             data.frame(annual_time_step = seq(min(annual_time_step), 
                                               max(annual_time_step),
                                               length = 300)))

p1 <- predict(m1$gam, newdata = pdat)
p2 <- predict(m2$gam, newdata = pdat)
lines(p1 ~ annual_time_step, data = pdat, col = "red")
lines(p2 ~ annual_time_step, data = pdat, col = "blue")
legend("topleft",
       legend = c("Uncorrelated Errors","AR(1) Errors"),
       bty = "n", col = c("red","blue"), lty = 1)

# *** Plot the derivatives and periods of change ----

m2.d <- Deriv(m2, n = 300)
plot(m2.d, sizer = TRUE, alpha = 0.01)

# Add periods of change to time series

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, type = "p", ylab = ylabel)
lines(p2 ~ annual_time_step, data = pdat)
CI <- confint(m2.d, alpha = 0.01)
S <- signifD(p2, m2.d$annual_time_step$deriv, 
             CI$annual_time_step$upper, 
             CI$annual_time_step$lower,
             eval = 0)
lines(S$incr ~ annual_time_step, data = pdat, lwd = 3, col = "blue")
lines(S$decr ~ annual_time_step, data = pdat, lwd = 3, col = "red")

# Recreate with ggplot

deriv_plot <- ggplot() +
  geom_point(data = gam_inputs,
             aes(x = annual_time_step, y = log(indicator_score)),
             alpha = 0.3) +
  geom_line(aes(x = gam_inputs$annual_time_step,
                y = p2)) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$incr), col = "blue", size = 2) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$decr), col = "red", size = 2) +
  labs(title = paste(gam_inputs$scenario[1],
                     gam_inputs$indicator[1], sep = " "))

deriv_plot

ggsave(file.path(analysis_plots_folder,
                 paste(gam_inputs$indicator[1],"_", gam_inputs$scenario[1], 
                       "derivative_plot.png", sep = "")),
       deriv_plot,  device = "png")

# *** Plot uncertainty ----

## simulate from posterior distribution of beta
Rbeta <- mvrnorm(n = 1000, coef(m2$gam), vcov(m2$gam))
Xp <- predict(m2$gam, newdata = pdat, type = "lpmatrix")
sim1 <- Xp %*% t(Rbeta)

# plot the observation and 25 of the 1000 trends
set.seed(321)
want <- sample(1000, 25)
# ylim <- range(sim1[,want], gam_inputs$annual_time_step)
ylim <- c(-1.2,1.2)
plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, ylim = ylim, ylab = ylabel)
matlines(pdat$annual_time_step, sim1[,want], col = "black", lty = 1, pch = NA)

# * LIVING PLANET INDEX ----

# ** Land use scenario ----

gam_inputs <- indicators_all_df %>% 
  filter(indicator == "LPI") %>% 
  filter(scenario == "100_Land_use") 

ylabel <- gam_inputs$indicator[1]
timesteps <- max(gam_inputs$annual_time_step)

row.names(gam_inputs) <- gam_inputs$annual_time_step
head(gam_inputs)

plot(indicator_score ~ annual_time_step, 
     data = gam_inputs, type = "o", ylab = ylabel)

# *** Fit a few different models ----

term <- 100
## Fit a smoother for Year to the data
m1 <- gamm(log(indicator_score+ 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs)
summary(m1$gam)

## look at autocorrelation in residuals:
acf(resid(m1$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m1$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best

## ...so fit the AR1
m2 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), 
           data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 1))
## ...and fit the AR2
m3 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), 
           data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 2))

# *** Model selection ----

anova(m1$lme, m2$lme)

# m2 looks best for our data too, lets have a look ...
selected_mod <- m2

plot(selected_mod$gam, residuals = TRUE, pch = 19, cex = 0.75)

gam.check(selected_mod$gam)

summary(selected_mod$gam)

# *** Plot the fitted trend ----

with(gam_inputs, tsDiagGamm(selected_mod, timevar = annual_time_step, 
                            observed = indicator_score))

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, 
     type = "p", ylab = ylabel)

pdat <- with(gam_inputs,
             data.frame(annual_time_step = seq(min(annual_time_step), 
                                               max(annual_time_step),
                                               length = 300)))

p1 <- predict(m1$gam, newdata = pdat)
p2 <- predict(selected_mod$gam, newdata = pdat)
lines(p1 ~ annual_time_step, data = pdat, col = "red")
lines(p2 ~ annual_time_step, data = pdat, col = "blue")
legend("topleft",
       legend = c("Uncorrelated Errors","AR(1) Errors"),
       bty = "n", col = c("red","blue"), lty = 1)

# *** Plot the derivatives and periods of change ----

selected_mod.d <- Deriv(selected_mod, n = 300)
plot(selected_mod.d, sizer = TRUE, alpha = 0.01)

# Add periods of change to time series

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, type = "p", ylab = ylabel)
lines(p2 ~ annual_time_step, data = pdat)
CI <- confint(selected_mod.d, alpha = 0.01)
S <- signifD(p2, selected_mod.d$annual_time_step$deriv, 
             CI$annual_time_step$upper, 
             CI$annual_time_step$lower,
             eval = 0)
lines(S$incr ~ annual_time_step, data = pdat, lwd = 3, col = "blue")
lines(S$decr ~ annual_time_step, data = pdat, lwd = 3, col = "red")

# Recreate with ggplot

deriv_plot <- ggplot() +
  geom_point(data = gam_inputs,
             aes(x = annual_time_step, y = log(indicator_score)),
             alpha = 0.3) +
  geom_line(aes(x = gam_inputs$annual_time_step,
                y = p2)) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$incr), col = "blue", size = 2) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$decr), col = "red", size = 2) +
  labs(title = paste(gam_inputs$scenario[1],
                     gam_inputs$indicator[1], sep = " "))

deriv_plot

ggsave(file.path(analysis_plots_folder,
                 paste(gam_inputs$indicator[1],"_", gam_inputs$scenario[1], 
                       "derivative_plot.png", sep = "")),
       deriv_plot,  device = "png")

# *** Plot uncertainty ----

## simulate from posterior distribution of beta
Rbeta <- mvrnorm(n = 1000, coef(selected_mod$gam), vcov(selected_mod$gam))
Xp <- predict(selected_mod$gam, newdata = pdat, type = "lpmatrix")
sim1 <- Xp %*% t(Rbeta)

# plot the observation and 25 of the 1000 trends
set.seed(321)
want <- sample(1000, 25)
# ylim <- range(sim1[,want], gam_inputs$annual_time_step)
ylim <- c(-2,0.5)
plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, 
     ylim = ylim, ylab = ylabel)
matlines(pdat$annual_time_step, sim1[,want], col = "black", lty = 1, pch = NA)

# ** Carnivore harvesting scenario ----
unique(indicators_all_df$indicator)

gam_inputs <- indicators_all_df %>% 
  filter(indicator == "LPI") %>% 
  filter(scenario == "200_Harvesting_carnivores") 

ylabel <- gam_inputs$indicator[1]
timesteps <- max(gam_inputs$annual_time_step)

row.names(gam_inputs) <- gam_inputs$annual_time_step
head(gam_inputs)

plot(indicator_score ~ annual_time_step, 
     data = gam_inputs, type = "o", ylab = ylabel)

# *** Fit a few different models ----

term <- 250
## Fit a smoother for Year to the data
m1 <- gamm(log(indicator_score+ 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs)
summary(m1$gam)

## look at autocorrelation in residuals:
acf(resid(m1$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m1$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best

## ...so fit the AR1
m2 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 1))
## ...and fit the AR2
m3 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 2))

m3 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 4))

# *** Model selection ----

anova(m1$lme, m2$lme, m3$lme)

# m2 looks best for our data too, lets have a look ...
m2 <- m2

plot(m2$gam, residuals = TRUE, pch = 19, cex = 0.75)

summary(m2$gam)

gam.check(m2$gam)

# *** Plot the fitted trend ----

with(gam_inputs, tsDiagGamm(m2, 
                            timevar = annual_time_step, 
                            observed = indicator_score))

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, 
     type = "p", ylab = ylabel)

pdat <- with(gam_inputs,
             data.frame(annual_time_step = seq(min(annual_time_step), 
                                               max(annual_time_step),
                                               length = 300)))

p1 <- predict(m1$gam, newdata = pdat)
p2 <- predict(m2$gam, newdata = pdat)
lines(p1 ~ annual_time_step, data = pdat, col = "red")
lines(p2 ~ annual_time_step, data = pdat, col = "blue")
legend("topleft",
       legend = c("Uncorrelated Errors","AR(1) Errors"),
       bty = "n", col = c("red","blue"), lty = 1)

# *** Plot the derivatives and periods of change ----

m2.d <- Deriv(m2, n = 300)
plot(m2.d, sizer = TRUE, alpha = 0.01)

# Add periods of change to time series

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, type = "p", ylab = ylabel)
lines(p2 ~ annual_time_step, data = pdat)
CI <- confint(m2.d, alpha = 0.01)
S <- signifD(p2, m2.d$annual_time_step$deriv, 
             CI$annual_time_step$upper, 
             CI$annual_time_step$lower,
             eval = 0)
lines(S$incr ~ annual_time_step, data = pdat, lwd = 3, col = "blue")
lines(S$decr ~ annual_time_step, data = pdat, lwd = 3, col = "red")

# Recreate with ggplot

deriv_plot <- ggplot() +
  geom_point(data = gam_inputs,
             aes(x = annual_time_step, y = log(indicator_score)),
             alpha = 0.3) +
  geom_line(aes(x = gam_inputs$annual_time_step,
                y = p2)) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$incr), col = "blue", size = 2) +
  geom_line(aes(x = gam_inputs$annual_time_step, 
                y = S$decr), col = "red", size = 2) +
  labs(title = paste(gam_inputs$scenario[1],
                     gam_inputs$indicator[1], sep = " "))

deriv_plot

ggsave(file.path(analysis_plots_folder,
                 paste(gam_inputs$indicator[1],"_", gam_inputs$scenario[1], 
                       "derivative_plot.png", sep = "")),
       deriv_plot,  device = "png")

# *** Plot uncertainty ----

## simulate from posterior distribution of beta
Rbeta <- mvrnorm(n = 1000, coef(m2$gam), vcov(m2$gam))
Xp <- predict(m2$gam, newdata = pdat, type = "lpmatrix")
sim1 <- Xp %*% t(Rbeta)

# plot the observation and 25 of the 1000 trends
set.seed(321)
want <- sample(1000, 25)
# ylim <- range(sim1[,want], gam_inputs$annual_time_step)
ylim <- c(-4,4)
plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, ylim = ylim, ylab = ylabel)
matlines(pdat$annual_time_step, sim1[,want], col = "black", lty = 1, pch = NA)

# ** Herbivore harvesting scenario ----

gam_inputs <- indicators_all_df %>% 
  filter(indicator == "LPI") %>% 
  filter(scenario == "300_Harvesting_herbivores") 

ylabel <- gam_inputs$indicator[1]
timesteps <- max(gam_inputs$annual_time_step)

row.names(gam_inputs) <- gam_inputs$annual_time_step
head(gam_inputs)

plot(indicator_score ~ annual_time_step, 
     data = gam_inputs, type = "o", ylab = ylabel)

# *** Fit a few different models ----

term <- 300
## Fit a smoother for Year to the data
m1 <- gamm(log(indicator_score+ 0.0001) ~ s(annual_time_step, k = term), data = gam_inputs)
summary(m1$gam)

## look at autocorrelation in residuals:
acf(resid(m1$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m1$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best

## ...so fit the AR1
m2 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), 
           data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 1))
## ...and fit the AR2
m3 <- gamm(log(indicator_score + 0.0001) ~ s(annual_time_step, k = term), 
           data = gam_inputs,
           correlation = corARMA(form = ~ annual_time_step, p = 2))

# *** Model selection ----

anova(m1$lme, m2$lme, m3$lme)

# m2 looks best for our data too, lets have a look ...

m2 <- m3

plot(m2$gam, residuals = TRUE, pch = 19, cex = 0.75)

summary(m2$gam)
gam.check(m2$gam)

# *** Plot the fitted trend ----

with(gam_inputs, tsDiagGamm(m2, timevar = annual_time_step, 
                            observed = indicator_score))

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, 
     type = "p", ylab = ylabel)

pdat <- with(gam_inputs,
             data.frame(annual_time_step = seq(min(annual_time_step), 
                                               max(annual_time_step),
                                               length = 300)))

p1 <- predict(m1$gam, newdata = pdat)
p2 <- predict(m2$gam, newdata = pdat)
lines(p1 ~ annual_time_step, data = pdat, col = "red")
lines(p2 ~ annual_time_step, data = pdat, col = "blue")
legend("topleft",
       legend = c("Uncorrelated Errors","AR(1) Errors"),
       bty = "n", col = c("red","blue"), lty = 1)

# *** Plot the derivatives and periods of change ----

m2.d <- Deriv(m2, n = 300)
plot(m2.d, sizer = TRUE, alpha = 0.01)

# Add periods of change to time series

plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, type = "p", ylab = ylabel)
lines(p2 ~ annual_time_step, data = pdat)
CI <- confint(m2.d, alpha = 0.01)
S <- signifD(p2, m2.d$annual_time_step$deriv, 
             CI$annual_time_step$upper, 
             CI$annual_time_step$lower,
             eval = 0)
lines(S$incr ~ annual_time_step, data = pdat, lwd = 3, col = "blue")
lines(S$decr ~ annual_time_step, data = pdat, lwd = 3, col = "red")

# Recreate with ggplot

deriv_plot <- ggplot() +
              geom_point(data = gam_inputs,
                         aes(x = annual_time_step, y = log(indicator_score)),
                             alpha = 0.3) +
              geom_line(aes(x = gam_inputs$annual_time_step,
                            y = p2)) +
              geom_line(aes(x = gam_inputs$annual_time_step, 
                               y = S$incr), col = "blue", size = 2) +
              geom_line(aes(x = gam_inputs$annual_time_step, 
                            y = S$decr), col = "red", size = 2) +
              labs(title = paste(gam_inputs$scenario[1],
                                 gam_inputs$indicator[1], sep = " "))

deriv_plot

ggsave(file.path(analysis_plots_folder,
                 paste(gam_inputs$indicator[1],"_", gam_inputs$scenario[1], 
                       "derivative_plot.png", sep = "")),
       deriv_plot,  device = "png")

# *** Plot uncertainty ----

## simulate from posterior distribution of beta
Rbeta <- mvrnorm(n = 1000, coef(m2$gam), vcov(m2$gam))
Xp <- predict(m2$gam, newdata = pdat, type = "lpmatrix")
sim1 <- Xp %*% t(Rbeta)

# plot the observation and 25 of the 1000 trends
set.seed(321)
want <- sample(1000, 25)
# ylim <- range(sim1[,want], gam_inputs$annual_time_step)
ylim <- c(-1.2,1.2)
plot(log(indicator_score) ~ annual_time_step, data = gam_inputs, ylim = ylim, ylab = ylabel)
matlines(pdat$annual_time_step, sim1[,want], col = "black", lty = 1, pch = NA)




