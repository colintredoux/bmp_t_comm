# Header------------------------------------------------------------------------
'########################################################
## This script creates some stat summaries that       ##
## explore the data from the BMP project;             ##
## specifically some patterns over T communities      ##
##  CG Tredoux April 2018                             ##
########################################################'

# Load libraries----------------------------------------------------------------
library(pacman)
p_load(tidyverse, magrittr, gmodels, lme4, lmerTest, stringr, 
       DHARMa, pbkrtest)

# Load data---------------------------------------------------------------------
bmp_tcom.dat <- readRDS("data/bmp_tcom.Rds")

# Explore data------------------------------------------------------------------
# How many data points per participant?  Histogram that
temp <- bmp_tcom.dat %$%
  as.data.frame(table(trackingid))
ggplot(temp, aes(x = Freq)) +
  geom_histogram(color = "red", fill = "white") +
  scale_x_log10()
# That seems a little widely dispersed, but what about the distribution
# of the total time spent per participant i.e. need sum first

bmp_tcom.dat %>%
  group_by(trackingid) %>%
  summarise(sum_mins = sum(mins, na.rm = T)) %>%
  ggplot(aes(x = sum_mins)) +
  geom_histogram(color = "red", fill = "white") +
  scale_x_log10()
# OK to that, when logged
# But one must ask what does the distribution look like over
# T Communities or per person


bmp_tcom.dat %>%
  group_by(trackingid, tcommunity) %>%
  summarise(sum_mins = sum(mins, na.rm = T)) %>%
  ggplot(aes(x = sum_mins)) +
  geom_histogram(color = "red", fill = "white") +
  scale_x_log10() +
  facet_wrap( ~ tcommunity)
# OK to that, when logged, but some small groupings
# EDIT: seems to have improved with Gemma's new data file

posdodge = position_dodge(0.9)
bmp_tcom.dat %>%
  filter(community != "Other") %>%
  group_by(trackingid, belonging, tcommunity, tmain, community) %>%
  summarise(
    sum_mins = sum(mins, na.rm = T),
    sd_mins  = sd(mins, na.rm = T),
    n_mins   = n()
  ) %>%
  group_by(belonging, community) %>%
  filter(community != "Other") %>%
  summarise(
    mean_sum_mins = mean(sum_mins, na.rm = T),
    sd_sum_mins   = sd(sum_mins, na.rm = T),
    n_sum_mins    = n(),
    ci            = 1.96 * (sd_sum_mins / sqrt(n_sum_mins))
  ) %>%
  ggplot(aes(x = belonging, y = mean_sum_mins, fill = community)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(
    aes(ymin = mean_sum_mins - ci, ymax = mean_sum_mins + ci),
    width = 0.2,
    position = posdodge
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.85, 0.75))
# graph shows clear pattern but quite a lot of error
# co-ords flipped below for different view

bmp_tcom.dat %>%
  filter(community != "Other") %>%
  group_by(trackingid, belonging, tcommunity, tmain, community) %>%
  summarise(
    sum_mins = sum(mins, na.rm = T),
    sd_mins  = sd(mins, na.rm = T),
    n_mins   = n()
  ) %>%
  group_by(tcommunity, community) %>%
  filter(community != "Other") %>%
  summarise(
    mean_sum_mins = mean(sum_mins, na.rm = T),
    sd_sum_mins   = sd(sum_mins, na.rm = T),
    n_sum_mins    = n(),
    ci            = 1.96 * sd_sum_mins / sqrt(n_sum_mins)
  ) %>%
  ggplot(aes(x = tcommunity, y = mean_sum_mins, fill = community)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(
    aes(ymin = mean_sum_mins - ci, ymax = mean_sum_mins + ci),
    width = 0.2,
    position = posdodge
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.85, 0.75)) +
  coord_flip()


# Model ------------------------------------------------------------------------
# First reduce file so that each individual has sum of time spent
# in areas, logged

bmp_tcom_sum.dat <-
  bmp_tcom.dat %>%
  mutate(mode = factor(mode)) %>%
  group_by(trackingid,
           mode,
           community,
           gender,
           belonging,
           tcommunity,
           tmain) %>%
  summarise(log_sum_mins = log(sum(mins, na.rm = T) + 1))

# Old code----------------------------------------------------------------------
# Ignore this section
lme_mod_0 <-
  bmp_tcom.dat %>%
  group_by(trackingid, belonging) %>%
  summarise(log_sum_mins = log(sum(mins, na.rm = T) + 1)) %>%
  spread(key = belonging, value = log_sum_mins) %>%
  replace_na(list(
    HomeT = 0,
    ingroup = 0,
    neutral = 0,
    `out group` = 0
  )) %>%
  gather(key = belonging, value = log_sum_mins, HomeT:`out group`) %>%
  lmer(log_sum_mins ~ belonging + (1 | trackingid),
       data = .,
       REML = FALSE)
summary(lme_mod_0)

resid1_fm0 <-
  HLMresid(lme_mod_0,
           level = 1,
           type = "LS",
           standardize = TRUE)
head(resid1_fm0)



# model time with predictors----------------------------------------------------
# Base model
lme_mod_0 <-
  lmer(log_sum_mins ~ (1 | trackingid),
       data =  bmp_tcom_sum.dat,
       REML = FALSE)
summary(lme_mod_0)

lme_mod_1 <- lmer(log_sum_mins ~ belonging +
                    (1 | trackingid),
                  data = bmp_tcom_sum.dat,
                  REML = FALSE)
summary(lme_mod_1)

# compare against null(base) model
anova(lme_mod_1, lme_mod_0)

# confints
confint(lme_mod_1, method = "boot", nsim = 100)

# residuals
qqnorm(resid(lme_mod_0))
plot(lme_mod_0)
# some heteroscedasticity in that

# explore some summ stats to amplify that
bmp_tcom_sum.dat %>%
  group_by(mode) %>%
  summarise(mean(log_sum_mins, na.rm = T))

bmp_tcom_sum.dat %>%
  group_by(belonging) %>%
  summarise(mean(log_sum_mins, na.rm = T))


# Conclude from all that that residuals don't look acceptable

# We need to re-think this, let's see whether removing zero counts improves
# things






# MODEL visit data--------------------------------------------------------------
# We can try to treat the data as count data - number of visits to an area

# Let's reduce the data to counts per person per region

x <- bmp_tcom.dat %>%
  group_by(trackingid, belonging, gender, mode) %>%
  summarise(visit_count = n()) %>%
  spread(key = belonging, value = visit_count, sep = "_") %>%
  replace_na(
    list(
      belonging_HomeT = 0,
      belonging_ingroup = 0,
      belonging_neutral = 0,
      `belonging_out group` = 0
    )
  ) %>%
  gather(key = belonging,
         value = visit_count,
         belonging_HomeT:`belonging_out group`)
y <- x$belonging
y <- str_replace(y, "belonging_", "")
x$belonging <- y
bmp_tcom_count.dat <- x
z <- 1:length(bmp_tcom_count.dat$trackingid)
bmp_tcom_count.dat <- x
bmp_tcom_count.dat$z <- z
rm(x, y, z)

# model the main effects, declaring just tracking id as random effect
glme_mod_1 <-
  glmer(
    visit_count ~ as.factor(gender) + as.factor(mode) +
      as.factor(belonging) +
      (1 | trackingid),
    family = "poisson",
    data = bmp_tcom_count.dat
  )
summary(glme_mod_1)

# Check residuals on model glme_mod_1 using DHARMa
simulationOutput <-
  simulateResiduals(fittedModel = glme_mod_1, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput)
plotResiduals(as.factor(bmp_tcom_count.dat$belonging),
              simulationOutput$scaledResiduals)
plotResiduals(as.factor(bmp_tcom_count.dat$gender),
              simulationOutput$scaledResiduals)
plotResiduals(as.factor(bmp_tcom_count.dat$mode),
              simulationOutput$scaledResiduals)
testOverdispersionParametric(glme_mod_1)
testZeroInflation(simulationOutput)
testUniformity(simulationOutput = simulationOutput)

# Those residuals are not good, need to take measures
# Perhaps best is to exclude zero count data, on notion
# that they reflect little activity

# First inspect zero counts by group space (here belonging)
bmp_tcom_count.dat %>% 
  select(belonging, visit_count) %>% 
  group_by(belonging) %>% 
  ggplot(aes(x = log(visit_count))) + 
    geom_histogram() +
    facet_wrap(~ belonging)
# That suggests a log of visit count ought to work, but would this mean
# we use Normal or Poisson?

# We try modeling with normal
lme_mod_visit_0 <-
  lmer(
    log(visit_count + 1) ~ as.factor(gender) + as.factor(mode) +
      as.factor(belonging) +
      (1 | trackingid),
    data = bmp_tcom_count.dat
  )
summary(lme_mod_visit_0)
# residuals
qqnorm(resid(lme_mod_visit_0))
plot(lme_mod_0)
# That looks moderately ok, except for band of residuals that likely
# correspond to the zero values

# We try modeling with normal again, removing zeros
lme_mod_visit_1 <-
  bmp_tcom_count.dat %>% 
  filter(visit_count != 0) %>% 
  lmer(
    log(visit_count) ~ as.factor(gender) + as.factor(mode) +
      as.factor(belonging) +
      (1 | trackingid),
    
    data = .
  )
# Throws an error, doesn't like the pipe, redo
bmp_tcom_count.dat %>% 
  filter(visit_count != 0) -> ty

lme_mod_visit_1 <-
  lmer(
    log(visit_count) ~ as.factor(gender) + as.factor(mode) +
      as.factor(belonging) +
      (1 | trackingid),
    data = ty
  )
summary(lme_mod_visit_1)
# residuals
qqnorm(resid(lme_mod_visit_1))
plot(lme_mod_visit_1)

# All seems ok, if marginal case of heteroscedasticity









