#########################################################
#           BELFAST MOBILITY PROJECT                    #
#           READ IN T COMMUNITY DATA                    #
#           Author CG TREDOUX                           #
#########################################################

# load some libraries-----------------------------------------------------------
library(pacman)
p_load(tidyverse, readxl, haven, magrittr)

# format numbers to avoid sci representation
options(scipen = 3)

# read in  data from Excel sheet------------------------------------------------
bmp_tcom.dat <- read_excel("data/UserTcommInfoCleaned.xlsx",
                           sheet = "Trackdata") %>%
  rename_all(tolower) %>%
  mutate(tcommunity = ifelse(tcommunity == "Neutral", "Shared", tcommunity)) %>%
  mutate(belonging = ifelse(belonging == "Neutral", "Shared", belonging)) %>%
  mutate(mode = case_when(mode == 1 ~ "Walk",
                          mode == 2 ~ "Unsure",
                          mode == 3 ~ "Vehicle",
                          TRUE ~ "Else"))

bmp_tcom.dat %>% 
  group_by(trackingid, belonging, tmain, mode) %>% 
  summarise(sum_mins = sum(mins, na.rm = T),
            sd_mins = sd(mins, na.rm = T),
            n_mins = n())
saveRDS(file = "data/bmp_tcom.Rds", bmp_tcom.dat )
