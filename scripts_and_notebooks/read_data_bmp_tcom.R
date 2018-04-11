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
bmp_tcom.dat <- read_excel("data/UserTcommInfo.xlsx",
                           sheet = "Trackdata") %>% 
  rename_all(tolower) 

bmp_tcom.dat %>% 
  group_by(trackingid, belonging, tmain, mode) %>% 
  summarise(sum_mins = sum(mins, na.rm = T),
            sd_mins = sd(mins, na.rm = T),
            n_mins = n())
save(file = "data/bmp_tcom.dat", bmp_tcom.dat )
