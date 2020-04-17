library(sqldf)
library(auk)
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)

# auk_set_ebd_path(path = "/Users/MarkRoth/Documents/Oregon State/Year 1/Research/eBird/Class Imbalance/", overwrite = TRUE)
auk_get_ebd_path()
f_in_weta <- "/Users/MarkRoth/Documents/Oregon State/Year 1/Research/ICB/data generation/ebd_weta_breeding_or_zf_no2020protocol.csv"
WETA <- read.delim(f_in_weta, header=TRUE, sep = ",")

############################################################
# auk functions expect data formatted in a specific way ####
#   that sometimes doesn't match the data we have ##########
#
# i have found that we need to append the columns ##########
#   "country_code" and "group_identifier" to the WETA data #
#   
# since we aren't concerned with these two columns in our ##
#   analysis and the filter_repeat_visits function doesn't #
#   use them, I wasn't overly concerned with what they are #
#   filled with. ###########################################
############################################################
WETA$country_code = "US"
WETA$group_identifier = "-"

WETA$observation_date
WETA$observation_date <- as.character(WETA$observation_date)
WETA$formatted_date <- mdy(WETA$observation_date)

# check to make sure it is in the %Y-%m-%d format
WETA$formatted_date

occ <- filter_repeat_visits(WETA, 
                            min_obs = 1, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "formatted_date",
                            site_vars = c("locality_id", "observer_id"))

occ_filtered <- filter_repeat_visits(WETA, 
                            min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "formatted_date",
                            site_vars = c("locality_id", "observer_id"))

occ_no_limits <- filter_repeat_visits(WETA, 
                                     min_obs = 1, max_obs = 1000000,
                                     annual_closure = TRUE,
                                     date_var = "formatted_date",
                                     site_vars = c("locality_id", "observer_id"))

summary(occ)
summary(occ_filtered)
summary(occ_no_limits)

############################################################
# to account for the output of auk's filter_repeat_visits, #
# if we want to look at site_specific information, we want #
# to filter out the duplicated counts ######################
############################################################
occ_no_dups <- subset(occ_no_limits, !duplicated(site))

# rename the site names as row numbers for display purposes
occ_no_dups$site_id <- as.numeric(rownames(occ_no_dups))
occ_first_100 <- occ_no_dups[1:100,]
ggplot(occ_first_100, aes(x = site_id, y = n_observations)) + geom_col()

num1s <- sqldf('SELECT * FROM occ_no_dups WHERE n_observations == 1')
occ_no_dups_ordered <- occ_no_dups[order(-occ_no_dups$n_observations),]
occ_first_100_ordered <- occ_no_dups_ordered[1:100,]
occ_first_100_ordered$site_id <- as.numeric(rownames(occ_first_100_ordered))

ggplot(occ_first_100_ordered, aes(x = site_id, y = n_observations)) + geom_col() + ggtitle("top 100 most visited sites")
ggplot(occ_no_dups_ordered, aes(x = site_id, y = n_observations)) + geom_col() + ggtitle("ordered visited sites")
