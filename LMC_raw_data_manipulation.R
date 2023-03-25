# Arianna Chiti
# Parasitoid wasps
# LMC experiment
# Created on 14.02.2023
# Last edited on 15.02.2023


# Libraries ####

library(here)
library(tidyverse)
library(readxl)
library(openxlsx)



# Data ####

LMC_data <- read_excel(here("../Raw data - 3.2.23.xlsx"))
Treatments <- read_excel(here("../LMC blocks 1-3 data.xlsx"))


# Data manipulation ####

# Check str of data objects

head(LMC_data)
str(LMC_data)

head(Treatments)
str(Treatments)


# Edit data format and type of treatment data set for later binding

Treatments <- 
  Treatments %>%
  mutate(Block = as.factor(Block),
         Tube_number = as.factor(ID),
         `Treatment (# females)` = as.factor(`Treatment (# females)`))


# Bind treatments to LMC raw data

LMC_with_treatment_data <- 
  LMC_data %>%
  mutate(Tube_number = as.factor(`Tube number`)) %>%
  left_join(Treatments %>% 
              select(Block, `Treatment (# females)`, `Date/time set up`, `Aphids added to plants`, Tube_number),
            by = "Tube_number") %>%
  select(Date, `Tub number`,Dish,Cup, Tube_number, `Treatment (# females)`, Block, Males, Females, `Male photo 1`, `Male photo 2`,`Male photo 3`, `Female photo 1`,`Female photo 2`, `Female photo 3`, `Collection method`, `Date/time set up`, `Aphids added to plants`) %>%
  mutate(`Date/time set up` = as.character(`Date/time set up`),
         `Aphids added to plants` = as.character(`Aphids added to plants`))


# Save to excel file

write.xlsx(LMC_with_treatment_data,
           here("../LMC_raw_data_with_treatment.xlsx"),
           colNames = TRUE,
           sheetName="Raw-with-treatment")






