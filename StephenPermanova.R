################################################################################
### PERMANOVAS

# Run PERMANOVAs to see if community structure is different between Flowcam (FC)
# and microscopy (QA) samples 


################################################################################
## Read in the relevant data

# Run code that has NMDS ordinations
# This is where most of the data gets prepared 
source("QuantitativeAssessment/figuresAndStats/nmdsQA.R")

library(readxl)
fcQaDf <- read_excel("max_strat4_abunds_Feb4.xlsx")
fcQaDf <- read_excel("five_thousand_strat4_abunds_Feb4.xlsx")
fcQaDf <- fcQaDf[fcQaDf$type != 'HM', ]
fcQaDf <- fcQaDf[fcQaDf$type != 'HI', ]
fcQaDf <- fcQaDf[fcQaDf$type != 'CI', ]
print(fcQaDf$type)

# Adjust dataframe, calculate rel abundances
ordPrepDf = fcQaDf %>%
  select(newName, FlowCamID, regionYear, abund, type) %>%
  # Get the relative abundance for each sample
  group_by(FlowCamID, type) %>%
  mutate(relAbund = abund / sum(abund)*100)
#mutate(relAbund = abund / max(abund)*100)
head(ordPrepDf)


# Adjust dataframe, calculate rel abundances
ordPrepDf = fcQaDf %>%
  select(newName, FlowCamID, regionYear, abund, type) %>%
  # Get the relative abundance for each sample
  group_by(FlowCamID, type) %>%
  mutate(relAbund = sqrt(abund) / sum(sqrt(abund))*100)
#mutate(relAbund = abund / max(abund)*100)
head(ordPrepDf)

# Convert dataframes to wide format

# # Make one for abundance in seawater
# ordPrepWideAbund = ordPrepDf %>% 
#   select(-relAbund) %>%
#   pivot_wider(names_from = newName, values_from = abund) %>%
#   mutate_all(~replace(., is.na(.), 0))   # replace NAs with 0 

# Make the other using relative abundance
ordPrepWideRelAbund = ordPrepDf %>% 
  select(-abund) %>% # have to remove this for some reason or data doesn't pivot properly
  pivot_wider(names_from = newName, values_from = relAbund) %>%
  mutate_all(~replace(., is.na(.), 0))   # replace NAs with 0 
head(ordPrepWideRelAbund)

# Separate out each of the datasets
permGulf2020Data = ordPrepWideRelAbund %>% 
  filter(regionYear == "Gulf")

permPac2021Data = ordPrepWideRelAbund %>% 
  filter(regionYear == "Pacific")

permNL2020Data = ordPrepWideRelAbund %>% 
  filter(regionYear == "NL 2020")

permNL2021Data = ordPrepWideRelAbund %>% 
  filter(regionYear == "NL 2021")

print(permNL2021Data)
################################################################################
## Run the tests

# Using Bray-Curtis dissimilarities on the relative abundance data
# Need to specify which data to look at: only include species data (starts at "Acartia spp., ends at the last column)

# GULF 2020
permGulf2020Result = adonis2(permGulf2020Data[,which(colnames(permGulf2020Data)== "Acartia spp."):ncol(permGulf2020Data)]~type, data = permGulf2020Data, method = "bray", permutations = 9999) %>%
  # Multiply R2 by 100 to convert to percentages
  mutate_at(vars(c(R2)), .funs = funs(.*100)) %>% 
  # Add mean sum of squares (MS) column. Note: this adds a cell in the "total" row. Will need to manually remove this.
  mutate(MS = SumOfSqs/Df, .before = R2) %>%
  # Round most columns to 2 decimal places
  mutate(across(c(SumOfSqs, MS, R2, F), round, 3)) %>%
  mutate(dataset = "Gulf 2020", .before = Df)

print(permGulf2020Result)

# PACIFIC 2021
permPac2021Result = adonis2(permPac2021Data[,which(colnames(permPac2021Data)== "Acartia spp."):ncol(permPac2021Data)]~type, data = permPac2021Data, method = "bray", permutations = 9999) %>%
  # Multiply R2 by 100 to convert to percentages
  mutate_at(vars(c(R2)), .funs = funs(.*100)) %>% 
  # Add mean sum of squares (MS) column. Note: this adds a cell in the "total" row. Will need to manually remove this.
  mutate(MS = SumOfSqs/Df, .before = R2) %>%
  # Round most columns to 2 decimal places
  mutate(across(c(SumOfSqs, MS, R2, F), round, 3))%>%
  mutate(dataset = "Pacific 2021", .before = Df)

# NEWFOUNDLAND 2020
permNL2020Result = adonis2(permNL2020Data[,which(colnames(permNL2020Data)== "Acartia spp."):ncol(permNL2020Data)]~type, data = permNL2020Data, method = "bray", permutations = 9999) %>%
  # Multiply R2 by 100 to convert to percentages
  mutate_at(vars(c(R2)), .funs = funs(.*100)) %>% 
  # Add mean sum of squares (MS) column. Note: this adds a cell in the "total" row. Will need to manually remove this.
  mutate(MS = SumOfSqs/Df, .before = R2) %>%
  # Round most columns to 2 decimal places
  mutate(across(c(SumOfSqs, MS, R2, F), round, 3))%>%
  mutate(dataset = "Newfoundland 2020", .before = Df)

# NEWFOUNDLAND 2021
permNL2021Result = adonis2(permNL2021Data[,which(colnames(permNL2021Data)== "Acartia spp."):ncol(permNL2021Data)]~type, data = permNL2021Data, method = "bray", permutations = 9999) %>%
  # Multiply R2 by 100 to convert to percentages
  mutate_at(vars(c(R2)), .funs = funs(.*100)) %>% 
  # Add mean sum of squares (MS) column. Note: this adds a cell in the "total" row. Will need to manually remove this.
  mutate(MS = SumOfSqs/Df, .before = R2) %>%
  # Round most columns to 2 decimal places
  mutate(across(c(SumOfSqs, MS, R2, F), round, 3))%>%
  mutate(dataset = "Newfoundland 2021", .before = Df)

# Put all results together to be exported
allPNresults = bind_rows(permGulf2020Result, permPac2021Result, permNL2020Result, permNL2021Result)
print(allPNresults)
# Write results to csv
write.csv(allPNresults, "trimmed_type_three_approaches_PNresults.csv")
