################################################################################
################################################################################
##### T-Tests

# Conduct t-tests (?) between relative abundance of taxa from flowcam vs microscopy

################################################################################

# Read in the script that creates bar charts for the different comparisons
# This has already created a dataframe of relative abundance for each sample
source("StephenBarCharts.R")

################################################################################
## Prep the data 

# Background: I need to be doing paired tests (wilcoxon or t-tests, TBD) to test for 
# differences in relative abundance between FlowCam (FC) and Microscopy (QA) methods
# Using paired tests (vs independent tests) because they are from the same sample and 
# therefore not independent.

# Paired tests require the same # of samples per group since there should be an equivalent
# value for each sample in FC/QA 
# Therefore, if a taxa was observed in one sample for QA/FC but not the other, I need to say the 
# relative abundance was 0 in the sample where it was not observed.
# Adding zeroes means that I need to use the complete() function to create all pairs of data by FlowCamID and taxa

# I also get errors when I'm comparing 0 relative abundance between 0 relative abundance (i.e., absent in both samples)
# So I need to remove those. BUT only when they are absent in NONE of the samples in that regionYear.


# From the bar chart figures code
# This tells me if taxa were observed in both of the regionYears, neither of the regionYears, FC Only, or MC only
# Note, this is for the overall regionYear, not by individual sample. e.g., if Acartia was found in 1 FC sample, but none
# of the others (and not in any of the microscopy samples), the sample would be ID'd as "FC Only"
sumWant = sampleSummaryAdj %>%
  select(newName, regionYear, type, presence) # only keep a few of the columns

# Get FlowCamID and regionYear (for filling in NAs later)
qaID.edit = qaID %>%
  select(FlowCamID, regionYear)

# Create the data frame described above with zeroes when there are no taxa
dfAllPairs = fcQaDf %>%
  # Get the relative abundance of each taxa for each sample
  group_by(FlowCamID, type) %>%
  mutate(relAbund = abund / sum(abund)*100) %>%
  # Need to ungroup after using group_by or functions won't work
  ungroup() %>%
  # Create all comparisons. When there's no data, put relAbund as zero
  complete(newName, FlowCamID, type, fill = list(relAbund =0)) %>%
  # Now, edit the df: I think there's a better way to do this. But above function will put NAs in some of the columns.
  # Need to fill in the regionYears that became NAs. Do this by joining with a dataframe of just FlowCamID/regionYear
  left_join(qaID.edit, by = c("FlowCamID")) %>%
  mutate(regionYear.x = coalesce(regionYear.x, regionYear.y)) %>%
  rename(regionYear = regionYear.x) %>%
  # Join with df that specifies whether taxa are observed in "both" FC/MC, just FC, or just MC for each regionYear
  full_join(sumWant) 


################################################################################
# Start running the stats

# Shapiro Wilk test
# This tests if the data are normally distributed: http://www.sthda.com/english/wiki/paired-samples-t-test-in-r
# Note: the data is normal if the p-value is ABOVE 0.05
shWilk = dfAllPairs %>%
  # Get an error message when the data from both pairs are all 0
  # Therefore, only conduct for taxa that are present in each regionYear for both C/QA 
  filter(presence == "Both") %>%
  # I want to conduct t-tests between all taxa for each regionYear 
  group_by(newName, regionYear) %>%
  # But need >2 observations total or I get an error message
  mutate(sample_size = n()) %>%
  filter(sample_size >2) %>%
  # Need to conduct on the DIFFERENCES between value (see point 5 just below References: https://statsandr.com/blog/wilcoxon-test-in-r-how-to-compare-2-groups-under-the-non-normality-assumption/#fnref5)
  group_by(newName, regionYear, FlowCamID) %>%
  summarize(dif = relAbund[type == "FC"] - relAbund[type =="QA"]) %>%
  # Conduct the test
  # summarize(p_value = shapiro.test(relAbund[type %in% c("FC", "QA")])$p.value, .groups = "drop") # idk why I made this so complicated lol
  group_by(newName, regionYear) %>%
  summarize(p_value = shapiro.test(dif)$p.value)

# Test: do I get the same results if I run it for just one taxa/regionYear
testSh = dfAllPairs %>%
  filter(newName == "Acartia spp." & regionYear == "Pac 21")
# Yes, it matches what's in my df above
shapiro.test(testSh$relAbund)


# Some of the data are NOT normally distributed
# Proceed with a Wilcoxon test instead of a t-test


# Wilcoxon test
wilcRes = dfAllPairs %>% 
  filter(presence == "Both") %>%
  
  # I want to conduct t-tests between all taxa for each regionYear 
  group_by(newName, regionYear) %>%
  
  # Run the t-tests for each grouping listed above
  do(t_test_result = tidy(wilcox.test(relAbund~type, data = ., paired = T))) %>%
  # If this isn't added, the t-test results get stored as lists instead of separate columns
  unnest_wider(t_test_result)


wilcox.test(testSh$relAbund~testSh$type, paired =T)
wilcox.test(testSh$relAbund~testSh$type, paired =T)

# T-TEST 
test = fcQaDf %>%
  # Get the relative abundance for each sample
  group_by(FlowCamID, type) %>%
  mutate(relAbund = countTot / sum(countTot)*100) %>%
  # Join with df that specifies whether taxa are observed in "both" FC/MC, just FC, or just MC for each regionYear
  left_join(sumWant) %>%
  # Need to ungroup or functions won't work
  ungroup() %>%
  # Create all comparisons. when there's no data, put relAbund as zero
  complete(newName, regionYear, FlowCamID, type, fill = list(relAbund =0)) %>%
  # I want to conduct t-tests between all taxa for each regionYear 
  group_by(newName, regionYear) %>%
  
  mutate(avg_relAbund = mean(relAbund, na.rm = T))

# Run the t-tests for each grouping listed above
do(t_test_result = tidy(t.test(relAbund~type, data = ., paired = T))) %>%
  # If this isn't added, the t-test results get stored as lists instead of separate columns
  unnest_wider(t_test_result) %>%
  mutate(avg_relAbund = mean(relAbund, na.rm = T))


################################################################################
###
# Test coin package for the paired Wilcoxon test
# This is to handle the case of "ties" https://library.virginia.edu/data/articles/the-wilcoxon-rank-sum-test

# I am having MAJOR issues with the R coin package
# Error in UseMethod("round_any") : no applicable method for 'round_any' applied to an object of class "c('ScalarIndependenceTestConfint', 'ScalarIndependenceTest', 'IndependenceTest')"
install.packages("coin")
library("coin")

wilcox.test(testSh$relAbund~testSh$type, paired =T)
x = wilcox_test(relAbund~as.factor(type), data = testSh, paired = T, distribution = "exact")
x = wilcox_test(relAbund~as.factor(type), data = testSh, paired = T)

pvalue(x)


# Wilcoxon test using coin package
wilcRes2 <- dfAllPairs %>%
  filter(presence == "Both") %>%
  group_by(newName, regionYear) %>%
  do(t_test_result = {
    dat <- data.frame(relAbund = .$relAbund, type = .$type)
    test <- wilcox_test(relAbund ~ as.factor(type), data = dat, distribution = "exact", paired = TRUE, conf.int = T)
    data.frame(
      est = statistic(test, type = "linear"),
      p.value = pvalue(test),
      CI.low = confint(test)$conf.int[1],
      CI.up = confint(test)$conf.int[2])
  }) %>%
  ungroup() %>%
  unnest_wider(t_test_result) %>%
  # Make adjustments to the dataset
  mutate(p.value = round(as.numeric(p.value),4),
         CI95 = paste0(round(CI.low,2), ", ", round(CI.up,2)))

# extracting the z score:
# https://stats.stackexchange.com/questions/330129/how-to-get-the-z-score-in-wilcox-test-in-r


## Example data
## Davis (1986, p. 140)
davis <- matrix(
  c(3, 6,
    2, 19),
  nrow = 2, byrow = TRUE
)

install.packages("coin")
library("coin")


diffusion = data.frame(
  pd = c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46,
         1.15, 0.88, 0.90, 0.74, 1.21),
  age = factor(rep(c("At term", "12-26 Weeks"), c(10, 5)))
)
## Exact Wilcoxon-Mann-Whitney test
## Hollander and Wolfe (1999, p. 111)
## (At term - 12-26 Weeks)
wt <- wilcox_test(pd ~ age, data = diffusion,
                  distribution = "exact", conf.int = TRUE)


pvalue(wt)

install.packages("coin")
library("coin")

#> 
#>  Wilcoxon rank sum exact test
#> 
#> data:  total_perMB by gender
#> W = 15, p-value = 0.6095
#> alternative hypothesis: true location shift is not equal to 0
## Expectation, variance, two-sided pvalue and confidence interval
expectation(wt)
covariance(wt)
pvalue(wt)
confint(wt)