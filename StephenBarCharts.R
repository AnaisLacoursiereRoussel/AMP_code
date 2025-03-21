################################################################################
################################################################################
##### Bar chart figures

# Create the code for different bar charts of FlowCam (FC) /microscopy 
# (QA = quantitative assessment) comparisons

# The first ones are comparisons of every sample between the FC/QA data
# Shows relative abundance and stacked bar charts with breakdown of plankton within each sample

# Next ones show side by side bar charts of each taxa (as average relative abundance), 
# and how those differ between regions.

################################################################################

# Read in the script that reads in the QA (microscopy) and Flowcam data
# after the adjustments to the taxa list have been made
#source("QuantitativeAssessment/dataProcessing/redistributeTaxa.R")

################################################################################

# If I want to use redistributed data:
# fcQaDf = fcQaDf.redist
library(readxl)
fcQaDf <- read_excel("just_FC_and_MC_5000_df.xlsx")

# Create a function that creates various types of bar charts
# A data with the data for each region/year (i.e., 10 samples) is passed in
barChart = function(regData, newSampleLabels, plotLabel){
  
  # Get a consistent legend!
  # If ALL taxa are shown, there will be too many classes. Instead, find the __
  # most abundant taxa. Label all others as "Other"
  datRanked = regData %>%
    # Want counts per taxa for the whole ("all") dataset
    group_by(newName) %>%
    # summarize(countAll = sum(countTot)) %>%
    summarize(abundAll = sum(abund)) %>%
    mutate(rank = rank(-abundAll),
           # Keep __ most abundant classes, make the rest "Other"
           classRanks = ifelse(rank <=11, newName, "Other")) %>%
    mutate(relAbund = abundAll/sum(abundAll)) # if i want the relative abundance
  # write_xlsx(datRanked, "nl20relabunds-taketwo.xlsx")
  # Add this these new classes as a column in the original data frame for plotting
  datPlot = datRanked %>%
    left_join(regData) %>%
    group_by(classRanks, FlowCamID, type) %>%
    # If you don't recompute counts, the "Other" class will have a bunch of black lines
    # if you set the outline colour to black in geom_bar
    summarise(abund = sum(abund)) %>%
    # mutate(type = if_else(type == "QA", "MC", "FC")) # Change microscopy abbreviation to "MC"
    mutate(type = type)
  # Change the labels of the FlowCam IDs otherwise the text within facet labels are too long
  datPlot$FlowCamID = as.factor(datPlot$FlowCamID)
  datPlot$Original = as.factor(datPlot$FlowCamID) # create, just in case I need to know which is which
  levels(datPlot$FlowCamID) = newSampleLabels
  
  # Relative abundance chart
  # This needs to be based off either abund or countTot due to counting of different proportions per taxa within each sample
  relAbund =
    ggplot() +
    geom_bar(datPlot, mapping = aes(x=type, y=abund, fill=classRanks), position = "fill", stat = "identity", col = "black", linewidth = 0.05) +
    facet_grid(.~FlowCamID)+
    scale_y_continuous(labels = scales::percent_format(), name = "Relative Abundance")+
    scale_x_discrete(name = "Sample")+
    ggtitle(plotLabel)+
    scale_fill_brewer(palette = "Set3", name = "")+
    theme_bw()+
    theme(
      axis.title.x = element_blank()
    )+
    guides(fill = guide_legend(ncol=2))
  
  # return(datRanked)
}


# Create new labels for samples
# If these aren't changed, the facet labels have too much text
gulfLabs = c("S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10")
pacLabs = c("S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20")
nl20Labs = c("S21", "S22", "S23", "S24", "S25", "S26", "S27", "S28", "S29", "S30")
nl21Labs = c("S31", "S32", "S33", "S34", "S35", "S36", "S37", "S38", "S39", "S40")

# Create the charts. Pass in combined flowcam and QA data, but separately for each region
# Also pass in the new labels for each of the samples 
gulf20Chart = barChart(fcQaDf %>% subset(regionYear == "Gulf 2020"), gulfLabs, "(A) Gulf 2020")
pac21Chart = barChart(fcQaDf %>% subset(regionYear == "Pac 21"), pacLabs, "(B) Pacific 2021")
nl20Chart = barChart(fcQaDf %>% subset(regionYear == "NL 2020"), nl20Labs, "(C) Newfoundland 2020")
nl21Chart = barChart(fcQaDf %>% subset(regionYear == "NL 2021"), nl21Labs, "(D) Newfoundland 2021")

ggarrange(gulf20Chart, pac21Chart, nl20Chart, nl21Chart, ncol = 1)



################################################################################
################################################################################
## Plots of taxa comparisons
# Create bar charts that show the average relative abundance for each taxa
# Have FlowCam data beside the Microscopy data

# Get summary stats for plotting
# I want the average relative abundance for each sample in each region
sampleSummary = fcQaDf %>%
  # Get the relative abundance for each sample
  group_by(FlowCamID, type) %>%
  # mutate(relAbund = countTot / sum(countTot)*100) %>%
  mutate(relAbund = abund / sum(abund)*100) %>%
  # Now get the AVERAGE relative abundance for each sample. and the standard deviation
#  group_by(newName, type, regionYear) %>%
#  mutate(meanAbund = mean(relAbund),
#         sdAbund = sd(relAbund)) %>%
#  select(newName, regionYear, meanAbund, sdAbund, type) %>%
#  distinct() 
write_xlsx(sampleSummary, "ALLrelabunds.xlsx")
head(sampleSummary)
# I want to know, for each regionYear, which taxa are found: in both flowcam and microsopy, flowcam only, and microscopy only
sampleSummaryAdj = sampleSummary %>%
  # First need to expand the list, so it has ALL comparisons 
  # Because previously, only taxa that are present are recorded for each regionYear/type
  # But if, e.g., Calanoida (unid) only present in FlowCam but not in microscopy, I need a new line of data for that too
  ungroup() %>%
  complete(newName, regionYear, type, fill = list(meanAbund =0, sdAbund = 0)) %>%
  {print(head(.)); .} %>%
  # Now do the comparisons.
  group_by(newName, regionYear) %>%
  mutate(presence = case_when(
    all(meanAbund == 0) ~ "neither", # Not present in either
    any(meanAbund >0) & all(meanAbund >0) ~ "Both", # present in both
    meanAbund[type == "MC"] == 0 & meanAbund[type == "FC"] > 0 ~ "FC Only", # present in flowcam only
    meanAbund[type == "MC"] > 0 & meanAbund[type == "FC"] == 0 ~ "MC Only" # present in microscopy only
  )) 


# Create a function to make bar chart comparisons for each region
barChartComparison = function(relData, regionYear, plotName) {
  
  # Don't include the data that are present in neither
  # That would include taxa that are present in other region years, but not these ones
  relData = relData %>%
    filter(presence != "neither")
  
  # Make a plot that has 3 panels: present in both, present in flowcam only, present in microscopy only
  ggplot(relData, aes(x = meanAbund, y = newName, fill = type))+
    
    geom_col(position = position_dodge2(preserve = "single"), width = 0.8)+
    #facet_wrap(~presence, scales = "free_y", ncol = 1)+
    # facet_grid seems better for getting bar width the same, but panel sizes to be different
    facet_grid(presence~., scales = "free", space = "free", switch = "y")+ 
    # Need an ifelse statement to remove error bars when the standard deviation is zero by making them NAs
    geom_errorbar(aes(xmin=ifelse(sdAbund == 0, NA, meanAbund), xmax=ifelse(sdAbund == 0, NA, meanAbund+sdAbund)), width=.2,
                  position=position_dodge(.8), col = "black", na.rm = T)+
    xlab("Average relative abundance (%)")+
    ggtitle(plotName)+
    scale_fill_discrete(labels=c('FlowCam (FC)', 'Microscopy (MC)'))+
    scale_color_discrete(labels = c("FlowCam (FC)", "Microscopy (MC)"))+
    theme_bw()+
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 12),
      legend.text = element_text(size = 11),
      legend.title = element_blank(),
      strip.text.y = element_text(size = 10)
    )
  
}

# Create the bar charts for each region separately
barChartComparison(sampleSummaryAdj %>% filter(regionYear == "Gulf 2020"), "Gulf 2020", "Gulf 2020") 
barChartComparison(sampleSummaryAdj %>% filter(regionYear == "Pac 21"), "Pac 21", "Pacific 2021") 
barChartComparison(sampleSummaryAdj %>% filter(regionYear == "NL 2020"), "NL 2020", "Newfoundland 2020")
barChartComparison(sampleSummaryAdj %>% filter(regionYear == "NL 2021"), "NL 2021", "Newfoundland 2021")


################################################################################

#### COMMENTING OUT OLD CODE. KEEPING FOR REFERENCE


# This is my OLD graph that shows multiple graphs (e.g., using raw counts, not just relative abundance)

# # Create a function that creates various types of bar charts
# # A data with the data for each region/year (i.e., 10 samples) is passed in 
# barChart = function(regData, newSampleLabels){
#   
#   # Get a consistent legend!
#   # If ALL taxa are shown, there will be too many classes. Instead, find the __
#   # most abundant taxa. Label all others as "Other"
#   datRanked = regData %>%
#     # Want counts per taxa for the whole ("all") dataset
#     group_by(newName) %>%
#     summarize(countAll = sum(countTot)) %>%
#     mutate(rank = rank(-countAll),
#            # Keep __ most abundant classes, make the rest "Other"
#            classRanks = ifelse(rank <=11, newName, "Other")) %>%
#     mutate(relAbund = countAll/sum(countAll)) # if i want the relative abundance
#   
#   # Add this these new classes as a column in the original data frame for plotting
#   datPlot = datRanked %>%
#     left_join(regData) %>% 
#     group_by(classRanks, FlowCamID, type) %>%
#     # If you don't recompute counts, the "Other" class will have a bunch of black lines
#     # if you set the outline colour to black in geom_bar
#     summarise(count = sum(count),
#               abund = sum(abund),
#               countTot = sum(countTot))
#   
#   # Change the labels of the FlowCam IDs otherwise the text within facet labels are too long
#   datPlot$FlowCamID = as.factor(datPlot$FlowCamID)
#   levels(datPlot$FlowCamID) = newSampleLabels
#   
#   
#   # Create stacked bar chart based on counts
#   # This is the number of individuals counted by the taxonomists
#   # In the Pacific, different proportions of taxa were ID'd based on if they were "common" or "rare"
#   # So this will NOT be reflective of the true relative abundances. Need to multiply by the split multiplier in the excel files
#   # This plot is just to show what was actually counted
#   stackedBar =
#     ggplot()+
#     geom_bar(datPlot, mapping = aes(x = type, y = count, fill = classRanks), col = "black", linewidth = 0.05, stat = "identity")+
#     facet_grid(.~FlowCamID)+ # this shows flowcam & microscopy data side by side
#     scale_y_continuous(name = "Count")+
#     scale_x_discrete(name = "Sample")+
#     scale_fill_brewer(palette = "Set3", name = "")+
#     ggtitle("Number of individuals ID'd by taxonomists")+
#     theme_bw()+
#     theme(
#       axis.title.x = element_blank()
#     )+
#     guides(fill = guide_legend(ncol=2))
#   
#   # Relative abundance chart
#   # This needs to be based off either abund or countTot due to counting of different proportions per taxa within each sample
#   relAbund =
#     ggplot() +
#     geom_bar(datPlot, mapping = aes(x=type, y=abund, fill=classRanks), position = "fill", stat = "identity", col = "black", linewidth = 0.05) +
#     facet_grid(.~FlowCamID)+
#     scale_y_continuous(labels = scales::percent_format(), name = "Relative Abundance")+
#     scale_x_discrete(name = "Sample")+
#     ggtitle("Relative abundance")+
#     scale_fill_brewer(palette = "Set3", name = "")+
#     theme_bw()+
#     theme(
#       axis.title.x = element_blank()
#     )+
#     guides(fill = guide_legend(ncol=2))
#   
#   # Abundance in seawater
#   stackedAbund = 
#     ggplot()+
#     geom_col(data = datPlot, aes(x=type, y=abund, fill=classRanks), col = "black", linewidth = 0.05) +
#     facet_grid(.~FlowCamID)+
#     scale_fill_brewer(palette = "Set3", name = "")+
#     ggtitle("Abundance in seawater")+
#     scale_x_discrete(name = "Sample")+
#     theme_bw()+
#     theme(
#       axis.title.x = element_blank()
#     )+
#     guides(fill = guide_legend(ncol=2))
#   
#   # Counts per subsample (1/4 of tow)
#   stackedCountTotal = 
#     ggplot()+
#     geom_col(data = datPlot, aes(x=type, y= countTot, fill=classRanks), col = "black", linewidth = 0.05) +
#     facet_grid(.~FlowCamID)+
#     ggtitle("Number of individuals in subsample (1/4 of tow)")+
#     scale_y_continuous(name = "Count")+
#     scale_fill_brewer(palette = "Set3", name = "")+
#     #scale_x_discrete(name = "Sample")+
#     theme_bw()+
#     theme(
#       axis.title.x = element_blank()
#     )+
#     guides(fill = guide_legend(ncol=2))
#   
#   
#   # Put everything together
#   # ggarrange makes plots line up even when legends are different sizes
#   ggarrange(stackedBar, stackedCountTotal, relAbund)
#   
# }




###############################################################################

# # Test doing log transforms of the data
# 
# # Get summary stats for plotting
# # I want the average relative abundance for each sample in each region
# sampleSummary = fcQaDf %>%
#   # Get the relative abundance for each sample
#   group_by(FlowCamID, type) %>%
#   # mutate(relAbund = countTot / sum(countTot)*100) %>%
#   mutate(relAbund = abund / sum(abund)*100) %>%
#   # Now get the AVERAGE relative abundance for each sample. and the standard deviation
#   group_by(newName, type, regionYear) %>%
#   mutate(meanAbund = mean(relAbund),
#          sdAbund = sd(relAbund)) %>%
#   select(newName, regionYear, meanAbund, sdAbund, type) %>%
#   distinct() 
# 
# # I want to know, for each regionYear, which taxa are found: in both flowcam and microsopy, flowcam only, and microscopy only
# # It is very annoying to do lol
# sampleSummaryAdj = sampleSummary %>%
#   # First need to expand the list, so it has ALL comparisons 
#   # Because previously, only taxa that are present are recorded for each regionYear/type
#   # But if, e.g., Calanoida (unid) only present in FlowCam but not in microscopy, I need a new line of data for that too
#   ungroup() %>%
#   complete(newName, regionYear, type, fill = list(meanAbund =0, sdAbund = 0)) %>%
#   
#   # Now do the comparisons.
#   group_by(newName, regionYear) %>%
#   mutate(presence = case_when(
#     all(meanAbund == 0) ~ "neither", # Not present in either
#     any(meanAbund >0) & all(meanAbund >0) ~ "Both", # present in both
#     meanAbund[type == "QA"] == 0 & meanAbund[type == "FC"] > 0 ~ "FC Only", # present in flowcam only
#     meanAbund[type == "QA"] > 0 & meanAbund[type == "FC"] == 0 ~ "MC Only" # present in microscopy only
#   )) 
# 
# 
# # Create a function to make bar chart comparisons for each region
# barChartComparison = function(relData, regionYear, plotName) {
#   
#   # Don't include the data that are present in neither
#   # That would include taxa that are present in other region years, but not these ones
#   relData = relData %>%
#     filter(presence != "neither")
#   
#   # Make a plot that has 3 panels: present in both, present in flowcam only, present in microscopy only
#   ggplot(relData, aes(x = log10(meanAbund*100+1), y = newName, fill = type))+
#     
#     geom_col(position = position_dodge2(preserve = "single"), width = 0.8)+
#     #facet_wrap(~presence, scales = "free_y", ncol = 1)+
#     # facet_grid seems better for getting bar width the same, but panel sizes to be different
#     facet_grid(presence~., scales = "free", space = "free", switch = "y")+ 
#     # Need an ifelse statement to remove error bars when the standard deviation is zero by making them NAs
#     # geom_errorbar(aes(xmin=ifelse(sdAbund == 0, NA, meanAbund), xmax=ifelse(sdAbund == 0, NA, meanAbund+sdAbund)), width=.2,
#     #               position=position_dodge(.8), col = "black", na.rm = T)+
#     xlab("Average relative abundance (%)")+
#     ggtitle(plotName)+
#     scale_fill_discrete(labels=c('FlowCam (FC)', 'Microscopy (MC)'))+
#     scale_color_discrete(labels = c("FlowCam (FC)", "Microscopy (MC)"))+
#     theme_bw()+
#     theme(
#       axis.title.y = element_blank(),
#       axis.title.x = element_text(size = 12),
#       legend.text = element_text(size = 11),
#       legend.title = element_blank(),
#       strip.text.y = element_text(size = 10)
#     )
#   
# }
# 
# # Create the bar charts for each region separately
# barChartComparison(sampleSummaryAdj %>% filter(regionYear == "Gulf 2020"), "Gulf 2020", "Gulf 2020") 
# barChartComparison(sampleSummaryAdj %>% filter(regionYear == "Pac 21"), "Pac 21", "Pacific 2021") 
# barChartComparison(sampleSummaryAdj %>% filter(regionYear == "NL 2020"), "NL 2020", "Newfoundland 2020")
# barChartComparison(sampleSummaryAdj %>% filter(regionYear == "NL 2021"), "NL 2021", "Newfoundland 2021")