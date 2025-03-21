################################################################################
################################################################################
#### NMDS FOR QUANTITATIVE ASSESSMENT DATA

# Note about transformations: arcsine square root transformed data, then Bray Curtis (ew)
# https://link.springer.com/article/10.1007/s10886-016-0784-x

# Square root transformed: https://peerj.com/articles/3437/

# square root: https://peerj.com/articles/3437/


################################################################################
# Prep the data

# Run the script that preps the QA and FlowCam data
#source("QuantitativeAssessment/dataProcessing/redistributeTaxa.R")

# If I want to use redistributed abundances (from redistributeTaxa.R), just set this:
#fcQaDf = fcQaDf.redist

# Function to load multiple packages
ipak = function(pkg){
  new.pkg = pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Choose necessary packages
packages = c("broom", "coin", "cowplot", "devtools", "dplyr", "egg", "ggiraph", "ggplot2", "ggrepel", "ggsignif", "ggsn", "ggnewscale", 
             "ggspatial", "ggthemes", "ggh4x", "gridExtra", "iNEXT", "jcolors", "leaflet", "lmodel2", "lsr", "mapr", "mapview",
             "patchwork", "pkgcond", "purrr", "readxl", "remotes", "reshape2", "rgdal", "rnaturalearth", "rnaturalearthdata", 
             "scales", "sf", "sp", "stringr", "tidyr", "tools", "useful", "vegan", "wbstats", "wpa", "writexl")
ipak(packages)


library(readxl)
fcQaDf <- read_excel("five_thousand_strat4_abunds_Feb4.xlsx")
fcQaDf <- read_excel("max_strat4_abunds_Feb4.xlsx")
fcQaDf <- read_excel("test_Jan22_max_strat4_abunds.xlsx")
fcQaDf <- read_excel("mixed_gulf_all_three_df.xlsx")
fcQaDf <- read_excel("relabunds.xlsx")
print(fcQaDf)
dim(fcQaDf)

#remove Bivalvia
fcQaDf <- fcQaDf[fcQaDf$newName != 'Bivalvia (larvae)', ]



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

# Create NMDS of both the FC and QA data
# Connect the samples from each method with a black line
# NMDS from each regionYear is made separately
nmdsMaker = function(df, regYearSelect, plotTitle){
  
  df = df %>%
    filter(regionYear == regYearSelect)
  
  # Select the beginning and end of the dataframe (with the species)
  beginNMDS = which(colnames(df)== "Acartia spp.")
  endNMDS = ncol(df)
  
  # Do NMDS ordination but only include species data
  ord = metaMDS(df[,c(beginNMDS:endNMDS)], distance = "bray", autotransform=FALSE)
  
  # Get NMDS coordinates from plot
  ordCoords = as.data.frame(scores(ord, display="sites")) %>%
    mutate(type = df$type,
           FlowCamID = df$FlowCamID)
  print(ordCoords)

  # Add NMDS stress
  # Note that round() includes UP TO 2 decimal places. Does not include 0s 
  ordStress = paste("2D Stress: ", format(round(ord$stress, digits=2), nsmall=2))
  
  # Make the ggplot item for each DFO region
  ggplot() + 
    geom_line(data = ordCoords, aes(x = NMDS1, y = NMDS2, group = FlowCamID), col = "gray20")+
    geom_point_interactive(data = ordCoords, aes(x=NMDS1, y=NMDS2, fill = type), pch = 21, size = 2, data_id = ordCoords$FlowCamID, tooltip = ordCoords$FlowCamID, onclick =ordCoords$FlowCamID)+ # Use pch=21 to get black outline circles
    ggtitle(plotTitle)+
    annotate("text", x = min(ordCoords$NMDS1), y=max(ordCoords$NMDS2), label = ordStress, size=4, hjust = -0.01)+
    scale_fill_discrete(labels=c('CI', 'HI', 'TM'))+
    theme_bw()+
    theme(axis.text = element_blank(),
          axis.title = element_text(size = 12), # don't want 
          axis.ticks = element_blank(),
          legend.text=element_text(size = 13),
          legend.title = element_blank(),
          #panel.border=element_rect(color="black", size=1), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          plot.margin=unit(c(0.3, 0.3, 0.3, 0.3),"cm"),
          plot.title = element_text(size=16))
  
}

#slightly different function that verifies that the labels are correct by printing the order
nmdsMaker = function(df, regYearSelect, plotTitle) {
  
  # Filter data for the selected region-year
  df = df %>%
    filter(regionYear == regYearSelect)
  
  # Select the beginning and end of the dataframe (with the species)
  beginNMDS = which(colnames(df) == "Acartia spp.")
  endNMDS = ncol(df)
  
  # Perform NMDS ordination using species data
  ord = metaMDS(df[, c(beginNMDS:endNMDS)], distance = "bray", autotransform = FALSE)
  
  # Extract NMDS coordinates and append relevant metadata
  ordCoords = as.data.frame(scores(ord, display = "sites")) %>%
    mutate(type = df$type,
           FlowCamID = df$FlowCamID)
  
  # Print ordCoords to verify mapping of 'type'
  print("NMDS Coordinates with Type Mapping:")
  print(ordCoords)
  file_path <- paste(regYearSelect, "max_nmds_coords.xlsx", sep = "_")
  write_xlsx(ordCoords, file_path)
  # Verify levels of the 'type' variable
  if (!is.factor(ordCoords$type)) {
    ordCoords$type <- factor(ordCoords$type)  # Convert to factor if not already
  }
  print("Levels of 'type' (used for legend mapping):")
  print(levels(ordCoords$type))
  
  # Add NMDS stress to plot title
  ordStress = paste("2D Stress:", format(round(ord$stress, digits = 2), nsmall = 2))
  
  # Generate the NMDS plot
  ggplot() + 
    geom_line(data = ordCoords, aes(x = NMDS1, y = NMDS2, group = FlowCamID), col = "gray20") +
    geom_point_interactive(
      data = ordCoords, 
      aes(x = NMDS1, y = NMDS2, fill = type), 
      pch = 21, size = 2, 
      data_id = ordCoords$FlowCamID, 
      tooltip = ordCoords$FlowCamID, 
      onclick = ordCoords$FlowCamID
    ) +
    ggtitle(plotTitle) +
    annotate("text", x = min(ordCoords$NMDS1), y = max(ordCoords$NMDS2), label = ordStress, size = 4, hjust = -0.01) +
    scale_fill_discrete(labels = c('CI', 'HI', 'HM')) +  # Customize labels if needed
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.title = element_text(size = 12),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 13),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"),
      plot.title = element_text(size = 16)
    )
}

# Using abundances
# gulf20.Ord = nmdsMaker(ordPrepWideAbund, "Gulf 2020", "(A) Gulf 2020")
# pac21.Ord = nmdsMaker(ordPrepWideAbund, "Pac 21", "(B) Pacific 2021")
# nl20.Ord = nmdsMaker(ordPrepWideAbund, "NL 2020", "(C) Newfoundland 2020")
# nl21.Ord = nmdsMaker(ordPrepWideAbund, "NL 2021", "(D) Newfoundland 2021")

# Using relative abundances
gulf20.OrdRel = nmdsMaker(ordPrepWideRelAbund, "Gulf", "(A) Gulf")
pac21.OrdRel = nmdsMaker(ordPrepWideRelAbund, "Pacific", "(B) Pacific")
nl20.OrdRel = nmdsMaker(ordPrepWideRelAbund, "NL 2020", "(C) Newfoundland 2020")
nl21.OrdRel = nmdsMaker(ordPrepWideRelAbund, "NL 2021", "(D) Newfoundland 2021")


# Put them all together
# # Abundance in seawater
# ggarrange(gulf20.Ord, pac21.Ord, nl20.Ord, nl21.Ord)

# Relative abundance
bigPlot = ggarrange(gulf20.OrdRel, pac21.OrdRel, nl20.OrdRel, nl21.OrdRel)
print(bigPlot)

################################################################################
# Make NMDS with all the data

nmdsMakerAllDat = function(df, plotTitle){
  
  beginNMDS = which(colnames(df)== "Acartia spp.")
  endNMDS = ncol(df)
  
  # Do NMDS ordination but only include species data
  ord = metaMDS(df[,c(beginNMDS:endNMDS)], distance = "bray", autotransform=FALSE)
  
  # Get NMDS coordinates from plot
  ordCoords = as.data.frame(scores(ord, display="sites")) %>%
    mutate(type = df$type,
           FlowCamID = df$FlowCamID,
           regionYear = df$regionYear) 
  print(ordCoords)
  write_xlsx(ordCoords, "nmds_coords.xlsx")
  # Add NMDS stress
  # Note that round() includes UP TO 2 decimal places. Does not include 0s 
  ordStress = paste("2D Stress: ", format(round(ord$stress, digits=2), nsmall=2))
  
  
  # Make the ggplot item for each DFO region
  ggplot() + 
    geom_line(data = ordCoords, aes(x = NMDS1, y = NMDS2, group = FlowCamID), col = "gray20")+
    geom_point_interactive(data = ordCoords, aes(x=NMDS1, y=NMDS2, fill = type, pch = regionYear), size = 3, data_id = ordCoords$FlowCamID, tooltip = ordCoords$FlowCamID, onclick =ordCoords$FlowCamID)+ # Use pch=21 to get black outline circles
    ggtitle(plotTitle)+
    scale_shape_manual(values=c(21:24), name = "Bay")+
    #stat_ellipse(data = ordCoords, aes(x = NMDS1, y = NMDS2, color = type), 
     #            level = 0.95, linetype = "dashed", size = 1) +  # Add ellipses for 95% confidence
    annotate("text", x = max(ordCoords$NMDS1), y=max(ordCoords$NMDS2), label = ordStress, size=4, hjust = 1)+
    scale_fill_discrete(labels=c('CI', 'HI', 'HM'))+
    theme_bw()+
    theme(axis.text = element_blank(),
          axis.title = element_text(size = 12), # don't want 
          axis.ticks = element_blank(),
          legend.text=element_text(size = 13),
          legend.title = element_blank(),
          #panel.border=element_rect(color="black", size=1), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          plot.margin=unit(c(0.3, 0.3, 0.3, 0.3),"cm"),
          plot.title = element_text(size=16))+
    guides(fill = guide_legend(override.aes = list(shape=21)))
  
}


nmdsMakerAllDat = function(df, plotTitle){
  
  beginNMDS = which(colnames(df)== "Acartia spp.")
  endNMDS = ncol(df)
  
  # Do NMDS ordination but only include species data
  ord = metaMDS(df[,c(beginNMDS:endNMDS)], distance = "bray", autotransform=FALSE)
  
  # Get NMDS coordinates from plot
  ordCoords = as.data.frame(scores(ord, display="sites")) %>%
    mutate(type = df$type,
           FlowCamID = df$FlowCamID,
           regionYear = df$regionYear) 
  print(ordCoords)
  #write_xlsx(ordCoords, "nmds_coords.xlsx")
  # Add NMDS stress
  ordStress = paste("2D Stress: ", format(round(ord$stress, digits=2), nsmall=2))
  
  # Make the ggplot item for each DFO region
  ggplot() + 
    geom_line(data = ordCoords, aes(x = NMDS1, y = NMDS2, group = FlowCamID), col = "gray20") +
    geom_point_interactive(data = ordCoords, aes(x=NMDS1, y=NMDS2, fill = type, pch = regionYear), 
                           size = 3, data_id = ordCoords$FlowCamID, tooltip = ordCoords$FlowCamID, 
                           onclick = ordCoords$FlowCamID) +
    ggtitle(plotTitle) +
    scale_shape_manual(values=c(21:24), name = "Bay") +
    annotate("text", x = max(ordCoords$NMDS1), y = max(ordCoords$NMDS2), 
             label = ordStress, size = 4, hjust = 1) +
    scale_fill_discrete(labels=c('CI', 'HI', 'HM')) +
    
    # Add an ellipse for each combination of regionYear and type
    stat_ellipse(data = ordCoords, aes(x = NMDS1, y = NMDS2, color = interaction(type, regionYear)), 
                 level = 0.95, linetype = "dashed", size = 1) +  # 95% confidence ellipses for each regionYear and type
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.title = element_text(size = 12), # don't want 
          axis.ticks = element_blank(),
          legend.text=element_text(size = 13),
          legend.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          plot.margin=unit(c(0.3, 0.3, 0.3, 0.3),"cm"),
          plot.title = element_text(size=16)) +
    guides(fill = guide_legend(override.aes = list(shape=21)))
}


# Run the function (for all data) using abundance in seawater
#nmdsAll = nmdsMakerAllDat(ordPrepWideAbund, "Abundance in seawater: all regions")
# Run using relative abundance

#allRegionPlot = ggarrange(nmdsAll.RelAbund, nmdsAll, ncol =1)

#png("sample_plot_max.png")
nmdsAllRelAbund = nmdsMakerAllDat(ordPrepWideRelAbund, "Relative abundance: all regions")
print(nmdsAllRelAbund)
#dev.off()
################################################################################
### EXPLORING THE DATA A BIT MORE

# Commented out since I am not using these

# # There are 3 Gulf samples with very high Bivalvia larvae that seem very different. Plot these separately
# gulf3Weird = nmdsMaker(ordPrepWideRelAbund %>%
#                          filter((FlowCamID %in% c("AMMP_Gulf_StPeters_1_20200903HT_250UM", 
#                                                   "AMMP_Gulf_StPeters_1_20200904HT_250UM", 
#                                                   "AMMP_Gulf_StPeters_3_20200903LT_250UM"))), "Gulf 2020", "Gulf 2020")
# 
# gulfOtherWeird = nmdsMaker(ordPrepWideRelAbund %>%
#                              filter(!(FlowCamID %in% c("AMMP_Gulf_StPeters_1_20200903HT_250UM", 
#                                                        "AMMP_Gulf_StPeters_1_20200904HT_250UM", 
#                                                        "AMMP_Gulf_StPeters_3_20200903LT_250UM"))), "Gulf 2020", "Gulf 2020")
# 
# # Plot these side by side
# ggarrange(gulf3Weird, gulfOtherWeird, ncol = 2)
# 
# # Try and see if I can figure out why 2 FC samples are overlapping?
# girafe(ggobj = gulf3Weird)
# 
# ################################################################################
# ## Make interactive plots to visualize "problem" samples (i.e., why does the Gulf figure look so weird???)
# 
# # Make it interactive so you can point & click and see the sampleCode
# girafe(ggobj = nmdsAll)
# girafe(ggobj = nmdsAll.RelAbund)
# 
# 
# 
# ### Track down some data
# # Two of the plots had distinct groupings
# # This is: Gulf 2020 (left hand side has 3 samples vs RHS has 7)
# # For Newfoundland, there is an Upper cluster (4 points) and Lower cluster (6 points) on the NMDS
# 
# # For Gulf, need to go back to the original data that specifies if they're inner/mid/outer etc
# # Too many sample labelling issues to get that from just the code 
# # For NL, just point and click on the interactive plot. Can get info from FlowCamCodes.
# 
# # This creates an interactive ggplot that I can click on to see the FlowCamCodes
# girafe(ggobj = gulf20.Ord)
# girafe(ggobj = nl21.Ord)
# 
# girafe(ggobj = gulf20.OrdRel)
# 
# # Gulf has 3 that form 1 group. They are:
# # THESE ARE ALL "OUTER" STATIONS
# # AMMP_Gulf_StPeters_1_20200903HT_250UM
# # AMMP_Gulf_StPeters_1_20200904HT_250UM
# # AMMP_Gulf_StPeters_3_20200903LT_250UM
# 
# # Need to go back to the original code from QAcodeMatches.R before i select for certain columns
# findWeirdGulfStations = qaID %>%
#   left_join(fcDataForQA, by = c("FlowCamID" = "flowcamCode")) %>%
#   # Only select the samples we're interested in
#   filter(selectForAnalysis == "Yes") %>%
#   left_join(fcTaxaChanges) %>%
#   filter(regionYear == "Gulf 2020") %>%
#   filter(!(FlowCamID %in% c("AMMP_Gulf_StPeters_1_20200903HT_250UM", 
#                           "AMMP_Gulf_StPeters_1_20200904HT_250UM", 
#                           "AMMP_Gulf_StPeters_3_20200903LT_250UM")))
# 
# # End result: 
# # For Gulf, LHS are Outer stations. RHS are Inner/Mid stations
# # For NL, upper is station 41 ("Mid-B"). Lower = station 17 ("Outer")
# 
# ################################################################################
# ### Find out what's behind these differences using SIMPER analysis
# 
# 
# testNL = ordPrepWide %>%
#   filter(regionYear == "NL 2021")
# 
# simTest = simper(sqrt(testNL[,which(colnames(testNL)== "Acartia spp."):ncol(testNL)]), 
#                  group=testNL$type)
# summary(simTest)


# Ok I think these three just have WAY more bivalvia larvae than the rest