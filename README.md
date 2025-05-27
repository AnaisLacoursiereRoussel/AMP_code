This is an explanation of the code in this repository as it relates to the research article titled ‘**Contrasting the efficiency of imaging systems for mesozooplankton indicators across Pacific and Atlantic coastal ecosystems**’  from the journal, Ecological Indicators. Code is summarized in the order that appears in the Results section of the article. If you have any questions about using the code, or reproducing results, contact lukecmclean@gmail.com.

# Preparation

* Much of the analysis in this paper depends on two dataframes called ml\_df\_5000 and ml\_df\_max. This stands for Machine Learning DataFrame and uses the models that were trained with a maximum training class size of 5000 and an unlimited maximum training class size (max). These dataframes are imported and produced in a file called ‘import\_ml\_df.ipynb’. The data is imported from ‘CyrilUnzipped/AMP FlowCam strategy predictions results and reports.’ The dataframes are stored as variables to be used by other files.  
* Other analysis depends on the relative abundance data, ‘five\_thousand\_strat4\_relabunds\_feb4.xlsx’ and ‘max\_strat4\_relabunds\_feb4.xlsx’. This data was compiled and organized in ‘combine\_all\_results.ipynb’ using ml\_df\_5000 and ml\_df\_max (described above) and ‘./FinnisData-20241007T133632Z-001/FinnisData/Method Paper/methodPaperDataForFigures.xlsx’

# 3\. Results

## 3.1 Evaluation of the ML Model

### 3.1.1 Effect of Training Parameters on Model Performance

The python file used for most of these calculations is called f1\_by\_strategy.ipynb

* It first loads in all the different training strategies and their corresponding weighted f1 scores from ‘CyrilUnzipped/AMP FlowCam strategy predictions results and reports’ and loads this into a dataframe with 256 rows, one for each strategy and region.  
* It then cleans the data by adding null to any column where the strategies were not actually applicable.  
* It then calculates and plots a bar graph showing the average effect of each parameter on all the models (line 378-379).  
* It then conducts some preliminary statistical tests demonstrating whether those differences were significant (not included in the paper due to redundancy / beyond scope).  
* It then computes graphs for each parameter. These were included in Appendix B, (except ‘Deep Features’ which were deemed inappropriate for the analysis due to being calibrated for the micro FlowCam rather than the macro FlowCam which was used in the study).  
* It then computes a few different figures. The clearest and most significant difference is that of increasing maximum training class size (Figure 2a).

The python file that calculates and compares the size of the training libraries is called ‘training\_data\_bargraph.ipynb’

* It first imports the Excel File ‘regional\_prediction\_strategies.xlsx’. This contains the library sizes.  
* It then organizes a bargraph ordered by the largest taxa libraries for each region (Figure 2b)

### 3.1.2. Performance of the Computer Model

The python file measures model performance as it relates to inside and outside of taxonomic class and order is called ‘error\_by\_phyla.ipynb’

* It first imports the ml\_df\_5000 and ml\_df\_max dataframes from ‘import\_ml\_df.ipynb’  
* It then defines phylogenetic trees for class and order as python dictionaries  
* It then calculates whether each ‘predicted’ and ‘actual’ image classification is in the same order or class and presents the amounts as bar graphs by region.  
* It then produces a boxplot comparing the results (Figure 3).

The python file that produces confusion matrices (Figure 4 and Appendix E) for the model performance is called ‘CM\_and\_prediction\_resampling.ipynb’

* It first imports the ml\_df\_5000 and ml\_df\_max dataframes from ‘import\_ml\_df.ipynb’  
* Then it defines the taxa order based off of phylogenetic trees that were produced in R (taxa\_trees.R)  
* Finally, it builds the confusion matrices (first by count, then by percentage)

The python file that runs regression of taxa F1 scores upon taxa training image count is called ‘f1\_by\_training\_image\_counts.ipynb’.

* It first imports the test image counts and f1 scores from ‘CyrilUnzipped/AMP FlowCam strategy predictions results and reports’  
* Then it imports the training image counts from ‘regional\_prediction\_strategies.xlsx’  
* Then it organizes the data and prepares a color scheme for the most abundant taxa  
* It calculates the regression statistics and creates scatterplots

### 3.1.3 Effect of Sample Characteristics on F1 Scores

The python file that explores sample characteristics and their effect on F1 scores is called ‘CreatureFeatures.ipynb’

* It is complex because it imports the various spreadsheets with sample features and each dataset has slightly different formatting for these features. It imports them one by one, each with its own cleaning to make them compatible  
* ml\_df\_5000 is imported from ‘import\_ml\_df.ipynb’  
* richness dictionaries are imported from ‘Richness\_BoxPlots.ipynb’  
* f1 scores are generated and organized from ml\_df\_5000  
* some Gulf sample names are corrected to be compatible across datasets  
* the above data is gathered together into a dataframe called final\_feature\_df  
* data is normalized and PCA is prepared  
* a heatmap (correlation table) is prepared for all features (Appendix F)  
* histograms are generated to inspect distributions  
* scatter plots and regression lines are computed for individual features  
* Multiple regression is conducted based on which features are input as X

## 3.2 A Comparative Analysis of Ecosystem Indicators

The relative taxa abundances and corresponding bar graph (Figure 7\) are produced in ‘StephenBarCharts\_one\_legend.R’

* The imported spreadsheets, such as ‘max\_strat4\_abunds\_Feb4.xlsx’ are produced in ‘combine\_all\_results.ipynb’. This python file combines relative abundance data from ‘./FinnisData-20241007T133632Z-001/FinnisData/Method Paper/methodPaperDataForFigures.xlsx’ with ml\_df imported from ‘import\_ml\_df.ipynb’  
* It can be adjusted to gather different results from the various models used. The paper uses strategy 4 (which represents the set of training parameters described in the paper) and either 5000 or unlimited maximum class size.  
* The R file’StephenBarCharts\_one\_legend.R’ then produces the bar graph.

The NMDS plot (Figure 8\) is also produced by multiple files.

* ‘StephenNMDS.R’ imports the same relative abundance data described above.  
* It produces NMDS plots and saves the coordinates as an excel file.  
* ‘NMDS\_ellipses.ipynb’ imports these coordinates and adds ellipses based on confidence intervals to the NMDS.

Sample dissimilarity (in relative abundances) (Table 2\) was computed and tested in ‘Sample\_Dissimilarity.ipynb’

* it imports the same relative abundance data described above  
* it adds any taxa that were represented in some but not all samples to have 0s rather than than absence  
* It then produces boxplots for the dissimilarity between samples  
* It then runs statistical tests

### 3.2.1 Differences in Relative Abundance

PERMANOVA analyses were run in ‘StephenPermanova.R’

* It imports the same relative abundance data described above  
* It organizes the data and runs the PERMANOVA

### 3.2.2 Biodiversity Indices

The biodiversity indices were computed and tested in ‘Richness\_BoxPlots.ipynb’

* It imports the raw counts from the microscopy and flowcam spreadsheets, ‘./FinnisData-20241007T133632Z-001/FinnisData/Method Paper/microscopyDataMethodPaper\_RawCounts.xlsx’ and ‘./FinnisData-20241007T133632Z-001/FinnisData/Method Paper/flowcamDataMethodPaper\_RawCounts.xlsx’  
* It tests the naming against the relative abundance spreadsheet, ‘./FinnisData-20241007T133632Z-001/FinnisData/Method Paper/methodPaperDataForFigures.xlsx’, as a double check on naming conventions  
* It then calculates the various indices, tests assumptions, runs tests, and produces boxplots