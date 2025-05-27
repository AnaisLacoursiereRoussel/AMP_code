library(tidyr)

richdf <- read_excel("richness_data.xlsx")
print(richdf)
gulf_df <- richdf[richdf$dataset == 'Gulf', ]
gulf_df
par(mfrow=c(1, 3))
hist(gulf_df$CI, prob=T)
hist(gulf_df$HI, prob=T)
hist(gulf_df$HM, prob=T)

# gulf_B<-matrix(nrow=30,ncol=4,c(rep(1:10,3),gulf_df$CI,
#                               fulg_df$HI,gulf_df$HM,rep(1,10),rep(2,10),rep(3,10)))

long_richness = read_excel("pivoted_richness_data.xlsx")
print(long_richness)

library(ez)

repeat1<ezANOVA(data=long_richness,dv=.(richness),wid=.(sample),within=.(sorting_process
),type=3)
repeat1

unique_datasets <- unique(long_richness$dataset)

# Iterate over each unique value
for (dataset_value in unique_datasets) {
  print(dataset_value)
  # Subset the data frame based on the current dataset_value
  region_df <- long_richness[long_richness$dataset == dataset_value, ]

  repeat1 <- ezANOVA(data = region_df, dv = .(richness), wid = .(sample), within = .(sorting_process), type = 3)
  print(repeat1)
}



long_pielou = read_excel("pivoted_pielou_data.xlsx")
print(long_pielou)

repeat1<ezANOVA(data=long_pielou,dv=.(pielou),wid=.(sample),within=.(sorting_process
),type=3)
repeat1

unique_datasets <- unique(long_pielou$dataset)

# Iterate over each unique value
for (dataset_value in unique_datasets) {
  print(dataset_value)
  # Subset the data frame based on the current dataset_value
  region_df <- long_pielou[long_pielou$dataset == dataset_value, ]
  
  repeat1 <- ezANOVA(data = region_df, dv = .(pielou), wid = .(sample), within = .(sorting_process), type = 3)
  print(repeat1)
}
