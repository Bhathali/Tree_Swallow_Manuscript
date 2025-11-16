```{r}
library("devtools") ## this needs to be installed and loaded prior to running
### the following:
# install_github("helixcn/seqRFLP")
# install_github("JCSzamosi/aftersl1p@*release")

## load packages:
library("phyloseq")
library("ShortRead") 
library("tidyverse") ## ggplot, dplyr, et c.
library("Polychrome") ## generate randomized discrete colour palettes
library("paletteer") ## generate colour palettes
library("readr") ## read_csv function
library("cluster")
library("data.table")
library("grid")
library("ape")
library("phangorn")
library("phytools")
library("vegan")
library("seqinr")
library("seqRFLP")
library("AfterSl1p")
```

```{r}
otufile <- ##Seqtab file

otu_df <- read.csv(otufile, row.names = 1)

seqs <- rownames(otu_df)

rownames(otu_df) = NULL

taxfile <- ##Taxa file

all(seqs == rownames(tax_df))

rownames(tax_df) = NULL

metfile <- ##Mapfile
map_df <- read.csv(metfile)

map_df$NAME <- gsub('-', '.', as.character(map_df$NAME))

map_df$NAME <- gsub(' ', '.', map_df$NAME)

map_df$NAME <- factor(map_df$NAME,
                      levels = c("EC1", "EC2","EC3",
                                 "EC4", "EC5", "EC6",
                                 "EC7", "EC8", "EC9", "EC10",  "EC11",  "EC12",  "EC13",  "EC14",  "EC15",  "EC16",  "EC17",  "EC18",  "EC19",  "EC20",  "EC21",  "EC22",  "EC23",  "EC24",  "EC25"))

rownames(map_df) <- map_df$NAME

rows <- sort(rownames(map_df)) 
cols <- sort(colnames(otu_df))

colnames(map_df) <- gsub(" ", ".", colnames(map_df))  # Replace spaces with dots
colnames(map_df) <- gsub("\\(|\\)", "", colnames(map_df))  # Remove parentheses
map_df$Weight.mg <- as.numeric(map_df$Weight..g.) # Ensuring weight is numbers
map_df$Wing.Chord.mm <- as.numeric(map_df$Wing.Chord..mm.) # Ensuring wing chord is numbers
map_df$Right.Tarsus.mm <- as.numeric(map_df$Right.Tarsus..mm.) # Ensuring right tarsus is numbers
map_df$X9th.Primary.mm <- as.numeric(map_df$X9th.Primary..mm.) # Ensuring 9th primary is numbers
map_df$Box.Temperature <- as.numeric(map_df$Box.Temperature) # Ensuring box temp is numbers
```

```{r}
colnames(map_df) <- gsub(" ", ".", colnames(map_df))  # Replace spaces with dots
colnames(map_df) <- gsub("\\(|\\)", "", colnames(map_df))  # Remove parentheses
map_df$Weight.mg <- as.numeric(map_df$Weight..g.) # Ensuring weight is numbers
map_df$Wing.Chord.mm <- as.numeric(map_df$Wing.Chord..mm.) # Ensuring wing chord is numbers
map_df$Right.Tarsus.mm <- as.numeric(map_df$Right.Tarsus..mm.) # Ensuring right tarsus is numbers
map_df$X9th.Primary.mm <- as.numeric(map_df$X9th.Primary..mm.) # Ensuring 9th primary is numbers
map_df$Box.Temperature <- as.numeric(map_df$Box.Temperature) # Ensuring box temp is numbers
```

```{r}
## Weight bar graph
ggplot(map_df, aes(x = Site, y = Weight.mg, fill = Site)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(title = "Average Weight by Site", x = "Site", y = "Weight (g)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Weight box plot — match BCI style (borders + no legend)
ggplot(map_df, aes(x = Site, y = Weight.mg, fill = Site)) +
  geom_boxplot(outlier.shape = NA, colour = "grey20") +
  labs(title = "Weight Distribution by Site",
       x = "Site", y = "Weight (g)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


## Manually assign Condition (MF & MD = Rural; all others = Waste)
map_df <- map_df %>%
  dplyr::mutate(
    Condition = ifelse(Site %in% c("MF","MD"), "Rural", "Waste"),
    Condition = factor(Condition, levels = c("Rural","Waste"))
  )

## Weight – compare by Condition (BCI-style, vertical y-axis title)
library(ggplot2)

ggplot(map_df, aes(x = Condition, y = Weight.mg, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, colour = "grey20") +
  geom_jitter(width = 0.10, alpha = 0.6, size = 2, colour = "grey40") +
  scale_fill_manual(values = c(Rural = "#2C6DA4",   # blue
                               Waste = "#8C1513")) + # deep red
  labs(title = "Nestling weight by site condition",
       x = NULL,
       y = "Weight (g)") +  # switch to (mg) if that’s your unit
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 0, hjust = 0.5),
    plot.title      = element_text(hjust = 0.5)
  )
```

```{r}
## Right Tarsus bar graph
ggplot(map_df, aes(x = Site, y = Right.Tarsus.mm, fill = Site)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(title = "Average Right Tarsus by Site", x = "Site", y = "Right Tarsus (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Right Tarsus box plot
ggplot(map_df, aes(x = Site, y = Right.Tarsus.mm, fill = Site)) +
  geom_boxplot() +
  labs(title = "Right Tarsus Distribution by Site", x = "Site", y = "Right Tarsus (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

```{r}
## 9th Primary bar graph
ggplot(map_df, aes(x = Site, y = X9th.Primary.mm, fill = Site)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(title = "Average 9th Primary by Site", x = "Site", y = "9th Primary (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## 9th Primary box plot
ggplot(map_df, aes(x = Site, y = X9th.Primary.mm, fill = Site)) +
  geom_boxplot() +
  labs(title = "9th Primary Distribution by Site", x = "Site", y = "9th Primary (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

```{r}
## Clutch Size bar graph
ggplot(map_df, aes(x = Site, y = Clutch.Size, fill = Site)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(title = "Average Clutch Size by Site", x = "Site", y = "Clutch Size") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Clutch Size box plot
ggplot(map_df, aes(x = Site, y = Clutch.Size, fill = Site)) +
  geom_boxplot() +
  labs(title = "Clutch Size by Site", x = "Site", y = "Clutch Size") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

```{r}
##Filtering
write.csv(seqs, file= ##ASV Sequences CSV File)

seqs.fasta <- dataframe2fas(seqs, file= ##ASV Sequences FASTA File)

dat <-  phyloseq(otu_table(otu_df, taxa_are_rows = TRUE), # or FALSE if false
               tax_table(as.matrix(tax_df)),
               sample_data(map_df))

samplesums <- sort(sample_sums(dat))

write.csv(samplesums, file=##Read counts CSV)

plot(samplesums)

dat_lessHOST <- subset_taxa(dat, Kingdom=="Bacteria", Family!="Mitochondria")
write.csv(otu_table(dat_lessHOST),file=## OTU table CSV)
write.csv(tax_table(dat_lessHOST),file=## Tax table CSV)

rel_abun_all <- transform_sample_counts(dat_lessHOST, function(x) x/sum(x))
rel_abun_all_prune <- prune_taxa(taxa_sums(rel_abun_all) > 0.001,
                                 rel_abun_all)                           
```

 
