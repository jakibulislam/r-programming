#packages----
install.packages('GEOquery')
install.packages('tidyverse')
install.packages('dplyr')
install.packages("ggplot2")
install.packages("BiocManager")
library(BiocGenerics)
library(Biobase)
BiocManager::install("DESeq2")
BiocManager::install("affy")
BiocManager::install("airway")
BiocManager::install("celldex")
BiocManager::install("SingleR")
BiocManager::install("Seurat")
BiocManager::install("pheatmap")
BiocManager::install("biomaRt")
BiocManager::install("annotables")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("EnsDb.Hsapiens.v86")
install.packages('survival')
install.packages('TCGAbiolinks')
install.packages('survminer')
install.packages('SummarizedExperiment')

#code for ".gz" file unzip----
#required 'GEOquery' package 
library(GEOquery)
library(dplyr)
library(tidyverse)

unzip("C:/Users/Asus/Downloads/GSE183947_fpkm.csv.gz")
aaa=gunzip("C:/Users/Asus/Downloads/GSE183947_fpkm.csv.gz")
aaa
setwd()
#-----
##Importing File-----
read.table(file.choose(), header =TRUE, sep = ',')


#direct import
dat=read.csv(file = "D:/Biochemistry and Biotechnology/Bioinformatics/
             GSE183947_fpkm.csv", header = TRUE)
dat
view(dat)

# get metadata --------
dim(geo) #check dimension 
library(GEOquery)
library(dplyr)
library(tidyverse)
gse = getGEO(GEO = 'GSE183947', GSEMatrix = TRUE) #GSE183947 is a geo_accession id of a paper
gse

#show more details about values from online [journal paper]
metadata <- pData(phenoData(gse[[2]]))
head(metadata)
view(metadata) #for easily visualize 


#filter/retrieve some data (column) using heading name or heading number
#required these functions "select, mutate, rename"
retrive= metadata %>%
  select(1,10,11,17) %>% #using "select" for selecting these column 
  head()
retrive

#rename the column name
  renamed=metadata %>%
  select (1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%  #characteristics_ch1 instead of tissue 
  rename(metastasis = characteristics_ch1.1) %>%  #characteristics_ch1 instead of metastasis
  head()
  renamed
  
#mute
muted=metadata %>%
  select (1,10,11,17) %>%
  mutate(tissue = gsub("tissue: ", "", characteristics_ch1)) %>%  #using function 'gsub' for update value
  mutate(metastasis = gsub("metastasis: ", "", characteristics_ch1.1))%>%
  head()
muted

#replacing/delete common value of a whole column 
muted=metadata %>%
  select (1,10,11,17) %>%
  mutate(tissue = gsub("tissue: ", "", characteristics_ch1)) %>%  #using function 'gsub' for update, delete of values
  mutate(metastasis = gsub("metastasis: ", " meta ", characteristics_ch1.1))%>%  #replacing meta to metastasis
  head()
muted

#select, mute and renamed  
  rename_muted=metadata %>%
  select (1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%  
  rename(metastasis = characteristics_ch1.1) %>%  
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))%>%
  head()
 view(rename_muted) 


# looking at gene expression data
  head(dat)
  
# reshaping data - from wide to long--------
  dat%>%
  rename(gene = Gene) %>%
    gather(key = 'samples', value = 'FPKM', -gene) %>% 
    #samples mean column name  (sample) and FPKM mean values
  head()
  
  dat.long <- dat %>% #changing structure dataframe
    rename(gene = Gene) %>%
    gather(key = 'samples', value = 'FPKM', -gene)
  
  view(dat.long)


# join dataframes = dat.long + rename_muted
  
  dat.long=dat.long %>%
    left_join(., rename_muted, by = c("samples" = "description")) 
  #(.) mean dat.long
head()
  
# explore data ------
# filter, group_by, summarize and arrange 
dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') 
head()


dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) 
head()

#summary
dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarize(FPKM) %>%
  head()  

#mean and median
dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarize(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM)) %>%
  arrange(-mean_FPKM)  #arrange for shorting and (-) for descending 
head()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    

#ggplot----
library(ggplot2)
library(tidyverse)

#filtering
dat.long%>%
filter(gene=='BRCA1')%>%
head()

dat.long%>%
  filter(gene=='BRCA1')%>%
#(.) mean using above data
ggplot(., aes(x= FPKM, fill=tissue))+geom_density()
#adjust obesity
dat.long%>%
  filter(gene=='BRCA1')%>%
#using "alpha" for adjusting obesity and "col" for outline color
  ggplot(., aes(x= FPKM, fill=tissue))+geom_density(alpha=0.5, col='red')

-------------------------------------------------------
#filtering
dat.long%>%
  filter(gene=='BRCA1' |gene=='BRCA2')%>%
  head()

dat.long%>%
  filter(gene=='BRCA1' |gene=='BRCA2')%>%
  spread(key=gene, value = FPKM) %>%
#row er differrent value ke colum hisabe show korbe [column er value hobe "FPKM" er man]
head()

#scatter_plot

abc= filter(dat.long, gene=='BRCA1' |gene=='BRCA2')
abc
dat.long%>%
  filter(gene=='BRCA1' |gene=='BRCA2')%>%
  spread(key=gene, value = FPKM)%>%
  ggplot(., aes(x= BRCA1, y= BRCA2, color= tissue))+
  geom_point()

dat.long%>%
  filter(gene=='BRCA1' |gene=='BRCA2')%>%
  spread(key=gene, value = FPKM)%>%
  ggplot(., aes(x= BRCA1, y= BRCA2, col= tissue))+
  geom_point()+ geom_smooth(method = lm, se=F)

dat.long%>%
  filter(gene=='BRCA1' |gene=='BRCA2')%>%
  spread(key=gene, value = FPKM)%>%
  ggplot(., aes(x= BRCA1, y= BRCA2, col= tissue))+
  geom_point()+ geom_smooth(method = lm, se=F)+ facet_grid()
-------------------------------------------------------
  dat.long%>%
  filter(gene=='BRCA1' |gene=='BRCA2')%>%
  spread(key=gene, value = FPKM)%>%
  ggplot(., aes(x= BRCA1, y= BRCA2, col= tissue, shape=metastasis))+
  geom_point()

scatter_line=dat.long%>%
  filter(gene=='BRCA1' |gene=='BRCA2')%>%
  spread(key=gene, value = FPKM)%>%
  ggplot(., aes(x= BRCA1, y= BRCA2, color= tissue))+
  geom_point()+geom_line()
scatter_line
-------------------------------------------------------
#heatmap
 genes.of.interest=c('BRCA1', 'BRCA2', 'TNMD', 'FGR','CFH', 'NFYA') 
dat.long%>%
  filter(gene %in% genes.of.interest )%>%
  ggplot(., aes(x= samples, y= gene, fill=FPKM))+geom_tile()
#color
genes.of.interest=c('BRCA1', 'BRCA2', 'TNMD', 'FGR','CFH', 'NFYA') 
dat.long%>%
  filter(gene %in% genes.of.interest )%>%
  ggplot(., aes(x= samples, y= gene, fill=FPKM))+geom_tile()+
  scale_fill_gradient(low = 'green', high = 2)
-------------------------------------------------------
#save file----
#scatter_line is plot name, scatter_line2.pdf is file name,  
#width and height is your choice [na dileo problem nai] 
#pdf, png, jpg jkono formate a save kora jabe

ggsave(scatter_line, filename='scatter_line2.jpg', width = 8, height = 8)

#7.DESeq2_workflow_Differential_Gene_Expression_Analysis----
install.packages("BiocManager")
BiocManager::install("DESeq2")
library(BiocGenerics)
library(Biobase)
BiocManager::install("airway")
library(DESeq2)
library(tidyverse)
library(airway)
#read.table(file.choose(), sep = ',', header = T)
counts_data <- read.csv(file = 'D:/R Programming/Bioinformatics_RNAseq
                        /Files/counts_data.csv', header = TRUE)

View(counts_data)

colData=read.table(file.choose(), sep = ',', header = T)
View(colData)
# Step 1: preparing count data 
# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData)) #%in% mean present in
#result= TRUE

# are they in the same order?
all(colnames(counts_data) == rownames(colData))
#result= TRUE

# Step 2: construct a DESeqDataSet object 
DESeqDataSet=DESeqDataSetFromMatrix(countData = counts_data,
                                    colData = colData,
                                    design = ~ dexamethasone)
#dexamethasone is a funcgtion
DESeqDataSet
#result=
#class: DESeqDataSet 
#dim: 64102 8 
#metadata(1): version
#assays(1): counts
#rownames(64102): ENSG00000000003 ENSG00000000005 ... LRG_98 LRG_99
#rowData names(0):
#  colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#colData names(2): cellLine dexamethasone
View(DESeqDataSet) #for more visualized run both functions


# keeping rows that have at least 10 reads total
read_more10 <- rowSums(counts(DESeqDataSet)) >= 10 #rowSums is a function
read_more10
View(read_more10)
#filtering: removing rows with low gene counts
dds <- DESeqDataSet[read_more10,]

dds
#result=
#class: DESeqDataSet 
#dim: 22369 8 [64102 rows]

# set the factor level
dds$dexamethasone= relevel(dds$dexamethasone, ref = 'untreated')
dds$dexamethasone
#result= Levels: untreated treated

# Step 3: Run DESeq ( collapse technical replicates)
dds= DESeq(dds)
#result= estimating size factors, estimating dispersions, gene-wise dispersion estimates
#mean-dispersion relationship, final dispersion estimates, fitting model and testing

# NOTE: collapse technical replicates
rst=results(dds)
rst

?results
#results extracts a result table from a DESeq analysis giving base means across samples, 
#log2 fold changes, standard errors, test statistics, p-values and adjusted p-values; 
#resultsNames returns the names of the estimated effects (coefficents) of the model;
#removeResults returns a DESeqDataSet object with results columns removed.

# Explore Results 
summary(rst)
#out of 22369 with nonzero total read count; adjusted p-value < 0.1
#LFC > 0 (up)       : 1884, 8.4%; LFC < 0 (down)     : 1502, 6.7%
#outliers [1]       : 51, 0.23%; low counts [2]     : 3903, 17%; #(mean count < 4)

rst0.01 <- results(dds, alpha = 0.01) #for changing pvalue
summary(rst0.01)
#out of 22369 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)[up regulator]  : 1030, 4.6%;  LFC < 0 (down)[down regulator]  : 708, 3.2%
#outliers [1]       : 51, 0.23%  low counts [2]     : 5200, 23%
#(mean count < 6)

View(dds)
resultsNames(dds)
# e.g.: treated_4hrs, treated_8hrs, untreated

results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))


#  dexamethasone is design factor and treated 4hr, untreated are compare level
results(dds, contrast = c('dexamethasone', 'treated 4hr', 'untreated'))










#8.How_to_read_and_normalize_microarray_data_in_R----
library(tidyverse)
library(affy)
library(GEOquery)

#get suplimentary file
getGEOSuppFiles("GSE148537")

#uzip tar file
untar('GSE148537/GSE148537_RAW.tar', exdir = 'data/') #exdir=extract direction

#reading .cel file
raw_data= ReadAffy(celfile.path = 'data/')

#performing RMA normalization
normalized_data= rma(raw_data)
normalized_data

#view/expression
view(normalized_data)
normalized_data_expression= exprs(normalized_data)
normalized_data_expression= as.data.frame(exprs(normalized_data))
normalized_data_expression


