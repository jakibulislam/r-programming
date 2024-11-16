
#packages----
install.packages("dplyr") #data manipulation
install.packages("plotrix") #plot manipulation
install.packages("GEOquery")#get data from NCBI Eene Expression Omnibus (GEO)
install.packages("tidyverse")
install.packages('car')
install.packages()
install.packages()
install.packages()
install.packages()
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

#save file----
#scatter_line is plot name, scatter_line2.pdf is file name,  
#width and height is your choice [na dileo problem nai] 
#pdf, png, jpg jkono formate a save kora jabe

ggsave(scatter_line, filename='scatter_line2.jpg', width = 8, height = 8)
#scatter_line [vector of code ]
#data assigning-----
x= c(5,6,9,8)
x
#5 6 9 8
y=x[-2]#2 no value is removed
y #5 9 8
v=c(1:15)#1,2,3....15
v#1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
v=c(6:15)
v[2] #find the specific value
v#7
v[c(8,10)] #show 8-10 no values
#13 15
v[-c(5:8)] #show without 8-13 no values
#6  7  8  9 **14 15
v[-5] #5 no value is removed  
#6  7  8  9 (10) 11 12 13 14 15
v[c(-2,-6)] #2 and 6 no values are removed 

v[-c(8:10)] #8-10 no values are removed

summary(v)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#6.00    8.25   10.50   10.50   12.75   15.00 
#matrix----
a=matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)
a
#      [,1] [,2]
#[1,]    1    4
#[2,]    2    5
#[3,]    3    6
b= matrix(c(16:35), nrow = 4, ncol = 5)
b#     [,1] [,2] [,3] [,4] [,5]
#[1,]   16   20   24   28   32
#[2,]   17   21   25   29   33
#[3,]   18   22   26   30   34
#[4,]   19   23   27   31   35

c=matrix(c("apple", 'mango', 'water malon', "banana", "cherry", "orange"), nrow = 2, ncol = 3)
c#     [,1]    [,2]          [,3]    
#[1,] "apple" "water malon" "cherry"
#[2,] "mango" "banana"      "orange"
c[2,] #2nd row  "mango"  "banana" "orange"
c[ ,1] #1st column   "apple" "mango"
c[2,3] #2nd row, 3rd column #specific Value  "orange"
c[5] #5 no specific Value   "cherry"



#Data frame #all data must be same length------ 

name= c('jakib', 'sadia', 'mim', 'jahid', 'bgd', 'sick')
fav.food=c("apple", 'mango', 'water malon', "banana", "cherry", "orange")
Pulse = c(100, 150, 120,90,110,120)
data.frame(name, fav.food, Pulse)
#   name    fav.food Pulse
1 jakib       apple   100
2 sadia       mango   150
3   mim water malon   120
4 jahid      banana    90
5   bgd      cherry   110
6  sick      orange   120

#another option
data.frame(
  name= c('jakib', 'sadia', 'mim', 'jahid', 'bgd', 'sick'),
  fav.food=c("apple", 'mango', 'water malon', "banana", "cherry", "orange"),
  Pulse = c(100, 150, 120,90,110,120)
  )
#   name    fav.food Pulse
1 jakib       apple   100
2 sadia       mango   150
3   mim  water malon  120
4 jahid      banana    90
5   bgd      cherry   110
6  sick      orange   120

Data_Frame <- data.frame (
  Training = c("Strength", "Stamina", "Other"),
  Pulse = c(100, 150, 120),
  Duration = c(60, 30, 45)
)

Data_Frame
#Training Pulse Duration
1 Strength   100       60
2  Stamina   150       30
3    Other   120       45

#summary (min, 1st quarter, median, mean, 3rd qurter, max)
summary(Data_Frame)
#Training             Pulse          Duration   
Length:3           Min.   :100.0   Min.   :30.0  
Class :character   1st Qu.:110.0   1st Qu.:37.5  
Mode  :character   Median :120.0   Median :45.0  
                   Mean   :123.3   Mean   :45.0  
                   3rd Qu.:135.0   3rd Qu.:52.5  
                   Max.   :150.0   Max.   :60.0  

#----
#array-----
Values <- c(1:24)

# An array with more than one dimension
multiarray <- array(Values, dim = c(4, 2,3 ))
multiarray

# ,,1  [,1] [,2]  , , 2 [,1] [,2]   , , 3  [,1] [,2] 
[1,]    1    5           9   13             17   21
[2,]    2    6          10   14             18   22
[3,]    3    7          11   15             19   23
[4,]    4    8          12   16             20   24


#Tutorial_1.6_Data_importing_and_checking----
#Direct csv file import
BG= read.csv(file = "D:/প্রয়োজনীয় কাগজপত্র/blood_group.csv")
BG
aa= read.csv(file = "D:/প্রয়োজনীয় কাগজপত্র/blood_group.csv", header = T)#header for heading
aa
#csv file import through choice 
bb= read.table(file.choose(), header = TRUE, sep = ',')
bb
View(bb)
#for text file
dat1=read.table(file.choose(), header = TRUE, sep = '\t')
dat1
dat2=read.delim(file.choose(), header = TRUE)
dat2

#checking/filter by some conditions
#checking variable name with 1st** 6 values
head(bb) #head(data name)
#checking variable name with last** 6 values
tail(bb) #tail(data name)

bb[2,] #2nd row
bb[ ,1] #1st column 
bb[2,3] #2nd row, 3rd column #specific Value
bb[5] #5 no specific Value
bb[-(2:4)] #show all columns without 2-4 no columns/delete 2-4 no columns
sub_b=bb[-c(2,3,4)] #show all columns without 2-4 no columns/delete 2-4 no columns
View(sub_b)
bb[ ,(1:5)] #show 1-5 no columns
bb[ ,c(1,2,3,4,5)] #show 1-5 no columns


#select a single column by $ sign [variable/table_name$column_name]
bb$Age #select a single column by $ sign [variable_name$column_name]

#check data type of each column
str(bb) #check data type of each column

#Tutorial_1.7_Data manipulation in R using dplyr----
#require 'dplyr' package 
#using inner data [starwars]
library(dplyr)
View(Starwars)
attach(starwars)
#subset [column contain all values]
stwr1=select(starwars,name, sex,mass, height, species)
View(stwr1)
starwares
stwr2=select(starwars,name, mass:species)#mass to species
View(stwr2)

stwr3=select(starwars, -mass, -height)#removed mass and height
View(stwr3)

#filtering [R is a case sensitive program so must be write column name is same to table ]
#column contain values by conditions
only_male=filter(starwars,sex=="male")
View(only_male)
only_Human=filter(stwr1,species=="Human")
View(only_Human)

#male_human
male_human=filter(starwars,
            sex=="male" &
              species=="Human")
View(male_human)

#female_human
female_human=filter(starwars,
                  sex=="female" &
                    species=="Human")
View(female_human)

species=filter(starwars,
            species=="Droid" &
              species=="Human")
View(species)
#OR
species=filter(starwars,
            species%in% c("Droid","Human")) 
View(species)

#male Droid and Human
male_species=filter(starwars,
            sex=="male" &
              species%in% c("Droid","Human"))
View(male_species)


#updating data by common operation
mass_height= mutate(starwars,
             height= height * 0.394,
             mass= mass * 2.205)
View(mass_height)

#ifelse
catagory= mutate(starwars,
             catagory  = ifelse(height180,
                                "Tall",
                                "Short"))
View(catagory)

#
height= mutate(starwars,
               height  = ifelse(height<75 | height200,
                                    'NA',
                                    height))
View(height)

height= mutate(starwars,
               height  = ifelse(height<75 | height200,
                                NA,
                                height))
View(height)
#
eye_color= mutate(starwars, #%in% = present in
                  Eye_color  = ifelse(eye_color%in% c('black', 'brown', 'blue'),
                                    "Noraml",
                                    "Different"))
View(eye_color)
#BMI
BMI=mutate(starwars, BMI=mass/(height/100)^2)
View(BMI)
#
Body=mutate(BMI, 
            Body  = ifelse(BMI18 | BMI<30,
                                    'Normal',
                                    'Thin/Thick'))
View(Body)

#identify the averages of height and mass
#average
me=female_human[c(1:7), ]
View(me)
mean(me$height) #NA values is disgusting

#average [gender wise mean of mass and height] (this is easier way)
average= group_by(only_Human, sex)
av_ht_wt= summarize(average, 
             mean_hieght=mean(height,na.rm=TRUE),
             mean_weight=mean(mass, na.rm=TRUE))
av_ht_wt

#store data(modified) from Rstudio----
write.csv(bb, file ='big_data.csv', row.names = FALSE )
write.csv(sub_b, "sub_lungcap.csv")

#Tutorial_1.8_Recoding_and_grouping_data_in_R_using_Dplyr-----
 human=na.omit(only_Human)
View(human)
#BMI
BMI=mutate(human, BMI=mass/(height/100)^2)
View(BMI)
#Body_condition
#
BMI$Body_condition=cut(BMI$BMI,
                       breaks = c(15,18.5,25,30,40),
                       right = FALSE)

BMI$Body_condition=cut(BMI$BMI,
                       breaks = c(15,18.5,25,30,40),
                       labels = c('Under weight', 'Normal', 'Over weight', 'Obese')
                       )
#
BMI$Body_condition=cut(BMI$BMI,
                       breaks = c(15,18.5,25,30,40),
                       labels = c('Under weight', 'Normal', 'Over weight', 'Obese'),
                       right = FALSE)
#grouping by conditions
BMI=within(BMI,{Height=NA #create a column "Height" and store in the column
Height[height<171]="Short"
Height[height=172 & height<=185]="Normal"
Height[height185]="tall"})

#Tutorial_1.9,10,11,12,13,14_Create_and_customize_histogram,bar_chart, 
    #boxplot, pie_chart, steam and leaf plot and "ggplot"----
install.packages(plotrix)
library(plotrix) #for 3d view
install.packages(ggplot2) #ggplot=grammar of graphical plot
library(ggplot2)


data1=read.table(file.choose(), header = TRUE, sep = ',')
lung=data1$LungCap
View(iris)#in built data
attach(iris)

 hist(lung)

 hist(lung, freq=FALSE) #using 'freq' for visualized in %
hist(lung, pro=T) #pro mean probability 
 hist(lung, pro=T, main="Histogram of lung_capacity") 
 #using 'main' for customizing the title


 ------------------------------------------------------- 
#ggplot_histogram
ggplot(data=iris,
       aes(Sepal.Width))+geom_histogram(bins =15, fill='green', col='red')

#fill=inner color, col=outer line, bins=number of column
ggplot(data=iris,
       aes(Sepal.Width))+geom_histogram(bins =15, fill='green', col='red')

#column filled by species [valvue of column is case sensitive]
ggplot(data=iris,
       aes(x=Sepal.Width, fill= Species))+geom_histogram(bins =15,col='red')

#data from lungcap
ggplot(data=data1,
       aes(x=LungCap))+geom_histogram(col='red')

#frequency polygon
ggplot(data=data1,
       aes(x=LungCap))+geom_freqpoly()
ggplot(data=data1,
       aes(x=LungCap, col=Gender))+geom_freqpoly(bins=30)

ggplot(data=data1,
       aes(x=LungCap, col=Smoke))+geom_freqpoly(bins=30)

-------------------------------------------------------
lines (density(lung))
 lines (density(lung),col=2,lwd=3) #lwd mean line_wide

table(data1$Gender) #frequency tab;e of gender
#data1$Gender mean select the gender column of data1 table
#easily selecting for any column 'attach' the table
attach(data1)
percent=table(Gender)/725*100
percent
bplot=barplot(percent) #value must be contain numerical value
 
barplot(percent, main="Barplot", xlab="Gender", ylab ="Percntage",names.arg=c("Female", "Male"))
#xlab mean X_axis and ylab mean Y_axis


ggplot(data=data1,
       aes(x=Smoke, fill= Gender))+geom_bar()
#data showed bt percentage 
ggplot(data=data1,
       aes(x=Smoke, fill= Gender))+geom_bar(position = 'fill')


-------------------------------------------------------
#scatter plot
table (Smoke, Gender) #cross table
table (Smoke, Gender, Caesarean) #cross table

#scatter plot [2ta catagory aksathe compare kore dekhar jonno use kora hoy]
table1<-table(Smoke, Gender) #visualized with cross table
barplot(table1)
barplot(table1, beside=T) #more visualization
barplot(table1, beside=T, legend.text=T)

barplot(table1, beside=T, legend.text=c("Non_somker", "Smoker"))
#using legend for indicating the bar with color
barplot(table1, beside=T,legend.text=c("Non-smoker", "Smoker"), 
        main="Gender and Smoking", xlab='Gender', col = c(2,4))
quantile(lung) #0,25,50,75 and 100% a average 
summary(data1) #work when values are numerical


-------------------------------------------------------
#boxplot
boxplot(lung) #(saw Min., 1st Qu., Median,3rd Qu., Max. values)
boxplot(lung, main="Boxplot", ylab="Lung Capacity")
boxplot(lung~data1$Gender) #scatter plot
boxplot(lung~data1$Gender, main="Boxplot", xlab="Gender", ylab="Lung Capaciy")
boxplot(lung~data1$Gender, main="Boxplot", xlab ="Gender",ylab="Lung Capaciy", col=c(2,4))


#ggplot with boxplot

ggplot(data=data1,
       aes(y=LungCap, col=Gender))+geom_boxplot()


ggplot(data=data1,
       aes(y=LungCap, x=Smoke, col=Gender ))+geom_boxplot()


-------------------------------------------------------
#pie_chart
smoke<-c(180,190,185,245)
pie(smoke)
#labels
pie (smoke, main="Simple Pie Chart")#"main" for heading
lbls=c("Chi", "Dha", "Khu", "Raj") #lbs=level=indicator
pie (smoke, labels=lbls, main="Simple Pie Chart")

#percentage
pct=round(smoke/sum(smoke)*100) #percentage 
lbls=paste (lbls," ",pct,"%", sep="") 
pie (smoke, labels=lbls, main="Simple Pie Chart")
#rainbow color
pie (smoke, labels=lbls, main="Simple Pie Chart",col=rainbow (4))
#control shape
pie (smoke, labels=lb1s, main="Simple Pie Chart",col=rainbow(4), radius = 1.5)

-------------------------------------------------------
#plotrix
library(plotrix)
#3d pie
pie3D (smoke, labels=lbls, main="Simple Pie Chart",col=rainbow (4))
pie3D (smoke, labels=lbls, main="Simple Pie Chart",col=rainbow (4),
       redious=2, explode=0.1) 
#"redious" for controlling size, 'explode' for inserting gap among several part

-------------------------------------------------------
#stem and leaf
stem(data1$LungCap) 
#left side express deca number and right side express single number


-------------------------------------------------------
#scatter plot
plot(Sepal.Length, Petal.Length)
#"pch" used for cganging shape of dot
plot(Sepal.Length, Petal.Length, main = 'Sepal Length vs Petal Length',
     ylab = 'Petal Length', xlab = 'Sepal Length', col="red", pch=2)
#ggplot
#aes mean axis, geom_point mean scatter plot
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length))+geom_point()

#colored by Species 
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length, col=Species))+geom_point()

#changing shape according to species
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length, col=Species, 
           shape=Species))+geom_point()

-------------------------------------------------------
#smooth plot
#show smooth line with  95% confidence intervel
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length))+geom_smooth()

#se=f means remove standard error (95%)
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length))+geom_smooth(se=F)

#visualized by Species
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length, col=Species))+geom_smooth(se=F)


#Scatter plot with smooth plot [loess method]
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length, col=Species))+
  geom_point()+ geom_smooth(method = lm, se=F)


#facteting
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length, col=Species))+
  geom_point()+ geom_smooth(method = lm, se=F)+ facet_grid()

#specieswise vizualization
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length, col=Species))+
  geom_point()+ geom_smooth(method = lm, se=F)+ facet_grid(Species)

-------------------------------------------------------
#theme----
obj1= ggplot(data=data1,
         aes(y=LungCap, col=Gender))+geom_boxplot()
#changing title
obj2= obj1+ labs(title = 'Boxplot of lung capacity', fill= Gender)
obj2
#changing background
obj3=obj2+ theme(panel.background = element_rect('lightgreen'))
obj3
#changing title shape and color 
#[plot.title=  title of plot, element_text= controlling text 
#(hjust=0.5 mean title ta maje thakbe), face= text type]
obj4=obj3+theme(plot.title = element_text(hjust = 0.5, face = 'bold', color = 'blue'))
obj4

#Tutorial_1.14_Data_visualization_in_R_with_ggplot2----
install.packages(ggplot2) #ggplot=grammar of graphical plot
library(ggplot2)
attach(data1)
attach(iris)

-------------------------------------------------------------

#ggplot_histogram
ggplot(data=iris,
       aes(Sepal.Width))+geom_histogram(bins =15, fill='green', col='red')
#fill=inner color, col=outer line, bins=number of column
ggplot(data=iris,
       aes(Sepal.Width))+geom_histogram(bins =15, fill='green', col='red')
#column filled by species [valvue of column is case sensitive]
ggplot(data=iris,
       aes(x=Sepal.Width, fill= Species))+geom_histogram(bins =15,col='red')

#data from lungcap
ggplot(data=data1,
       aes(x=LungCap))+geom_histogram(col='red')

#frequency polygon
ggplot(data=data1,
       aes(x=LungCap))+geom_freqpoly()
ggplot(data=data1,
       aes(x=LungCap, col=Gender))+geom_freqpoly(bins=30)

ggplot(data=data1,
       aes(x=LungCap, col=Smoke))+geom_freqpoly(bins=30)



ggplot(data=data1,
       aes(x=Smoke, fill= Gender))+geom_bar()

#data showed bt percentage 
ggplot(data=data1,
       aes(x=Smoke, fill= Gender))+geom_bar(position = 'fill')



-------------------------------------
#ggplot with boxplot

ggplot(data=data1,
       aes(y=LungCap, col=Gender))+geom_boxplot()


ggplot(data=data1,
       aes(y=LungCap, x=Smoke, col=Gender))+geom_boxplot()


-------------------------------------------------------
#scatter plot
plot(Sepal.Length, Petal.Length)
#"pch" used for cganging shape of dot
plot(Sepal.Length, Petal.Length, main = 'Sepal Length vs Petal Length',
     ylab = 'Petal Length', xlab = 'Sepal Length', col="red", pch=2)
#ggplot
#aes mean axis, geom_point mean scatter plot
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length))+geom_point()
#colored by Species 
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length, col=Species))+geom_point()
#changing shape according to species
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length, col=Species, 
           shape=Species))+geom_point()

-------------------------------------------------------
#smooth plot
#show smooth line with  95% confidence intervel
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length))+geom_smooth()
#se=f means remove standard error (95%)
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length))+geom_smooth(se=F)
#visualized by Species
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length, col=Species))+geom_smooth(se=F)

-------------------------------------------------------
#Scatter plot with smooth plot [loess method]
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length, col=Species))+
  geom_point()+ geom_smooth(method = lm, se=F)

-------------------------------------------------------
#faceting
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length, col=Species))+
  geom_point()+ geom_smooth(method = lm, se=F)+ facet_grid()

-------------------------------------------------------
#specieswise vizualization
ggplot(data = iris,
       aes(x=Sepal.Length, Petal.Length, col=Species))+
  geom_point()+ geom_smooth(method = lm, se=F)+ facet_grid(Species)

-------------------------------------------------------
#theme
obj1= ggplot(data=data1,
             aes(y=LungCap, col=Gender))+geom_boxplot()
#changing title
obj2= obj1+ labs(title = 'Boxplot of lung capacity', fill= Gender)
obj2
#changing background
obj3=obj2+ theme(panel.background = element_rect('lightgreen'))
obj3
#changing title shape and color 
#[plot.title=  title of plot, element_text= controlling text 
#(hjust=0.5 mean title ta maje thakbe), face= text type]
obj4=obj3+theme(plot.title = element_text(hjust = 0.5, face = 'bold', color = 'blue'))
obj4

#Tutorial_1.15_Create_and_Customize_Scatterplot----
View(cars)#in build data
attach(cars)
plot(speed, dist) # scatter plot

plot(speed, dist, ylab='distance') #change text in axis

#Changing colour in the plot using 'col'
plot(speed, dist, main="scatterplot")
plot(speed, dist, main="scatterplot", col=5)
#Changing colour in the plot using 'col'
plot(speed, dist, main="scatterplot")
plot (speed, dist, main="scatterplot", col=5)
plot(speed,dist, main="scatterplot",col=5)
plot(speed, dist, main="scatterplot",col=5,col.main=- 4)
plot(speed, dist, main="scatterplot")
plot(speed, dist, main="scatterplot",col=5)
plot(speed, dist, main="scatterplot",col=5,col.main= 4)
plot(speed, dist, main="scatterplot",col=5,col.main= 4,col.lab=2)
plot(speed, dist, main="scatterplot",col=5,col.main= 4,col.lab=2,col.axis=3)
plot(speed, dist, main="scatterplot")
plot(speed, dist, main="scatterplot", pch=2)
plot(speed, dist, main="scatterplot", pch=3)  plot(speed, dist, main="scatterplot", pch=4)
plot(speed,dist, main="scatterplot", pch=2)
plot(speed, dist, main="scatterplot")
plot(speed, dist, main="scatterplot", pch=2)
plot(speed, dist, main="scatterplot", pch=3)
plot(speed, dist, main="scatterplot", pch=4)
plot(speed, dist, main="scatterplot", pch="m")
plot (speed, dist, main="scatterplot")
abline (lm(dist-speed))
abline(lm(dist-speed),col=6)
abline (lm(dist-speed),col=6,lwd=3)
plot (speed,dist, main="scatterplot")
plot(speed, dist, main="scatterplot", axes=F)
axis (side=1,at=c(10,15,20))
axis(side=2,at=c(20,70,100))
box()
#Tutorial_1.16_Line_plot_and_Time_series_plot----
#install.packages("dplyr")
library(dplyr)
library(ggplot2)
library(lubridate)

rainfall=read.table(file.choose(), header = TRUE, sep = ',')
#data_monthly_rainfall.csv
View(rainfall)
Barisal<-filter (rainfall, Station=="Barisal")
View(Barisal)
Barisal_16<-filter (rainfall, Station=="Barisal", Year=='2016' )
View(Barisal_16)
attach (Barisal)
plot (Rainfall)
plot(Rainfall, type="l")#type=plot type, l=line
plot (Rainfall,type="l",lwd=2,col="red")

#insert a table using value of 2/more column 
Bsl_dt<-mutate(Barisal,date=make_datetime (Year,Month))
View(Bsl_dt)
plot(Bsl_dt$date, Rainfall,type="l",lwd=2,col="red")
plot(Bsl_dt$date, Rainfall,type="l",lwd=3,col="red",
     xlab="Time", ylab="Monthly Rainfall", main="Barisal")

#ggplot
ggplot(Bsl_dt, aes (date, Rainfall))+geom_line()
ggplot(Bsl_dt, aes (date, Rainfall))+geom_line(colour="red")
plot1<-ggplot(Bsl_dt,aes(date, Rainfall))+geom_line(colour="red")
plot1+ggtitle('Monthly Rainfall at Barishal')+ labs(x='Year')

Stations<-filter (rainfall, Station=="Chittagong"| Station=="Dhaka"|
 Station=="Khulna"| Station=="Rajshahi"| Station=="Sylhet" |Station=='Barishal')
View(Stations)
st2=filter(Stations, Year>='2000')
View(st2)

#insert a table using value of 2/more column 
st2_dt= mutate(st2, date=make_date(Year, Month))
View(st2_dt)

ggplot(st2_dt, aes (date, Rainfall))+geom_line(lwd=1) 
ggplot(st2_dt, aes (date, Rainfall,color=Station))+geom_line()
plot2<-ggplot(st2_dt, aes (date, Rainfall,color=Station))+geom_line(lwd=1)
rf_st=plot2+ggtitle("Monthly Rainfall")+labs (x="Date")
rf_st

#save file
ggsave(rf_st, filename='rainfall of several devision.jpg', width = 20, height = 8)
ggsave(rf_st, filename='rainfall of several devision.pdf', width = 20, height = 8)











































































#Tutorial_2.1_Descriptive_statistics__mean_median_variance_sd----
lung=read.table(file.choose(), header = TRUE, sep = ',')
View(lung)
attach(lung)
min(LungCap)#minimum
mean(LungCap)#mean/average
max(LungCap)#maximum
range(LungCap)
sum(LungCap)
sd(LungCap)#standard deviation
quantile(LungCap)#0,25,50,75,100 % a value
var(LungCap)#variation
median(LungCap)#je value ta besibar ase
mad(LungCap)
summary(LungCap)
summary(Smoke)#value character a ase tai read korte  pareni
smoke=factor(Smoke)#character value ke read korte factor use kora hoy
summary(smoke)
levels(smoke)#smoke column a je doroner value royese
percentage_of_smoker=summary(smoke)/725*100#725 ta sample ase
percentage_of_smoker
names(lung)#column gulir nam dekhte chaile
class(Gender)#show data_type

#Tutorial_2.2_Single_mean_test_and_equality_of_two_means_test,t-test----
#t test(assumption test) and leven's test(homojenesity test)

install.packages('car')
library(car)
View(data1)
attach(data1)
lung=data1$LungCap
#t_test [t test akta assumption, kasakasi valuer idea dey but fixed na]
#HO:mu=8 and HA:mu<8 [HO=hypothesis, HA= null hypothesis] [mu= mean ]
t.test(lung,mu=8, alternative = "less", conf.level = 0.95)
#t = -1.3842, df [degree of freedom] 724, p-value = 0.08336 [>0.05 (accepted)]
t.test(lung,mu=8, alt = "greater", conf = 0.95)#[conf=confidence, alt=alternative]
#p-value = 0.9166 [>>0.05(accepted)]
#two sided
t.test(lung,mu=8, alt = "two.sided", conf = 0.95)
#p-value = 0.1667 [>0.05(accepted)]

#jsob value amra check korte pari
tst=t.test(lung,mu=8, alt = "two.sided", conf = 0.95)
attributes(tst)
#$names
#[1] "statistic"   "parameter"   "p.value"     "conf.int"    "estimate"   
#[6] "null.value"  "stderr"      "alternative" "method"      "data.name"  
#$class [1] "htest"
#check
tst$conf.int  #7.669052 8.057243 attr(,"conf.level") 0.95

#equality of two mean test
#lungCap of male=female (mu=0 mean mu er difference 0) but variance male=/female 
t.test(LungCap~Gender , mu=0, alt = 'two.sided', conf.level = 0.95, 
       var.eq = F ) #var.equal = F [variance male=/female ]
#result
#LungCap by Gender
#t = -4.6362, df = 722.7, p-value = 4.21e-06 so alternative hypothesis: true difference 
#mean in group female(7.405746 )   mean in group male (  8.309332 )
          
t.test(LungCap~Gender) #****short form


#leven's test(homogenesity test)
gender=factor(Gender)
leveneTest(LungCap~gender)
#result: Df F value Pr(>F)
#group   1  0.6757 0.4114 [variance are equal]

bp=read.table(file.choose(),sep = ',', header = T)
bp
attach(bp)
View(bp)
#medication is effective for lowering SBP
#HO: mean diff of SBP is 0; one sided test (greater)
t.test(Before, After, mu=0, alternative = 'greater', paired = T, conf.level = 0.99)
#p-value = 0.0003493 [<.01 (rejected)]
t.test(Before, After, mu=0, alternative = 'two.sided', paired = T, conf.level = 0.99)
#t = 3.8882, df = 24, p-value = 0.0006986[rejected]
#99 percent confidence interval: 2.245279 13.754721 & mean difference is 8 
***
t.test(Before, After, mu=8, alternative = 'two.sided', paired = T, conf.level = 0.99)
#p-value = 1

#Tutorial_2.3_Equality_of_multiple_means_test_in_RStudio__ANOVA & Kruskal_Wallis----
wt=read.table(file.choose(), sep = ',', header = T)
View(wt)
attach(wt)
#HO: weightloss is same for all Diets
#HA: weightloss is not same for all Diets
#Anova test
diet=factor(Diet)
class(diet)
aov(WeightLoss~diet)
#diet Residuals
#Sum of Squares   97.32983 296.98667
#Deg. of Freedom         3        56
av1=aov(WeightLoss~diet)
summary(av1)
#             Df Sum Sq   Mean Sq   F value  Pr(>F)   
#diet         3  97.33    32.44     6.118   0.00113 **
#Residuals   56  296.99    5.30 
attributes(av1)
#"coefficients"  "residuals" "effects" "rank" "fitted.values" "assign" "qr" 
#"df.residual"  "contrasts" "xlevels" "call" "terms"  "model"        

av1$coefficients
#(Intercept)       dietB       dietC       dietD  
#9.1800000  -0.2733333   2.9333333   1.3600000 [A er sathe B, C, D er difference er value ]
#[dietA ke reference hisebe dore niye baki value ke dietA er sathe tulona korese ]














