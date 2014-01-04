# Set local working directory
setwd("/users/stuart/R_Files/Aerolink")

# The analysis contained in this aerolink.R (R-file) supports an effort to identify and model 'anomalous'
# Aerolink activity. Though Aerolink has several file commands available, we selected the 'fetch'
# command (a file extraction task) upon which to base a user's behavior, and ultimately, describe
# or identify anomalous changes in behavior.' This effort is an exploration as to which
# of the various data analytical techniques are most appropriate for this problem. 

# Analytical techniques used:
#  General plotting: to get a feel for data, distribution types, etc., Scattered throughout major sections below
#  Null Hypothesis Significance Testing (NHST): to determine with a given confidence a change is significant
#  Frequency Analysis (Fourier Transforms): to identify patterns, fluctuations
#  Clustering: to identify groups of Aerolink users sharing common useage patterns

# Load required libraries
library(zoo)
library(ggplot2)
library(psych)
library(lsr)
library(car)
library(plyr)
library(xts)
library(cluster)
library(fpc)
library(mclust)

# Read in data files:
fetch <- read.table("aerolink_fetch_multi.dat", header = T, sep=',')
fetch2 <- read.table("aerolink_fetch_multi_noheads.dat", sep = ',')
fetch2010 <- read.table("aerolink_fetch_access_20092010.dat", sep = ',')
fetch2011 <- read.table("aerolink_fetch_access_20102011.dat", sep = ',')
fetch2012 <- read.table("aerolink_fetch_access_20112012.dat", sep = ',')
fetch2013 <- read.table("aerolink_fetch_access_20122013.dat", sep = ',')
fetch2014 <- read.table("aerolink_fetch_access_20132014.dat", sep = ',')

fetch2010wh <- read.table("aerolink_fetch_access_wh_20092010.dat", header = T, sep=',')
fetch2011wh <- read.table("aerolink_fetch_access_wh_20102011.dat", header = T, sep=',')
fetch2012wh <- read.table("aerolink_fetch_access_wh_20112012.dat", header = T, sep=',')
fetch2013wh <- read.table("aerolink_fetch_access_wh_20122013.dat", header = T, sep=',')
fetch2014wh <- read.table("aerolink_fetch_access_wh_20132014.dat", header = T, sep=',')

kerr_fetch_times <- read.table("stu_fetch_sorted.dat")
dash_fetch_times <- read.table("dashofy_fetch_sort.dat")

daily_2010 <- read.table("fetch_daily_2009_2010", header = F, sep=',')

# Data over the web interface:
x <- read.table("http://tcsnc3:8080/yearly/2010/25876",header=FALSE,sep=",")

# Transformations and subsets (Note: most transformations happen in each section)

stu10 <- subset(fetch2010, V1 == 25876)
stu11 <- subset(fetch2011, V1 == 25876)
stu12 <- subset(fetch2012, V1 == 25876)
stu13 <- subset(fetch2013, V1 == 25876)

stu_hrs10 <- stu10[2:25]
stu_hrs11 <- stu11[2:25]
stu_hrs12 <- stu12[2:25]
stu_hrs13 <- stu13[2:25]

stu_hrs10_numeric <- as.numeric(stu_hrs10)
stu_hrs11_numeric <- as.numeric(stu_hrs11)
stu_hrs12_numeric <- as.numeric(stu_hrs12)
stu_hrs13_numeric <- as.numeric(stu_hrs13)

dash10 <- subset(fetch2010, V1 == 28006)
dash11 <- subset(fetch2011, V1 == 28006)
dash12 <- subset(fetch2012, V1 == 28006)
dash13 <- subset(fetch2013, V1 == 28006)

dash_hrs10 <- dash10[2:25]
dash_hrs11 <- dash11[2:25]
dash_hrs12 <- dash12[2:25]
dash_hrs13 <- dash13[2:25]

dash_hrs10_numeric <- as.numeric(dash_hrs10)
dash_hrs11_numeric <- as.numeric(dash_hrs11)
dash_hrs12_numeric <- as.numeric(dash_hrs12)
dash_hrs13_numeric <- as.numeric(dash_hrs13)

# Capture badges to be used later as column names
badges <- fetch2$V1
badges2010 <- fetch2010$V1
badges2011 <- fetch2011$V1
badges2012 <- fetch2012$V1
badges2013 <- fetch2013$V1
badges2014 <- fetch2014$V1

######### TIME SERIES ###############################################
# Note: Did not make a lot of headway here. Know restrictions on ts()

stu_hrs10_numeric_ts <- ts(stu_hrs10_numeric)
stu_hrs_df <- rbind(stu_hrs10,stu_hrs11, stu_hrs12, stu_hrs13)
col_names <- c('h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12','h13','h14','h15','h16','h17','h18','h19','h20','h21','h22','h23','h24')
colnames(stu_hrs_df) <- col_names
rownames(stu_hrs_df) <- c('2010','2011','2012','2013')
stu_hrs_df
str(stu_hrs_df)
#stu_hrs_tot <- rbind(stu_hrs10_numeric,stu_hrs11_numeric, stu_hrs12_numeric, stu_hrs13_numeric)
#stu_mat <- as.matrix(stu_hrs_tot)

stu_ts <- ts(stu_hrs_df,start = c(2010,1), end = c(2013,24),frequency = 24)
stu_ts
ts.plot(stu_ts)
start(stu_ts)
end(stu_ts)
frequency(stu_ts)
plot(stu_ts)
decompose(stu_ts)


######### FOURIER ANALYSIS ###################################
# Let's try FFT of a) fetches and b) the time differences between fetches
# Reference: http://homepage.univie.ac.at/erhard.reschenhofer/pdf/zr/Spec.pdf

#stu_mat_ts_agg <- aggregate(stu_mat_ts)/24

plot(abs(fft(stu_hrs10_numeric)))
plot( Im(fft(stu_hrs10_numeric))/Re(fft(stu_hrs10_numeric)),type='l' )

# Estimate spectral density via a periodogram ######

# get data that contains fetch times

kft <- kerr_fetch_times$V1
eft <- dash_fetch_times$V1
result <- head(kft) - head(eft)

# calculate time differences between fetches
kerr_fetch_times_diff <- diff(kft)
dash_fetch_times_diff <- diff(eft)

par(mfrow=c(2,1))
spectrum(kerr_fetch_times_diff, method = c('pgram'))
spectrum(dash_fetch_times_diff, method = c('pgram'))

par(mfrow=c(2,2))
spectrum(stu_hrs10_numeric,method = c('pgram'))
spectrum(stu_hrs11_numeric,method = c('pgram'))
spectrum(stu_hrs12_numeric,method = c('pgram'))
spectrum(stu_hrs13_numeric,method = c('pgram'))

par(mfrow=c(2,2))
spectrum(dash_hrs10_numeric,method = c('pgram'))
spectrum(dash_hrs11_numeric,method = c('pgram'))
spectrum(dash_hrs12_numeric,method = c('pgram'))
spectrum(dash_hrs13_numeric,method = c('pgram'))

# Let's do spectrum on all 4 sequential years of hourly fetch data 
EmpA_all <- c(stu_hrs10_numeric, stu_hrs11_numeric, stu_hrs12_numeric,stu_hrs13_numeric)
EmpB_all <- c(dash_hrs10_numeric, dash_hrs11_numeric, dash_hrs12_numeric,dash_hrs13_numeric)

par(mfrow=c(2,1))
plot(EmpA_all,type="l", main = "Aggregate Hourly Fetches for 4 Years")
spectrum(EmpA_all, spans = c(3,3))

plot(EmpB_all,type="l", main = "Aggregate Hourly Fetches for 4 Years")
spectrum(EmpB_all, spans = c(3,3))

spectrum(EmpA_all)
spectrum(EmpA_all, spans = c(3,3))

# Let's do spectrum on 24 hour time period, where 4 years of hourly stacked fetches are stacked
stu_rbind <- rbind(stu_hrs10,stu_hrs11,stu_hrs12,stu_hrs13)
EmpA_rbind_summed <- as.numeric(colSums(stu_rbind))

par(mfrow=c(2,1))
plot(EmpA_rbind_summed, type='l', main = "Aggregate Fetches per hour for 4 Years", ylab = c('EmpA Fetches'), xlab = c('24 hours'),ylim=c(0,700))
spectrum(EmpA_rbind_summed, method=c("pgram"))

dash_rbind <- rbind(dash_hrs10, dash_hrs11, dash_hrs12, dash_hrs13)
EmpB_rbind_summed <- as.numeric(colSums(dash_rbind))

par(mfrow=c(2,1))
plot(EmpB_rbind_summed, type='l', main = "Aggregate Fetches per hour for 4 Years", ylab = c('EmpB Fetches'), xlab = c('24 hours'), ylim=c(0,700))
spectrum(EmpB_rbind_summed, method = c('pgram'))

######## AREA TO FIGURE OUT HOW TO PLOT FETCHES BY BADGE NUMBER ###############

# Order badge numbers so that we can plot newest/oldest or oldest/newest employees
arranged <- arrange(fetch,desc(Badge))
arranged2010 <- arrange(fetch2010wh,desc(Badge))
arranged2011 <- arrange(fetch2011wh,desc(Badge))
arranged2012 <- arrange(fetch2012wh,desc(Badge))
arranged2013 <- arrange(fetch2013wh,desc(Badge))
arranged2014 <- arrange(fetch2014wh,desc(Badge))

##########  Compare badges for differences between years #############
# Get some feel and check the data
head(arranged2011, n = 30L)
arranged2010[1:10,1]
str(arranged2011[1:10, 1])

# Total number of fetches per employee for a given year
rowsums2010 <- rowSums(arranged2010[1:nrow(arranged2010),2:25])
rowsums2011 <- rowSums(arranged2011[1:nrow(arranged2011),2:25])
rowsums2012 <- rowSums(arranged2012[1:nrow(arranged2012),2:25])
rowsums2013 <- rowSums(arranged2013[1:nrow(arranged2013),2:25])
rowsums2014 <- rowSums(arranged2014[1:nrow(arranged2014),2:25])

par(mfrow=c(1,1))
# Let's look at some stats for contiguous badges in groups of 1K (descending)
# Stats for newest 1K employees
summary(rowsums2010[1:1000])
summary(rowsums2011[1:1000])
summary(rowsums2012[1:1000])
summary(rowsums2013[1:1000])

# Stats for next 1K employees
summary(rowsums2010[1001:2000])
summary(rowsums2011[1001:2000])
summary(rowsums2012[1001:2000])
summary(rowsums2013[1001:2000])

# Stats for next 1K employees
summary(rowsums2010[2001:3000])
summary(rowsums2011[2001:3000])
summary(rowsums2012[2001:3000])
summary(rowsums2013[2001:3000])

# Stats for next 1K employees
summary(rowsums2010[3001:4000])
summary(rowsums2011[3001:4000])
summary(rowsums2012[3001:4000])
summary(rowsums2013[3001:4000])

boxplot(log10(rowsums2010[1001:2000]))  ### do a box plot with a vertical rug
log10(rowsums2010)

par(mfrow=c(2,2))
plot(rowsums2010[1:1000], xlim = c(0,1000), ylim = c(0,8000))
plot(rowsums2011[1:1000], xlim = c(0,1000), ylim = c(0,8000))
plot(rowsums2012[1:1000], xlim = c(0,1000), ylim = c(0,8000))
plot(rowsums2013[1:1000], xlim = c(0,1000), ylim = c(0,8000))

plot(rowsums2010[1001:2000], xlim = c(0,1000), ylim = c(0,8000))
plot(rowsums2011[1001:2000], xlim = c(0,1000), ylim = c(0,8000))
plot(rowsums2012[1001:2000], xlim = c(0,1000), ylim = c(0,8000))
plot(rowsums2013[1001:2000], xlim = c(0,1000), ylim = c(0,8000))

plot(rowsums2010[2001:3000], xlim = c(0,1000), ylim = c(0,8000))
plot(rowsums2011[2001:3000], xlim = c(0,1000), ylim = c(0,8000))
plot(rowsums2012[2001:3000], xlim = c(0,1000), ylim = c(0,8000))
plot(rowsums2013[2001:3000], xlim = c(0,1000), ylim = c(0,8000))

rowsums2011[2001:3000]
########### Get differences between fetch times - any pattern?################
EmpA_fetch_times_diff <- diff(kerr_fetch_times$V1)
EmpB_fetch_times_diff <- diff(dash_fetch_times$V1)

# Look at the files
head(kerr_fetch_times_diff)
str(kerr_fetch_times_diff)
describe(EmpA_fetch_times_diff)
describe(EmpB_fetch_times_diff)

# Plot time differences using histograms of density and actual numbers
par(mfrow=c(2,1))
hist(log10(EmpA_fetch_times_diff), ylim = c(0,800), xlim = c(0,7))
hist(log10(EmpB_fetch_times_diff), ,ylim = c(0,800), xlim = c(0,7))

par(mfrow=c(2,1))
hist(log10(EmpA_fetch_times_diff),freq=FALSE, ylim = c(0,.4))
hist(log10(EmpB_fetch_times_diff),freq=FALSE, ylim = c(0,0.4))

# How many individual badges are we tracking?
nrow(fetch2010)

### Now get numerical matrix, remove badges and transpose
hours <- subset(fetch2,select = c(V2:V25) )

# Get pure numerical df
hours2010 <- subset(fetch2010, select = c(V2:V25) )
hours2011 <- subset(fetch2011, select = c(V2:V25) )
hours2012 <- subset(fetch2012, select = c(V2:V25) )
hours2013 <- subset(fetch2013, select = c(V2:V25) )
hours2014 <- subset(fetch2014, select = c(V2:V25) )

# transpose df so that 24 hours are observations (rows)
hours_t <- data.frame(t(hours))
hours2010_t <- data.frame(t(hours2010))
hours2011_t <- data.frame(t(hours2011))
hours2012_t <- data.frame(t(hours2012))
hours2013_t <- data.frame(t(hours2013))
hours2014_t <- data.frame(t(hours2014))
hours2014_t[,1]  # test to see if transpose occured

# rowbind badges so that it is the top "row" in the dataframes - to be made into col names
hours_wb <- rbind(badges,hours_t)
head(hours_wb)
hours2010_rb <- rbind(badges2010,hours2010_t)
hours2011_rb <- rbind(badges2011,hours2011_t)
hours2012_rb <- rbind(badges2012,hours2012_t)
hours2013_rb <- rbind(badges2013,hours2013_t)
hours2014_rb <- rbind(badges2014,hours2014_t)

# assign badges as the col names
colnames(hours_wb) <- hours_wb[1,]
head(hours_wb) 
colnames(hours2010_rb) <- hours2010_rb[1,]
colnames(hours2011_rb) <- hours2011_rb[1,]
colnames(hours2012_rb) <- hours2012_rb[1,]
colnames(hours2013_rb) <- hours2013_rb[1,]
colnames(hours2014_rb) <- hours2014_rb[1,]

# so remove first row of integer badge numbers, and we're done 
hours_mat <- hours_wb[-1,]
str(hours_mat)
str(hours2010_mat)
hours2010_mat[1,]

hours2010_mat <- hours2010_rb[-1,]
hours2011_mat <- hours2011_rb[-1,]
hours2012_mat <- hours2012_rb[-1,]
hours2013_mat <- hours2013_rb[-1,]
hours2014_mat <- hours2014_rb[-1,]

####  Test and verify######################################
# For individual
par(mfrow=c(3,2))
plot(hours2010_mat$'28006')
plot(hours2011_mat$'28006')
plot(hours2012_mat$'28006')
plot(hours2013_mat$'28006')
plot(hours2014_mat$'28006')
hours2010_mat$'25876'

# Create a dataframe for an individual user for years 2010-2013
# rbind the years, transpose to get years as columns, then give columns names
# Kerr and Dashofy
# Because plots show "names" -- changing the badge tags of Kerr and Dashofy to EmpA and EmpB respectively

user25876 <- rbind(hours2010_mat$'25876',hours2011_mat$'25876',hours2012_mat$'25876',hours2013_mat$'25876')
EmpA_t <- data.frame(t(user25876))
colnames(EmpA_t) <- c(2010,2011,2012,2013)

user28006 <- rbind(hours2010_mat$'28006',hours2011_mat$'28006',hours2012_mat$'28006',hours2013_mat$'28006')
EmpB_t <- data.frame(t(user28006))
colnames(EmpB_t) <- c(2010,2011,2012,2013)

colSums(user25876_t)
plot(colSums(EmpA_t))
plot(colSums(EmpA_t))

###### ANALYSIS:  NHST (t-test and Wilcoxon) ################################################
# Do a Dependent t-test (same person sampled during 2 sequential years)
# Because t-test assumes "homogeniety of variance" - we must test for this

var.test(EmpA_t$'2012',EmpA_t$'2013') 

# p-value is large. Can i can assume that homogeniety of variance exists in order to proceed with t-test.
# Note: If homogeniety of variance does not exist, use a Welch two sample test or consider Wilcoxan
# Note: When performing > 2 dependent t-test, use "Repeated ANOVA"

t.test(EmpA_t$'2010',EmpA_t$'2011', paired=T)
t.test(EmpA_t$'2011',EmpA_t$'2012', paired=T)
t.test(EmpA_t$'2012',EmpA_t$'2013', paired=T)

# But wait, the distributions are not normal, so switch to Wilcoxon

wilcox.test(EmpA_t$'2010',EmpA_t$'2011', paired = T)
wilcox.test(EmpA_t$'2011',EmpA_t$'2012', paired = T)
wilcox.test(EmpA_t$'2012',EmpA_t$'2013', paired = T)

# Plots to support the above wilcoxan tests
par(mfrow = c(2,2))
plot(as.numeric(stu_hrs10),type='l', ylab =c('EmpA Fetches'),xlab = c('Clock Hour'), main = c('Fetch vs Hour 2010'),ylim = c(0,70))
plot(as.numeric(stu_hrs11),type='l', ylab =c('EmpA Fetches'),xlab = c('Clock Hour'), main = c('Fetch vs Hour 2011'),ylim = c(0,70))
plot(as.numeric(stu_hrs12),type='l', ylab =c('EmpA Fetches'),xlab = c('Clock Hour'), main = c('Fetch vs Hour 2012'),ylim = c(0,70))
plot(as.numeric(stu_hrs13),type='l', ylab =c('EmpA Fetches'),xlab = c('Clock Hour'), main = c('Fetch vs Hour 2013'),ylim = c(0,70))

# NOTE: There exists a significant difference between 2012 and 2013, but not between other years

# Summary stats for stu's hours

sum(hours_mat$'25876')
summary(hours_mat$'25876')
hours_mat$'25876'

# For group
plot(colSums(hours_mat),ylim = c(0,2500))
par(mfrow=c(1,1))
plot(colSums(hours2010_mat), ylim = c(0,2500))
plot(colSums(hours2011_mat), ylim = c(0,2500))
plot(colSums(hours2012_mat), ylim = c(0,2500))
plot(colSums(hours2013_mat), ylim = c(0,2500))
plot(colSums(hours2014_mat), ylim = c(0,2500))

# Plot total fetches by badge number
badges_2010 <- hours2010_mat[1,]
col_sums2010 <- colSums(hours2010_mat)
plot(badges_2010,col_sums2010, type = 'l', ylim = c(0,100), xlim = c(10000,10500))

col_sums <- colSums(hours_mat)
col_sums_0 <- subset(col_sums, col_sums == 0)
col_sums_20 <- subset(col_sums, col_sums <= 20)
col_sums_40 <- subset(col_sums, col_sums <= 40)
col_sums_60 <- subset(col_sums, col_sums <= 60)
col_sums_80 <- subset(col_sums, col_sums <= 80)
col_sums_100 <- subset(col_sums, col_sums <= 100)

length(col_sums_0)
length(col_sums_20)
length(col_sums_40)
length(col_sums_60)
length(col_sums_80)
length(col_sums_100)

###Histogram and other descriptive stats of Badge numbers

depth <- ggplot(fetch,aes(Badge)) + xlim(3000,28000)
depth + geom_histogram(binwidth = 100)

badge_20000 <- subset(fetch, fetch$Badge <= 20000)
nrow(badge_20000)

############## General Fetch useage Plots #########################
colsumsbig <- colSums(fetch)
colsumsbighours <- colsumsbig[2:25]
plot(colsumsbighours, type = 'l',ylab = 'Number of Fetchs (2010-2014)', xlab = '24-hour Clock (PST)')

# Now do above for individual years
colsums2010 <- colSums(fetch2010)
colsums2010_hours <- colsums2010[2:25]

colsums2011 <- colSums(fetch2011)
colsums2011_hours <- colsums2011[2:25]

colsums2012 <- colSums(fetch2012)
colsums2012_hours <- colsums2012[2:25]

colsums2013 <- colSums(fetch2013)
colsums2013_hours <- colsums2013[2:25]

colsums2014 <- colSums(fetch2014)
colsums2014_hours <- colsums2014[2:25]

par(mfrow=c(2,3))
plot(colsumsbighours, type = 'l',ylab = 'Number of Fetchs (2010-2014)', xlab = '24-hour Clock (PST)')
plot(colsums2010_hours, type = 'l',ylim = c(0,120000), ylab = 'Number of Fetchs (2010)', xlab = '24-hour Clock (PST)')
plot(colsums2011_hours, type = 'l', ylim = c(0,120000),ylab = 'Number of Fetchs (2011)', xlab = '24-hour Clock (PST)')
plot(colsums2012_hours, type = 'l', ylim = c(0,120000),ylab = 'Number of Fetchs (2012)', xlab = '24-hour Clock (PST)')
plot(colsums2013_hours, type = 'l',ylim = c(0,120000), ylab = 'Number of Fetchs (2013)', xlab = '24-hour Clock (PST)')
plot(colsums2014_hours, type = 'l',ylim = c(0,120000),ylab = 'Number of Fetchs (2014)', xlab = '24-hour Clock (PST)')


########### Cluster Analysis ##########################################
# Generally speaking, cluster analysis methods are of either of two types: 
# 1. 'Partitioning' methods: Algorithms that divide the data jet into k clusters, where the integer k needs 
# to be specified by the user. Typically, the user runs the algorithm for a range of k-values. For each k, 
# the algorithm carries out the clustering and also yields a “quality index”, which allows the user to select a value of k afterwards.
# 2. 'Hierarchical' methods: Algorithms yielding an entire hierarchy of clusterings of the data set. Agglomerative methods 
# start with the situation where each object in the data set forms its own little cluster, and then successively merge clusters until 
# only one large cluster remains which is the whole data set. Divisive methods start by considering the whole data set as one cluster,
# and then split up clusters until each object is separate. 

# Cluster analysis involves the partitioning of data into meaningful subgroups, when the
# number of subgroups and other information about their composition may be unknown. It is designed
# to identify groups of people that share certain common characteristics (Aerolink useage)
# Cluster analysis is frequently employed as a classification tool or to represent the structure in a data set. 
# It is a statistical method, where unlike other methods, makes no prior assumptions about a population. 
# Here, it is used to segment the Aerolink user community. 

# Remove badge number column for analysis, will add it back later
test2010 <- fetch2010wh[2:25]     
test2011 <- fetch2011wh[2:25]
test2012 <- fetch2012wh[2:25]
test2013 <- fetch2013wh[2:25]

# Let's scale the data sets (each column has mean=0, sd=1)
test2010_s <- scale(test2010)
test2011_s <- scale(test2011)
test2012_s <- scale(test2012)
test2013_s <- scale(test2013)

# K-Means is an unsupervised learning algorithm.
# The aim of K-means clustering is to divide M points, in N dimensions (24) into K clusters (determined apriori) so that the within-cluster
# sum of squares is minimized. We seek local optima solutions so that no movement of a point from one cluster to 
# another will reduce the within-cluster sum of squares. The algorithm minimizes the sum of the squared distances to the cluster centers.
# Unfortunately there is no general theoretical solution to find the optimal number of clusters for any given data set. 
# A simple approach is to compare the results of multiple runs with different k classes and choose the best one according to a 
# given criterion; for instance the Schwarz Criterion. Note that the Schwarz Criterion is related to the Bayes Information Criterion (BIC) 

# This is a 24-dimensional space. Selecting 5 clusters (as a guess)
fit2010 <- kmeans(test2010,5)
fit2011 <- kmeans(test2011,5) 
fit2012 <- kmeans(test2012,5) 
fit2013 <- kmeans(test2013,5) 

# Doing same as above, but for scaled data
fit2010_s <- kmeans(test2010_s,5)
fit2011_s <- kmeans(test2011_s,5)
fit2012_s <- kmeans(test2012_s,5)
fit2013_s <- kmeans(test2013_s,5)

# Cluster centers for kmeans k = 5 (by hour)
fit2010$centers 
fit2011$centers
fit2012$centers
fit2013$centers

# Same as above for scaled data
fit2010_s$centers 
fit2011_s$centers 
fit2012_s$centers 
fit2013_s$centers 

# Let's try to analytically figure out an appropriate number of clusters (k)
# Naturally, the within groups sum of squares decreases as we increase the number of clusters. However, 
# there is a trend of diminishing marginal returns as we increase the number of clusters. I select the number
# of clusters based on the point at which the marginal return of adding one more cluster is less than was the 
# marginal return for adding the clusters prior to that, or simply eyeballig the graph looking for a 'knee'

ssPlot <- function(data,maxCluster = 9) {
  # Initialize within sum of squares
  SSa <- (nrow(data) - 1) * sum(apply(data,2,var))
  SSa <- vector()
  for (i in 2:maxCluster) {
    SSa[i] <- sum(kmeans(data, centers = i)$withinss)
  }
  plot(1:maxCluster,SSa, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares")
  title("AeroLink kMeans Cluster Data")
}
# call function and k = 4 or 5 seems to be a good number (year 2013 is odd)
ssPlot(test2010_s)
ssPlot(test2011_s)
ssPlot(test2012_s)
ssPlot(test2013_s)

### Let's try another analytical approach to determine the number of clusters
pamk(test2010,2:7, usepam=FALSE, criterion="multiasw")

# Reshape original dataframe to include badge, k-means cluster assignment, and hourly fetch data
test2010_clustered <- cbind(fetch2010wh[1],fit2010$cluster,test2010)  
test2011_clustered <- cbind(fetch2011wh[1],fit2011$cluster,test2011)
test2012_clustered <- cbind(fetch2012wh[1],fit2012$cluster,test2012)
test2013_clustered <- cbind(fetch2013wh[1],fit2013$cluster,test2013)

# Same as above for scaled data
test2010_clustered_s <- cbind(fetch2010wh[1],fit2010_s$cluster,test2010) 
test2011_clustered_s <- cbind(fetch2011wh[1],fit2011_s$cluster,test2011)
test2012_clustered_s <- cbind(fetch2012wh[1],fit2012_s$cluster,test2012)
test2013_clustered_s <- cbind(fetch2013wh[1],fit2013_s$cluster,test2013)

# Here we use another (preferred) clustering technique, partitioning around medoids (pam). The function pam is based on the search 
# for k representative objects, called medoids; among the actual objects of the data set (Kaufman and Rousseeuw, 1987). 
# These medoids are computed such that the total 'dissimilarity' of all objects to their nearest medoid is minimal
# Dissimilarity matrix is the distinguishing feature of pam versus k-means.
# In the k-means method the center of each cluster is defined as the mean (average) 
# of all objects of the cluster. Thus, kmeans requires the input of an n x p data matrix, 
# unlike pam. The goal of kmeans is to minimize a sum of squared euclidean distances, 
# implicitly assuming that each cluster has a spherical normal distribution. The function 
# pam is more robust because it minimizes a sum of unsquared 'dissimilarities,' reducing distortion caused
# by outliers. Moreover pam does not need initial guesses for the cluster centers, contrary to kmeans.

fit2010_pam <- pam(test2010,5, diss=FALSE)
fit2011_pam <- pam(test2011,5, diss=FALSE)
fit2012_pam <- pam(test2012,5, diss=FALSE)
fit2013_pam <- pam(test2013,5, diss=FALSE)

# pam for scaled data
fit2010_pam_s <- pam(test2010_s,5, diss=FALSE)
fit2011_pam_s <- pam(test2011_s,5, diss=FALSE)
fit2012_pam_s <- pam(test2012_s,5, diss=FALSE)
fit2013_pam_s <- pam(test2013_s,5, diss=FALSE)

# try same above but without setting diss=FALSE (dissimilarity matrix is used) - no difference
fit2010_pam_s_d <- pam(test2010_s,5)

# Reshape original dataframe to include badge, pam cluster assignment, and hourly fetch data
test2010_pam_clustered <- cbind(fetch2010wh[1],fit2010_pam$clustering,test2010)  
test2011_pam_clustered <- cbind(fetch2011wh[1],fit2011_pam$clustering,test2011)
test2012_pam_clustered <- cbind(fetch2012wh[1],fit2012_pam$clustering,test2012)
test2013_pam_clustered <- cbind(fetch2013wh[1],fit2013_pam$clustering,test2013)

# Reshape using scaled data
test2010_pam_clustered_s <- cbind(fetch2010wh[1],fit2010_pam_s$clustering,test2010)
test2011_pam_clustered_s <- cbind(fetch2011wh[1],fit2011_pam_s$clustering,test2011)
test2012_pam_clustered_s <- cbind(fetch2012wh[1],fit2012_pam_s$clustering,test2012)
test2013_pam_clustered_s <- cbind(fetch2013wh[1],fit2013_pam_s$clustering,test2013)

### Graphing Clusters
# Clusplot Creates a bivariate (2-dim) plot visualizing a partition (clustering) of the data. 
# All observation are represented by points in the plot.
# Clusplot uses principal components analysis PCA or multidimensional scaling in order to reduce 24-space to 2-space 
# and uses the first 2 principal components. Principal components are a sequence of projections of the data, 
# mutually uncorrelated (orthogonal) and ordered by variance. For these projections, the two most significant eigenvectors are used.
# An ellipse is drawn around each cluster. 

# Graphing k-means
library(cluster)
par(mfrow=c(1,1))
clusplot(test2010,fit2010$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)
clusplot(test2011,fit2011$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)
clusplot(test2012,fit2012$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)
clusplot(test2013,fit2013$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)

# Same as above, but using scaled data
clusplot(test2010_s,fit2010_s$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)
clusplot(test2011_s,fit2011_s$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)
clusplot(test2012_s,fit2012_s$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)
clusplot(test2013_s,fit2013_s$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)

# Graphing pam
clusplot(test2010,fit2010_pam$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)
clusplot(test2011,fit2011_pam$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)
clusplot(test2012,fit2012_pam$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)
clusplot(test2013,fit2013_pam$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)

# Graphing pam: Same as above for scaled data
clusplot(test2010_s,fit2010_pam_s$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)
clusplot(test2011_s,fit2011_pam_s$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)
clusplot(test2012_s,fit2012_pam_s$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)
clusplot(test2013_s,fit2013_pam_s$cluster,color=TRUE, shade=TRUE,labels=2,lines=0)

# Get summary information (medoids, clustering, diss, etc., )
# From this data, identifying the heavy medoids/cluster is used in creating subsets for heavy users and can
# be used in cluster descriptive statistics (e.g., % of 'light' users)
head(fit2010_pam)
head(fit2011_pam)
head(fit2012_pam)
head(fit2013_pam)

# Same as for scaled data:
head(fit2010_pam_s)
head(fit2011_pam_s)
head(fit2012_pam_s)
head(fit2013_pam_s)

# Focusing on just the pam clusters of scaled data
fit2010_pam_s$clusinfo
fit2011_pam_s$clusinfo
fit2012_pam_s$clusinfo
fit2013_pam_s$clusinfo

# Subset heavy users per year using pam clustered data
heavy1_2010_pam <- subset(test2010_pam_clustered,test2010_pam_clustered[2] == 5 )
heavy2_2010_pam <- subset(test2010_pam_clustered,test2010_pam_clustered[2] == 4 )

heavy1_2011_pam <- subset(test2011_pam_clustered,test2011_pam_clustered[2] == 5 )
heavy2_2011_pam <- subset(test2011_pam_clustered,test2011_pam_clustered[2] == 4 )

heavy1_2012_pam <- subset(test2012_pam_clustered,test2012_pam_clustered[2] == 5 )
heavy2_2012_pam <- subset(test2012_pam_clustered,test2012_pam_clustered[2] == 4 )

heavy1_2013_pam <- subset(test2013_pam_clustered,test2013_pam_clustered[2] == 5 )
heavy2_2013_pam <- subset(test2013_pam_clustered,test2013_pam_clustered[2] == 3 )

# Same as above but for scaled kmeans clusters
heavy1_2010_s <- subset(test2010_clustered_s,test2010_clustered_s[2] == 5 )
heavy2_2010_s <- subset(test2010_clustered_s,test2010_clustered_s[2] == 4 )
heavy1_2010_s 
heavy2_2010_s

# Compare with non-scaled kmeans clusters
heavy1_2010 <- subset(test2010_clustered,test2010_clustered[2] == 1 )
heavy2_2010 <- subset(test2010_clustered,test2010_clustered[2] == 5 )
heavy1_2010
heavy2_2010


# Model Based Clustering (uses library mclust)
fit2010_mclust <- Mclust(test2010)
fit2011_mclust <- Mclust(test2011)
fit2012_mclust <- Mclust(test2012)
fit2013_mclust <- Mclust(test2013)

plot(fit2) # plot results 

summary(fit) # display the best model
