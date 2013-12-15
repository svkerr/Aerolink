# Set local working directory
setwd("/Users/stuart/R_Files/Aerolink")

# Load required libraries
library(zoo)
library(ggplot2)
library(psych)
library(lsr)
library(car)
library(plyr)
library(xts)

# Badge References:
Dashofy: 28006
B. Wong: 21462
Razouk: 19402

### Read in data files (prefer to have all reads from external files here)

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

# Example of reading data over the web interface:
x <- read.table("http://tcsnc3:8080/yearly/2010/25876",header=FALSE,sep=",")

### Transformations and subsets of the above datasets required to do all 
### analysis in this file

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

# Capture the badges to be used later as column names
badges <- fetch2$V1
badges2010 <- fetch2010$V1
badges2011 <- fetch2011$V1
badges2012 <- fetch2012$V1
badges2013 <- fetch2013$V1
badges2014 <- fetch2014$V1

### TIME SERIES ###########

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


### FOURIER ANALYSIS #################
# Let's try FFT of a) fetches and b) the time differences between fetches
# Reference: http://homepage.univie.ac.at/erhard.reschenhofer/pdf/zr/Spec.pdf

#stu_mat_ts_agg <- aggregate(stu_mat_ts)/24

plot(abs(fft(stu_hrs10_numeric)))
plot( Im(fft(stu_hrs10_numeric))/Re(fft(stu_hrs10_numeric)) )
stu_hrs10_numeric

#### Estimate spectral density via a periodogram ##########

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
############################################################################################
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

#############################################################################################
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
# Order badge numbers so that we can plot newest/oldest to oldest/newest employees
arranged <- arrange(fetch,desc(Badge))
arranged2010 <- arrange(fetch2010wh,desc(Badge))
arranged2011 <- arrange(fetch2011wh,desc(Badge))
arranged2012 <- arrange(fetch2012wh,desc(Badge))
arranged2013 <- arrange(fetch2013wh,desc(Badge))
arranged2014 <- arrange(fetch2014wh,desc(Badge))

##########  Compare badges for differences between years #############
head(arranged2011, n = 30L)
arranged2010[1:10,1]
str(arranged2011[1:10, 1])

rowsums2010 <- rowSums(arranged2010[1:nrow(arranged2010),2:25])
rowsums2011 <- rowSums(arranged2011[1:nrow(arranged2011),2:25])
rowsums2012 <- rowSums(arranged2012[1:nrow(arranged2012),2:25])
rowsums2013 <- rowSums(arranged2013[1:nrow(arranged2013),2:25])
rowsums2014 <- rowSums(arranged2014[1:nrow(arranged2014),2:25])

head(rowsums2013, n = 100L)
par(mfrow=c(1,1))
summary(rowsums2010[1:1000])
summary(rowsums2011[1:1000])
summary(rowsums2012[1:1000])
summary(rowsums2013[1:1000])

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

## create a dataframe for an individual user for years 2010-2013
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

###### ANALYSIS:  t-test and Wilcoxon  ################################################
# Do a Dependent t-test (same person sampled during 2 sequential years)
# Because t-test assumes "homogeniety of variance" - we must test for this

var.test(EmpA_t$'2012',EmpA_t$'2013') 

# Since p-value is large, i can assume that homogeniety of variance is valid thus proceed with t-test.
# Note: If homogeniety of variance does not exist, use a Welch two sample test
# Note: When performing > 2 dependent t-test, use "Repeated ANOVA"

t.test(EmpA_t$'2010',EmpA_t$'2011', paired=T)
t.test(EmpA_t$'2011',EmpA_t$'2012', paired=T)
t.test(EmpA_t$'2012',EmpA_t$'2013', paired=T)

# But wait, the distributions are not normal, hence i should defer to Wilcoxon

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

############## General Fetch useage plots #########################
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
