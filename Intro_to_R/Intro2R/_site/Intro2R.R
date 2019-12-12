## ------------------------------------------------------------------------
# assign number 150 to variable a.
a <- 150
a
# assign a number in scientific format to variable b.
b <- 3e-2
b

## ------------------------------------------------------------------------
# assign a string "Professor" to variable title
title <- "Professor"
title
# assign a string "Hello World" to variable hello
hello <- "Hello World"
hello

## ------------------------------------------------------------------------
# assign logical value "TRUE" to variable is_female
is_female <- TRUE
is_female
# assign logical value "FALSE" to variable is_male
is_male <- FALSE
is_male
# assign logical value to a variable by logical operation
age <- 20
is_adult <- age > 18
is_adult

## ------------------------------------------------------------------------
class(is_female)
# To check whether the variable is a specific type
is.numeric(hello)
is.numeric(a)
is.character(hello)

## ------------------------------------------------------------------------
as.numeric(is_female)
as.numeric(is_male)

## ------------------------------------------------------------------------
b
as.character(b)

## ------------------------------------------------------------------------
friend_ages <- c(21, 27, 26, 32)
friend_ages

friend_names <- c("Mina", "Ella", "Anna", "Cora")
friend_names

## ------------------------------------------------------------------------
# assign names to a vector by specifying them
names(friend_ages) <- c("Mina", "Ella", "Anna", "Carla")
friend_ages

# assign names to a vector using another vector
names(friend_ages) <- friend_names
friend_ages

## ------------------------------------------------------------------------
friend_ages <- c(Mina=21, Ella=27, Anna=26, Cora=32)
friend_ages

## ------------------------------------------------------------------------
length(friend_ages)

## ------------------------------------------------------------------------
friend_ages[2]
friend_ages["Ella"]

friend_ages[c(1,3)]
friend_ages[c("Mina", "Anna")]

# selecting elements of a vector by excluding some of them.
friend_ages[-3]

## ------------------------------------------------------------------------
my_friends <- c("Mina", "Ella", "Anna", "Cora")
my_friends
has_child <- c("TRUE", "TRUE", "FALSE", "TRUE")
has_child

my_friends[has_child == "TRUE"]

## ------------------------------------------------------------------------
col1 <- c(1,3,8,9)
col2 <- c(2,18,27,10)
col3 <- c(8,37,267,19)

my_matrix <- cbind(col1, col2, col3)
my_matrix

rownames(my_matrix) <- c("row1", "row2", "row3", "row4")
my_matrix

t(my_matrix)

## ------------------------------------------------------------------------
ncol(my_matrix)
nrow(my_matrix)
dim(my_matrix)

## ------------------------------------------------------------------------
my_matrix[1,3]
my_matrix["row1", "col3"]

## ------------------------------------------------------------------------
my_matrix[1,]
my_matrix[,3]

## ------------------------------------------------------------------------
my_matrix[col3 > 20,]

## ------------------------------------------------------------------------
my_matrix * 3
log10(my_matrix)

## ------------------------------------------------------------------------
rowSums(my_matrix)

## ------------------------------------------------------------------------
colSums(my_matrix)

## ------------------------------------------------------------------------
apply(my_matrix, 1, mean)

## ------------------------------------------------------------------------
apply(my_matrix, 1, median)

## ------------------------------------------------------------------------
friend_groups <- factor(c(1,2,1,2))
friend_groups

## ------------------------------------------------------------------------
levels(friend_groups)

## ------------------------------------------------------------------------
levels(friend_groups) <- c("best_friend", "not_best_friend")
friend_groups

## ------------------------------------------------------------------------
levels(friend_groups) <- c("not_best_friend", "best_friend")
friend_groups

## ------------------------------------------------------------------------
friend_groups <- factor(c("not_best_friend", "best_friend", "not_best_friend", "best_friend"))
friend_groups

## ------------------------------------------------------------------------
friend_groups <- factor(c("not_best_friend", "best_friend", "not_best_friend", "best_friend"), levels=c("not_best_friend", "best_friend"))
friend_groups

## ------------------------------------------------------------------------
summary(friend_groups)

## ------------------------------------------------------------------------
table(friend_groups)

## ------------------------------------------------------------------------
# creating a data frame using previously defined vectors
friends <- data.frame(name=friend_names, age=friend_ages, child=has_child)
friends

## ------------------------------------------------------------------------
is.data.frame(friends)

## ------------------------------------------------------------------------
is.data.frame(my_matrix)

## ------------------------------------------------------------------------
class(my_matrix)

## ------------------------------------------------------------------------
my_data <- as.data.frame(my_matrix)
class(my_data)

## ------------------------------------------------------------------------
my_data

## ------------------------------------------------------------------------
t(my_data)

## ------------------------------------------------------------------------
friends["Mina",]

## ------------------------------------------------------------------------
# The columns of a data frame can be referred to by the names of the columns
friends
friends$age
friends[friends$age > 26,]

## ------------------------------------------------------------------------
friends[friends$child == "TRUE",]

## ------------------------------------------------------------------------
# select friends that are older than 26
subset(friends, age > 26)

# select the information of the ages of friends
subset(friends, select=age)

## ------------------------------------------------------------------------
# add a column that has the information on the marrital status of friends
friends$married <- c("YES", "YES", "NO", "YES")
friends

## ------------------------------------------------------------------------
# add a column that has the information on the salaries of friends
cbind(friends, salary=c(4000, 8000, 2000, 6000))

## ------------------------------------------------------------------------
my_list <- list(mother="Sophia", father="John", sisters=c("Anna", "Emma"), sister_age=c(5, 10))
my_list

## ------------------------------------------------------------------------
# names of elements in the list
names(my_list)

## ------------------------------------------------------------------------
# number of elements in the list
length(my_list)

## ------------------------------------------------------------------------
my_list$mother

## ------------------------------------------------------------------------
my_list[["mother"]]

## ------------------------------------------------------------------------
my_list[[1]]

## ------------------------------------------------------------------------
my_list[[3]]

## ------------------------------------------------------------------------
my_list[[3]][2]

## ------------------------------------------------------------------------
# If you have downloaded the raw_counts.txt file to your working directory, you may use the following command to read it in.
data <- read.table(file="raw_counts.txt", sep="\t", header=T, stringsAsFactors=F)

# There is a more convenient way to read files from the internet.
data <- read.table(file="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2017-August-Variant-Analysis-Workshop/master/friday/Intro2R/raw_counts.txt", sep="\t", header=T, stringsAsFactors=F)


## ------------------------------------------------------------------------
head(data)

## ------------------------------------------------------------------------
# If you have downloaded the raw_counts.csv file to your working directory, you may use the following command to read it in.
data2 <- read.csv(file="raw_counts.csv", stringsAsFactors=F)

# Otherwise, you may read the file by providing the url to the read.csv() function.
data2 <- read.csv(file="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2017-August-Variant-Analysis-Workshop/master/friday/Intro2R/raw_counts.csv", stringsAsFactors=F)

# To look at the file:
head(data2)

## ------------------------------------------------------------------------
# To write to a file called "output.txt" in your current working directory.
write.table(data2[1:20,], file="output.txt", sep="\t", quote=F, row.names=T, col.names=T)

## ----echo=FALSE, results='asis'------------------------------------------
cat("\\newpage")

## ------------------------------------------------------------------------
# to find out the current working directory
getwd()

# to set a different working directory, use setwd
setwd("/Users/jli/Desktop")

# to list all variables in the environment
ls()

# to create a vector from 2 to 3, usin increment of 0.1
seq(2, 3, by=0.1)

# to create a vector with repeated elements
rep(1:3, times=3)
rep(1:3, each=3)

# to get help information on a function in R: ?function.name()
?seq()
?sort()
?rep()


## ------------------------------------------------------------------------
str(data2)

## ------------------------------------------------------------------------
#?lapply()

data <- as.data.frame(matrix(rnorm(49), ncol=7), stringsAsFactors=F)
dim(data)
lapply(1:dim(data)[1], function(x){sum(data[x,])})
apply(data, MARGIN=1, sum)
lapply(1:dim(data)[1], function(x){log10(sum(data[x,]))})

## ------------------------------------------------------------------------
# To check the syntax of using sapply():
#?sapply()

sapply(1:dim(data)[1], function(x){log10(sum(data[x,]))})

## ------------------------------------------------------------------------
sapply(1:dim(data)[1], function(x){log10(sum(data[x,]))}, simplify=FALSE)

## ----echo=FALSE, results= 'asis'-----------------------------------------
library(knitr)
kable(data.frame(Description=c("Mean", "Standard deviation", "Variance", "Minimum", "Maximum", "Median", "Range of values: minimum and maximum", "Sample quantiles", "Generic function", "Interquartile range"), R_function=c("mean()", "sd()", "var()", "min()", "max()", "median()", "range()", "quantile()", "summary()", "IQR()"), stringsAsFactors=F), align='c')

## ------------------------------------------------------------------------
apply(data, 2, mean)

## ------------------------------------------------------------------------
apply(data, 2, range)

## ------------------------------------------------------------------------
apply(data, 2, quantile)

## ------------------------------------------------------------------------
x <- c(1:50)
y <- 1 + sqrt(x)/2
plot(x,y)

plot(x,y, type="l")

# plot both the points and lines
## first plot points
plot(x,y)
lines(x,y, type="l")
## lines() can only be used to add information to a graph, while it cannot produce a graph on its own.


## ------------------------------------------------------------------------
boxplot(data, xlab="Sample ID", ylab="Raw Counts")

## ------------------------------------------------------------------------
x <- rnorm(1000)
boxplot(x)

## ------------------------------------------------------------------------
hist(x)

# use user defined break points
hist(x, breaks=seq(range(x)[1]-1, range(x)[2]+1, by=0.5))

## ------------------------------------------------------------------------
# clear plotting device/area
dev.off()

## ------------------------------------------------------------------------
source("http://bioconductor.org/biocLite.R")
## install core packages
#biocLite()
## install specific packages
#biocLite("RCircos")
#biocLite(c("IdeoViz", "devtools"))

## ------------------------------------------------------------------------
#install.packages("ggplot2", repos="http://cran.us.r-project.org")

## ------------------------------------------------------------------------
library(devtools)
install_github("stephenturner/qqman")

## ------------------------------------------------------------------------
#biocLite("BiocUpgrade")

## ------------------------------------------------------------------------
savehistory(file="Sept6.history")

#loadhistory(file="Sept6.history")

## ------------------------------------------------------------------------
save(list=c("x", "data"), file="Sept6.RData")

#load("Sept6.RData")
