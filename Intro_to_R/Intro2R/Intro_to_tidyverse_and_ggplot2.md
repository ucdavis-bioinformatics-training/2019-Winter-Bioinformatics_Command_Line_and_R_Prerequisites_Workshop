---
title: "An introduction to tidyverse and GGPLOT2"
author: "Bioinformatics Core"
output: 
  html_document:
    keep_md: TRUE
---



### Getting started with Tidyverse

Tidyverse is a set of package for doing data science.  [R for Data Science](https://r4ds.had.co.nz/index.html) by Garrett Grolemund, Hadley Wickham model Data Science in the following way.

![](./Intro_to_tidyverse_and_ggplot2_images/R-data-science.png)


We will start learning about Tidyverse tools by starting at the first step in this process, importing data.

*** 

### Step 1: Import data with the [readr](https://readr.tidyverse.org/) package ![](./Intro_to_tidyverse_and_ggplot2_images/readr.png){width=100px}


> “The goal of 'readr' is to provide a fast and friendly way to read rectangular data (like 'csv', 'tsv', and 'fwf'). It is designed to flexibly parse many types of data found in the wild, while still cleanly failing when data unexpectedly changes.”

The readr package gets loaded automatically when you use library(tidyverse), or you can load it directly.

```r
library(readr)
```

***


#### readr supports a number of file formats with different read_* functions including:

* read_csv(): comma separated (CSV) files
* read_tsv(): tab separated files
* read_delim(): general delimited files (you must supply delimiter!)
* read_fwf(): fixed width files           
* read_table(): tabular files where columns are separated by white-space.
* read_log(): web log files

#### Readr also has functions write data in a number of formats with various write_* functions:

* write_csv(): comma separated (CSV) files
* write_tsv(): tab separated files
* write_delim(): general delimited files
* write_excel_csv(): comma separated files for Excel

*** 

#### Get some data and try out these functions:

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-Winter-Bioinformatics_Command_Line_and_R_Prerequisites_Workshop/master/intro2R/mpg.tsv", "mpg.tsv")
```


The file has a ".tsv" extension, so it is probably a tab separated values file. Lets check this assumption (note that the output from the system() function does not appear in the markdown document, and this approach may not work on windows computers):

```r
getwd()
```

```
## [1] "/bio/CoreWork/workshops/2019-Winter-Bioinformatics_Command_Line_and_R_Prerequisites_Workshop/Intro_to_R/Intro2R"
```

```r
dir(pattern="*.tsv")
```

```
## [1] "break_readr.tsv" "mpg.tsv"
```

```r
system('head mpg.tsv')

system('wc -l mpg.tsv')
```

Alternatively we could use another of the readr functions, **read_lines**, to look at the first few lines of the file:

```r
read_lines('mpg.tsv', n_max = 5)
```

```
## [1] "manufacturer\tmodel\tdispl\tyear\tcyl\ttrans\tdrv\tcty\thwy\tfl\tclass"
## [2] "audi\ta4\t1.8\t1999\t4\tauto(l5)\tf\t18\t29\tp\tcompact"               
## [3] "audi\ta4\t1.8\t1999\t4\tmanual(m5)\tf\t21\t29\tp\tcompact"             
## [4] "audi\ta4\t2\t2008\t4\tmanual(m6)\tf\t20\t31\tp\tcompact"               
## [5] "audi\ta4\t2\t2008\t4\tauto(av)\tf\t21\t30\tp\tcompact"
```

*How many lines does the file have?*


*What is the first line of the file?*


*What separates the values in the file?*

Read the file and store it in an object:

```r
mpg <- read_tsv('mpg.tsv')
```

```
## Parsed with column specification:
## cols(
##   manufacturer = col_character(),
##   model = col_character(),
##   displ = col_double(),
##   year = col_double(),
##   cyl = col_double(),
##   trans = col_character(),
##   drv = col_character(),
##   cty = col_double(),
##   hwy = col_double(),
##   fl = col_character(),
##   class = col_character()
## )
```

#### Column specifications

In order to read a file, readr skims the [first 1000 lines](https://readr.tidyverse.org/articles/readr.html) of the file, investigating the values it finds there, and building a guess at the format of the file.

*** 



d = read_csv(readr_example("challenge.csv"), guess_max = 2000)


### Step 1.5: Detour for [Tibbles](https://tibble.tidyverse.org/) ![](./Intro_to_tidyverse_and_ggplot2_images/hex-tibble.png){width=100}

Tibbles are 


TODO



#### Readr and Tibble Exercises 

1) In the R code chunk below try to generate a *parsing failure* in readr. 
Make a trecherous tibble. Write it out. Read it in again.



2) Take a look at the documentation for read_delim. Enter R code below that successfully loads the trecherous tibble you created in the last exercise.




### Step 2: Tidying data with [tidyr](https://tidyr.tidyverse.org/) package![](./Intro_to_tidyverse_and_ggplot2_images/hex-tidyr.png){width=100}

### tidyr


### dplyr

### stringr


### ggplot2


### Real data example


