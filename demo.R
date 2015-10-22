library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ccaPP)
library(parallel)

## Pipelines
## Flow of data
## r vs python (FP vs OOP for data programming)


## Shapes of data
##  - data entry
##  - modelling
##  - plotting
##  - tabulating


## EHR data entry
## Problem - How has the way that doctors code diagnoses of medical conditions 
## changed over time?
##  - New regulations
##  - New definitions
##  - "drift"

## Why?
## - improve research
## - make better comparisons between studies
## - understand the effects of policy changes

## data:

## 1. All clinical events in a domain
load("data/hf_incidence.rda")
View(hf_incidence)

## merge with clinical entity names
master <- read.delim("data/medical.txt")
a <- inner_join(hf_incidence, master) # same as base::merge
View(a)

## 2. Frequencies of all clinical codes used per year
load("data/hf_frequencies.rda")
View(hf_incidence_freqs)


## Calculate, evenness, richness and diversity of clinical codes:

conditions <- c("hf", "learning_disability")

## load data into new environment.  This keeps your top-level namespace clean
## For loops are fine in this situation because you have side effects (i.e. you are performing I/O)
## If there are no side effects, use lapply or similar
frequencies <- new.env()
for(condition in conditions){
    load(paste0("data/", condition, "_frequencies.rda"), envir = frequencies)
}

## forget  <- do.call("rbind", lapply(frequencies, function(x) x))
## Use bind_rows - it is cleaner and quicker
all_freqs <- bind_rows(lapply(frequencies, function(x) x))



#' calculates the Shannon index of diversity/entropy
#' I think I stole this from the "vegan" package
shannon <- function(x) 
{
    x <- drop(as.matrix(x))
    if (length(dim(x)) > 1) {
        total <- apply(x, 1, sum)
        x <- sweep(x, 1, total, "/")
    }
    else {
        x <- x/sum(x)
    }
    x <- -x * log(x)
    
    if (length(dim(x)) > 1) 
        H <- apply(x, sum, na.rm = TRUE)
    else H <- sum(x, na.rm = TRUE)
    H
}

# get indices for all years and conditions:
H_all <- all_freqs %>%
    group_by(year, condition) %>% # set these as groups
    summarise(H = shannon(count),   # calculate your statistics for all groups
              richness = n(),
              evenness = H / log(richness)) %>%
    ungroup() %>%
    gather(index, measure, H:evenness) # equivalent to reshape2::melt, need to ungroup first







## this code is exactly equivalent to (without the pipe operator):
## Notice how the previous way, you are reading the flow of comands in a linear way.
## here you need to read the code from the inside out in order to understand what is going on.
gather(
    ungroup(
      summarise(
          group_by(
              all_freqs, year, condition), 
          H = shannon(count),
          richness = n(),
          evenness = H / log(richness))), 
    index, measure, H:evenness)
  


## Best alternative in base R I could come up with!
## cf with line 87 code!

H_all2 <- aggregate(all_freqs$count, 
                    list(year = all_freqs$year, condition = all_freqs$condition), 
                    FUN = shannon)
names(H_all2)[names(H_all2) == "x"] <- "H"
H_all2$richness <- aggregate(all_freqs$count, 
                             list(all_freqs$year, 
                                  all_freqs$condition), 
                             FUN = length)$x
H_all2$evenness <- H_all2$H / log(H_all2$richness)
H_all2 <- reshape(H_all2, v.names = "measure",
                  varying = c("H", "richness", "evenness"),
                  direction = "long")
H_all2$index <- "H"
H_all2$index[ H_all2$time == 2] <- "richness"
H_all2$index[ H_all2$time == 3] <- "evenness"
H_all2$index <- factor(H_all2$index, levels = c("H", "richness", "evenness"))
H_all2$time <- NULL
H_all2$id <- NULL
H_all2 <- H_all2[order( H_all2$index, H_all2$year), 
                 c("year", "condition", "index", "measure")]

all(H_all == H_all2) # return the same values




## compare the two methods using microbenchmark:
library(microbenchmark)

microbenchmark(dplyr = {H_all <- all_freqs %>%
                    group_by(year, condition) %>%
                    summarise(H = shannon(count),
                              richness = n(),
                              evenness = H / log(richness)) %>%
                    ungroup() %>%
                    gather(index, measure, H:evenness)
}, base = {H_all2 <- aggregate(all_freqs$count, 
                        list(year = all_freqs$year, condition = all_freqs$condition), 
                        FUN = shannon)
    names(H_all2)[names(H_all2) == "x"] <- "H"
    H_all2$richness <- aggregate(all_freqs$count, 
                                 list(all_freqs$year, all_freqs$condition), 
                                 FUN = length)$x
    H_all2$evenness <- H_all2$H / log(H_all2$richness)
    H_all2 <- reshape(H_all2, v.names = "measure",
                      varying = c("H", "richness", "evenness"),
                      direction = "long")
    H_all2$index <- "H"
    H_all2$index[ H_all2$time == 2] <- "richness"
    H_all2$index[ H_all2$time == 3] <- "evenness"
    H_all2$index <- factor(H_all2$index, levels = c("H", "richness", "evenness"))
    H_all2$time <- NULL
    H_all2$id <- NULL
    H_all2 <- H_all2[order( H_all2$index, H_all2$year), 
                     c("year", "condition", "index", "measure")]}, times = 500)

# dplyr is 10 x faster:
# Unit: milliseconds
# expr       min        lq    median        uq       max neval
# dplyr  2.396384  2.528917  2.658345  2.837134  5.146062   500
# base 20.728356 21.163055 21.548402 23.728467 75.431676   500


## plot in ggplot2
p <- ggplot(H_all, aes(x = year, y = measure))
p + geom_line() +   # line plot
    facet_grid(index ~ condition, scale = "free_y") +  # indices on the vertical, conditions on the horizontal facets
    expand_limits(y = 0) +    # include 0 in all y axes
    geom_vline(x = 2004, colour = "red", linetype = "longdash") + # pay for performance introduced in 2004
    theme_bw() +   # more minimal theme than the default
    labs(title  = "Clinical code usage") 




#' get robust cannononical correlation for two dataframes using ccaPP package
get_cc <- function(mat1, mat2, ...){ ## dots allow you to pass optional arguments to ccaProj
    x <- as.matrix(mat1)
    y <- as.matrix(mat2)
    ccaProj(x, y,standardize = FALSE, ...) ## calculate correlation
}


for(condition in conditions){
    ## get a database of medcode frequencies by year and practice
    load(paste0("data/", condition, "_incidence.rda"), envir = frequencies)
    
    ## Build a matrix of frequencies n x m
    ## n is the different GP practices
    ## m is the different clinical codes
    freq_all <- get(paste0(condition, "_incidence"), envir = frequencies) %>%
        mutate(practid = str_sub(as.character(patid), -3, -1), # extract practice number from patient id
               medcode = factor(medcode)) %>%
        group_by(medcode, year, practid) %>%
        summarise(count = n()) %>% # count rows for each code,year,practice combination
        ungroup() %>%
        group_by(year, practid) %>% # regroup
        mutate(total = sum(count), 
               freq = count / total) %>% # calculate frequencies for each year and practice
        select(-count, -total) %>%
        ungroup() %>%
        spread(medcode, freq) %>% # put into wide format
        mutate_each(funs(replace(., which(is.na(.)), 0))) # replace NAs with zeros
    
    ## split into a list of matrices, one for each year
    freq_mats <- lapply(2000:2013, function(x){
        freq_all %>%
            filter(year == x) %>% 
            select(-year) # remove year column
    })
    
    ## include only those practices that are common across years
    ## Reduce makes a function that only acccepts 2 arguments into a general vectorised function over a list
    common_practices <- Reduce(intersect, lapply(freq_mats, function(x) x$practid))
    freq_mats <- lapply(freq_mats, function(x){
        x %>%
            filter(practid %in% common_practices) %>%
            select(-practid)
    })
    names(freq_mats) <- 2000:2013
    
    ## combinations for comparisons of matrices:
    combinations <- lapply(list(
        # year-on-year
        c(2000, 2001),
        c(2001, 2002),
        c(2002, 2003),
        c(2003, 2004),
        c(2004, 2005),
        c(2005, 2006),
        c(2006, 2007),
        c(2007, 2008),
        c(2008, 2009),
        c(2009, 2010),
        c(2010, 2011),
        c(2011, 2012),
        c(2012, 2013)), as.character)
    
    ## get correlations for all combinations of frequency matrices
    assign(paste0(condition, "_year_cors"),
           mclapply(combinations, function(this_year){
               mat1 <- freq_mats[[this_year[1]]] # get the matrices for that year
               mat2 <- freq_mats[[this_year[2]]]
               get_cc(mat1, mat2, 
                      method = "spearman")$cor # calculate the correlation
               
           }, mc.cores = 2))
}

## extract the data and bind for all years
all_corrs <- bind_rows(lapply(conditions, function(condition){
    data_frame(year = 2001:2013, # for each condition, get the data for each year 
               condition = condition,
               index = "correlation",
               measure = unlist(get(paste0(condition, "_year_cors"))))
}))

## repeat of above plot with extra correlation data
all_indices <- bind_rows(H_all, all_corrs) # stick the dataframes on top of each other
p <- ggplot(all_indices, aes(x = year, y = measure))
p + geom_line() +
    facet_grid(index ~ condition, scale = "free_y") +
    expand_limits(y = 0) +
    geom_vline(x = 2004, colour = "red", linetype = "longdash") +
    theme_bw() +
    labs(title  = "Clinical code usage") 



