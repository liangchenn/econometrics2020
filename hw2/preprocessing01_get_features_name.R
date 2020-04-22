# Preprocessing for hw2.R
# get the group data variable names

# input string for var. names
string <- "age
male
black
Asian
hisp
race other
both par
less hs
more hs
momedu miss
welfare
momjob miss
prof
job other
sport
white
yr school
gpa
overage"


feature.to.use <- "age, male, black, Asian, hisp, race other, less hs, more hs, momedu miss, welfare, momjob miss, prof, job other" %>%
  strsplit(., split = ', ') %>%
  `[[`(., 1)
# split the string, extract element from output list object
feature.names <- `[[`(strsplit(string, split = '\n'), 1)

# keep feature.names
rm(string)

