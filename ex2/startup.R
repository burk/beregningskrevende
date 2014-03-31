library(spam) # load the data
str(Oral) # see structure of data
load('Rmat.Rdata')
#'data.frame': 544 obs. of 3 variables:
# $ Y : int 18 62 44 12 18 27 20 29 39 21 . . .
# $ E : num 16.4 45.9 44.7 16.3 26.9 . . .
# $ SMR: num 1.101 1.351 0.985 0.735 0.668 . . .
attach(Oral) # allow direct referencing to Y and E
# generate some plots
library(fields, warn.conflict=FALSE)
library(colorspace)
col <- diverge_hcl(8) # blue - red
# alternative colors
# col <- rev(gray(0:8 / 8)) # gray scales
# col <- rev(heat_hcl(64))
# use the function provided by spam
map.landkreis(log(Oral$Y),col=col)
map.landkreis(Oral$Y/Oral$E,col=col)