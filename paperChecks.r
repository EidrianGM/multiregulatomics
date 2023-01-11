
library(openxlsx)
anamatiadata <- "/home/eidriangm/Downloads/mmc2.xlsx"
anamatiaDF <-  read.xlsx(anamatiadata)


View(anamatiaDF)

sum(duplicated(anamatiaDF[4:9]))

View(anamatiaDF[duplicated(anamatiaDF[4:9]),])
