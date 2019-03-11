###try some real data
##for this exercise, use first only behavior observed. use LU as forst variable

Raw<-read.csv("~/Google Drive/FuturAves/Occupancy_Model/data/JZBirdSurveys_2017_rep.csv")

#remove those out of range
Raw<-Raw[Raw$Distance!="+",]

#find most observed species - RNWR, rufous-naped wren
Raw<-Raw[Raw$Species=="RNWR",]

#repeat rows by number of observations
summary()