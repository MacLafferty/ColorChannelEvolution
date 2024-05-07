#adaptive walk exploration
#first we should grab our adaptive walk functions we use for testing
source("adaptive_walk_functions.R")

#here's a little log of the various wandering we've done
#runs 1-3 were done on the 25 spectra strausberg set from FReD
#v1 was 500 moves; v2 was 2000; v3 scaled the number of walks from 5 to 10
#v4 I made some big model changes, shifting to one channel mutating at a time, and increased to 3k moves
#also scaled up to the 50 spectra strausberg set
#I screwed up directionality here, fixed for v5
#also realized the channel addition moving to tri-chromat was being funny so I addressed that before the 3D switch



#load up our current data set
setwd("/Users/davidmorris/Biology/ColorModeling/StimulusSets/StrausbergLarge_FReD/ModifiedSpectra/")
large_strauss_raw<-getspec(where=getwd(),ext="csv",lim=c(300,700),sep = ",")
large_strauss<-procspec(as.rspec(large_strauss_raw), opt="smooth", fixneg="zero",span=.05)
  
#first, let's save some model output from our adaptive walk approach
w1<-medianCalc(large_strauss,Channels = c(350,450),Ratios = c(1,1),expNoise = TRUE)  
w2<-medianCalc(large_strauss,Channels = c(370,430),Ratios = c(1,1),expNoise = TRUE)
w3<-medianCalc(large_strauss,Channels = c(345,650),Ratios = c(1,1),expNoise = TRUE)

#now reload medianCalc from discrimSpace
d1<-medianCalc(large_strauss,Channels = c(350,450),Ratios = c(1,1),expNoise = TRUE) 
d2<-medianCalc(large_strauss,Channels = c(370,430),Ratios = c(1,1),expNoise = TRUE) 
d3<-medianCalc(large_strauss,Channels = c(345,650),Ratios = c(1,1),expNoise = TRUE)

#let's try loading in our current dataset landscape, and then spot check
setwd("/Users/davidmorris/Biology/ColorModeling/landscapes/")
straus_large_landscape<-read.csv("Large_Strausberg_expNoise",header = FALSE)
str(straus_large_landscape)

#so this spot should show the same value as w1/d1, but doesn't; diagonals are consistent
l1<-straus_large_landscape[30,130]
l1.2<-straus_large_landscape[130,30]
#next check
l2<-straus_large_landscape[50,110]
l3<-straus_large_landscape[25,330]

#the values aren't scaled by a consistent factor, I guess we try rerunning this landscape
#see where the weirdness creeps in?

  
#I wanna double check the landscape v walk models and see that output is the same
#test that stsuff here

#we wanna have our wandering function be able to point at previous simulations and use those files to start
#so we output 2D files first
goWandering(myData=large_strauss,numChannels=2,generations=3000,numWalks=10,baseString="test_v5")

#if we're doing higher chromacy, we need to build off the previous files
goWandering(myData=StrausbergLarge,numChannels=3,generations=3000,numWalks=10,baseString="test_v5")



#grab our information from our output files for plotting how we're doing
round1.1<-plateauFraming("test_v1_2D_",5)
round1.2<-plateauFraming("test_v1_3D_",5)
round2.1<-plateauFraming("test_v2_2D_",5)
round2.2<-plateauFraming("test_v2_3D_",5)
round3.1<-plateauFraming("test_v3_2D_",10)
round3.2<-plateauFraming("test_v3_3D_",10)
round4.1<-plateauFraming("test_v4_2D_",10)
round4.2<-plateauFraming("test_v4_3D_",10)
round5.1<-plateauFraming("test_v5_2D_",10)

#set things up so we can do our walk overlap more easily? anyway, grab landscape and do the thing
#rework a bunch of this into the adaptive walk functions
setwd("/Users/davidmorris/Biology/ColorModeling/landscapes/")
Straus_scape_Large<-as.data.frame(read.csv("Large_Strausberg_expNoise",header = FALSE))
#fix labels so we're working with proper wavelengths, from 300-699nm
colnames(Straus_scape_Large)<-320:680
rownames(Straus_scape_Large)<-320:680
#melt to prepare for heatmap
#double check to see whether different landscape values default to strings or not
Straus_melt<-Straus_scape_Large %>% tibble::rownames_to_column() %>% gather(colname, value, -rowname)
Straus_melt$rowname<-as.numeric(Straus_melt$rowname)
Straus_melt$colname<-as.numeric(Straus_melt$colname)

#load up our file of interest; consider function-ifying our heat map approach?
#work in the new 3D channel stuff here

t4.1<-landscape_overlay(Straus_melt,"t1",walk_file="test_v4_2D_1")
t4.2<-landscape_overlay(Straus_melt,"t2",walk_file="test_v4_2D_2")
t4.3<-landscape_overlay(Straus_melt,"t3",walk_file="test_v4_2D_3")
t4.4<-landscape_overlay(Straus_melt,"t4",walk_file="test_v4_2D_4")
t5.1<-landscape_overlay(Straus_melt,"t5.1",walk_file = "test_v5_2D_1")
t5.2<-landscape_overlay(Straus_melt,"t5.2",walk_file = "test_v5_2D_2")
t5.3<-landscape_overlay(Straus_melt,"t5.3",walk_file = "test_v5_2D_3")
t5.4<-landscape_overlay(Straus_melt,"t5.4",walk_file = "test_v5_2D_4")
t5.5<-landscape_overlay(Straus_melt,"t5.5",walk_file = "test_v5_2D_5")
t5.10<-landscape_overlay(Straus_melt,"t5.10",walk_file = "test_v5_2D_10")

B1Plot<-ggplot(data=Straus_melt, aes(x=rowname,y=colname,fill=value)) +
  geom_tile(color="white") + scale_fill_gradient(low="cadetblue2",high="red") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  geom_point(data=blue1,aes(x=V1,y=V2,color=myMedians),inherit.aes=F) +
  labs(title = "v4_test")
