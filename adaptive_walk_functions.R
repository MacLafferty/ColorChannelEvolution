#functions for exploring discrimination spaces with adaptive walks
#by David Morris

#necessary libraries for receptor noise limited modeling, data organization, and graphical visualization
library(pavo)
library(tidyverse)
library(ggplot2)


################################################################################################
###################Receptor Noise Limited Modeling + adaptive walk functions####################
################################################################################################


#function to generate color discriminability values (delta S) given photoreceptor parameters
#myData=Pavo rspec dataset ; Channels = list of peak sensitivities e.g. c(450,550)
#Ratios = the ratio between different photoreceptor types e.g. c(1,1), c(1,1,3) etc.
#expNoise = exponential noise term that causes photoreceptors to be noisier if they're sensitive to longer wavelengths

medianCalc <- function(myData,Channels,Ratios,expNoise=TRUE) {
  #store the upper and lower bands of our data set, lets us dynamically scale our govardovskii templates below
  lowerWave<-myData$wl[1]
  upperWave<-myData$wl[length(myData$wl)]
  #check the length of our channels, adjust based on chromacy level
  numChannels<-length(Channels)
  numRatios<-length(Ratios)
  #if these aren't the same we should report to user, then maybe fix the ratios to match the right length?
  
  #we use the Channels list to make curves from lambda max values w/ govardovskii templates
  sensitivityCurves<-sensmodel(Channels,range = c(lowerWave,upperWave))
  
  #specify parameters to calculate quantum catch:
  #D65 (standard daylight) as illuminant
  #Always using Fechner quantum catch, don't see the benefit of varying given Qi isn't physiologically relevant
  #We don't want to deal with color constancy stuff, hence no VK or background data necessary
  #note that, per pavo advice, the scaling factor is useful when dealing w/ relative values like their D65
  #they recommend 500 for dim light or up to 10k for very bright. I erred on the lower side here
  quantumCatches<-vismodel(rspecdata = myData,visual = sensitivityCurves,illum = "D65",
                           qcatch ="fi",vonkries = F,scale = 1000, relative = F)
  
  #we want to be able to use different noise calculations, so we'll conditionally switch based on parameter
  if (expNoise==TRUE){
    #exponential noise function to penalize longer wavelength photoreceptors due to increased sensitivity
    # noise  = A * B^(X) + C where
    #A= intercept, B= how the curve is shaped, C= base noise
    #X value = i = peak sensitivity of each photoreceptor
    noiseVals<-c()
    for (i in 1:numChannels){
      X<-Channels[i]
      Aval<-.03
      Bval<-1.01
      Cval<-.05
      noiseValue<-Aval * (Bval^(X-550)) + Cval
      noiseVals<-append(noiseVals,noiseValue)
    } 
  } else {
    noiseVals<-.05
  }
  
  #next up we put it all together to measure our standard deviation of noise, or delta S
  #note: NOT JNDs because of no associated behavioral data
  deltaS  <-  coldist(modeldata = quantumCatches,noise = "quantum",n=Ratios,weber=noiseVals)
  
  return(median(deltaS$dS,na.rm = TRUE))
}

#function to propose lambda max shifts for color channels
#iterates through vector of lambda max values
#equally likely to propose positive, negative, or no shift
#eventually needs arguments for distribution information?
#calibrate this stuff from the Hagen review paper
#for single site substitutions + wavelength sensitivity shifts:
#12 instances b/w 1-5nm, 7 shifts b/w 6-10nm, 3 shifts above 11nm
#gamma distribution, adding one to the result such that the lowest bound =1
#alpha = 1 and scale beta = 4 makes shorter shifts more likely, with longer shifts still possible
#i.e. a 1nm shift is ~22% likely while an 8nm shift is ~4%, and so forth
#make our slightly more likely event (.334 instead of .333) the do-nothing case

#we can shift our move proposal such that we only change one channel at a time?
#also we want to include some sort of 'conserved channel' where one doesn't shift, as much or at all?
#maybe just passed a conserved Channel index or something

moveProposal<- function(Channels,Conserved){
  #randomly choose a channel to perturb?
  #maybe have a 
  
  for (i in 1:numChannels){
    Proposal<-runif(1,0,1)
    #no change/negative/positive
    if (Proposal < .333){
      Mod<- -1*floor(rgamma(n=1,shape=1,scale=4) + 1)
    } else if (Proposal < .666){
      Mod<- floor(rgamma(n=1,shape=1,scale=4) + 1)
    } else {
      Mod<-0
    }
    if (Mod > 20){
      Mod <- 20
    }
    #add the proposed sensitivity shift to the previous lambda max
    newMax<-Channels[i]+Mod
    #check to make sure proposed sensitivity isn't outside empirical bounds (i.e. edges of landscape)
    if (newMax > 320 & newMax < 680){
      newChannels[i]<-newMax
    } else {
      newChannels[i]<-Channels[i]
    }
  }
  return(newChannels)
}

moveProposal<- function(Channels,Conserved=0){
  #conditionals to select a channel we want to mutate
  if (Conserved==0){
    #we aren't selecting a conserved channel; all channels are fair game
    whichChannel<-ceiling(runif(1,0,length(Channels)))
    oldMax<-Channels[whichChannel]
  } else if (Conserved > length(Channels)) {
    #we can't conserve a channel we don't have
    print("ya done fucked up, try again")
    whichChannel<-0
    oldMax<-0
  } else{
    #we set aside a conserved channel we won't change
    #first drop the channel we don't intend to mutate
    subset<-Channels[-c(Conserved)]
    #select which of the remaining ends up mutating
    whichChannel<-ceiling(runif(1,0,length(subset)))
    oldMax<-subset[whichChannel]
    #readjust this variable so it'll let us index against the original list of Channels
    if (whichChannel >=Conserved){
      whichChannel<-whichChannel+1
    }
  }
  #decide whether we're increasing or decreasing sensitivity
  Proposal<-runif(1,0,1)
  if (Proposal <= .5){
    Mod<- 1*floor(rgamma(n=1,shape=1,scale=4) + 1)
  } else {
    Mod<- -1*floor(rgamma(n=1,shape=1,scale=4) + 1)
  }
  #hard limit at 20 to reduce any sort of CRAZY big jump
  if (Mod > 20){
    Mod <- 20
  }
  #modify the old lambda max with our move proposal
  newMax<-oldMax+Mod
  #double check we haven't gone outside our boundaries
  if (newMax >= 320 & newMax <= 680){
    Channels[whichChannel]<-newMax
  }
  return(Channels)
}


#function for performing our MCMC-like adaptive walk procedure
#myData= pavo-formatted rspec data
#startLocations = list of starting wavelengths for each channel e.g. c(450,450), c(390,545) etc.
#currently auto-generates a uniform receptor ratio but we could change this?
#numGens= the number of move proposals the simulation will undergo prior to termination

theWanderer<-function(myData,startLocations,numGens,Conserved=0){
  #number of channels determined by length of startLocations vector
  numChannels<-length(startLocations)
  ratios<-rep(1, numChannels)
  #need to make a matrix for the walk locations
  #rows = number of generations + 1 (for initial starting locations)
  #columns = number of channels, putting medians in another vector
  myMatrix<-matrix(0L, nrow = numGens+1, ncol = numChannels)
  #save starting locations to the matrix
  for (location in 1:numChannels){
    myMatrix[1,location]<-startLocations[location]
  }
  #create vector to store the medians
  myMedians<-vector(mode="double",length = numGens+1)
  #determine starting values given starting location
  currentMedian<-medianCalc(myData=myData,Channels=myMatrix[1,],Ratios = ratios,expNoise = TRUE)
  myMedians[1]<-currentMedian
  
  #iterate through the number of generations listed while moving channels
  #starting on index 2 since index 1 are starting values; means we add 1 to numGens
  numGens<-numGens+1
  for (i in 2:numGens){
    #generate new proposed lambda maxes
    newPeaks<-moveProposal(myMatrix[i-1,],Conserved)
    #generate new median from new lambda maxes
    newMedian<-medianCalc(myData=myData,Channels = newPeaks,Ratios=ratios,expNoise = TRUE)
    #compare proposed median to the current best median
    #if it's better, replace it and make the move
    if (newMedian > currentMedian){
      #set new best median and fill matrix with data
      currentMedian<-newMedian
      myMatrix[i,]<-newPeaks
    } else {
      #function for bad move penalty; is it steep enough?
      target<-newMedian^3/ (2*currentMedian^3)
      #sometimes this gets so small it NaNs which messes up our comparison below
      #check for NaN and set to 0 so our comparisons can continue
      if (is.nan(target)){
        target<-0
      }
      
      #if our random draw exceeds our target frequency value, we still accept the move
      #there's some weird error here where very occasionally our first comparison fails
      #how is that possible? I dunno, hopefully we'll figure it out
      roll<-runif(1,0,1)
      if (roll < target){
        #update best median and matrix
        currentMedian<-newMedian
        myMatrix[i,]<-newPeaks
      }else{
        #current median remains the same; no update needed
        #coordinate for current generation remain the same as previous
        myMatrix[i,]<-myMatrix[i-1,]
      }
    }
    myMedians[i]<-currentMedian
  }
  #stick the location matrix and medians together, then return them
  results<-cbind(myMatrix,myMedians)
  results<-as.data.frame(results)
  return(results)
}

#need some sort of wrapper to our wanderer above so we can manage multiple simulations
#needs to be able to read in files/write out files
#check out RCPP


goWandering<-function(myData,numChannels,generations,numWalks,baseString,Conserved=0){
  #so we need to know how many channels we'll be working with; meaning we need at least SOME parameter
  
  #if we're working with tri+ chromacy, then we have some options
  if (numChannels > 2){
    #tri+ chromacy means we need to look at previous output and confirm we can build off that data
    prevNum<-numChannels-1
    checkFile<-paste(baseString,"_",prevNum,"D_",numWalks,sep="")
    #check to see if file is present; if yes, proceed
    if (file.exists(checkFile)){
      print("yay we have the data we need")
      #now we loop
      for (i in 1:numWalks){
        #need to open the correct file and get our starting data
        prevFile<-paste(baseString,"_",prevNum,"D_",i,sep="")
        prevData<-as.data.frame(read.csv(prevFile,header=TRUE))
        lastGen<-length(prevData$myMedians)
        #reformat our starting coordinates for this series of runs
        startingCoords<-as.numeric(prevData[lastGen,-length(prevData)])
        R1<-prevData$V1[lastGen]
        R2<-prevData$V2[lastGen]
        #here's where we need to add in our channel
        if (Conserved==0){
          #we aren't selecting a conserved channel; all channels are fair game
          whichChannel<-ceiling(runif(1,0,numChannels-1))
          updatedCoords<-append(startingCoords,startingCoords[whichChannel],after = whichChannel)
        } else if (Conserved > numChannels) {
          #we can't conserve a channel we don't have
          #maybe change some values and catch any errors here?
          print("ya done fucked up, try again")
          #whichChannel<-0
          #oldMax<-0
        } else{
          #we set aside a conserved channel we won't change
          subset<-startingCoords[-c(Conserved)]
          #select which of the remaining channels ends up duplicating
          whichChannel<-ceiling(runif(1,0,length(subset)))
          if (whichChannel >=Conserved){
            whichChannel<-whichChannel+1
          }
          #stick the value back into the starting coordinates for our run
          updatedCoords<-append(startingCoords,startingCoords[whichChannel],after=whichChannel)
          #readjust this variable so it'll let us index against the original list of Channels
        }
        #we need to pick a random channel, or check conserved, to make a choice about how to change stuff here
        walkOutput<-theWanderer(myData,startLocations = updatedCoords,numGens=generations,Conserved)
        fileName<-paste(baseString,"_",numChannels,"D_",i,sep="")
        write.csv(walkOutput,fileName,row.names=FALSE)
      }
    } else{
      print("you don't have data from lower chromacy, go fix that")
    }
  } else if (numChannels == 2){
    for (i in 1:numWalks){
      #if it's dichromatic space then we dont have previous data to draw from, so we pick random start
      seed<-floor(runif(1,350,651))
      walkOutput<-theWanderer(myData,startLocations = c(seed,seed),numGens = generations,Conserved)
      #write data to file using provided data
      fileName<-paste(baseString,"_",numChannels,"D_",i,sep="")
      write.csv(walkOutput,fileName,row.names=FALSE)
    }
  } else {
    print("numChannels needs to be 2 or greater")
  }
  #once we're done do we put summary stuff in this function? or another?
}

################################################################################################
############################Adaptive Walk Visualization Functions###############################
################################################################################################

#interesting that when starting at green we never make it to the longer wavelengths
#how do we define our plateau?? probably some sort of median index cutoff
findPlateau<-function(walk_data,trimPerc=0){
  keepPerc<-1-trimPerc
  stopIndex<-floor(length(walk_data$myMedians) * keepPerc)
  keptMedians<-walk_data$myMedians[1:stopIndex]
  medMed<-median(keptMedians)
  medBool<-keptMedians > medMed
  firstMed<-match(TRUE,medBool)
  return(firstMed)
}

plateauPermute<-function(myWalk){
  walkTrim<-c()
  for (i in 0:9){
    perc<-i*.1
    walkTrim<-append(walkTrim,findPlateau(myWalk,perc))
  }
  return(rev(walkTrim))
}

#need some visualization functions
walkPanel<-function()

#reworking this so it auto-updates the title of the graph would be nice
plateauFraming<-function(baseString,runNumber){
  plateauFrame<-c()
  
  #make the multi-plot frame; add some more labeling and such; how to save?
  par(mfrow=c(3,4))
  for (i in 1:runNumber){
    fileString<-paste(baseString,i,sep="")
    walkInfo<-read.csv(fileString)
    plot(walkInfo$myMedians,type="l",main=as.character(i))
    walkPlateau<-plateauPermute(walkInfo)
    plateauFrame<-rbind(plateauFrame,walkPlateau)
  }
  #do some relabeling and reformatting
  percentLabel<-c("10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")
  percentOrdered<-ordered(percentLabel, levels = percentLabel)
  plateauFrame<-as.data.frame(plateauFrame)
  colnames(plateauFrame)<-percentOrdered
  plateauMelt<-plateauFrame %>% tibble::rownames_to_column() %>% gather(colname, value, -rowname)
  plateauMelt$colname<-ordered(plateauMelt$colname,levels(percentOrdered))
  #now graph it!
  plateauMap<-ggplot(data=plateauMelt, aes(x=colname, y=rowname, fill = value, label = value)) +
    geom_tile() +
    geom_text(col = "black") +
    theme_minimal() +
    scale_fill_gradient2(low = "white", mid = "yellow", high = "red") +
    scale_color_gradient2(low = "white", mid = "yellow", high = "red") +
    labs(title = "Adaptive Walks - 2000 Moves", x= "Percent Cut",y = "Adaptive Walks",fill="Threshold")
  return(plateauMap)
}

#we're gonna need: (melted) landscape, file w/ walk information, plot output name
landscape_overlay<-function(melt_frame,plot_title,walk_file=""){
  
  walk_data<-read.csv(walk_file)
  
  ourPlot<-ggplot(data=melt_frame, aes(x=rowname,y=colname,fill=value)) +
    geom_tile(color="white") + scale_fill_gradient(low="cadetblue2",high="red") +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    geom_point(data=walk_data,aes(x=V1,y=V2,color=myMedians),inherit.aes=F) +
    labs(title = plot_title)
  
  return(ourPlot)
}