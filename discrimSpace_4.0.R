#exploration of color channel evolution
#by David Morris

#load libraries
library(pavo)
library(lattice)

#the most fundamental procedure we do 
#is calculate median discriminability of environmental stimuli with a given parameter set
#parameters of interest are: stimuli, background, number of channels, receptor ratios
#maybe also have a flag for dark noise?

##############################################
################STEP ONE######################
##############################################

medianCalc <- function(myData,Channels,Ratios,expNoise=TRUE) {
  #store the upper and lower bands of our data set, lets us dynamically scale our govardovskii templates below
  lowerWave<-myData$wl[1]
  upperWave<-myData$wl[length(myData$wl)]
  #check the length of our channels and such, allows the method to work as we change the number of channels
  numChannels<-length(Channels)
  numRatios<-length(Ratios)
  #if these aren't the same we should throw an error
  
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

#function seems to do its deal which is nice
test1<-medianCalc(myData = Silvretta,Channels = c(400,550),Ratios = c(1,1),expNoise = TRUE)
test2<-medianCalc(myData = Silvretta,Channels = c(400,550),Ratios = c(1,1),expNoise = FALSE)

#next up, we can use our median values to:
#describe entire discrimination landscapes
#walk through those landscapes

#we need specified data as well as the range of wavelengths we want to bound the landscape
generateLandscape<- function (spectralData,bounds=c(320,680),Ratios=c(1,1),Noise=TRUE){
  #create an array to store data
  numWaves<-(bounds[2]-bounds[1])+1
  storageArray<-array(bounds[1]:bounds[2],c(numWaves,numWaves))
  #I know, loops in R....
  for (i in 1:numWaves){
    for (j in i:numWaves){
      #specify the lambda max for each channel
      channels<-c(i,j)
      #we can  only get away with  our  trick of using one value across the matrix diagonal if symmetrical
      #so we'll do  two calculations  for each wavelength pairing and store each
      #make sure to include receptor stuff via parameter soon
      medianDeltaS1<-medianCalc(spectralData,c(i,j),Ratios,Noise)
      medianDeltaS2<-medianCalc(spectralData,c(j,i),Ratios,Noise)
      #make sure the ordering works
      storageArray[i,j]=medianDeltaS1
      storageArray[j,i]=medianDeltaS2
      
    }
  }
  return(storageArray)
}

#with assymetric matrices we'll rarely want only half the matrix, but it  could be useful?
#comment these later
DiagonalExtraction<-function(myMatrix){
  receptorOne<-c()
  receptorTwo<-c()
  discrimValues<-c()
  
  myColumns<-ncol(myMatrix)
  myRows<-nrow(myMatrix)
  counts<-0
  for (i in 1:myRows){
    i
    for (j in i:myColumns){
      receptorOne[counts]<-i+299
      receptorTwo[counts]<-j+299
      discrimValues[counts]<-myMatrix[i,j]
      counts<-counts+1
    }
  }
  arrayBind<-cbind(receptorOne,receptorTwo,discrimValues)
  arrayBind<-as.data.frame(arrayBind)
  colnames(arrayBind)<-c('X','Y','Z')
  return(arrayBind)
}

FullMatrixExtraction<-function(myMatrix,startWL){
  receptorOne<-c()
  receptorTwo<-c()
  discrimValues<-c()
  
  myColumns<-ncol(myMatrix)
  myRows<-nrow(myMatrix)
  counts<-0
  for (i in 1:myRows){
    i
    for (j in i:myColumns){
      receptorOne[counts]<-i+startWL
      receptorTwo[counts]<-j+startWL
      discrimValues[counts]<-myMatrix[i,j]
      receptorOne[counts+1]<-j+startWL
      receptorTwo[counts+1]<-i+startWL
      discrimValues[counts+1]<-myMatrix[j,i]
      counts<-counts+2
    }
  }
  arrayBind<-cbind(receptorOne,receptorTwo,discrimValues)
  arrayBind<-as.data.frame(arrayBind)
  colnames(arrayBind)<-c('X','Y','Z')
  return(arrayBind)
}

#load some data
setwd("/Users/davidmorris/Biology/ColorModeling/StimulusSets/Silvretta_Austria_FReD/modifiedSpectra")
Silvretta_raw<-getspec(where=getwd(),ext="csv",lim=c(300,700),sep = ",")
Silvretta<-procspec(as.rspec(Silvretta_raw), opt="smooth", fixneg="zero",span=.05)
#run landscape
#first let's do symmetrical landscape w/ noise
Silvretta_Matrix<-generateLandscape(Silvretta)
write(Silvretta_Matrix,file="Silvretta_forest",ncolumns = 381,sep=",")
#rename
Silvretta_Forest<-Silvretta_Matrix
Silvretta_Forest_extract<-FullMatrixExtraction(Silvretta_Forest,300)
wireframe(Silvretta_Forest_extract$Z~Silvretta_Forest_extract$X*Silvretta_Forest_extract$Y,shade=T,xlab='Receptor1',ylab='Receptor2',zlab='deltaS',main='Silvretta_Forest')
#next up, let's do bluesky
Silvretta_Blue<-generateLandscape(Silvretta)
write(Silvretta_Blue,file="Silvretta_blue",ncolumns = 381,sep=",")
Silvretta_Blue_extract<-FullMatrixExtraction(Silvretta_Blue,300)
wireframe(Silvretta_Blue_extract$Z~Silvretta_Blue_extract$X*Silvretta_Blue_extract$Y,shade=T,xlab='Receptor1',ylab='Receptor2',zlab='deltaS',main='Silvretta_Bluesky')
#last for the regular ratio is D65
Silvretta_day <-generateLandscape(Silvretta)
write(Silvretta_day,file="Silvretta_day",ncolumns=381,sep=",")
Silvretta_day_extract<-FullMatrixExtraction(Silvretta_day,300)
wireframe(Silvretta_day_extract$Z~Silvretta_day_extract$X*Silvretta_day_extract$Y,shade=T,xlab='Receptor1',ylab='Receptor2',zlab='deltaS',main='Silvretta_D65')

#for kicks we should shift the background too and confirm effects
Silvretta_whiteBG<-generateLandscape(Silvretta)
write(Silvretta_whiteBG,file="Silvretta_whiteBG",ncolumns = 381, sep=",")
Silvretta_whiteBG_extract<-FullMatrixExtraction(Silvretta_whiteBG,300)
wireframe(Silvretta_whiteBG_extract$Z~Silvretta_whiteBG_extract$X*Silvretta_whiteBG_extract$Y,shade=T,xlab='Receptor1',ylab='Receptor2',zlab='deltaS',main='Silvretta_whiteBG')
#now it's time to unbalance the receptors
Silvretta_ratio<-generateLandscape(Silvretta)
write(Silvretta_ratio,file="Silvretta_1:2",ncolumns = 381, sep=",")
SilvrettaOneToTwo<-FullMatrixExtraction(Silvretta_ratio,300)
wireframe(SilvrettaOneToTwo$Z~SilvrettaOneToTwo$X*SilvrettaOneToTwo$Y,shade=T,xlab='Receptor1',ylab='Receptor2',zlab='deltaS',main='Silvretta_1:2')
#let's push it higher
Silvretta_ratio2<-generateLandscape(Silvretta)
write(Silvretta_ratio2,file="Silvretta_1:3",ncolumns = 381,sep=",")


MallesMatrix<-generateLandscape(Malles)
MallesMatrix_forestshade<-generateLandscape(Malles)
MallesFrame<-FullMatrixExtraction(MallesMatrix,300)
write(MallesMatrix,file="MallesData",ncolumns = 401,sep=",")
write(MallesNoiseMatrix,file="MallesNoise_forEvo",ncolumns = 401,sep=",")


