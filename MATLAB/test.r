##-------------------------------------------------------------------------------------------------------------------##
##-------------------------------------------------------------------------------------------------------------------##
##-------------------------------------------------------------------------------------------------------------------##
##-------------------------------##
##---teasing out tim's program---##
##-------------------------------##

## The variable Stored contains the formant values, the VT start and
## end values and the vowel label, will append new values on previous
## ones - but make sure you press the save formantvalue button.


##-------------------------------------------------------------------------------------------------------------------##

### function to join up tubes - I think this will work - haven't tested it yet.
"all.tubes"<-function(tdata,first= TRUE,olddata=NULL)
{
#takes in the tubeMA1Selection variable from Tim's work and appends things.
#find length
len<-length(tdata$x)
num<-nrow(olddata$index)
#sort out new index
if(first)
{ st= 1} else {st= olddata$index[num,2]+1}

en = st+len-1
if(first)
{list(x=tdata$x,y=tdata$y,length=len,index=matrix(c(st,en),ncol=2))}
else
{list(x=c(olddata$x,tdata$x),y=c(olddata$y,tdata$y),length=c(olddata$length,len),index=rbind(olddata$index,c(st,en)))}
}


##-------------------------------------------------------------------------------------------------------------------##

### Now write a plotting program
"plot.vt"<-function(tubes,ind,ver=4,ttl= "")
{

for(i in 1:ver)
{

	plot(tubes$x[tubes$index[ind,1]:tubes$index[ind,2]]+1,tubes$y[tubes$index[ind,1]:tubes$index[ind,2]],type="l",ylab="cross-sectional area(cm-sq)",xlab="from lips (cm)",xlim=c(-0.5,20.5),ylim=c(0,16.5),col=i)
	if(i!=ver)
	{
		par(new= T)
		ind=ind+1
	}
}
title(ttl)
}


##-------------------------------------------------------------------------------------------------------------------##
##-------------------------------------------------------------------------------------------------------------------##
##-------------------------------------------------------------------------------------------------------------------##
##------------------------##
##---acoustic.tubes.new---##
##------------------------##

## after a whole lot of testing I think there was something wrong in
## the both the rc2lpc and lpc2spec functions. I have rewritten them
## and will include them here, along with the files that are needed
## prior to this, and along with the results. I am keeping the old
## versions just incase. But from now on I think I should use these
## functions.

## MARCH 9TH 2007


##-------------------------------------------------------------------------------------------------------------------##

calc.reflection.coef <- function(areadata,diameter=FALSE)
{
# calculates reflection coefficents from the cross-sectional
# areas. The cross-sectional areas are ordered from the lips to the
# glottis.

# Need to remember the last tube is the area at the glottis. Which is
# zero for a closed tube

# If vocal tract areas are diameters of sections ( diameter=TRUE),
# convert to cross-sectional areas

if(diameter)
	{
	areadata <- ((areadata/2)**2)*pi
	}

# The glottis is modelled as a tube of area 0. This provides a lossless tube model
areadata <- c(areadata, 0)

#make vector for reflection coeff.
refl.coef <- vector(mode="numeric",length=length(areadata))

# first example(command)need to get reflection coefficent between lips and outside, so
# set areadata[0]to be area outside lips and make this infinity. This
# give the first reflection coefficent as 1 ( substitute into formulae
# below to check)

refl.coef[1] <-1

# calculate reflection coefficent between lips and glottis
for(i in 2:length(areadata))
      {
      refl.coef[i] <- (areadata[i-1]-areadata[i])/(areadata[i-1]+areadata[i])
      }

return(refl.coef)
}


##-------------------------------------------------------------------------------------------------------------------##

rc2lpc.II <-function(reflcoef)
{
# this function converts reflection coefficients to linear prediction
# coefficents using a recursive algorithm. This version uses the algorithm
# in JMH's book p216
# modified march 2nd 2007.

# note this function expects that the first reflection coeffcient is 1
if(reflcoef[1]!=1)
{
cat("reflcoef coef not in expected format, first coef must be 1")
return()
}

# store values in matrix - recursive process which requires values
# from previous iteration
lpc.order <- length(reflcoef)
lin.pred.coef <- matrix(nrow=lpc.order, ncol=lpc.order, byrow=TRUE)


#recursive process
lin.pred.coef[1,1] <- reflcoef[1]

#the order p is one less than the length of rc vector
for(j in 2:lpc.order)
{
lin.pred.coef[j,j] <- reflcoef[j]	
	 for(i in 1:j-1)
	 {
	 lin.pred.coef[j,i] <- lin.pred.coef[j-1,i]-(reflcoef[j]*lin.pred.coef[j-1,j-i])
	 }
}

#only return last row of matrix
return(lin.pred.coef[lpc.order,])
}


##-------------------------------------------------------------------------------------------------------------------##

lpc2spec.II <-function(filtcoef,n=128)
{
# calculates spectrum from lpc coef. Could just zero extend and put
# into fft, but wanted to evalue my self
# modified 2/3/07 -
# filtercoef = number of LPC coef
# n = number of steps around ½ the unit circle


# make a vector of the negative powers of z
order <-1:length(filtcoef)

# make a row vector of the n steps around 1/2 the unit circle. Note: R
# vector default to be column vector, want this to be a row vector so
# must transpose.
tn <-seq(0,0.5,length=(n+1))
tn<-t(tn[1:n])
# calculate the power of the exponential(note how you make it imaginary)
# NOTE: MATRIX MULTIPLICATION
# "j times 2 pi times (0:0.5 in steps of n)" remember this has been normalized 
zpow <- 2i*pi*tn

# now solve the formulae H(z)= 1/A(z), where
# A(z) = 1 - SUM(i=1 to order)a[i]*z^-i
# R coerses the filtercoef vector so it can be multiplied by the z vector
hz <- 1/(1-(filtcoef%*%exp(-order%*%zpow)))

# return the absolute value of this
return(Mod(hz))
}


##-------------------------------------------------------------------------------------------------------------------##

lpc2spec.III <-function(filtcoef,n=128)
{

# calculates spectrum from lpc coef. Could just zero extend and put
# into fft, but wanted to evalue my self
# modified 2/3/07 -
# filtercoef = number of LPC coef
# n = number of steps around ½ the unit circle


# make a vector of the negative powers of z
order <-1:length(filtcoef)

# make a row vector of the n steps around 1/2 the unit circle. Note: R
# vector default to be column vector, want this to be a row vector so
# must transpose.
tn <-seq(0,0.5,length=(n+1))
tn<-t(tn[1:n])
# calculate the power of the exponential(note how you make it imaginary)
# NOTE: MATRIX MULTIPLICATION
# "j times 2 pi times (0:0.5 in steps of n)" remember this has been normalized 
zpow <- 2i*pi*tn

# now solve the formulae H(z)= 1/A(z), where
# A(z) = 1 - SUM(i=1 to order)a[i]*z^-i
# R coerses the filtercoef vector so it can be multiplied by the z vector
hz <- (1-exp(-zpow))/(1-(filtcoef%*%exp(-order%*%zpow)))

# return the absolute value of this
return(Mod(hz))
}


##-------------------------------------------------------------------------------------------------------------------##

"area2spec"<-function(areadata,diameter=F,vtlength,doPlot=T,n=512,limit=F,lim=NULL,ylim=NULL, col = "black")
{
#note for lossless tube, the area of the glottis is 0, i.e  Area[p]=0
#calculate reflection coef from area data
	rc<-calc.reflection.coef(areadata,diameter)

# need an even number of areacoef, other wize don't have pole pairs,
# let user know if the number of areacoef is odd)
  if(length(rc)%%2)
 	cat("WARNING: Odd number of area coefficients \n")

# calculate the lpc coef
	coef<-rc2lpc.II(rc)

# calculate the spectrum from the linear prediction coef. See notes in
# acoustic.tubes.new.txt reason for lpc2spec.III

  spec<-lpc2spec.III(coef,n)

# calculate sampling frequency from 
# fs=(speed of sound in cm^2/s)*(No. tubes)/(2*vtlength)

  fs <-34000*length(areadata)/(2*vtlength)
	
  if(doPlot)
  {
	if(limit)
	{
		#find the index of the limit,this rounds down
		nlim<-as.integer(lim*2*n/fs)
		if(lim%%2) nlim=nlim+1
		x <- (1:nlim)*fs/(2*n)
		plot(x,log(spec[1:nlim]),type="l",axes=T,xlab="",ylab="",ylim=ylim, col = col)
		#axis(1,c(1,nlim/2,nlim),c(0,ceiling(nlim*fs/(4*n)),ceiling(nlim*fs/(2*n))))
		#axis(2)
	}
	else
	{
		x <- (1:n)*fs/(2*n)
		plot(x,log(spec),type="l",axes=T,xlab="",ylab="",ylim=ylim, col = col)
		#axis(1,c(1,n/2,n),c(0,fs/4,fs/2))
		#axis(2)		
	}
  }
  return(list(rc=rc,lpc=coef,spec=spec,freqs=x,fs=fs))
}


##-------------------------------------------------------------------------------------------------------------------##

## I check ed out the above function on the same two wakita vowels -
## looked the same. Then I did it on /i/ in the Fant's data listed
## p255 in Lass Book.

# Things are looking good - although am a bit worried by the NaN in
# the first value- but think that would be caused by 0/0 (i.e the
# effect of mulitplying 1/(A(z) by (1-z^-1)

# Need to now start thinking about what data to use - Fant's Wakitas
# etc ?, and the number of tubes to use - Fant had more than Wakita -
# but I can always interpolate Fant's tubes.  Do this on Friday


##-------------------------------------------------------------------------------------------------------------------##
##-------------------------------------------------------------------------------------------------------------------##
##-------------------------------------------------------------------------------------------------------------------##
##-----------------##
##---Summer Code---##
##-----------------##

# Written by Duncan Eason - Summer 08/09
# Based on code written by Tim Judkins - Summer 07/08
# R based tcltk GUI for analysis of pharangometer files

# When storing formant data, the tract start and ends, formant values and 
# labels are accessible in the variable 'Stored' under 'Stored$FormantValues',
# 'Stored$Labels', 'Stored$HVD', 'Stored$VTStart' and 'Stored$VTEnd'


##-------------------------------------------------------------------------------------------------------------------##
##---Loading Libraries---##

require(tcltk)
require(tkrplot)
require(emu)
require(sound)
addTclPath(path="C:\\Tcl\\lib\\bwidget1.8")  #put by CIW 15/4/08
tclRequire("BWidget")


##-------------------------------------------------------------------------------------------------------------------##
##---Initialising---##

imgFull <- NULL
imgSelection <- NULL
imgTubes <- NULL
imgSpectrum <- NULL

# Data variables
D <- NULL
DAve <- NULL
MA <- NULL
MAAve <- NULL
MAAdjPos <- NULL
MAAdjNeg <- NULL
StdDv <- NULL
StdDvAve <- NULL
tubesListInitial <- list(NULL, NULL, NULL, NULL, NULL)
tubesListTcling <- list(NULL, NULL, NULL, NULL, NULL)
vtlengthList <- list(NULL, NULL, NULL, NULL, NULL)
Stored <- NULL

# Misc variables
setNumber <- 0
loaded <<- FALSE
beenModified <- c(F, F, F, F, F)
initialSpeechSampleImpulse <- list(NULL, NULL, NULL, NULL, NULL)
initialSpeechStoredImpulse <- c(F, F, F, F, F)
initialSpeechSampleGlottus <- list(NULL, NULL, NULL, NULL, NULL)
initialSpeechStoredGlottus <- c(F, F, F, F, F)
showInitial <<- TRUE

# The tract will be slipt later into a number of 'tubes'
# tubeNumber defines how many parts to INITIALLY split it into and needs to be odd 
# to work with area2spec.
tubeNumber <<- 11

##-------------------------------------------------------------------------------------------------------------------##
##---Necessary functions---##

# make greatest common divisor and lowest common multiple methods (R doesn't have?)
gcd <- function(a,b) ifelse (b==0, a, gcd(b, a %% b))
lcm <- function(a,b) ((a/gcd(a,b))*b)

# Function to reset bindings to do nothing
doNothing <- function() {}


##-------------------------------------------------------------------------------------------------------------------##
##---Data loading---##

# Loads the data from the file from the pharangometer.
# 'D' the distances, 'MA' the areas and 'StdDv' the standard deviations.  
# Each are aranged as a matrix with the 4 different columns as the 4 data sets.
LoadData <- function() 
{
	fileName <- tclvalue(tkgetOpenFile()) 
	if (!nchar(fileName)) 
		{
       	tkmessageBox(title = "Error", message = "No file was selected.")
		return()
        	}
	a <- read.table(fileName, nrows = 84, skip = 1)
	n <<- scan(fileName, what=character(), strip.white=T, nlines = 4, skip = 94, sep="\n", quiet =T)

	# Store patient name
	PatientNameTemp <<- scan(fileName, what=character(), strip.white=T, nlines = 1, skip = 87, sep="\n", quiet =T)
	PatientNameTemp <<- substr(PatientNameTemp, 16, nchar(PatientNameTemp))
	tclvalue(PatientName) <<- PatientNameTemp

	# Pull information from file
	Names <<- list(substr(n[1],13,nchar(n[1])), substr(n[2],13,nchar(n[2])),
	               substr(n[3],13,nchar(n[3])), substr(n[4],13,nchar(n[4])),
			   substr(n[4],13,nchar(n[4])))
	offset = 1
	
	D <<- list(a[offset:length(a[,1]),1], a[offset:length(a[,1]),4], 
	           a[offset:length(a[,1]),7], a[offset:length(a[,1]),10])

	MA <<- list(a[offset:length(a[,1]),2], a[offset:length(a[,1]),5],
	            a[offset:length(a[,1]),8], a[offset:length(a[,1]),11])

	StdDv <<- list(a[offset:length(a[,1]),3], a[offset:length(a[,1]),6],
	               a[offset:length(a[,1]),9], a[offset:length(a[,1]),12])

	# Calculate averages area of 4 samples
	for(i in 1:length(D[[1]]))
      	{
		DAve[i] <<- ((D[[1]][i] + D[[2]][i] + D[[3]][i] + D[[4]][i])/4)
		MAAve[i] <<- ((MA[[1]][i] + MA[[2]][i] + MA[[3]][i] + MA[[4]][i])/4)
		StdDvAve[i] <<- ((StdDv[[1]][i] + StdDv[[2]][i] + StdDv[[3]][i] + StdDv[[4]][i])/4)
      	}

	# Append averages onto lists
	D[[5]] <<- DAve
	MA[[5]] <<- MAAve
	StdDv[[5]] <<- StdDvAve

	# Adjust data according to the standard deviation
	MAAdjPos <<- list(MA[[1]] + StdDv[[1]], MA[[2]] + StdDv[[2]],
	                  MA[[3]] + StdDv[[3]], MA[[4]] + StdDv[[4]],
				MA[[5]] + StdDv[[5]])

	MAAdjNeg <<- list(MA[[1]] - StdDv[[1]], MA[[2]] - StdDv[[2]],
	                  MA[[3]] - StdDv[[3]], MA[[4]] - StdDv[[4]],
				MA[[5]] - StdDv[[5]])
}

# Gets the value from the dropdown selector to decide which dataset will be plotted,
# loads it's values to be used later and then estimates teact start and end.
ChooseDataSet <- function()
{
	setNumber <<- as.numeric(tclvalue(tcl(DataSet.comboBox,"getvalue")))+1
	if(setNumber == 0)
		setNumber <<- 1
	D1 <<- D[[setNumber]]
	MA1 <<- MA[[setNumber]]
	StdDv1 <<- StdDv[[setNumber]]
	MA1AdjPos <<- MAAdjPos[[setNumber]]
	MA1AdjNeg <<- MAAdjNeg[[setNumber]]

	# Some information to be displayed in the GUI is adjusted
	tclvalue(CurrentVowelSound) <- Names[[setNumber]]
	tclvalue(selectedStartPosition) <<- -1
	tclvalue(selectedGlottisPosition) <<- D1[which.min(MA1[51:62])+50]
}


##-------------------------------------------------------------------------------------------------------------------##
##---Full Plot---##

# Function running half of the plotting of the entire area data.  Has to call the other for use in tkrplot.
# Splines out all the values just for a nicer plot
PlotFull <- function() 
{
	splineNumFull <- lcm(tubeNumber, length(D1))
	tubeMA1 <<- spline(D1, MA1, n = splineNumFull)
	tubeMA1AdjPos <<- spline(D1, MA1AdjPos, n = splineNumFull)
	tubeMA1AdjNeg <<- spline(D1, MA1AdjNeg, n = splineNumFull)
	if( is.null(imgFull) )
      	{
		imgFull <<- tkrplot(tt,fun=PlotFullFunction,hscale=1.1,vscale=1.1)
		tkgrid(imgFull, row = 1, rowspan = 17, column = 7, columnspan = 1)
      	}
	else
      	{
		tkrreplot(imgFull)
      	}
}

# Function half of of the plotting of the entire area data.  A red line is drawn
# at the selectedStartPosition of the tract and a red at the selectedGlottisPosition.
PlotFullFunction <- function()
{
	plot(tubeMA1AdjNeg$x, tubeMA1AdjNeg$y, type="l", main="Full data plot", xlab="Distance along with the vocal tract (cm)",
	     ylab="Cross-sectional area of the vocal tract (cm^2)", col="blue", ylim=c(0, max(MA1AdjPos)/10+max(MA1AdjPos)))
	par(new=T)
	plot(tubeMA1AdjPos$x, tubeMA1AdjPos$y, type="l", col="black",ylim=c(0, max(MA1AdjPos)/10+max(MA1AdjPos)),xlab="",ylab="")
	par(new=T)
	plot(tubeMA1$x, tubeMA1$y, type="l", col="red",ylim=c(0, max(MA1AdjPos)/10+max(MA1AdjPos)),xlab="",ylab="")
	abline(v = as.numeric(tclvalue(selectedStartPosition)), col = "red")
	abline(v = as.numeric(tclvalue(selectedGlottisPosition)), col = "red")

# Information about the graph grabbed to be used for the interactive start/end selection
	parPlotSize <<- par("plt")
	usrCoords   <<- par("usr")
}


##-------------------------------------------------------------------------------------------------------------------##
##---Selection Choosing---##

# Binds the left mouse click with the full plot image so when clicked StartPosFunction 
# is run.
SelectStartPos <- function()
{
	tkbind(imgFull, "<Button-1>", StartPosFunction)
	tkconfigure(imgFull, cursor="crosshair")
}

# Takes the position of the cursor on the full plot and asks if the user wants
# to move the vocal tract start line there.  Replots after point chosen.
# NOTE - doesn't take into consideration canceling, once user starts they must
#        choose a point before it unbinds.  Could maybe use a 'cancel' button.
StartPosFunction <- function(x,y) 
{
	xClick <- x
	require(tcltk)
	width	<- as.numeric(tclvalue(tkwinfo("reqwidth",imgFull)))

	xMin <- parPlotSize[1] * width
	xMax <- parPlotSize[2] * width

	rangeX <- usrCoords[2] - usrCoords[1]
	imgXcoords <- (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
	xClick <- as.numeric(xClick)+0.5
	xPlotCoord <<- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
	
	msg <- paste("Choose vocal tract start position to be at: \n",
				"x =",format(xPlotCoord,digits=3),"?")
	mbval <- tkmessageBox(title="Mark tract start position at selected position",
					message=msg,type="yesno",icon="question")
	if (tclvalue(mbval)=="yes")
      	{
		tclvalue(selectedStartPosition) <<- round(xPlotCoord, 2)
		tkbind(imgFull, "<Button-1>", doNothing)
		tkconfigure(imgFull, cursor="arrow")
		DoAllPlots()
      	}
}

# Binds the left mouse click with the full plot image so when clicked GlottisPosFunction 
# is run.
SelectGlottisPos <- function()
{
	tkbind(imgFull, "<Button-1>", GlottisPosFunction)
	tkconfigure(imgFull, cursor="crosshair")
}

# Takes the position of the cursor on the full plot and asks if the user wants
# to move the vocal tract end line there.  Replots after point chosen.
# NOTE - Again, canceling not taken into consideration
GlottisPosFunction <- function(x,y) 
{
	xClick <- x
	require(tcltk)
	width	<- as.numeric(tclvalue(tkwinfo("reqwidth",imgFull)))

	xMin <- parPlotSize[1] * width
	xMax <- parPlotSize[2] * width

	rangeX <- usrCoords[2] - usrCoords[1]
	imgXcoords <- (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
	xClick <- as.numeric(xClick)+0.5
	xPlotCoord <<- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
	
	msg <- paste("Choose glottis position to be at: \n",
				"x =",format(xPlotCoord,digits=3),"?")
	mbval <- tkmessageBox(title="Mark glottis position at selected position",
					message=msg,type="yesno",icon="question")
	if (tclvalue(mbval)=="yes")
      	{
		tclvalue(selectedGlottisPosition) <<- round(xPlotCoord, 2)
		tkbind(imgFull, "<Button-1>", doNothing)
		tkconfigure(imgFull, cursor="arrow")
		DoAllPlots()
      	}
}


##-------------------------------------------------------------------------------------------------------------------##
##---Selection Plot---##

# Function running half of the selection plot.  Grabs out the data between the start and end position
# calculates the vocal tract length, splines it out for splitting into tubes later and tkplots.
PlotSelection <- function()
{
	start <- min(which(D1 >= as.numeric(tclvalue(selectedStartPosition))))
	end <- max(which(D1 <= as.numeric(tclvalue(selectedGlottisPosition))))
	D1Selection <<- D1[start:end]
	MA1Selection <<- MA1[start:end]
	MA1SelectionAdjPos <<- MA1AdjPos[start:end]
	MA1SelectionAdjNeg <<- MA1AdjNeg[start:end]	
	vtlength <<- as.numeric(tclvalue(selectedGlottisPosition))-as.numeric(tclvalue(selectedStartPosition))
	
	splineNumSelected <<- lcm(tubeNumber, length(D1Selection))
	tubeMA1Selection <<- spline(D1Selection, MA1Selection, n = splineNumSelected)
	tubeMA1SelectionAdjPos <<- spline(D1Selection, MA1SelectionAdjPos, n = splineNumSelected)
	tubeMA1SelectionAdjNeg <<- spline(D1Selection, MA1SelectionAdjNeg, n = splineNumSelected)
	if( is.null(imgSelection) )
      	{
		imgSelection <<- tkrplot(tt,fun=PlotSelectionFunction,hscale=1.1,vscale=1.1)
		tkgrid(imgSelection, row = 1, rowspan=17, column = 8)
      	}
	else
      	{
		tkrreplot(imgSelection)
      	}
}

# Function half of the selection plot.  Plots out the section of the graph between the 
# start and end points selected.
PlotSelectionFunction <- function()
{
	plot(tubeMA1SelectionAdjNeg$x, tubeMA1SelectionAdjNeg$y, type="l", main="Vocal tract plot", xlab="Distance along with the vocal tract (cm)",
	     ylab="Cross-sectional area of the vocal tract (cm^2)", col="blue", ylim=c(0, 1.1*max(MA1SelectionAdjPos)))
	par(new=T)
	plot(tubeMA1SelectionAdjPos$x, tubeMA1SelectionAdjPos$y, type="l", col="black",ylim=c(0, 1.1*max(MA1SelectionAdjPos)),xlab="",ylab="")
	par(new=T)
	plot(tubeMA1Selection$x, tubeMA1Selection$y, type="l", col="red",ylim=c(0, 1.1*max(MA1SelectionAdjPos)),xlab="",ylab="")
}


##-------------------------------------------------------------------------------------------------------------------##
##---Tubes Plot---##

# Function running half of the tubes plot, splits the selection into the selected number of pieces
# and averages the area in them. These values are then "tcled" so that they may be changed in the GUI.
# tkplot and PlotTubesFunction used to plot them.
PlotTubes <- function(x) 
{
	if (x == "initial")
      	{
		tubeWidth = splineNumSelected / tubeNumber;
		startPt <- 1;
		finishPt <- tubeWidth;
		tubes <<- array(dim = tubeNumber);
		radiusInitial <<- array(dim = tubeNumber);

		for(i in 1:tubeNumber) 
            	{
			tubes[i] <<- mean(tubeMA1Selection$y[startPt : finishPt]);
			startPt <- finishPt;
			finishPt <- finishPt + tubeWidth;
            	}

		# "tcl" tube radii
		tclvalue(radiusTcl1) <- sqrt(tubes[1]/pi)
		tclvalue(radiusTcl2) <- sqrt(tubes[2]/pi)
		tclvalue(radiusTcl3) <- sqrt(tubes[3]/pi)
		tclvalue(radiusTcl4) <- sqrt(tubes[4]/pi)
		tclvalue(radiusTcl5) <- sqrt(tubes[5]/pi)
		tclvalue(radiusTcl6) <- sqrt(tubes[6]/pi)
		tclvalue(radiusTcl7) <- sqrt(tubes[7]/pi)
		tclvalue(radiusTcl8) <- sqrt(tubes[8]/pi)
		tclvalue(radiusTcl9) <- sqrt(tubes[9]/pi)
		tclvalue(radiusTcl10) <- sqrt(tubes[10]/pi)
		tclvalue(radiusTcl11) <- sqrt(tubes[11]/pi)

		# the initial radii are stored for various later uses
		radiusInitial <<- sqrt(tubes/pi)
      	}

	radius <<- array(dim = tubeNumber);
	radius <<- NULL;

	# "tcled" radius values assigned to regular array so that they may be plotted
	radius[1] <<- as.numeric(tclvalue(radiusTcl1));
	radius[2] <<- as.numeric(tclvalue(radiusTcl2));
	radius[3] <<- as.numeric(tclvalue(radiusTcl3));
	radius[4] <<- as.numeric(tclvalue(radiusTcl4));
	radius[5] <<- as.numeric(tclvalue(radiusTcl5));
	radius[6] <<- as.numeric(tclvalue(radiusTcl6));
	radius[7] <<- as.numeric(tclvalue(radiusTcl7));
	radius[8] <<- as.numeric(tclvalue(radiusTcl8));
	radius[9] <<- as.numeric(tclvalue(radiusTcl9));
	radius[10] <<- as.numeric(tclvalue(radiusTcl10));
	radius[11] <<- as.numeric(tclvalue(radiusTcl11));

	if ( is.null(imgTubes) )
      	{
		imgTubes <<- tkrplot(tt,fun=PlotTubesFunction,hscale=1.1,vscale=1.1)
		tkgrid(imgTubes, row = 18, rowspan=13,  column = 7)
      	}
	else
      	{
		tkrreplot(imgTubes)
      	}
}

# Function half of the tubes plot.  The tubes are a model of the vocal tract split into
# a number of tubes, with an average area across their size.  To show this in a plot, a
# bar graph with half above and half below the zero line is used, with the height being
# what the radius would be. if showInitial is true then the initial radii will also be
# plotted
PlotTubesFunction <- function()
{
	barplot(radius, names.arg = c(1:11), space=FALSE, main="Tubes plot", xlab="Tube Number", ylab="Tube radius (cm)",
				     				  col="red", ylim=c(-1.25*max(radiusInitial), 1.25*max(radiusInitial)))
	par(new=T)
	barplot(-radius, space=FALSE, col="red", ylim=c(-1.25*max(radiusInitial), 1.25*max(radiusInitial)))

	# plotting initial radii
	if (showInitial)
		{
		par(new=T)
		barplot(radiusInitial, space=FALSE, col="blue", ylim=c(-1.25*max(radiusInitial), 1.25*max(radiusInitial)), density = 10)
		par(new=T);
		barplot(-radiusInitial, space=FALSE, col="blue", ylim=c(-1.25*max(radiusInitial), 1.25*max(radiusInitial)), density = 10)
		}	
}


##-------------------------------------------------------------------------------------------------------------------##
##---Spectrum Plot---##

# Function running half of spectrum plot, plots the spectrum and binds the graph left mouse button
# and the graph to FindSpectPointFunction, so points can be found when clicking the graph.
PlotSpectrum <- function(x) 
{
	# Varible that determines whether just the initial or the initial and modified spectrum
	# will be plotted. Assigned by method of function call
	whatPlotting <<- x

	spectTubes <<- array(dim = tubeNumber)

	for(i in 1:tubeNumber) 
      	{
		spectTubes[i] <<- (radius[i]^2)*pi
      	}

	if (is.null(imgSpectrum))
      	{
		imgSpectrum <<- tkrplot(tt,fun=PlotSpectrumFunction,hscale=1.1,vscale=1.1)
		tkgrid(imgSpectrum, row = 18, rowspan=13, column = 8)
		tkbind(imgSpectrum, "<Button-1>", FindSpectPointFunction )
		tkconfigure(imgSpectrum, cursor="crosshair")
      	}
	else
      	{
		tkrreplot(imgSpectrum)		
      	}
}

# Function half of the spectrum plot, uses the area2spec method to do the spectrum calculation
# and plotting. Initial and modified spectrum may be plotted
PlotSpectrumFunction <- function() 
{
	# Plot just initial spectrum
	if (whatPlotting == "underlay")
		{
		spectrumData <<- area2spec(tubes, F, vtlength , ylim=c(-2.5,6))
		title("Spectrum Analysis", xlab="Frequency (Hz)", ylab="Magnitude")

		FindSpectrumPeaks()
		}

	# Plot both initial and modified spectrum	
	if (whatPlotting == "overlay")
		{
		spectrumData <<- area2spec(tubes, F, vtlength , ylim=c(-2.5,6))
		title("Spectrum Analysis", xlab="Frequency (Hz)", ylab="Magnitude")
		par(new = T)
		spectrumData <<- area2spec(spectTubes, F, vtlength , ylim=c(-2.5,6), col = "red")

		FindSpectrumPeaks()
		}

# Grab graph data to be used to find points on the graphs by clicking.
	parSpectPlotSize <<- par("plt")
	usrSpectCoords   <<- par("usr")

# Store this datasets values used for R spectogram plotting. Current and modified radii stored
	tubesListInitial[[setNumber]] <<- tubes
	tubesListTcling[[setNumber]] <<- spectTubes
	vtlengthList[[setNumber]] <<- vtlength

}


##-------------------------------------------------------------------------------------------------------------------##
##---Summary Funtions---##

# Calls all 4 plot functions, updates labels, and sets beenModified to FALSE for currently selected set
DoAllPlots <- function()
{
	PlotFull()
	PlotSelection()
	PlotTubes(x = "initial")
	PlotSpectrum(x = "underlay")

	UpdateInfoLabels(x = "initial")

	beenModified[setNumber] <<- FALSE
}

# Plots modified bar graph and spectrum, updates labels, and sets beenModified to TRUE for currently 
# selected set. Steps have been taken to ensure that nothing is plotted if nothing has been modified 
PlotTcling <- function(...)
{
	if (loaded)
		{
		PlotTubes(x = "tcling")

		count <- 0	
		for (i in 1:11)
			{
			if (round(radius[i], 7) == round(radiusInitial[i], 7))
				count <- count + 1	
			}

		if (count == 11) # i.e. if current radius is the same as initial radius for all tubes
			{
			if (!beenModified[setNumber])
				tkmessageBox(title = "Warning", message = "Tubes unmodified from initial radii. Will not continue.")
			else
				{
				tkmessageBox(title = "Warning", message = "Tubes have been changed back to initial radii.")
				DoAllPlots()
				}
			}
		else
			{
			PlotSpectrum(x = "overlay")
			UpdateInfoLabels(x = "tcling")
			beenModified[setNumber] <<- TRUE
			}
		}
}

# Loads and plots new data. Also wipes any speech currently stored and changes "stored to not"
# variables to FALSE
DoEverything <- function()
{
	LoadData()
	ChooseDataSet()
	DoAllPlots()

	loaded <<- TRUE

	for (i in 1:5)
		{
		beenModified[i] <<- F
		initialSpeechSampleImpulse[[i]] <<- NULL
		initialSpeechStoredImpulse[i] <<- F
		initialSpeechSampleGlottus[[i]] <<- NULL
		initialSpeechStoredGlottus[i] <<- F
		}

	tclvalue(initialSpeechStoredTclImpulse) <- "Initial speech\nstored: No"
	tclvalue(initialSpeechStoredTclGlottus) <- "Initial speech\nstored: No"
}

# Function called when a new set is selected. Also checks whether or not speech is stored for that
# dataset, and changes "tcl" variables accordingly 
ChooseAndPlotDataSet <- function()
{
	ChooseDataSet()
	DoAllPlots()

	if (initialSpeechStoredImpulse[setNumber])
		tclvalue(initialSpeechStoredTclImpulse) <- "Initial speech\nstored: Yes"
	else
		tclvalue(initialSpeechStoredTclImpulse) <- "Initial speech\nstored: No"

	if (initialSpeechStoredGlottus[setNumber])
		tclvalue(initialSpeechStoredTclGlottus) <- "Initial speech\nstored: Yes"
	else
		tclvalue(initialSpeechStoredTclGlottus) <- "Initial speech\nstored: No"
}


##-------------------------------------------------------------------------------------------------------------------##
##---Spectrum Point Finder---##

# Finds the first 3 maximium peaks in the spectrum
FindSpectrumPeaks <- function()
{
		FirstFormant <<- which.min((spectrumData $spec[2:512]-spectrumData $spec[1:511])>0)
		nextMin <- FirstFormant + which.min((spectrumData$spec[FirstFormant+2:512]-spectrumData$spec[FirstFormant+1:511])<0)
		SecondFormant <<- nextMin + which.min((spectrumData $spec[nextMin+2:512]-spectrumData $spec[nextMin+1:511])>0)
		nextMin <- SecondFormant + which.min((spectrumData$spec[SecondFormant +2:512]-spectrumData$spec[SecondFormant +1:511])<0)
		ThirdFormant <<- nextMin + which.min((spectrumData $spec[nextMin+2:512]-spectrumData $spec[nextMin+1:511])>0)

		FirstFormant <<- spectrumData$freq[FirstFormant]
		SecondFormant <<- spectrumData$freq[SecondFormant]
		ThirdFormant <<- spectrumData$freq[ThirdFormant]
}

# Finds the position the cursor is at when clicked on the spectrum plot and displays it.
FindSpectPointFunction <- function(x,y)
{
  xClick <- x
  yClick <- y
  width  <- as.numeric(tclvalue(tkwinfo("reqwidth",imgSpectrum)))
  height <- as.numeric(tclvalue(tkwinfo("reqheight",imgSpectrum)))

  xMin <- parSpectPlotSize[1] * width
  xMax <- parSpectPlotSize[2] * width
  yMin <- parSpectPlotSize[3] * height
  yMax <- parSpectPlotSize[4] * height

  rangeX <- usrSpectCoords[2] - usrSpectCoords[1]
  rangeY <- usrSpectCoords[4] - usrSpectCoords[3]

  imgXcoords <- (xCoords-usrSpectCoords[1])*(xMax-xMin)/rangeX + xMin
  imgYcoords <- (yCoords-usrSpectCoords[3])*(yMax-yMin)/rangeY + yMin

  xClick <- as.numeric(xClick)+0.5
  yClick <- as.numeric(yClick)+0.5
  yClick <- height - yClick

  xPlotCoord <- usrSpectCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
  yPlotCoord <- usrSpectCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)

  msg <- paste("Point selected is approximately: \n",
               "x =",format(xPlotCoord,digits=5),",y =",format(yPlotCoord,digits=2))
  mbval<- tkmessageBox(title="Spectrum point",
                       message=msg,type="ok",icon="question")
}


##-------------------------------------------------------------------------------------------------------------------##
##---Upadates The Info Labels--##

# Updates information labels in a different way depending on circumstance.
UpdateInfoLabels <- function(x)
{
	if (x == "initial")
		{	
		tclvalue(FormantValuesInitial) <<- paste("Initial Formant Values:\nF1:\t", round(FirstFormant),
	                                          	    			    " Hz\nF2:\t", round(SecondFormant), 
	                                          		    		    " Hz\nF3:\t", round(ThirdFormant), " Hz")

		tclvalue(FormantValuesTcling) <<- ("Modified Formant Values:\nF1:\nF2:\nF3:")
		}

	if (x == "tcling")
		{
		tclvalue(FormantValuesTcling) <<- paste("Modified Formant Values:\nF1:\t", round(FirstFormant),
	      	                                    	    		    " Hz\nF2:\t", round(SecondFormant), 
	            	                              	    		    " Hz\nF3:\t", round(ThirdFormant), " Hz")
		}

	tclvalue(Tract) <<- paste("Vocal tract length:\n", vtlength, " cm")
	tclvalue(SampFreq) <<- paste("Sampling Frequency:\n", round(spectrumData$fs), " Hz")
}

##-------------------------------------------------------------------------------------------------------------------##
##---Store Formant Data---##

# Adds F1, F2, F3 values, names, HVD, whether the data has been modified or not, and the vocal tract
# starts and end to a one table which contains all previous stored ones too.
StoreFormants <- function(x)
{
	mbval <- tkmessageBox(title="Add Formant Data",
			message="Are you sure you want to add these formant values to the stored values?",type="yesno",icon="question")
	if (tclvalue(mbval)=="yes")
      	{
		Stored$FormantsValues <<- rbind(Stored$FormantsValues, c(FirstFormant, SecondFormant, ThirdFormant))
		Stored$Labels <<- rbind(Stored$Labels, paste(tclvalue(PatientName), "- dataset", setNumber, sep = " "))
		Stored$HVD <<- rbind(Stored$HVD, Names[[setNumber]])
		Stored$Modified <<- rbind(Stored$Modified, beenModified[setNumber])
		Stored$VTStart <<- rbind(Stored$VTStart, as.numeric(tclvalue(selectedStartPosition)))
		Stored$VTEnd <<- rbind(Stored$VTEnd, as.numeric(tclvalue(selectedGlottisPosition)))
		tkmessageBox(title="Add Formant Data", message="Formant values added.")
      	}
}

# Method to clear the stored formant data
ClearStoredFormants <- function()
{
	mbval <- tkmessageBox(title="Clear Formant Data",
				message="Are you sure you want to clear all the stored formant values?" ,type="yesno",icon="question")
	if (tclvalue(mbval)=="yes")
      	{
		Stored$FormantsValues <<- NULL
		Stored$Labels <<- NULL
		Stored$HVD <<- NULL
		Stored$Modified <<- NULL
		Stored$VTStart <<- NULL
		Stored$VTEnd <<- NULL
		tkmessageBox(title="Clear Formant Data", message="Stored formant values cleared.")
      	}
}


##-------------------------------------------------------------------------------------------------------------------##
##---R Plotting Formant Plot---##

# Formant plots
MakeFormantPlots <- function()
{
	if (length(Stored) == 0)
		tkmessageBox(title="Error", message="Must store formant values.")

	par(mfrow=c(1,1))
	eplot(Stored$FormantsValues[,1:2], Stored$HVD, formant = T, dopoints = T)
}


##-------------------------------------------------------------------------------------------------------------------##
##---R Plotting All 4 Spectograms---##

# Spectrum Plots. Also uses area2spec and able to plot initial and modified spectrum
FourSetSpectRPlot <- function()
{
	if (length(tubesListInitial[[2]]) == 0 || length(tubesListInitial[[3]]) == 0 || length(tubesListInitial[[4]]) == 0 )
		tkmessageBox(title="Warning", message="Spectrogram not complete. Look at all datasets to ensure they are plotted.")
		
	par(mfrow=c(2,2))
	for(i in 1:4)
		{
		area2spec(tubesListInitial[[i]],F,vtlengthList[[i]], ylim=c(-2,6.5))
		title(main=Names[[i]])
		if (beenModified[i])
			{
			par(new=T)
			area2spec(tubesListTcling[[i]],F,vtlengthList[[i]], ylim=c(-2,6.5), col = "red")
			}
		}
}


##-------------------------------------------------------------------------------------------------------------------##
##---Misc Functions---##

ToggleInitial <- function()
{
	showInitial <<- xor(showInitial, 1)
	PlotTubes(x = "tcling")
}


##-------------------------------------------------------------------------------------------------------------------##
##---Speech Generation Functions---#

# Creates and impulse train with user inputable duration and frequeny values. Sampling freq is
# taken into account
CreateImpulseTrain <- function()
{
	impulseTrain <<- array(dim = round((as.numeric(tclvalue(impTrainDuration)))*spectrumData$fs))
	temp <- round(round(spectrumData$fs) / (as.numeric(tclvalue(impTrainDesiredFreq))))

	for (i in 1:length(impulseTrain))
		{
		if ((i %% temp) == 0)
			impulseTrain[i] <<- 1
		else
			impulseTrain[i] <<- 0
		}

	# Able to display freq of created train
	tclvalue(impTrainActualFreq) <- round(spectrumData$fs) / temp
}

# Crude method of creating a 150Hz glottal pulse. Sampling freq must be 16050. Should try importing
# from a wav
CreateGlottalPulse <- function()
{
	glottalPulse <<- 0

	for(i in 1:150)
		{
		glottalPulse <<- c(glottalPulse, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
							   0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20,
							   0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40,
							   0.43, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.70,
							   0.73, 0.76, 0.79, 0.82, 0.85, 0.88, 0.91, 0.94, 0.97, 0.98,
							   0.99, 1.00, 0.99, 0.97, 0.90, 0.70, 0.50, 0.40, 0.30, 0.25,
							   0.20, 0.15, 0.10, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
							   0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00)
						   

#		glottalPulse <<- c(glottalPulse, 0.1414, 0.1414, 0.1420, 0.1433, 0.1452, 0.1474, 0.1504, 0.1539, 0.1580, 0.1628,
#  							   0.1712, 0.1810, 0.1921, 0.2043, 0.2176, 0.2329, 0.2562, 0.2817, 0.3104, 0.3510,
# 							   0.3971, 0.4589, 0.5369, 0.6373, 0.7677, 0.9361, 1.1418, 1.3727, 1.6270, 1.9026,
#  							   2.1990, 2.5142, 2.8490, 3.2004, 3.5641, 3.9429, 4.3325, 4.7277, 5.1227, 5.4637,
# 							   5.7554, 5.9995, 6.1937, 6.3148, 6.4251, 6.5022, 6.5215, 6.5263, 6.5146, 6.4332,
#							   6.3430, 6.2476, 6.1427, 6.0217, 5.8519, 5.6689, 5.4884, 5.3070, 5.1182, 4.9248,
# 							   4.7167, 4.5118, 4.3003, 4.0843, 3.8682, 3.6487, 3.4274, 3.2143, 3.0067, 2.8073,
# 							   2.6100, 2.4152, 2.2186, 2.0354, 1.8669, 1.7065, 1.5559, 1.4076, 1.2749, 1.1524,
#							   1.0410, 0.9335, 0.8326, 0.7482, 0.6666, 0.5913, 0.5268, 0.4687, 0.4153, 0.3683,
# 							   0.3284, 0.2923, 0.2631, 0.2394, 0.2178, 0.2037, 0.1916, 0.1806, 0.1709, 0.1623,
# 							   0.1567, 0.1528, 0.1494, 0.1468, 0.1445, 0.1429, 0.1419)
		}
}

# Filters the impulse train or the glottal pulse using a lattice filter method
FilterInput <- function(inputType)
{
	print("Lattice filtering going down!")

	# Chooses input
	if(inputType == "impulse")
		input <- c(impulseTrain)
	if(inputType == "glottus")
		input <- c(glottalPulse)

	# Calculates reflection coefficients. Final coefficient is ignored as this is how it is..
	rcTemp <<- calc.reflection.coef(spectTubes,F)
	rc <<- rcTemp[1:11]
	print("Reflection Coefficients:")
	print(rc)

	order <- length(rc)
	speechLength <- length(input)

	# Creates necessary matrices
	speech  <<- matrix(nrow = 1, ncol = speechLength)
	forward  <- matrix(nrow = speechLength, ncol = order + 1)
	backward <- matrix(nrow = speechLength, ncol = order)
	
	foo <- vector(length = speechLength)

	# Lattice filter algorithm
	n <- 1
	for(i in order:1)
		backward[n,i] <- -(rc[i]*input[n])

	speech[1,n] <<- input[n]

	for (n in 2:speechLength)
		{
		forward[n,order+1] <- input[n]
		for (i in order:2)
			{
			forward[n,i]  <- forward[n,i+1] + rc[i]*backward[n-1,i-1]
			backward[n,i] <- -(rc[i]*forward[n,i]) + backward[n-1,i-1]
			}
		forward[n,1]  <- forward[n,2] + rc[1]*speech[1,n-1]
		backward[n,1] <- -(rc[1]*forward[n,1]) + speech[1,n-1]

		speech[1,n] <<- forward[n,1]
		}

	# Differentiates speech to account for earlier pre emphasis
	for (i in 2:length(speech[1,]))
 		foo[i] <- speech[1,i]-speech[1,i-1]

	speech[1,] <<- foo

	# If an unmodified vocal tract is doing the filtering, the speech created for that particular
	# set is stored and the appropriate index in initialSpeechStored... is set to TRUE. Will do
	# this for either input source depending on what is beind used
	if(inputType == "impulse")
		{
		speechSample <<- as.Sample(speech[1,], round(spectrumData$fs),  16)
		if (!beenModified[setNumber])
			{
			initialSpeechSampleImpulse[[setNumber]] <<- speechSample
			initialSpeechStoredImpulse[[setNumber]] <<- T
			tclvalue(initialSpeechStoredTclImpulse) <- "Initial speech\nstored: Yes"
			}
		}	
	if(inputType == "glottus")
		{
		speechSample <<- as.Sample(speech[1,], 12000, 16)
		if (!beenModified[setNumber])
			{
			initialSpeechSampleGlottus[[setNumber]] <<- speechSample
			initialSpeechStoredGlottus[[setNumber]] <<- T
			tclvalue(initialSpeechStoredTclGlottus) <- "Initial speech\nstored: Yes"
			}
		}

#	plot(input[1:320], typ = "l", ylim = c(min(speech[1,]), max(speech[1,])), title = "Output speech", xlab = "Sample number", ylab = "Amplitude")
#	par(new = T)
#	plot(speech[1,1:320], typ = "l", col = "red", ylim = c(min(speech[1,]), max(speech[1,])))
}

# Pretty self explainitory function..
CreateFilterPlayImpulse <- function()
{
	CreateImpulseTrain()
	FilterInput("impulse")
	play(speechSample)
}

# Pretty self explainitory function..
CreateFilterPlayGlottus <- function()
{
	FilterInput("glottus")
	play(speechSample)
}

# Pretty self explainitory function..
PlayInitialSpeechImpulse <- function()
{
	if (!initialSpeechStoredImpulse[setNumber])
		{
		tkmessageBox(title="Error", message="Speech filtered by initial vocal tract has not been created.")
		return()
		}

	play(initialSpeechSampleImpulse[[setNumber]])
}

# Pretty self explainitory function..
PlayInitialSpeechGlottus <- function()
{
	if (!initialSpeechStoredGlottus[setNumber])
		{
		tkmessageBox(title="Error", message="Speech filtered by initial vocal tract has not been created.")
		return()
		}

	play(initialSpeechSampleGlottus[[setNumber]])
}


##-------------------------------------------------------------------------------------------------------------------##
##---Tcltk Window Creation Etc.---##

# Create the tcltk window.
tt <- tktoplevel()

# Create the GUI items for use in the interface.
Load.but <- tkbutton(tt, text = "Load new file", command = DoEverything, width = 14)

PatientName <- tclVar()
PatientNameEntry.label <- tklabel(tt, text="Patient Name:")
PatientNameEntry.entry <- tkentry(tt, text= tclvalue(PatientName), width = 6)
tkconfigure(PatientNameEntry.entry, textvariable=PatientName)

CurrentVowelSound <- tclVar()
CurrentVowelSound.label <- tklabel(tt,text="Current HVD: ")
CurrentVowelSound.entry <- tklabel(tt,text=tclvalue(CurrentVowelSound))
tkconfigure(CurrentVowelSound.entry, textvariable=CurrentVowelSound)

DataSet.label <- tklabel(tt,text="Dataset:        ")
DataSet <- c("1","2","3","4","Ave")
DataSet.comboBox <- tkwidget(tt,"ComboBox",editable=FALSE, values=DataSet, width=4, modifycmd = ChooseAndPlotDataSet)

StartSel.but <- tkbutton(tt, text = "Choose tract start\nby hand", command = SelectStartPos, width = 14)
GlottisSel.but <- tkbutton(tt, text = "Choose tract end\nby hand", command = SelectGlottisPos, width = 14)

selectedStartPosition <- tclVar()
StartEntry.label <- tklabel(tt, text="Tract start:    ")
StartEntry.entry <- tkentry(tt, text= tclvalue(selectedStartPosition), width = 14)
tkconfigure(StartEntry.entry, textvariable=selectedStartPosition)

selectedGlottisPosition <- tclVar()
EndEntry.label <- tklabel(tt,text="Tract end:     ")
EndEntry.entry <- tkentry(tt, textvariable = tclvalue(selectedGlottisPosition), width = 14)
tkconfigure(EndEntry.entry, textvariable = selectedGlottisPosition)
SetStartAndEnd.but <-tkbutton(tt, text = "Set entered\ntract values", command = DoAllPlots, width = 14)

Tract <- tclVar("Vocal tract length:\n")
Tract.label <- tklabel(tt, text=tclvalue(Tract), justify="left")
tkconfigure(Tract.label, textvariable=Tract)

SampFreq <- tclVar("Sampling Frequency:\n")
SampFreq.label <- tklabel(tt, text=tclvalue(SampFreq), justify="left")
tkconfigure(SampFreq.label, textvariable=SampFreq)

sliderMax <- 3

radiusTcl1 <- tclVar()
Tube1Entry.label <- tklabel(tt, text="Tube 1 radius: ")
Tube1Entry.slider <- tkscale(tt, from = 0, to = sliderMax, showvalue = F, variable = radiusTcl1, resolution = 0.001, orient = "horiz", width = 14, command = PlotTcling)
Tube1Entry.entry <- tkentry(tt, text= tclvalue(radiusTcl1), width = 14)
tkconfigure(Tube1Entry.entry, textvariable=radiusTcl1)

radiusTcl2 <- tclVar()
Tube2Entry.label <- tklabel(tt, text="Tube 2 radius: ")
Tube2Entry.slider <- tkscale(tt, from = 0, to = sliderMax, showvalue = F, variable = radiusTcl2, resolution = 0.001, orient = "horiz", width = 14, command = PlotTcling)
Tube2Entry.entry <- tkentry(tt, text= tclvalue(radiusTcl2), width = 14)
tkconfigure(Tube2Entry.entry, textvariable=radiusTcl2)

radiusTcl3 <- tclVar()
Tube3Entry.label <- tklabel(tt, text="Tube 3 radius: ")
Tube3Entry.slider <- tkscale(tt, from = 0, to = sliderMax, showvalue = F, variable = radiusTcl3, resolution = 0.001, orient = "horiz", width = 14, command = PlotTcling)
Tube3Entry.entry <- tkentry(tt, text= tclvalue(radiusTcl3), width = 14)
tkconfigure(Tube3Entry.entry, textvariable=radiusTcl3)

radiusTcl4 <- tclVar()
Tube4Entry.label <- tklabel(tt, text="Tube 4 radius: ")
Tube4Entry.slider <- tkscale(tt, from = 0, to = sliderMax, showvalue = F, variable = radiusTcl4, resolution = 0.001, orient = "horiz", width = 14, command = PlotTcling)
Tube4Entry.entry <- tkentry(tt, text= tclvalue(radiusTcl4), width = 14)
tkconfigure(Tube4Entry.entry, textvariable=radiusTcl4)

radiusTcl5 <- tclVar()
Tube5Entry.label <- tklabel(tt, text="Tube 5 radius: ")
Tube5Entry.slider <- tkscale(tt, from = 0, to = sliderMax, showvalue = F, variable = radiusTcl5, resolution = 0.001, orient = "horiz", width = 14, command = PlotTcling)
Tube5Entry.entry <- tkentry(tt, text= tclvalue(radiusTcl5), width = 14)
tkconfigure(Tube5Entry.entry, textvariable=radiusTcl5)

radiusTcl6 <- tclVar()
Tube6Entry.label <- tklabel(tt, text="Tube 6 radius: ")
Tube6Entry.slider <- tkscale(tt, from = 0, to = sliderMax, showvalue = F, variable = radiusTcl6, resolution = 0.001, orient = "horiz", width = 14, command = PlotTcling)
Tube6Entry.entry <- tkentry(tt, text= tclvalue(radiusTcl6), width = 14)
tkconfigure(Tube6Entry.entry, textvariable=radiusTcl6)

radiusTcl7 <- tclVar()
Tube7Entry.label <- tklabel(tt, text="Tube 7 radius: ")
Tube7Entry.slider <- tkscale(tt, from = 0, to = sliderMax, showvalue = F, variable = radiusTcl7, resolution = 0.001, orient = "horiz", width = 14, command = PlotTcling)
Tube7Entry.entry <- tkentry(tt, text= tclvalue(radiusTcl7), width = 14)
tkconfigure(Tube7Entry.entry, textvariable=radiusTcl7)

radiusTcl8 <- tclVar()
Tube8Entry.label <- tklabel(tt, text="Tube 8 radius: ")
Tube8Entry.slider <- tkscale(tt, from = 0, to = sliderMax, showvalue = F, variable = radiusTcl8, resolution = 0.001, orient = "horiz", width = 14, command = PlotTcling)
Tube8Entry.entry <- tkentry(tt, text= tclvalue(radiusTcl8), width = 14)
tkconfigure(Tube8Entry.entry, textvariable=radiusTcl8)

radiusTcl9 <- tclVar()
Tube9Entry.label <- tklabel(tt, text="Tube 9 radius: ")
Tube9Entry.slider <- tkscale(tt, from = 0, to = sliderMax, showvalue = F, variable = radiusTcl9, resolution = 0.001, orient = "horiz", width = 14, command = PlotTcling)
Tube9Entry.entry <- tkentry(tt, text= tclvalue(radiusTcl9), width = 14)
tkconfigure(Tube9Entry.entry, textvariable=radiusTcl9)

radiusTcl10 <- tclVar()
Tube10Entry.label <- tklabel(tt, text="Tube 10 radius: ")
Tube10Entry.slider <- tkscale(tt, from = 0, to = sliderMax, showvalue = F, variable = radiusTcl10, resolution = 0.001, orient = "horiz", width = 14, command = PlotTcling)
Tube10Entry.entry <- tkentry(tt, text= tclvalue(radiusTcl10), width = 14)
tkconfigure(Tube10Entry.entry, textvariable=radiusTcl10)

radiusTcl11 <- tclVar()
Tube11Entry.label <- tklabel(tt, text="Tube 11 radius: ")
Tube11Entry.slider <- tkscale(tt, from = 0, to = sliderMax, showvalue = F, variable = radiusTcl11, resolution = 0.001, orient = "horiz", width = 14, command = PlotTcling)
Tube11Entry.entry <- tkentry(tt, text= tclvalue(radiusTcl11), width = 14)
tkconfigure(Tube11Entry.entry, textvariable=radiusTcl11)

SetTubeRadii.but <-tkbutton(tt, text = "Set entered\ntube radii", command = PlotTcling, width = 14)
RestoreTubeRadii.but <-tkbutton(tt, text = "Restore initial\ntube radii", command = DoAllPlots, width = 14)
ToggleInitial.but <-tkbutton(tt, text = "Toggle showing\nof initial tubes", command = ToggleInitial, width = 14)

FormantValuesInitial <- tclVar("Initial Formant Values:\nF1:\nF2:\nF3:")
FormantValuesInitial.label <- tklabel(tt, text=tclvalue(FormantValuesInitial), justify="left")
tkconfigure(FormantValuesInitial.label, textvariable=FormantValuesInitial)

FormantValuesTcling <- tclVar("Modified Formant Values:\nF1:\nF2:\nF3:")
FormantValuesTcling.label <- tklabel(tt, text=tclvalue(FormantValuesTcling), justify="left")
tkconfigure(FormantValuesTcling.label, textvariable=FormantValuesTcling)

StoreFormants.but <- tkbutton(tt,text = "Store current\nformant values", command = StoreFormants, width = 14)
ClearStoredFormants.but <- tkbutton(tt, text = "Clear stored\nformant values", command = ClearStoredFormants, width = 14)
MakeEmuFormantPlot.but <- tkbutton(tt, text = "Plot formants", command = MakeFormantPlots, width = 14)
FourRSpect.but <- tkbutton(tt, text = "Plot spectograms", command = FourSetSpectRPlot, width = 14)

impTrainHeader <- tclVar("Input: IMPULSE TRAIN")
impTrainHeader.label <- tklabel(tt, text=tclvalue(impTrainHeader))

impTrainDuration <- tclVar(1)
impTrainDuration.label <- tklabel(tt, text="Time duration of\nimpulse train (s): ")
impTrainDuration.entry <- tkentry(tt, text= tclvalue(impTrainDuration), width = 14)
tkconfigure(impTrainDuration.entry, textvariable=impTrainDuration)

impTrainDesiredFreq <- tclVar(150)
impTrainDesiredFreq.label <- tklabel(tt, text="Desired freq of\nimpulse train (Hz): ")
impTrainDesiredFreq.entry <- tkentry(tt, text= tclvalue(impTrainDesiredFreq), width = 14)
tkconfigure(impTrainDesiredFreq.entry, textvariable=impTrainDesiredFreq)

impTrainActualFreq <- tclVar()
impTrainActualFreq.label <- tklabel(tt, text="Actual freq of\nimpulse train (Hz): ")
impTrainActualFreq.entry <- tkentry(tt, text= tclvalue(impTrainActualFreq), width = 14)
tkconfigure(impTrainActualFreq.entry, textvariable=impTrainActualFreq)

impTrainCreateFilterPlay.but <-tkbutton(tt, text = "Filter then play\nwith current radii", command = CreateFilterPlayImpulse, width = 14)

initialSpeechStoredTclImpulse <- tclVar("Initial speech\nstored: No")
initialSpeechStoredImpulse.label <- tklabel(tt, text = tclvalue(initialSpeechStoredTclImpulse))
tkconfigure(initialSpeechStoredImpulse.label, textvariable=initialSpeechStoredTclImpulse)

playInitialSpeechImpulse.but <-tkbutton(tt, text = "Play initial\nspeech", command = PlayInitialSpeechImpulse, width = 14)

pulseHeader <- tclVar("Input: GLOTTAL PULSE")
pulseHeader.label <- tklabel(tt, text=tclvalue(pulseHeader))

pulseDurationOne.label <- tklabel(tt, text="Time duration of\nglottal pulse (s): ")
pulseDurationTwo.label <- tklabel(tt, text="1                          ")

pulseFreqOne.label <- tklabel(tt, text="Actual freq of\nglottal pulse (Hz): ")
pulseFreqTwo.label <- tklabel(tt, text="150                      ")

pulseCreateFilterPlay.but <-tkbutton(tt, text = "Filter then play\nwith current radii", command = CreateFilterPlayGlottus, width = 14)

initialSpeechStoredTclGlottus <- tclVar("Initial speech\nstored: No")
initialSpeechStoredGlottus.label <- tklabel(tt, text = tclvalue(initialSpeechStoredTclGlottus))
tkconfigure(initialSpeechStoredGlottus.label, textvariable=initialSpeechStoredTclGlottus)

playInitialSpeechGlottus.but <-tkbutton(tt, text = "Play initial\nspeech", command = PlayInitialSpeechGlottus, width = 14)

ColumnSkip <- tklabel(tt, text = "")

# Place the starting GUI items on the window and focus to it.
tkgrid(Load.but, columnspan = 6)
tkgrid(PatientNameEntry.label, PatientNameEntry.entry, columnspan = 3)
tkgrid(CurrentVowelSound.label, CurrentVowelSound.entry, columnspan = 3)
tkgrid(DataSet.label, DataSet.comboBox, columnspan = 3)
tkgrid(StartSel.but, GlottisSel.but, columnspan = 3)
tkgrid(StartEntry.label, StartEntry.entry, columnspan = 2)
tkgrid(EndEntry.label, EndEntry.entry, columnspan = 2)
tkgrid(SetStartAndEnd.but, column = 4, columnspan = 2, row = 5, rowspan = 2)
tkgrid(Tract.label, SampFreq.label, columnspan = 3)
tkgrid(Tube1Entry.label, Tube1Entry.slider, Tube1Entry.entry, columnspan = 2)
tkgrid(Tube2Entry.label, Tube2Entry.slider, Tube2Entry.entry, columnspan = 2)
tkgrid(Tube3Entry.label, Tube3Entry.slider, Tube3Entry.entry, columnspan = 2)
tkgrid(Tube4Entry.label, Tube4Entry.slider, Tube4Entry.entry, columnspan = 2)
tkgrid(Tube5Entry.label, Tube5Entry.slider, Tube5Entry.entry, columnspan = 2)
tkgrid(Tube6Entry.label, Tube6Entry.slider, Tube6Entry.entry, columnspan = 2)
tkgrid(Tube7Entry.label, Tube7Entry.slider, Tube7Entry.entry, columnspan = 2)
tkgrid(Tube8Entry.label, Tube8Entry.slider, Tube8Entry.entry, columnspan = 2)
tkgrid(Tube9Entry.label, Tube9Entry.slider, Tube9Entry.entry, columnspan = 2)
tkgrid(Tube10Entry.label, Tube10Entry.slider, Tube10Entry.entry, columnspan = 2)
tkgrid(Tube11Entry.label, Tube11Entry.slider, Tube11Entry.entry, columnspan = 2)
tkgrid(ToggleInitial.but, RestoreTubeRadii.but, SetTubeRadii.but, columnspan = 2)
tkgrid(FormantValuesInitial.label, FormantValuesTcling.label, columnspan = 3)
tkgrid(StoreFormants.but, ClearStoredFormants.but, columnspan = 3)
tkgrid(MakeEmuFormantPlot.but, FourRSpect.but, columnspan = 3)
tkgrid(impTrainHeader.label, columnspan = 6)
tkgrid(impTrainDuration.label, impTrainDuration.entry, impTrainCreateFilterPlay.but, columnspan = 2)
tkgrid(impTrainDesiredFreq.label, impTrainDesiredFreq.entry, initialSpeechStoredImpulse.label, columnspan = 2)
tkgrid(impTrainActualFreq.label, impTrainActualFreq.entry, playInitialSpeechImpulse.but, columnspan = 2)
tkgrid(pulseHeader.label, columnspan = 6)
tkgrid(pulseDurationOne.label, pulseDurationTwo.label, pulseCreateFilterPlay.but, columnspan = 2)
tkgrid(pulseFreqOne.label, pulseFreqTwo.label, initialSpeechStoredGlottus.label, columnspan = 2)
tkgrid(ColumnSkip, ColumnSkip, playInitialSpeechGlottus.but, columnspan = 2)

DoEverything()
CreateImpulseTrain()
CreateGlottalPulse()
tkfocus(tt)
