library(ChainLadder)

## reading in the triangles from the CSV file
all.data <- read.csv("F:\Data\\RunOff Data.csv", header=F)

Tri.arr <- array(0, dim=c(10,10,50)) 

for(i in 1:50){
	Tri.arr[,,i] <- as.matrix(all.data[((i*11)-10):((i*11)-1), 1:10])
}

## formating each of the triangles as seperate to read into the formula
triangle_1 <- as.triangle(Tri.arr[,,1])
triangle_2 <- as.triangle(Tri.arr[,,2])
triangle_3 <- as.triangle(Tri.arr[,,3])
triangle_4 <- as.triangle(Tri.arr[,,4])
triangle_5 <- as.triangle(Tri.arr[,,5])
triangle_6 <- as.triangle(Tri.arr[,,6])
triangle_7 <- as.triangle(Tri.arr[,,7])
triangle_8 <- as.triangle(Tri.arr[,,8])
triangle_9 <- as.triangle(Tri.arr[,,9])
triangle_10 <- as.triangle(Tri.arr[,,10])
triangle_11 <- as.triangle(Tri.arr[,,11])
triangle_12 <- as.triangle(Tri.arr[,,12])
triangle_13 <- as.triangle(Tri.arr[,,13])
triangle_14 <- as.triangle(Tri.arr[,,14])
triangle_15 <- as.triangle(Tri.arr[,,15])
triangle_16 <- as.triangle(Tri.arr[,,16])
triangle_17 <- as.triangle(Tri.arr[,,17])
triangle_18 <- as.triangle(Tri.arr[,,18])
triangle_19 <- as.triangle(Tri.arr[,,19])
triangle_20 <- as.triangle(Tri.arr[,,20])
triangle_21 <- as.triangle(Tri.arr[,,21])
triangle_22 <- as.triangle(Tri.arr[,,22])
triangle_23 <- as.triangle(Tri.arr[,,23])
triangle_24 <- as.triangle(Tri.arr[,,24])
triangle_25 <- as.triangle(Tri.arr[,,25])
triangle_26 <- as.triangle(Tri.arr[,,26])
triangle_27 <- as.triangle(Tri.arr[,,27])
triangle_28 <- as.triangle(Tri.arr[,,28])
triangle_29 <- as.triangle(Tri.arr[,,29])
triangle_30 <- as.triangle(Tri.arr[,,30])
triangle_31 <- as.triangle(Tri.arr[,,31])
triangle_32 <- as.triangle(Tri.arr[,,32])
triangle_33 <- as.triangle(Tri.arr[,,33])
triangle_34 <- as.triangle(Tri.arr[,,34])
triangle_35 <- as.triangle(Tri.arr[,,35])
triangle_36 <- as.triangle(Tri.arr[,,36])
triangle_37 <- as.triangle(Tri.arr[,,37])
triangle_38 <- as.triangle(Tri.arr[,,38])
triangle_39 <- as.triangle(Tri.arr[,,39])
triangle_40 <- as.triangle(Tri.arr[,,40])
triangle_41 <- as.triangle(Tri.arr[,,41])
triangle_42 <- as.triangle(Tri.arr[,,42])
triangle_43 <- as.triangle(Tri.arr[,,43])
triangle_44 <- as.triangle(Tri.arr[,,44])
triangle_45 <- as.triangle(Tri.arr[,,45])
triangle_46 <- as.triangle(Tri.arr[,,46])
triangle_47 <- as.triangle(Tri.arr[,,47])
triangle_48 <- as.triangle(Tri.arr[,,48])
triangle_49 <- as.triangle(Tri.arr[,,49])
triangle_50 <- as.triangle(Tri.arr[,,50])

# ------------------------------------------------------------------------ #

## Bootstrapping (Attempt!) 
## Over dispersed Poisson 

## Creating a vector for the reserve estimates for one triangle 
boot.output <- BootChainLadder(triangle_2, 1000, process.distr="od.pois") 
summary(boot.output) 

boot.total.reserves <- numeric(11) 

for(i in 1:10){ 
	## Bootstrapped Reserve for year i 
	boot.total.reserves[i] <- summary(boot.output)$ByOrigin[i,2] 

	## Bootstrapped Total Reserve 
	boot.total.reserves[11] <-summary(boot.output)$Totals[2,] 
}

boot.total.reserves

## Creating a matrix with the Reserve estimates for each of the triangles
## Rows correspond to each triangle
## Columns 1 - 10 are for each year and 11 is for the total reserves

boot.total.reserves.mat <- matrix(0,50,11) 

for(j in 1:50){ 
	boot.output <- BootChainLadder(Tri.arr[,,j], 1000, process.distr="od.pois") 
	boot.total.reserves <- numeric(11) 
	
	for(i in 1:10){ 
		## Bootstrapped Reserve for year i 
		boot.total.reserves[i] <- summary(boot.output)$ByOrigin[i,2] 
		## Bootstrapped Total Reserve 
		boot.total.reserves[11] <-summary(boot.output)$Totals[2,] 
	} 
	boot.total.reserves.mat[j,] <- boot.total.reserves 
} 

boot.total.reserves.mat[2,]

# ------------------------------------------------------ #

## Scaling the triangles 
Tri.arr.scale <- array(0, dim=c(10,10,50)) 

for(i in 1:50){
	if(i > 0) Tri.arr.scale[,,i] <- Tri.arr[,,i]*1000
	if(i > 15) Tri.arr.scale[,,i] <- Tri.arr[,,i]*5000
	if(i > 30) Tri.arr.scale[,,i] <- Tri.arr[,,i]*10000
}

## Check
Tri.arr.scale[,,49]/10000

Tri.arr[,,49]

Tri.arr.scale[,,2]

# ------------------------------------------------------------------- #

Tri.arr[,,43]

boot.output <- BootChainLadder(Tri.arr.scale[,,1], 1000, process.distr="od.pois") 
summary(boot.output) 

boot.output <- BootChainLadder(triangle_50, 1000, process.distr="od.pois") 
summary(boot.output) 

names(boot.output)

hist(boot.output$IBNR.Totals)

x11()

# Total Reserves
hist(boot.output$IBNR.Totals)

boot.output$IBNR.ByOrigin[,,2]
boot.output$IBNR.Triangles[,,2]

boot.reserves <- array(0, dim=c(10,1000,50))

for( j in 1:50){
	boot.output <- BootChainLadder(Tri.arr.scale[,,j], 1000, process.distr="od.pois")
	for( i in 1:1000){
		boot.reserves[,i,j] <- boot.output$IBNR.ByOrigin[,,i]
	}
}

boot.output$ChainLadder.Residuals

boot.reserves[,,1]

hist(boot.reserves[2,,1])

boot.reserves[3,100,]
length(subset(boot.reserves[,,1], boot.reserves[,,1]<0))

mean(boot.reserves[10,,43])

# --------------------------------------------------------------------------------- #
## Changing to a gamma process distribution to see about the negative increments

boot.reserves.gamma <- array(0, dim=c(10,100,50))

for( j in 1:50){
	boot.output <- BootChainLadder(Tri.arr.scale[,,j], 100, process.distr="gamma")
	for( i in 1:100){
		boot.reserves.gamma[,i,j] <- boot.output$IBNR.ByOrigin[,,i]
	}
}

hist(boot.reserves.gamma[10,,43])

boot.reserves.gamma[,,]

boot.output.gamma <- BootChainLadder(Tri.arr.scale[,,43], 1000, process.distr="gamma") 
summary(boot.output.gamma) 

# ----------------------------------------------------------------------------#

boot.output.TA <- BootChainLadder(GenIns, 1000, process.distr="od.pois") 
summary(boot.output.TA)

names(boot.output.TA)

hist(boot.output.TA$IBNR.Totals)

boot.output.TA$ChainLadder.Residuals

boot.output.RD <- ?BootChainLadder(Tri.arr.scale[,,43], 1000, process.distr="od.pois") 
summary(boot.output.RD)


boot.output.RD$ChainLadder.Residuals

boot.output.TA$f
boot.output.RD$f

boot.output.TA$NYCost.ByOrigin
boot.output.RD$f


# ------------------------------------------------------------------------ #

## Comparing the amount of negative values in each of the acc years

negatives <- matrix(0,50,11)

for( i in 1:50){
	for ( j in 1:10 ){
	negatives[i,j] <- length(subset(boot.reserves[j,,i], boot.reserves[j,,i]<0))
	}
	negatives[i,11] <- sum(negatives[i,])
}

length(subset(boot.reserves[3,,50], boot.reserves[3,,50]<0))


negatives.TA <- matrix(0,1,11)

for( i in 1:10 ){
	negatives.TA[1,i] <- length(subset(boot.output.TA$IBNR.ByOrigin[i,1,], boot.output.TA$IBNR.ByOrigin[i,1,]<0))
}	
negatives.TA[1,11] <- sum(negatives.TA[1,])


boot.output.TA$IBNR.ByOrigin

length(subset(boot.output.TA$IBNR.ByOrgin[2,1,], boot.output.TA$IBNR.ByOrgin[2,1,]<0))


# ------------------------------------------------------------------------ #

## 


glmfit <- glmReserve(triangle_43,  mse.method = c("bootstrap"), nsim = 1000)
summary(glmfit)

glmfit <- glmReserve(GenIns,  mse.method = c("bootstrap"), nsim = 1000)
summary(glmfit)


# ------------------------------------------------------------------------ #

## 






