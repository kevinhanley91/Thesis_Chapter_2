library(ChainLadder)

## reading in the triangles from the CSV file
all.data <- read.csv("F:Data\\RunOff Data.csv", header=F)

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




fit43 <- glmReserve(triangle_43, mse.method="bootstrap", nsim=1000)
fit43$FullTriangle

names(fit43)


## Looking at glmReserve for our triangles: 
## Getting the GLM results for all triangles 

#glm.results.all <- list() 
#for(i in 1:50){ 
# 	glm.results.all[[i]] <- glmReserve(as.triangle(Tri.arr[,,i]), cum=T, mse.method="bootstrap", nsim=1000) 
#} 

## Looking at getting S.D / reserve values for all triangles along with summaries as above 

p.e <- matrix(0, 50, 10) 
claim.sim.m <- array(0,dim=c(100,9,50)) 

glm.results.all.1 <- list() 

for(i in 1:50){ 
	glm.results.all.1[[i]] <- glmReserve(as.triangle(Tri.arr[,,i]), cum=T, mse.method="bootstrap", nsim=100) 
	p.e[i,] <- glm.results.all.1[[i]]$summary[,5]/glm.results.all.1[[i]]$summary[,4] 
	claim.sim.m[,,i] <- glm.results.all.1[[i]][[8]] 
} 

## Naming the output of our loops 
colnames(p.e) <- c("Year 2", "Year 3", "Year 4", "Year 5", "Year 6", "Year 7", "Year 8", "Year 9", "Year 10", "Totals") 
colnames(claim.sim.m) <- c("Year 2", "Year 3", "Year 4", "Year 5", "Year 6", "Year 7", "Year 8", "Year 9", "Year 10") 

## testing (the above) 
#glm.results.all.1[[1]] 
#p.e[1,] 
## Looking at histograms 
## Getting the simulated total reserves for each triangle (i.e. 100 values for each triangle): 

## Individual 
total.sim.res <- numeric(100) 
for(i in 1:100){ 
	total.sim.res[i] <- sum(claim.sim.m[i,,2]) 
} 

## Now for all triangles 
total.sim.res.m <- matrix(0, 50, 100) 
for(k in 1:50){ 
	for(i in 1:100){ 
		total.sim.res.m[k,i] <- sum(claim.sim.m[i,,k]) 
	} 
} 

hist(total.sim.res.m[9,], breaks=10) 

## Test 
#total.sim.res.m[2,]

















