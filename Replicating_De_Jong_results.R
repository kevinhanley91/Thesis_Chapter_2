library(ChainLadder)
library(MASS)

## Entering Data from De Jong


a<-rep(0,10);b<-rep(0,10);osr<-rep(0,1000)
tab <- array(dim=c(10,10))
tab[1,1]=5012
tab[1,2]=8269
tab[1,3]=10907
tab[1,4]=11805
tab[1,5]=13539
tab[1,6]=16181
tab[1,7]=18009
tab[1,8]=18608
tab[1,9]=18662
tab[1,10]=18834
tab[2,1]=106
tab[2,2]=4285
tab[2,3]=5396
tab[2,4]=10666
tab[2,5]=13782
tab[2,6]=15599
tab[2,7]=15496
tab[2,8]=16169
tab[2,9]=16704
tab[3,1]=3410
tab[3,2]=8992
tab[3,3]=13873
tab[3,4]=16141
tab[3,5]=18735
tab[3,6]=22214
tab[3,7]=22863
tab[3,8]=23466
tab[4,1]=5655
tab[4,2]=11555
tab[4,3]=15766
tab[4,4]=21266
tab[4,5]=23425
tab[4,6]=26083
tab[4,7]=27067
tab[5,1]=1092
tab[5,2]=9565
tab[5,3]=15836
tab[5,4]=22169
tab[5,5]=25955
tab[5,6]=26180
tab[6,1]=1513
tab[6,2]=6445
tab[6,3]=11702
tab[6,4]=12935
tab[6,5]=15852
tab[7,1]=557
tab[7,2]=4020
tab[7,3]=10946
tab[7,4]=12314
tab[8,1]=1351
tab[8,2]=6947
tab[8,3]=13112
tab[9,1]=3133
tab[9,2]=5395
tab[10,1]=2063

Y=as.vector(tab)
D=rep(1:10,each=10)
A=rep(1:10,10)
base=data.frame(Y,D,A)

triangle <- as.triangle(base, origin="A", dev="D", value="Y")


## Fitting the Mack Model

Mack <- MackChainLadder(tab, est.sigma="Mack")
Mack


## Fitting the OD Poisson GLM

PoisGLM <- glmReserve(triangle)

# this method doesn't work for this data because there is a negative increment

summary(PoisGLM)


PoisGLM1 <- glm(Y ~ as.factor(D) + as.factor(A), data=base, family=poisson(link="log"))

PoisGLM$coef

predict.value <- predict(PoisGLM, typ="response", newdata = base)



