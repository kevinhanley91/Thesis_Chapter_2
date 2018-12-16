library(ChainLadder)

## reading in the triangles from the CSV file
all.data <- read.csv("F:\Data\\RunOff Data.csv", header=F)

Tri.arr <- array(0, dim=c(10,10,50)) 

for(i in 1:50){
	Tri.arr[,,i] <- as.matrix(all.data[((i*11)-10):((i*11)-1), 1:10])
}



## Looking at Bootstrapping the Mack Model: 

tab<-array(dim=c(10,10));t<-array(dim=c(10,10));res<-rep(0,1000);rs2<-rep(0,1000);rs3<-rep(0,1000);rs4<-rep(0,1000);rs5<-rep(0,1000);rs6<-rep(0,1000);rs7<-rep(0,1000);rs8<-rep(0,1000);rs9<-rep(0,1000);rs10<-rep(0,1000)
mack.boot.ta.res.arr <-array(0, dim=c(10, 1000, 50))
for(w in 1:50){
	tab <- cum2incr(Tri.arr[,,w])
	Y=as.vector(tab)
	D=rep(1:10,each=10)
	A=rep(2001:2010,10)
	base=data.frame(Y,D,A)

 	for (j in 1:10){
		for (i in 1:11-j){
			if (i==1) {t[j,1]<-tab[j,1]} 
			else {t[j,i]<-t[j,i-1]+tab[j,i]
		}
	}
}

f<-rep(0,10)
f[1]=(t[1,2]+t[2,2]+t[3,2]+t[4,2]+t[5,2]+t[6,2]+t[7,2]+t[8,2]+t[9,2]) / (t[1,1]+t[2,1]+t[3,1]+t[4,1]+t[5,1]+t[6,1]+t[7,1]+t[8,1]+t[9,1])
f[2]=(t[1,3]+t[2,3]+t[3,3]+t[4,3]+t[5,3]+t[6,3]+t[7,3]+t[8,3]) / (t[1,2]+t[2,2]+t[3,2]+t[4,2]+t[5,2]+t[6,2]+t[7,2]+t[8,2])
f[3]=(t[1,4]+t[2,4]+t[3,4]+t[4,4]+t[5,4]+t[6,4]+t[7,4]) / (t[1,3]+t[2,3]+t[3,3]+t[4,3]+t[5,3]+t[6,3]+t[7,3])
f[4]=(t[1,5]+t[2,5]+t[3,5]+t[4,5]+t[5,5]+t[6,5]) / (t[1,4]+t[2,4]+t[3,4]+t[4,4]+t[5,4]+t[6,4])
f[5]=(t[1,6]+t[2,6]+t[3,6]+t[4,6]+t[5,6]) / (t[1,5]+t[2,5]+t[3,5]+t[4,5]+t[5,5])
f[6]=(t[1,7]+t[2,7]+t[3,7]+t[4,7]) / (t[1,6]+t[2,6]+t[3,6]+t[4,6])
f[7]=(t[1,8]+t[2,8]+t[3,8]) / (t[1,7]+t[2,7]+t[3,7])
f[8]=(t[1,9]+t[2,9]) / (t[1,8]+t[2,8])
f[9]=(t[1,10] / t[1,9])
f[10]=1


tfc<-array(dim=c(10,10))
tfi<-array(dim=c(10,10))
tot=0;
for (i in 2:10){
	for (j in (12-i):10){
		if ((i+j)==12) {tfc[i,j]=t[i,j-1]*f[j-1]} 
			else {tfc[i,j]=tfc[i,j-1]*f[j-1]}
		if ((i+j)==12) {tfi[i,j]=tfc[i,j]-t[i,j-1]} 
			else {tfi[i,j]=tfc[i,j]-tfc[i,j-1]}
		tot=tot+tfi[i,j]
	}
}

aa<-rep(0,9)
for (j in 1:9){
	a=0;
	for (i in 1:(10-j)){
		a=a+(sqrt(t[i,j])*(t[i,j+1]/t[i,j]-f[j]))^2
	}
	aa[j]=a/(10-j-1)
}

aa[9]=446.61655

q=0
r<-rep(0,44)
for (j in 1:9){
	for (i in 1:(10-j)){
	q<-q+1
	r[q]=sqrt((10-j)/(9-j))*(sqrt(t[i,j])*(t[i,j+1]/t[i,j]-f[j]))/sqrt(aa[j])
	}
}

r<-r[!is.na(r)]

for (k in 1:1000){
	rr=sample(r,replace=T,size=45)-mean(r)
	a=0
	s<-array(dim=c(10,10))
	for (i in 1:9){
		for (j in 2:(11-i)){
			a<-a+1;s[i,j]=rr[a]
		}
	}
#a=0;s<-array(dim=c(10,10))
#for (i in 1:9){
#for (j in 2:(11-i)){
#a<-a+1;s[i,j]=r[a]-mean(r)
#}
#}

#s[9,2]=r[44]-mean(r)

ss<-array(dim=c(10,10))
for (i in 1:9){
	for (j in 2:(11-i)){
		ss[i,j]=s[i,j]*sqrt(aa[j-1]/t[i,j-1])+f[j-1]
	}
}


f1<-rep(0,10)
f1[10]=1
f1[9]=ss[1,10]
f1[8]=(ss[1,9]*t[1,8]+ss[2,9]*t[2,8])/(t[1,8]+t[2,8])
f1[7]=(ss[1,8]*t[1,7]+ss[2,8]*t[2,7]+ss[3,8]*t[3,7])/(t[1,7]+t[2,7]+t[3,7])
f1[6]=(ss[1,7]*t[1,6]+ss[2,7]*t[2,6]+ss[3,7]*t[3,6]+ss[4,7]*t[4,6])/	(t[1,6]+t[2,6]+t[3,6]+t[4,6])
f1[5]=(ss[1,6]*t[1,5]+ss[2,6]*t[2,5]+ss[3,6]*t[3,5]+ss[4,6]*t[4,5]+ss[5,6]*t[5,5]) /(t[1,5]+t[2,5]+t[3,5]+t[4,5]+t[5,5])
f1[4]=(ss[1,5]*t[1,4]+ss[2,5]*t[2,4]+ss[3,5]*t[3,4]+ss[4,5]*t[4,4]+ss[5,5]*t[5,4]+ss[6,5]*t[6,4])/(t[1,4]+t[2,4]+t[3,4]+t[4,4]+t[5,4]+t[6,4])
f1[3]=(ss[1,4]*t[1,3]+ss[2,4]*t[2,3]+ss[3,4]*t[3,3]+ss[4,4]*t[4,3]+ss[5,4]*t[5,3]+ss[6,4]*t[6,3]+ss[7,4]*t[7,3])/(t[1,3]+t[2,3]+t[3,3]+t[4,3]+t[5,3]+t[6,3]+t[7,3])
f1[2]=(ss[1,3]*t[1,2]+ss[2,3]*t[2,2]+ss[3,3]*t[3,2]+ss[4,3]*t[4,2]+ss[5,3]*t[5,2]+ss[6,3]*t[6,2]+ss[7,3]*t[7,2]+ss[8,3]*t[8,2])/(t[1,2]+t[2,2]+t[3,2]+t[4,2]+t[5,2]+t[6,2]+t[7,2]+t[8,2])
f1[1]=(ss[1,2]*t[1,1]+ss[2,2]*t[2,1]+ss[3,2]*t[3,1]+ss[4,2]*t[4,1]+ss[5,2]*t[5,1]+ss[6,2]*t[6,1]+ss[7,2]*t[7,1]+ss[8,2]*t[8,1]+ss[9,2]*t[9,1])/(t[1,1]+t[2,1]+t[3,1]+t[4,1]+t[5,1]+t[6,1]+t[7,1]+t[8,1]+t[9,1])

 
tfc1<-array(dim=c(10,10))
tfi1<-array(dim=c(10,10))
tot1=0
for (i in 2:10){
	for (j in (12-i):10){
		if ((i+j)==12){tfc1[i,j]=rnorm(1,t[i,j-1]*f1[j-1],sqrt(aa[j-1]*t[i,j-1]))} 
		else {tfc1[i,j]=rnorm(1,tfc1[i,j-1]*f1[j-1],sqrt(aa[j-1]*abs(tfc1[i,j-1])))}
		#if ((i+j)==12) {tfc1[i,j]=t[i,j-1]*f1[j-1]} 
		#else {tfc1[i,j]=tfc1[i,j-1]*f1[j-1]}
		if ((i+j)==12) 
{tfi1[i,j]=tfc1[i,j]-t[i,j-1]} 
		else {tfi1[i,j]=tfc1[i,j]-tfc1[i,j-1]}
	tot1=tot1+tfi1[i,j]
	}
}

res[k]=tot1
rs2[k]=tfi1[2,10]
rs3[k]=tfi1[3,10]+tfi1[3,9]
rs4[k]=tfi1[4,10]+tfi1[4,9]+tfi1[4,8]
rs5[k]=tfi1[5,10]+tfi1[5,9]+tfi1[5,8]+tfi1[5,7]
rs6[k]=tfi1[6,10]+tfi1[6,9]+tfi1[6,8]+tfi1[6,7]+tfi1[6,6]
rs7[k]=tfi1[7,10]+tfi1[7,9]+tfi1[7,8]+tfi1[7,7]+tfi1[7,6]+tfi1[7,5]
rs8[k]=tfi1[8,10]+tfi1[8,9]+tfi1[8,8]+tfi1[8,7]+tfi1[8,6]+tfi1[8,5]+tfi1[8,4]
rs9[k]=tfi1[9,10]+tfi1[9,9]+tfi1[9,8]+tfi1[9,7]+tfi1[9,6]+tfi1[9,5]+tfi1[9,4]+tfi1[9,3]
rs10[k]=tfi1[10,10]+tfi1[10,9]+tfi1[10,8]+tfi1[10,7]+tfi1[10,6]+tfi1[10,5]+tfi1[10,4]+tfi1[10,3]+tfi1[10,2]
}

 

mean(res)
sd(res)

sd(rs2)
sd(rs3)
sd(rs4)
sd(rs5)
sd(rs6)
sd(rs7)
sd(rs8)
sd(rs9)
sd(rs10)

mack.boot.ta.res.arr[1,,w] <- rs2
mack.boot.ta.res.arr[2,,w] <- rs3
mack.boot.ta.res.arr[3,,w] <- rs4
mack.boot.ta.res.arr[4,,w] <- rs5
mack.boot.ta.res.arr[5,,w] <- rs6
mack.boot.ta.res.arr[6,,w] <- rs7
mack.boot.ta.res.arr[7,,w] <- rs8
mack.boot.ta.res.arr[8,,w] <- rs9
mack.boot.ta.res.arr[9,,w] <- rs10
mack.boot.ta.res.arr[10,,w] <- res
}


# mack.boot.ta.res.arr[,,4]

## Histogram of total reserve estimates for triangle 4
hist(mack.boot.ta.res.arr[10,,4])

## Mean & SD of total reserve estimates for triangle 4
mean(mack.boot.ta.res.arr[10,,4])
sd(mack.boot.ta.res.arr[10,,4])

## Prediction Error Estimate
sd(mack.boot.ta.res.arr[10,,4])/mean(mack.boot.ta.res.arr[10,,4])

## Prediction errors for all triangles for total reserves (NOT! individual years)
p.e.mack.boot.totals <- numeric(50)
for(i in 1:50){
	p.e.mack.boot.totals[i] <- sd(mack.boot.ta.res.arr[10,,i])/mean(mack.boot.ta.res.arr[10,,i])
}

## Prediction errors for all triangles for each year as well as total reserves
p.e.mack.boot <- matrix(0, 50, 10)

for(j in 1:10){
	for(i in 1:50){
		p.e.mack.boot[i,] <- sd(mack.boot.ta.res.arr[j,,i])/	mean(mack.boot.ta.res.arr[j,,i])
	}
}

colnames(p.e.mack.boot) <- c("Year 2", "Year 3", "Year 4", "Year 5", "Year 6", "Year 7", "Year 8", "Year 9", "Year 10", "Totals")
p.e.mack.boot
















