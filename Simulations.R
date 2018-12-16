library(ChainLadder)

## Lognormal Orignial Run

o.1=0 # o counts the number of times the 99 percentile is breached
o.5=0 # o counts the number of times the 95 percentile is breached
o.10=0 # o counts the number of times the 90 percentile is breached
o.20=0 # o counts the number of times the 80 percentile is breached
o.30=0 # o counts the number of times the 70 percentile is breached
o.50=0 # o counts the number of times the 50 percentile is breached
o.70=0 # o counts the number of times the 30 percentile is breached
o.80=0 # o counts the number of times the 20 percentile is breached
o.90=0 # o counts the number of times the 10 percentile is breached
o.95=0 # o counts the number of times the 5 percentile is breached
o.99=0 # o counts the number of times the 1 percentile is breached

oo=0 # oo totals the BCL estimate so you can average it

tab <- matrix(NA, 10,10)
for (tr in 1:1000){


f<-rep(0,9);h<-rep(0,9);EV=0;SD=0;V=0;c<-array(dim=c(10,10));ch<-array(dim=c(10,10));R<-rep(0,10)
f[1]=4.289;f[2]=2.064;f[3]=1.502;f[4]=1.268;f[5]=1.150
f[6]=1.085;f[7]=1.048;f[8]=1.027;f[9]=1.015

ac1=0
for (i in 1:10){
tab[i,1]=rlnorm(1,-0.5*log(2),sqrt(log(2)))
for (j in 2:10){
b=sqrt(log(1+(1/(tab[i,j-1]*(f[j-1]-1)^2))))
a=log(tab[i,j-1]*(f[j-1]-1))-.5*b*b
tab[i,j]=rlnorm(1,a,b)+tab[i,j-1];
}
}

c<-array(dim=c(10,10))
for (i in 1:10){
for (j in 1:(11-i)){
c[i,j]=tab[i,j]
}
}
ac1=c[10,1]+c[9,2]+c[8,3]+c[7,4]+c[6,5]+c[5,6]+c[4,7]+c[3,8]+c[2,9]+c[1,10]
ac2=sum(tab[,10])-ac1


f[9]=c[1,10]/c[1,9]
f[8]=(c[1,9]+c[2,9])/(c[1,8]+c[2,8])
f[7]=(c[1,8]+c[2,8]+c[3,8])/(c[1,7]+c[2,7]+c[3,7])
f[6]=(c[1,7]+c[2,7]+c[3,7]+c[4,7])/(c[1,6]+c[2,6]+c[3,6]+c[4,6])
f[5]=(c[1,6]+c[2,6]+c[3,6]+c[4,6]+c[5,6])/(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5])
f[4]=(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5]+c[6,5])/(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4])
f[3]=(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4]+c[7,4])/(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3])
f[2]=(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3]+c[8,3])/(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2])
f[1]=(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2]+c[9,2])/(c[1,1]+c[2,1]+c[3,1]+c[4,1]+c[5,1]+c[6,1]+c[7,1]+c[8,1]+c[9,1])

for (i in 1:10){
for (j in 1:10){
if (i+j<12) {ch[i,j]<-c[i,j]} else {ch[i,j]<-ch[i,j-1]*f[j-1]}
}
}


for (i in 1:8){
a=0;for (j in 1:(10-i)){
a<-a+c[j,i]*(c[j,i+1]/c[j,i]-f[i])^2
}
h[i]=a/(9-i)
}
h[9]=min(h[8],h[7])

for (i in 2:10){
a=0;for (k in (11-i):9){
a1=0;for (j in 1:(10-k)){
a1<-a1+c[j,k]}
a<-a+((h[k]/(f[k])^2))*(1/ch[i,k]+1/a1)
}
R[i]=sqrt(a*(ch[i,10])^2)
}

V=0;for (i in 2:10){
a1=R[i]^2
b1=0
if (i<10) {for (j in (1+i):10){b1<-b1+ch[j,10]};a2=b1*ch[i,10]} else {a2=ch[i,10]}
b1=0;for (k in (11-i):9){
b2=0;for (n in 1:(10-k)){b2<-b2+c[n,k]}
b1<-b1+(2*h[k]/(f[k]^2))/b2
}
V=V+a1+a2*b1
}

EV=c[10,1]*f[1]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[9,2]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[8,3]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[7,4]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[6,5]*f[5]*f[6]*f[7]*f[8]*f[9]+c[5,6]*f[6]*f[7]*f[8]*f[9]+c[4,7]*f[7]*f[8]*f[9]+c[3,8]*f[8]*f[9]+c[2,9]*f[9]+c[1,10]-ac1
SD=sqrt(V)

sigma=sqrt(log((1+V/(EV^2))))
mu=log(EV)-0.5*sigma^2
p=plnorm(ac2,mu,sigma)
if (p>0.99) {o.1 <- o.1+1}
if (p>0.95) {o.5 <- o.5+1}
if (p>0.90) {o.10 <- o.10+1}
if (p>0.80) {o.20 <- o.20+1}
if (p>0.70) {o.30 <- o.30+1}
if (p>0.50) {o.50 <- o.50+1}
if (p>0.30) {o.70 <- o.70+1}
if (p>0.20) {o.80 <- o.80+1}
if (p>0.10) {o.90 <- o.90+1}
if (p>0.05) {o.95 <- o.95+1}
if (p>0.01) {o.99 <- o.99+1}
oo=oo+EV
}


o.1/1000
o.5/1000
o.10/1000
o.20/1000
o.30/1000
o.50/1000
o.70/1000
o.80/1000
o.90/1000
o.95/1000
o.99/1000

oo

# ----------------------------------------------------------------------------------- #

## Now looking at changing the standard deviation to 2 instead of 1:

x.1=0 # o counts the number of times the 99 percentile is breached
x.5=0 # o counts the number of times the 95 percentile is breached
x.10=0 # o counts the number of times the 90 percentile is breached
x.20=0 # o counts the number of times the 80 percentile is breached
x.30=0 # o counts the number of times the 70 percentile is breached
x.50=0 # o counts the number of times the 50 percentile is breached
x.70=0 # o counts the number of times the 30 percentile is breached
x.80=0 # o counts the number of times the 20 percentile is breached
x.90=0 # o counts the number of times the 10 percentile is breached
x.95=0 # o counts the number of times the 5 percentile is breached
x.99=0 # o counts the number of times the 1 percentile is breached
xx=0 # oo totals the BCL estimate so you can average it

tab <- matrix(NA,10,10)
for (tr in 1:1000){


f<-rep(0,9);h<-rep(0,9);EV=0;SD=0;V=0;c<-array(dim=c(10,10));ch<-array(dim=c(10,10));R<-rep(0,10)
f[1]=4.289;f[2]=2.064;f[3]=1.502;f[4]=1.268;f[5]=1.150
f[6]=1.085;f[7]=1.048;f[8]=1.027;f[9]=1.015

## generating a triangle
ac1=0
for (i in 1:10){
tab[i,1]=rlnorm(1,-0.5*log(5),sqrt(log(5)))
for (j in 2:10){
b=sqrt(log(1+(4/(tab[i,j-1]*(f[j-1]-1)^2))))
a=log(tab[i,j-1]*(f[j-1]-1))-.5*b*b
tab[i,j]=rlnorm(1,a,b)+tab[i,j-1];
}
}

c<-array(dim=c(10,10))
for (i in 1:10){
for (j in 1:(11-i)){
c[i,j]=tab[i,j]
}
}
ac1=c[10,1]+c[9,2]+c[8,3]+c[7,4]+c[6,5]+c[5,6]+c[4,7]+c[3,8]+c[2,9]+c[1,10]
ac2=sum(tab[,10])-ac1


f[9]=c[1,10]/c[1,9]
f[8]=(c[1,9]+c[2,9])/(c[1,8]+c[2,8])
f[7]=(c[1,8]+c[2,8]+c[3,8])/(c[1,7]+c[2,7]+c[3,7])
f[6]=(c[1,7]+c[2,7]+c[3,7]+c[4,7])/(c[1,6]+c[2,6]+c[3,6]+c[4,6])
f[5]=(c[1,6]+c[2,6]+c[3,6]+c[4,6]+c[5,6])/(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5])
f[4]=(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5]+c[6,5])/(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4])
f[3]=(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4]+c[7,4])/(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3])
f[2]=(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3]+c[8,3])/(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2])
f[1]=(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2]+c[9,2])/(c[1,1]+c[2,1]+c[3,1]+c[4,1]+c[5,1]+c[6,1]+c[7,1]+c[8,1]+c[9,1])

for (i in 1:10){
for (j in 1:10){
if (i+j<12) {ch[i,j]<-c[i,j]} else {ch[i,j]<-ch[i,j-1]*f[j-1]}
}
}


for (i in 1:8){
a=0;for (j in 1:(10-i)){
a<-a+c[j,i]*(c[j,i+1]/c[j,i]-f[i])^2
}
h[i]=a/(9-i)
}
h[9]=min(h[8],h[7])

for (i in 2:10){
a=0;for (k in (11-i):9){
a1=0;for (j in 1:(10-k)){
a1<-a1+c[j,k]}
a<-a+((h[k]/(f[k])^2))*(1/ch[i,k]+1/a1)
}
R[i]=sqrt(a*(ch[i,10])^2)
}

V=0;for (i in 2:10){
a1=R[i]^2
b1=0
if (i<10) {for (j in (1+i):10){b1<-b1+ch[j,10]};a2=b1*ch[i,10]} else {a2=ch[i,10]}
b1=0;for (k in (11-i):9){
b2=0;for (n in 1:(10-k)){b2<-b2+c[n,k]}
b1<-b1+(2*h[k]/(f[k]^2))/b2
}
V=V+a1+a2*b1
}

EV=c[10,1]*f[1]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[9,2]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[8,3]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[7,4]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[6,5]*f[5]*f[6]*f[7]*f[8]*f[9]+c[5,6]*f[6]*f[7]*f[8]*f[9]+c[4,7]*f[7]*f[8]*f[9]+c[3,8]*f[8]*f[9]+c[2,9]*f[9]+c[1,10]-ac1
SD=sqrt(V)

sigma=sqrt(log((1+V/(EV^2))))
mu=log(EV)-0.5*sigma^2
p=plnorm(ac2,mu,sigma)
if (p>0.99) {x.1 <- x.1+1}
if (p>0.95) {x.5 <- x.5+1}
if (p>0.90) {x.10 <- x.10+1}
if (p>0.80) {x.20 <- x.20+1}
if (p>0.70) {x.30 <- x.30+1}
if (p>0.50) {x.50 <- x.50+1}
if (p>0.30) {x.70 <- x.70+1}
if (p>0.20) {x.80 <- x.80+1}
if (p>0.10) {x.90 <- x.90+1}
if (p>0.05) {x.95 <- x.95+1}
if (p>0.01) {x.99 <- x.99+1}
xx=xx+EV
}

x.1/1000
x.5/1000
x.10/1000
x.20/1000
x.30/1000
x.50/1000
x.70/1000
x.80/1000
x.90/1000
x.95/1000
x.99/1000

xx

# ----------------------------------------------------------------- #

## Now looking at changing the standard deviation and mean to 5 instead of 1:

v.1=0 # o counts the number of times the 99 percentile is breached
v.5=0 # o counts the number of times the 95 percentile is breached
v.10=0 # o counts the number of times the 90 percentile is breached
v.20=0 # o counts the number of times the 80 percentile is breached
v.30=0 # o counts the number of times the 70 percentile is breached
v.50=0 # o counts the number of times the 50 percentile is breached
v.70=0 # o counts the number of times the 30 percentile is breached
v.80=0 # o counts the number of times the 20 percentile is breached
v.90=0 # o counts the number of times the 10 percentile is breached
v.95=0 # o counts the number of times the 5 percentile is breached
v.99=0 # o counts the number of times the 1 percentile is breached
vv=0 # oo totals the BCL estimate so you can average it

tab <- matrix(NA,10,10)
for (tr in 1:1000){


f<-rep(0,9);h<-rep(0,9);EV=0;SD=0;V=0;c<-array(dim=c(10,10));ch<-array(dim=c(10,10));R<-rep(0,10)
f[1]=4.289;f[2]=2.064;f[3]=1.502;f[4]=1.268;f[5]=1.150
f[6]=1.085;f[7]=1.048;f[8]=1.027;f[9]=1.015

## generating a triangle
ac1=0
for (i in 1:10){
tab[i,1]=rlnorm(1,(log(5)-0.5*log(2)),sqrt(log(2)))
for (j in 2:10){
b=sqrt(log(1+(25/(5 + tab[i,j-1]*(f[j-1]-1)^2))))
a=log(5 + (tab[i,j-1]*(f[j-1]-1)))-.5*b*b
tab[i,j]=rlnorm(1,a,b)+tab[i,j-1];
}
}

c<-array(dim=c(10,10))
for (i in 1:10){
for (j in 1:(11-i)){
c[i,j]=tab[i,j]
}
}
ac1=c[10,1]+c[9,2]+c[8,3]+c[7,4]+c[6,5]+c[5,6]+c[4,7]+c[3,8]+c[2,9]+c[1,10]
ac2=sum(tab[,10])-ac1


f[9]=c[1,10]/c[1,9]
f[8]=(c[1,9]+c[2,9])/(c[1,8]+c[2,8])
f[7]=(c[1,8]+c[2,8]+c[3,8])/(c[1,7]+c[2,7]+c[3,7])
f[6]=(c[1,7]+c[2,7]+c[3,7]+c[4,7])/(c[1,6]+c[2,6]+c[3,6]+c[4,6])
f[5]=(c[1,6]+c[2,6]+c[3,6]+c[4,6]+c[5,6])/(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5])
f[4]=(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5]+c[6,5])/(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4])
f[3]=(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4]+c[7,4])/(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3])
f[2]=(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3]+c[8,3])/(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2])
f[1]=(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2]+c[9,2])/(c[1,1]+c[2,1]+c[3,1]+c[4,1]+c[5,1]+c[6,1]+c[7,1]+c[8,1]+c[9,1])

for (i in 1:10){
for (j in 1:10){
if (i+j<12) {ch[i,j]<-c[i,j]} else {ch[i,j]<-ch[i,j-1]*f[j-1]}
}
}


for (i in 1:8){
a=0;for (j in 1:(10-i)){
a<-a+c[j,i]*(c[j,i+1]/c[j,i]-f[i])^2
}
h[i]=a/(9-i)
}
h[9]=min(h[8],h[7])

for (i in 2:10){
a=0;for (k in (11-i):9){
a1=0;for (j in 1:(10-k)){
a1<-a1+c[j,k]}
a<-a+((h[k]/(f[k])^2))*(1/ch[i,k]+1/a1)
}
R[i]=sqrt(a*(ch[i,10])^2)
}

V=0;for (i in 2:10){
a1=R[i]^2
b1=0
if (i<10) {for (j in (1+i):10){b1<-b1+ch[j,10]};a2=b1*ch[i,10]} else {a2=ch[i,10]}
b1=0;for (k in (11-i):9){
b2=0;for (n in 1:(10-k)){b2<-b2+c[n,k]}
b1<-b1+(2*h[k]/(f[k]^2))/b2
}
V=V+a1+a2*b1
}

EV=c[10,1]*f[1]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[9,2]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[8,3]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[7,4]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[6,5]*f[5]*f[6]*f[7]*f[8]*f[9]+c[5,6]*f[6]*f[7]*f[8]*f[9]+c[4,7]*f[7]*f[8]*f[9]+c[3,8]*f[8]*f[9]+c[2,9]*f[9]+c[1,10]-ac1
SD=sqrt(V)

sigma=sqrt(log((1+V/(EV^2))))
mu=log(EV)-0.5*sigma^2
p=plnorm(ac2,mu,sigma)
if (p>0.99) {v.1 <- v.1+1}
if (p>0.95) {v.5 <- v.5+1}
if (p>0.90) {v.10 <- v.10+1}
if (p>0.80) {v.20 <- v.20+1}
if (p>0.70) {v.30 <- v.30+1}
if (p>0.50) {v.50 <- v.50+1}
if (p>0.30) {v.70 <- v.70+1}
if (p>0.20) {v.80 <- v.80+1}
if (p>0.10) {v.90 <- v.90+1}
if (p>0.05) {v.95 <- v.95+1}
if (p>0.01) {v.99 <- v.99+1}
vv=vv+EV
}

v.1/1000
v.5/1000
v.10/1000
v.20/1000
v.30/1000
v.50/1000
v.70/1000
v.80/1000
v.90/1000
v.95/1000
v.99/1000

vv

# ------------------------------------------------------------------------------------ #

## Now looking at changing the standard deviation to 0.5 instead of 1:

u.1=0 # o counts the number of times the 99 percentile is breached
u.5=0 # o counts the number of times the 95 percentile is breached
u.10=0 # o counts the number of times the 90 percentile is breached
u.20=0 # o counts the number of times the 80 percentile is breached
u.30=0 # o counts the number of times the 70 percentile is breached
u.50=0 # o counts the number of times the 50 percentile is breached
u.70=0 # o counts the number of times the 30 percentile is breached
u.80=0 # o counts the number of times the 20 percentile is breached
u.90=0 # o counts the number of times the 10 percentile is breached
u.95=0 # o counts the number of times the 5 percentile is breached
u.99=0 # o counts the number of times the 1 percentile is breached
uu=0 # oo totals the BCL estimate so you can average it

tab <- matrix(NA,10,10)
for (tr in 1:1000){


f<-rep(0,9);h<-rep(0,9);EV=0;SD=0;V=0;c<-array(dim=c(10,10));ch<-array(dim=c(10,10));R<-rep(0,10)
f[1]=4.289;f[2]=2.064;f[3]=1.502;f[4]=1.268;f[5]=1.150
f[6]=1.085;f[7]=1.048;f[8]=1.027;f[9]=1.015

## generating a triangle
ac1=0
for (i in 1:10){
tab[i,1]=rlnorm(1,-0.5*log(1.5),sqrt(log(1.5)))
for (j in 2:10){
b=sqrt(log(1+(0.25/(tab[i,j-1]*(f[j-1]-1)^2))))
a=log(tab[i,j-1]*(f[j-1]-1))-.5*b*b
tab[i,j]=rlnorm(1,a,b)+tab[i,j-1];
}
}

c<-array(dim=c(10,10))
for (i in 1:10){
for (j in 1:(11-i)){
c[i,j]=tab[i,j]
}
}
ac1=c[10,1]+c[9,2]+c[8,3]+c[7,4]+c[6,5]+c[5,6]+c[4,7]+c[3,8]+c[2,9]+c[1,10]
ac2=sum(tab[,10])-ac1


f[9]=c[1,10]/c[1,9]
f[8]=(c[1,9]+c[2,9])/(c[1,8]+c[2,8])
f[7]=(c[1,8]+c[2,8]+c[3,8])/(c[1,7]+c[2,7]+c[3,7])
f[6]=(c[1,7]+c[2,7]+c[3,7]+c[4,7])/(c[1,6]+c[2,6]+c[3,6]+c[4,6])
f[5]=(c[1,6]+c[2,6]+c[3,6]+c[4,6]+c[5,6])/(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5])
f[4]=(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5]+c[6,5])/(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4])
f[3]=(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4]+c[7,4])/(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3])
f[2]=(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3]+c[8,3])/(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2])
f[1]=(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2]+c[9,2])/(c[1,1]+c[2,1]+c[3,1]+c[4,1]+c[5,1]+c[6,1]+c[7,1]+c[8,1]+c[9,1])

for (i in 1:10){
for (j in 1:10){
if (i+j<12) {ch[i,j]<-c[i,j]} else {ch[i,j]<-ch[i,j-1]*f[j-1]}
}
}


for (i in 1:8){
a=0;for (j in 1:(10-i)){
a<-a+c[j,i]*(c[j,i+1]/c[j,i]-f[i])^2
}
h[i]=a/(9-i)
}
h[9]=min(h[8],h[7])

for (i in 2:10){
a=0;for (k in (11-i):9){
a1=0;for (j in 1:(10-k)){
a1<-a1+c[j,k]}
a<-a+((h[k]/(f[k])^2))*(1/ch[i,k]+1/a1)
}
R[i]=sqrt(a*(ch[i,10])^2)
}

V=0;for (i in 2:10){
a1=R[i]^2
b1=0
if (i<10) {for (j in (1+i):10){b1<-b1+ch[j,10]};a2=b1*ch[i,10]} else {a2=ch[i,10]}
b1=0;for (k in (11-i):9){
b2=0;for (n in 1:(10-k)){b2<-b2+c[n,k]}
b1<-b1+(2*h[k]/(f[k]^2))/b2
}
V=V+a1+a2*b1
}

EV=c[10,1]*f[1]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[9,2]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[8,3]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[7,4]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[6,5]*f[5]*f[6]*f[7]*f[8]*f[9]+c[5,6]*f[6]*f[7]*f[8]*f[9]+c[4,7]*f[7]*f[8]*f[9]+c[3,8]*f[8]*f[9]+c[2,9]*f[9]+c[1,10]-ac1
SD=sqrt(V)

sigma=sqrt(log((1+V/(EV^2))))
mu=log(EV)-0.5*sigma^2
p=plnorm(ac2,mu,sigma)
if (p>0.99) {u.1 <- u.1+1}
if (p>0.95) {u.5 <- u.5+1}
if (p>0.90) {u.10 <- u.10+1}
if (p>0.80) {u.20 <- u.20+1}
if (p>0.70) {u.30 <- u.30+1}
if (p>0.50) {u.50 <- u.50+1}
if (p>0.30) {u.70 <- u.70+1}
if (p>0.20) {u.80 <- u.80+1}
if (p>0.10) {u.90 <- u.90+1}
if (p>0.05) {u.95 <- u.95+1}
if (p>0.01) {u.99 <- u.99+1}
uu=uu+EV
}

u.1/1000
u.5/1000
u.10/1000
u.20/1000
u.30/1000
u.50/1000
u.70/1000
u.80/1000
u.90/1000
u.95/1000
u.99/1000

uu

# ------------------------------------------------------------------------------------ #


## Faster Run Off

y.1=0 # o counts the number of times the 99 percentile is breached
y.5=0 # o counts the number of times the 95 percentile is breached
y.10=0 # o counts the number of times the 90 percentile is breached
y.20=0 # o counts the number of times the 80 percentile is breached
y.30=0 # o counts the number of times the 70 percentile is breached
y.50=0 # o counts the number of times the 50 percentile is breached
y.70=0 # o counts the number of times the 30 percentile is breached
y.80=0 # o counts the number of times the 20 percentile is breached
y.90=0 # o counts the number of times the 10 percentile is breached
y.95=0 # o counts the number of times the 5 percentile is breached
y.99=0 # o counts the number of times the 1 percentile is breached
yy=0 # oo totals the BCL estimate so you can average it


tab <- matrix(NA,10,10)
for (tr in 1:1000){


f<-rep(0,9);h<-rep(0,9);EV=0;SD=0;V=0;c<-array(dim=c(10,10));ch<-array(dim=c(10,10));R<-rep(0,10)
f[1]=6.029;f[2]=4.564;f[3]=2.301;f[4]=1.1;f[5]=1.05
f[6]=1.005;f[7]=1.048;f[8]=1.027;f[9]=1.015

ac1=0
for (i in 1:10){
tab[i,1]=rlnorm(1,-0.5*log(2),sqrt(log(2)))
for (j in 2:10){
b=sqrt(log(1+(1/(tab[i,j-1]*(f[j-1]-1)^2))))
a=log(tab[i,j-1]*(f[j-1]-1))-.5*b*b
tab[i,j]=rlnorm(1,a,b)+tab[i,j-1];
}
}

c<-array(dim=c(10,10))
for (i in 1:10){
for (j in 1:(11-i)){
c[i,j]=tab[i,j]
}
}
ac1=c[10,1]+c[9,2]+c[8,3]+c[7,4]+c[6,5]+c[5,6]+c[4,7]+c[3,8]+c[2,9]+c[1,10]
ac2=sum(tab[,10])-ac1


f[9]=c[1,10]/c[1,9]
f[8]=(c[1,9]+c[2,9])/(c[1,8]+c[2,8])
f[7]=(c[1,8]+c[2,8]+c[3,8])/(c[1,7]+c[2,7]+c[3,7])
f[6]=(c[1,7]+c[2,7]+c[3,7]+c[4,7])/(c[1,6]+c[2,6]+c[3,6]+c[4,6])
f[5]=(c[1,6]+c[2,6]+c[3,6]+c[4,6]+c[5,6])/(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5])
f[4]=(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5]+c[6,5])/(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4])
f[3]=(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4]+c[7,4])/(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3])
f[2]=(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3]+c[8,3])/(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2])
f[1]=(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2]+c[9,2])/(c[1,1]+c[2,1]+c[3,1]+c[4,1]+c[5,1]+c[6,1]+c[7,1]+c[8,1]+c[9,1])

for (i in 1:10){
for (j in 1:10){
if (i+j<12) {ch[i,j]<-c[i,j]} else {ch[i,j]<-ch[i,j-1]*f[j-1]}
}
}


for (i in 1:8){
a=0;for (j in 1:(10-i)){
a<-a+c[j,i]*(c[j,i+1]/c[j,i]-f[i])^2
}
h[i]=a/(9-i)
}
h[9]=min(h[8],h[7])

for (i in 2:10){
a=0;for (k in (11-i):9){
a1=0;for (j in 1:(10-k)){
a1<-a1+c[j,k]}
a<-a+((h[k]/(f[k])^2))*(1/ch[i,k]+1/a1)
}
R[i]=sqrt(a*(ch[i,10])^2)
}

V=0;for (i in 2:10){
a1=R[i]^2
b1=0
if (i<10) {for (j in (1+i):10){b1<-b1+ch[j,10]};a2=b1*ch[i,10]} else {a2=ch[i,10]}
b1=0;for (k in (11-i):9){
b2=0;for (n in 1:(10-k)){b2<-b2+c[n,k]}
b1<-b1+(2*h[k]/(f[k]^2))/b2
}
V=V+a1+a2*b1
}

EV=c[10,1]*f[1]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[9,2]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[8,3]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[7,4]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[6,5]*f[5]*f[6]*f[7]*f[8]*f[9]+c[5,6]*f[6]*f[7]*f[8]*f[9]+c[4,7]*f[7]*f[8]*f[9]+c[3,8]*f[8]*f[9]+c[2,9]*f[9]+c[1,10]-ac1
SD=sqrt(V)

sigma=sqrt(log((1+V/(EV^2))))
mu=log(EV)-0.5*sigma^2
p=plnorm(ac2,mu,sigma)
if (p>0.99) {y.1 <- y.1+1}
if (p>0.95) {y.5 <- y.5+1}
if (p>0.90) {y.10 <- y.10+1}
if (p>0.80) {y.20 <- y.20+1}
if (p>0.70) {y.30 <- y.30+1}
if (p>0.50) {y.50 <- y.50+1}
if (p>0.30) {y.70 <- y.70+1}
if (p>0.20) {y.80 <- y.80+1}
if (p>0.10) {y.90 <- y.90+1}
if (p>0.05) {y.95 <- y.95+1}
if (p>0.01) {y.99 <- y.99+1}
yy=yy+EV
}

y.1/1000
y.5/1000
y.10/1000
y.20/1000
y.30/1000
y.50/1000
y.70/1000
y.80/1000
y.90/1000
y.95/1000
y.99/1000

yy
	
# ----------------------------------------------------------------------------------- #

## Slower Run Off

z.1=0 # o counts the number of times the 99 percentile is breached
z.5=0 # o counts the number of times the 95 percentile is breached
z.10=0 # o counts the number of times the 90 percentile is breached
z.20=0 # o counts the number of times the 80 percentile is breached
z.30=0 # o counts the number of times the 70 percentile is breached
z.50=0 # o counts the number of times the 50 percentile is breached
z.70=0 # o counts the number of times the 30 percentile is breached
z.80=0 # o counts the number of times the 20 percentile is breached
z.90=0 # o counts the number of times the 10 percentile is breached
z.95=0 # o counts the number of times the 5 percentile is breached
z.99=0 # o counts the number of times the 1 percentile is breached
zz=0 # oo totals the BCL estimate so you can average it

tab <- matrix(NA,10,10)
for (tr in 1:1000){


f<-rep(0,9);h<-rep(0,9);EV=0;SD=0;V=0;c<-array(dim=c(10,10));ch<-array(dim=c(10,10));R<-rep(0,10)
f[1]=2.525;f[2]=2.464;f[3]=2.268;f[4]=2.189;f[5]=2.08
f[6]=1.958;f[7]=1.748;f[8]=1.427;f[9]=1.015

ac1=0
for (i in 1:10){
tab[i,1]=rlnorm(1,-0.5*log(2),sqrt(log(2)))
for (j in 2:10){
b=sqrt(log(1+(1/(tab[i,j-1]*(f[j-1]-1)^2))))
a=log(tab[i,j-1]*(f[j-1]-1))-.5*b*b
tab[i,j]=rlnorm(1,a,b)+tab[i,j-1];
}
}

c<-array(dim=c(10,10))
for (i in 1:10){
for (j in 1:(11-i)){
c[i,j]=tab[i,j]
}
}
ac1=c[10,1]+c[9,2]+c[8,3]+c[7,4]+c[6,5]+c[5,6]+c[4,7]+c[3,8]+c[2,9]+c[1,10]
ac2=sum(tab[,10])-ac1


f[9]=c[1,10]/c[1,9]
f[8]=(c[1,9]+c[2,9])/(c[1,8]+c[2,8])
f[7]=(c[1,8]+c[2,8]+c[3,8])/(c[1,7]+c[2,7]+c[3,7])
f[6]=(c[1,7]+c[2,7]+c[3,7]+c[4,7])/(c[1,6]+c[2,6]+c[3,6]+c[4,6])
f[5]=(c[1,6]+c[2,6]+c[3,6]+c[4,6]+c[5,6])/(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5])
f[4]=(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5]+c[6,5])/(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4])
f[3]=(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4]+c[7,4])/(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3])
f[2]=(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3]+c[8,3])/(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2])
f[1]=(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2]+c[9,2])/(c[1,1]+c[2,1]+c[3,1]+c[4,1]+c[5,1]+c[6,1]+c[7,1]+c[8,1]+c[9,1])

for (i in 1:10){
for (j in 1:10){
if (i+j<12) {ch[i,j]<-c[i,j]} else {ch[i,j]<-ch[i,j-1]*f[j-1]}
}
}


for (i in 1:8){
a=0;for (j in 1:(10-i)){
a<-a+c[j,i]*(c[j,i+1]/c[j,i]-f[i])^2
}
h[i]=a/(9-i)
}
h[9]=min(h[8],h[7])

for (i in 2:10){
a=0;for (k in (11-i):9){
a1=0;for (j in 1:(10-k)){
a1<-a1+c[j,k]}
a<-a+((h[k]/(f[k])^2))*(1/ch[i,k]+1/a1)
}
R[i]=sqrt(a*(ch[i,10])^2)
}

V=0;for (i in 2:10){
a1=R[i]^2
b1=0
if (i<10) {for (j in (1+i):10){b1<-b1+ch[j,10]};a2=b1*ch[i,10]} else {a2=ch[i,10]}
b1=0;for (k in (11-i):9){
b2=0;for (n in 1:(10-k)){b2<-b2+c[n,k]}
b1<-b1+(2*h[k]/(f[k]^2))/b2
}
V=V+a1+a2*b1
}

EV=c[10,1]*f[1]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[9,2]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[8,3]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[7,4]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]+c[6,5]*f[5]*f[6]*f[7]*f[8]*f[9]+c[5,6]*f[6]*f[7]*f[8]*f[9]+c[4,7]*f[7]*f[8]*f[9]+c[3,8]*f[8]*f[9]+c[2,9]*f[9]+c[1,10]-ac1
SD=sqrt(V)

sigma=sqrt(log((1+V/(EV^2))))
mu=log(EV)-0.5*sigma^2
p=plnorm(ac2,mu,sigma)
if (p>0.99) {z.1 <- z.1+1}
if (p>0.95) {z.5 <- z.5+1}
if (p>0.90) {z.10 <- z.10+1}
if (p>0.80) {z.20 <- z.20+1}
if (p>0.70) {z.30 <- z.30+1}
if (p>0.50) {z.50 <- z.50+1}
if (p>0.30) {z.70 <- z.70+1}
if (p>0.20) {z.80 <- z.80+1}
if (p>0.10) {z.90 <- z.90+1}
if (p>0.05) {z.95 <- z.95+1}
if (p>0.01) {z.99 <- z.99+1}
zz=zz+EV
}

z.1/1000
z.5/1000
z.10/1000
z.20/1000
z.30/1000
z.50/1000
z.70/1000
z.80/1000
z.90/1000
z.95/1000
z.99/1000

zz

# ----------------------------------------------------------------------------------- #

## Smaller Triangle 5x5

z.5.1=0 # o counts the number of times the 99 percentile is breached
z.5.5=0 # o counts the number of times the 95 percentile is breached
z.5.10=0 # o counts the number of times the 90 percentile is breached
z.5.20=0 # o counts the number of times the 80 percentile is breached
z.5.30=0 # o counts the number of times the 70 percentile is breached
z.5.50=0 # o counts the number of times the 50 percentile is breached
z.5.70=0 # o counts the number of times the 30 percentile is breached
z.5.80=0 # o counts the number of times the 20 percentile is breached
z.5.90=0 # o counts the number of times the 10 percentile is breached
z.5.95=0 # o counts the number of times the 5 percentile is breached
z.5.99=0 # o counts the number of times the 1 percentile is breached
zz.5=0 # oo totals the BCL estimate so you can average it

tab <- matrix(NA,5,5)
for (tr in 1:1000){


f<-rep(0,4);h<-rep(0,4);EV=0;SD=0;V=0;c<-array(dim=c(5,5));ch<-array(dim=c(5,5));R<-rep(0,5)
f[1]=3.289;f[2]=2.064;f[3]=1.268;f[4]=1.015

ac1=0
for (i in 1:5){
tab[i,1]=rlnorm(1,-0.5*log(2),sqrt(log(2)))
for (j in 2:5){
b=sqrt(log(1+(1/(tab[i,j-1]*(f[j-1]-1)^2))))
a=log(tab[i,j-1]*(f[j-1]-1))-.5*b*b
tab[i,j]=rlnorm(1,a,b)+tab[i,j-1];
}
}

c<-array(dim=c(5,5))
for (i in 1:5){
for (j in 1:(6-i)){
c[i,j]=tab[i,j]
}
}
ac1=c[5,1]+c[4,2]+c[3,3]+c[2,4]+c[1,5]
ac2=sum(tab[,5])-ac1


f[4]=c[1,5]/c[1,4]
f[3]=(c[1,4]+c[2,4])/(c[1,3]+c[2,3])
f[2]=(c[1,3]+c[2,3]+c[3,3])/(c[1,2]+c[2,2]+c[3,2])
f[1]=(c[1,2]+c[2,2]+c[3,2]+c[4,2])/(c[1,1]+c[2,1]+c[3,1]+c[4,1])

for (i in 1:5){
for (j in 1:5){
if (i+j<7) {ch[i,j]<-c[i,j]} else {ch[i,j]<-ch[i,j-1]*f[j-1]}
}
}


for (i in 1:3){
a=0;for (j in 1:(5-i)){
a<-a+c[j,i]*(c[j,i+1]/c[j,i]-f[i])^2
}
h[i]=a/(4-i)
}
h[4]=min(h[3],h[2])

for (i in 2:5){
a=0;for (k in (6-i):4){
a1=0;for (j in 1:(5-k)){
a1<-a1+c[j,k]}
a<-a+((h[k]/(f[k])^2))*(1/ch[i,k]+1/a1)
}
R[i]=sqrt(a*(ch[i,5])^2)
}

V=0;for (i in 2:5){
a1=R[i]^2
b1=0
if (i<5) {for (j in (1+i):5){b1<-b1+ch[j,5]};a2=b1*ch[i,5]} else {a2=ch[i,5]}
b1=0;for (k in (6-i):4){
b2=0;for (n in 1:(5-k)){b2<-b2+c[n,k]}
b1<-b1+(2*h[k]/(f[k]^2))/b2
}
V=V+a1+a2*b1
}

EV=c[5,1]*f[1]*f[2]*f[3]*f[4] + c[4,2]*f[2]*f[3]*f[4] + c[3,3]*f[3]*f[4] + c[2,4]*f[4] + c[1,5] - ac1
SD=sqrt(V)

sigma=sqrt(log((1+V/(EV^2))))
mu=log(EV)-0.5*sigma^2
p=plnorm(ac2,mu,sigma)
if (p>0.99) {z.5.1 <- z.5.1+1}
if (p>0.95) {z.5.5 <- z.5.5+1}
if (p>0.90) {z.5.10 <- z.5.10+1}
if (p>0.80) {z.5.20 <- z.5.20+1}
if (p>0.70) {z.5.30 <- z.5.30+1}
if (p>0.50) {z.5.50 <- z.5.50+1}
if (p>0.30) {z.5.70 <- z.5.70+1}
if (p>0.20) {z.5.80 <- z.5.80+1}
if (p>0.10) {z.5.90 <- z.5.90+1}
if (p>0.05) {z.5.95 <- z.5.95+1}
if (p>0.01) {z.5.99 <- z.5.99+1} 
zz.5=zz.5+EV
}

z.5.1/1000
z.5.5/1000
z.5.10/1000
z.5.20/1000
z.5.30/1000
z.5.50/1000
z.5.70/1000
z.5.80/1000
z.5.90/1000
z.5.95/1000
z.5.99/1000

zz.5

# ----------------------------------------------------------------------------------- #
## Running Simulations for a 15x15 Triangle

z.15.1=0 # o counts the number of times the 99 percentile is breached
z.15.5=0 # o counts the number of times the 95 percentile is breached
z.15.10=0 # o counts the number of times the 90 percentile is breached
z.15.20=0 # o counts the number of times the 80 percentile is breached
z.15.30=0 # o counts the number of times the 70 percentile is breached
z.15.50=0 # o counts the number of times the 50 percentile is breached
z.15.70=0 # o counts the number of times the 30 percentile is breached
z.15.80=0 # o counts the number of times the 20 percentile is breached
z.15.90=0 # o counts the number of times the 10 percentile is breached
z.15.95=0 # o counts the number of times the 5 percentile is breached
z.15.99=0 # o counts the number of times the 1 percentile is breached
zz.15=0 # oo totals the BCL estimate so you can average it


tab <- matrix(NA,15,15)
for (tr in 1:1000){


f<-rep(0,14);h<-rep(0,14);EV=0;SD=0;V=0;c<-array(dim=c(15,15));ch<-array(dim=c(15,15));R<-rep(0,15)
f[1]=4.64;f[2]=3.156;f[3]=2.525;f[4]=2.064;f[5]=1.748
f[6]=1.502;f[7]=1.268;f[8]=1.15;f[9]=1.085;f[10]=1.048;
f[11]=1.027;f[12]=1.019;f[13]=1.008;f[14]=1.001

ac1=0
for (i in 1:15){
tab[i,1]=rlnorm(1,-0.5*log(2),sqrt(log(2)))
for (j in 2:15){
b=sqrt(log(1+(1/(tab[i,j-1]*(f[j-1]-1)^2))))
a=log(tab[i,j-1]*(f[j-1]-1))-.5*b*b
tab[i,j]=rlnorm(1,a,b)+tab[i,j-1];
}
}

c<-array(dim=c(15,15))
for (i in 1:15){
for (j in 1:(16-i)){
c[i,j]=tab[i,j]
}
}
ac1=c[15,1]+c[14,2]+c[13,3]+c[12,4]+c[11,5]+c[10,6]+c[9,7]+c[8,8]+c[7,9]+c[6,10]+c[5,11]+c[4,12]+c[3,13]+c[2,14]+c[1,15]
ac2=sum(tab[,15])-ac1

f[14]=c[1,15]/c[1,14]
f[13]=(c[1,14]+c[2,14])/(c[1,13]+c[2,13])
f[12]=(c[1,13]+c[2,13]+c[3,13])/(c[1,12]+c[2,12]+c[3,12])
f[11]=(c[1,12]+c[2,12]+c[3,12]+c[4,12])/(c[1,11]+c[2,11]+c[3,11]+c[4,11])
f[10]=(c[1,11]+c[2,11]+c[3,11]+c[4,11]+c[5,11])/(c[1,10]+c[2,10]+c[3,10]+c[4,10]+c[5,10])
f[9]=(c[1,10]+c[2,10]+c[3,10]+c[4,10]+c[5,10]+c[6,10])/(c[1,9]+c[2,9]+c[3,9]+c[4,9]+c[5,9]+c[6,9])
f[8]=(c[1,9]+c[2,9]+c[3,9]+c[4,9]+c[5,9]+c[6,9]+c[7,9])/(c[1,8]+c[2,8]+c[3,8]+c[4,8]+c[5,8]+c[6,8]+c[7,8])
f[7]=(c[1,8]+c[2,8]+c[3,8]+c[4,8]+c[5,8]+c[6,8]+c[7,8]+c[8,8])/(c[1,7]+c[2,7]+c[3,7]+c[4,7]+c[5,7]+c[6,7]+c[7,7]+c[8,7])
f[6]=(c[1,7]+c[2,7]+c[3,7]+c[4,7]+c[5,7]+c[6,7]+c[7,7]+c[8,7]+c[9,7])/(c[1,6]+c[2,6]+c[3,6]+c[4,6]+c[5,6]+c[6,6]+c[7,6]+c[8,6]+c[9,6])
f[5]=(c[1,6]+c[2,6]+c[3,6]+c[4,6]+c[5,6]+c[6,6]+c[7,6]+c[8,6]+c[9,6]+c[10,6])/(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5]+c[6,5]+c[7,5]+c[8,5]+c[9,5]+c[10,5])
f[4]=(c[1,5]+c[2,5]+c[3,5]+c[4,5]+c[5,5]+c[6,5]+c[7,5]+c[8,5]+c[9,5]+c[10,5]+c[11,5])/(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4]+c[7,4]+c[8,4]+c[9,4]+c[10,4]+c[11,4])
f[3]=(c[1,4]+c[2,4]+c[3,4]+c[4,4]+c[5,4]+c[6,4]+c[7,4]+c[8,4]+c[9,4]+c[10,4]+c[11,4]+c[12,4])/(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3]+c[8,3]+c[9,3]+c[10,3]+c[11,3]+c[12,3])
f[2]=(c[1,3]+c[2,3]+c[3,3]+c[4,3]+c[5,3]+c[6,3]+c[7,3]+c[8,3]+c[9,3]+c[10,3]+c[11,3]+c[12,3]+c[13,3])/(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2]+c[9,2]+c[10,2]+c[11,2]+c[12,2]+c[13,2])
f[1]=(c[1,2]+c[2,2]+c[3,2]+c[4,2]+c[5,2]+c[6,2]+c[7,2]+c[8,2]+c[9,2]+c[10,2]+c[11,2]+c[12,2]+c[13,2]+c[14,2])/(c[1,1]+c[2,1]+c[3,1]+c[4,1]+c[5,1]+c[6,1]+c[7,1]+c[8,1]+c[9,1]+c[10,1]+c[11,1]+c[12,1]+c[13,1]+c[14,1])


for (i in 1:15){
for (j in 1:15){
if (i+j<17) {ch[i,j]<-c[i,j]} else {ch[i,j]<-ch[i,j-1]*f[j-1]}
}
}


for (i in 1:13){
a=0;for (j in 1:(15-i)){
a<-a+c[j,i]*(c[j,i+1]/c[j,i]-f[i])^2
}
h[i]=a/(14-i)
}
h[14]=min(h[13],h[12])

for (i in 2:15){
a=0;for (k in (16-i):14){
a1=0;for (j in 1:(15-k)){
a1<-a1+c[j,k]}
a<-a+((h[k]/(f[k])^2))*(1/ch[i,k]+1/a1)
}
R[i]=sqrt(a*(ch[i,15])^2)
}

V=0;for (i in 2:15){
a1=R[i]^2
b1=0
if (i<15) {for (j in (1+i):15){b1<-b1+ch[j,15]};a2=b1*ch[i,15]} else {a2=ch[i,15]}
b1=0;for (k in (16-i):4){
b2=0;for (n in 1:(15-k)){b2<-b2+c[n,k]}
b1<-b1+(2*h[k]/(f[k]^2))/b2
}
V=V+a1+a2*b1
}

EV=( c[15,1]*f[1]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]*f[10]*f[11]*f[12]*f[13]*f[14]
	+ c[14,2]*f[2]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]*f[10]*f[11]*f[12]*f[13]*f[14]
	+ c[13,3]*f[3]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]*f[10]*f[11]*f[12]*f[13]*f[14]   
	+ c[12,4]*f[4]*f[5]*f[6]*f[7]*f[8]*f[9]*f[10]*f[11]*f[12]*f[13]*f[14]
	+ c[11,5]*f[5]*f[6]*f[7]*f[8]*f[9]*f[10]*f[11]*f[12]*f[13]*f[14]
	+ c[10,6]*f[6]*f[7]*f[8]*f[9]*f[10]*f[11]*f[12]*f[13]*f[14]
	+ c[9,7]*f[7]*f[8]*f[9]*f[10]*f[11]*f[12]*f[13]*f[14]
	+ c[8,8]*f[8]*f[9]*f[10]*f[11]*f[12]*f[13]*f[14]
	+ c[7,9]*f[9]*f[10]*f[11]*f[12]*f[13]*f[14]
	+ c[6,10]*f[10]*f[11]*f[12]*f[13]*f[14]
	+ c[5,11]*f[11]*f[12]*f[13]*f[14]
	+ c[4,12]*f[12]*f[13]*f[14]
	+ c[3,13]*f[13]*f[14]
	+ c[2,14]*f[14]
	+ c[1,15]
	-ac1 )
SD=sqrt(V)

sigma=sqrt(log((1+V/(EV^2))))
mu=log(EV)-0.5*sigma^2
p=plnorm(ac2,mu,sigma)
if (p>0.99) {z.15.1 <- z.15.1+1}
if (p>0.95) {z.15.5 <- z.15.5+1}
if (p>0.90) {z.15.10 <- z.15.10+1}
if (p>0.80) {z.15.20 <- z.15.20+1}
if (p>0.70) {z.15.30 <- z.15.30+1}
if (p>0.50) {z.15.50 <- z.15.50+1}
if (p>0.30) {z.15.70 <- z.15.70+1}
if (p>0.20) {z.15.80 <- z.15.80+1}
if (p>0.10) {z.15.90 <- z.15.90+1}
if (p>0.05) {z.15.95 <- z.15.95+1}
if (p>0.01) {z.15.99 <- z.15.99+1}
zz.15=zz.15+EV
}


z.15.1/1000
z.15.5/1000
z.15.10/1000
z.15.20/1000
z.15.30/1000
z.15.50/1000
z.15.70/1000
z.15.80/1000
z.15.90/1000
z.15.95/1000
z.15.99/1000

zz.15

# --------------------------------------------------------- 


## Mack Bootstrap Simulations













