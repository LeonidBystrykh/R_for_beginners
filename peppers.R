## data from Rdatasets
#here we exersise with basic plots using base graphics
#install.packages("AER")
library(AER) #it contains data you see in Rdatasets

#load data 
data("PepperPrice", package = "AER")

head(PepperPrice, 5)
class(PepperPrice)
summary(PepperPrice)
#original plot
plot(PepperPrice, plot.type = "single", col = c("blue","red"), 
     main="Prices for pepper in years")
grid()
#for other purposes it is better to create a data frame
#convert to df and extract time scale
df<-data.frame(PepperPrice)
#recover time scale
time<-data.frame(PepperPrice, date=time(PepperPrice))[,3]
time
typeof(time)
class(time)
#df$time<-time
#descriptive statistics
#HISTOGRAMS
hist(as.matrix(df), main="Prices for both peppers",
     xlab="Prices", breaks=30)
# Change the plot region color
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "#f7f7f7") # Color
grid()
#repeat hist in add=T mode, so the background remains
hist(as.matrix(df), main="Prices for both peppers",
     xlab="Prices", breaks=30, add=T)

#do it separately
hist(df$white, breaks=20, main="Price for White pepper",
     xlab="White", col=rgb(1,0.8,0.8,0.5))
hist(df$black, breaks=20, main="Prices", add=T, col=rgb(0.8,0.8,1,0.5))

#DENSITY PLOTS
#for all
densityPlot(as.matrix(df))
#separately 
densityPlot(df$black, col="blue")
d2<-densityPlot(df$white)
#restore 1st plot for black, add 2nd using lines
lines(d2$x, d2$y, col="red", lwd=2)
#the same in limma
library(limma)
plotDensities(df, legend ="topright", col=c("blue","red"))

#BOX plots and JITTER
boxplot(df, main="Prices", col=c("blue","red"))
stripchart(df, vert = T, 
           method = "jitter", 
           add = TRUE, 
           pch = 19, col=c("red","blue"), cex=0.2)

#VIOLIN PLOT
library(vioplot)

vioplot(df, col=c("lightblue","tomato"), main="Prices for white and black")
grid()
histoplot(df, col=c(rgb(0.8,0.8,1),rgb(1,0.8,0.8)))

#DESCRIPTIVE STATS
summary(df)
m1<-mean(df[,1])
m2<-mean(df[,2])
abline(h=m1, col="blue")
abline(h=m2, col="red")
s1<-sd(df[,1])
s2<-sd(df[,2])
lines(c(1,1),c(m1-s1,m1+s1),lwd=5,col="blue")
lines(c(2,2),c(m2-s2,m2+s2),lwd=5,col="red")

#plot in time
plot(as.numeric(time), df$white, "l", col="blue", ylab="Price")
lines(as.numeric(time), df$black, "l", col="red")
legend("topleft", colnames(df), col=c("blue","red"), lty=1)
title("Prices for black and white pepper")

#usless for this case
dotchart(as.matrix(df), labels="")
barplot(df) #wrong
barplot(colMeans(df), main="Means of Prices for pepper") #right
barplot(df$white[1:50], col="#FF0000")
barplot(df$black[1:50], col="#0000FF",add=T)

#analytical part

#are prices significantly different for black and white p.?
test<-t.test(df$black, df$white, paied=T)
test
test<-wilcox.test(df$black, df$white,paired=T)
test
fc<-test$estimate[2]/test$estimate[1]
fc
#do they correlate?
plot(df$black, df$white, xlab="black", 
     ylab="white", main="Do prices correlate?")
#we can use correlation test
cor.test(df$white, df$black, data=df)

#or we can use lm
fit<-lm(white~0+black, data=df)
fit
summary(fit)
predicted<-predict.lm(fit)
lines(df$black, predicted, "l", col="red", lwd=2, lty=2)

