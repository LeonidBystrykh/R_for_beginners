#try ggplot2 graphs with peppers data
# useful discussion in 
#https://stackoverflow.com/questions/9531904/plot-multiple-columns-on-the-same-graph-in-r
## data from Rdatasets
#install.packages("AER")
library(AER) #it contains data you see in Rdatasets
library(ggplot2) #major lib for ggplot
library(reshape2) #for melt
library(ggpubr) #for compare means and boxplot 
#library(devtools) #to install lib below 
#install_github("kassambara/easyGgplot2")
library(easyGgplot2) #another way of the boxplot
#load data 
data("PepperPrice", package = "AER")

head(PepperPrice, 5)
class(PepperPrice)
summary(PepperPrice)

#for other purposes it is better to create a data frame
#convert to df and extract time scale
df<-data.frame(PepperPrice)
#recover time scale
time<-data.frame(PepperPrice, date=time(PepperPrice))[,3]
df$time<-time

#you can draw two lines from df directly
ggplot(data=df)+
  geom_line(mapping=aes(time, black, col="black"))+
  geom_line(mapping=aes(time,white, col="white"))+
  scale_color_manual(values=c("blue", "red"))+
    ggtitle("Prices for pepper") +
    xlab("Time") +
    ylab("Price") 

#however, the commont practice is to melt table and make it "long"    
#melt data before any use
long<-melt(df, id.vars="time")
class(long)
#make proper column names
colnames(long)
#[1] "time"     "variable" "value" 
colnames(long)<-c("time","pepper","price")
#Time plot
ggplot(data=long, aes(time, price, col=pepper)) + 
   geom_line() + 
   scale_color_manual(values=c("blue", "red"))+
   ggtitle("Prices for pepper") +
   xlab("Time") +
   ylab("Price") 

#histograms
ggplot(long, aes(price, col=pepper))+
  geom_histogram(fill="#DDDDDD")+
  scale_color_manual(values=c("blue", "red"))+
  ggtitle("Prices for pepper") 

#density plot
ggplot(long, aes(price, col=pepper))+
  geom_density()+
  scale_color_manual(values=c("blue", "red"))+
  ggtitle("Prices for pepper") 

#boxplot
ggplot(long, aes(pepper, price, col=pepper) ) + 
  geom_boxplot(fill=c("#DDDDFF","#FFDDDD")) +
  scale_color_manual(values=c("blue", "red"))+
  geom_jitter(size=0.1) +
  ggtitle("Prices for pepper") +
  stat_compare_means(method="wilcox")
  ylab("Price") 

#another option
#library(ggpubr)
ggboxplot(long, x="pepper", y = "price",
          color = "pepper", palette = c("blue","red"),
          fill="#DDEEDD"
          ) + #it cannot make two colors
  ggtitle("Prices for pepper")+
  #stat_compare_means()
  theme_gray()

#and another one

ggplot2.boxplot(data=long, xName='pepper',yName='price',groupName='pepper',
                groupColors = c("blue", "red"))+
  geom_jitter(size=0.1, colour="brown")+ #only one color!
  ggtitle("Prices for pepper")+
  stat_compare_means(method="t.test")
  
#there was one more: qplot. Now it is deprecated

#violin plot
ggplot(long, aes(pepper, price, col=pepper)) + 
  scale_color_manual(values=c("blue", "red"))+
  geom_violin(draw_quantiles = c(0,0.5,1)) +
  geom_jitter(size=0.1) +
  ggtitle("Prices for pepper") +
  xlab("Time") +
  stat_compare_means(method="t.test") +
  ylab("Price") 

#alternative from ggpubr
ggviolin(long, x="pepper", y="price", fill="pepper",
         palette = c("#CCCCFF","#FFCCCC"), add="boxplot",
         add.params = list(fill="white"))+
  geom_jitter(size=0.1)+
  stat_compare_means(method="t.test")+
  theme_linedraw()
         
#one more from easyGgplot2
#let's try to add mean and stDev values on the plot

sd1<-sd(df$black)
sd2<-sd(df$white)
m1<-mean(df$black)
m2<-mean(df$white)

ggplot2.violinplot(data=long, xName='pepper',yName='price',
                   groupName='pepper', meanPointShape = 19, 
                   meanPointSize = 5,
                   addMean=T,
                   groupColors=c("#CCCCFF","#FFCCCC"))+
#  geom_jitter(size=0.1)+
  ggtitle("Two guitars")+
  geom_segment(aes(x=1, y=m1-sd1, xend=1,yend=m1+sd1))+
  geom_segment(aes(x=2, y=m2-sd2, xend=2,yend=m2+sd2))+
  stat_compare_means(method="t.test")

#try correlation plot. Note, we are back to df table!
p1<-ggplot(df, aes(x=black, y=white))+
  geom_point()+
  stat_smooth(method="lm", se=T, formula= y ~ x)+
  stat_cor(method = "pearson")
p1
summary(p1)

fit<-lm(white~ black, data=df) 
summary(fit)

