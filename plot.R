#try line and scatter

###one line, one plot with plot()
Time<-c(1,2,3,4,5,6,7)
Effect<-c(2,3,5,3,1,0,0)
plot(Time, Effect)#, main="Effect in time", type="b")
plot(Time, Effect, main="Effect in time", type="b", 
     col="red", lwd=2, lty=2)

#Note the next plot overrides the previous one

#One plot, two lines
Eff2<-runif(length(Effect),min(Effect),max(Effect))
plot(Time,Effect, col="red", "b", ylim=c(0,9), ylab="Effect")
#but you can use option
par(new=T)
plot(Time,Eff2, col="blue", "b", main="Two effects",
     ylim=c(0,9), ylab="Effect")
#note: y scale is ruined, use ylim or below

#Better solution:
plot(Time, Time, main="Time and Effect",pch=19, 
     col=rgb(0,0,1), xlab="", ylab="Time and Effect")
points(Time, Effect, col="red", pch=19, cex=2)
lines(Time, Eff2, col="orange", lwd="2", lty=2)
text(4,6, "dots and lines")

#line or curve through
scatter.smooth(Time, Effect, degree=0, col="blue",
               main="smooth curve through the points")


#can I draw a line through?
fit<-lm(Effect~Time)
fit
predicted<-predict.lm(fit)
lines(Time, predicted, col="red")


#draw curve by function
x<-1:5
eq<-function(x){x^2}
curve(eq(x), to=5, main="Curve by function")
points(x, x^2)


#what if i want multiple lines from data frame?
df<-data.frame(Time, Effect, Eff2)
plot(df) #it is completely different

#this will work
matplot(1:7, df, "b", lty=1, xlab="Time", ylab="Effect")
title(main="Multiple lines" )
legend("topleft", legend=colnames(df), col=1:3, fill=1:3)



#Appendix
grid(col="lightgray")
par(bg="#FFFFDD")
#draw step line through data (this is cumulative)
plot.stepfun(Time,Effect)
