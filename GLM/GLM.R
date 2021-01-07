#install.packages("titanic")
library(titanic)

#Dataset
data.raw <- titanic_train

#convert missing values to NA
data.raw[data.raw==""] <- NA

#Check for missing values
sapply(data.raw,function(x) sum(is.na(x)))

data <- subset(data.raw,select=c(2,3,5,7,8,9,10,12))

#glm binomial // log(success/failure)=log(odds)=-2.5137
model <- glm(Survived ~ Sex, family = binomial, data = data)
summary(model)
#confint(model)

attach(data)


#Fitting with all the features
M1 <- glm(Survived ~ Sex + Fare + Pclass + Embarked + SibSp + Parch , family = binomial, data = data)
#names(M1)
summary(M1)
plot(M1)

#Significant features
Model <- glm(Survived ~ Sex + Pclass, family = binomial, data = data)
#names(Model)
summary(Model)
#plot(Model$residuals)


#Different families comparison
M3 <- glm(Survived ~ Sex + Pclass , family = gaussian, data = data)
summary(M3)
plot(M3) # This looks good
pchisq(133.25,888)
M4 <- glm(Survived ~ Sex + Pclass ,family = poisson, data = data)
summary(M4)
plot(M4) #This doesn't look good
pchisq(465.2,888)
M5 <- glm(Survived ~ Sex + Pclass ,family = quasi, data = data)
summary(M5)
pchisq(133.25,888)
M6 <- glm(Survived ~ Sex + Pclass ,family = quasibinomial, data = data)
summary(M6)
pchisq(827.2,888)
M7 <- glm(Survived ~ Sex + Pclass ,family = quasipoisson, data = data)
summary(M6)
pchisq(465.2,888)

#Hence, using either gaussian or quasi family we get the smallest deviances and 
# the smallest pchisq, also the AIC for the gaussian is 843.52 and for quasi is not applicable.

