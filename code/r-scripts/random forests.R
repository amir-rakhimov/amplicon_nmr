library(randomForest)
# install.packages("randomForest")
# Airquality ####
# view structure of airquality dataset
str(airquality)
head(airquality)
# find num of rows with missing values
sum(!complete.cases(airquality))

airq.df<-airquality
#replace NAs with column medians
for(i in 1:ncol(airq.df)) {
  airq.df[ , i][is.na(airq.df[ , i])] <- median(airq.df[ , i], na.rm=TRUE)
}

head(airq.df)
head(airquality)

#make this example reproducible
set.seed(1)

# fit the random forest model
rf.model<-randomForest(
  formula = Ozone ~ .,
  data = airq.df
)
rf.model

#find number of trees that produce lowest test MSE
which.min(rf.model$mse)

#find RMSE of best model
sqrt(rf.model$mse[which.min(rf.model$mse)])

#plot the test MSE by number of trees
plot(rf.model)

#plot the test MSE by number of trees
varImpPlot(rf.model)

model_tuned<- tuneRF(
  x=airq.df[,-1],
  y=airq.df$Ozone,
  ntreeTry = 500,
  mtryStart = 4,
  stepFactor = 1.5,
  improve = 0.01,
  trace = FALSE
)

# Use the final model to make predictions
new<-data.frame(Solar.R=150, Wind=8, Temp=70,
                Month=5, Day= 5)
predict(rf.model,newdata = new)

# readingSkills dataset ####
# install.packages("party")
library(party)
print(head(readingSkills))
output.forest<-randomForest(nativeSpeaker~ age+shoeSize+score,
                            data=readingSkills)
print(output.forest)
print(importance(output.forest,type=2))

