#####
#Load packages
#####
#load libraries 
library(R.matlab)
library(caret)
library(dplyr)
library(ggplot2)
library(ggcorrplot)
library(GGally)
library(rpart)
library(DMwR)
library(ROSE)
library(car)
library(caret)
library(tidyr)
library(rattle)
library(pROC)
library(purrr)
library(MASS)
library(corrplot)
library(ROCR)
library(ggthemr)
ggthemr("light")
library(scales)

#Function to visualize CM
draw_confusion_matrix <- function(cm) {
  
  total <- sum(cm$table)
  res <- as.numeric(cm$table)
  
  # Generate color gradients. Palettes come from RColorBrewer.
  greenPalette <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
  redPalette <- c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
  getColor <- function (greenOrRed = "green", amount = 0) {
    if (amount == 0)
      return("#FFFFFF")
    palette <- greenPalette
    if (greenOrRed == "red")
      palette <- redPalette
    colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
  }
  
  # set the basic layout
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  classes = colnames(cm$table)
  rect(150, 430, 240, 370, col=getColor("green", res[1]))
  text(195, 435, classes[1], cex=1.2)
  rect(250, 430, 340, 370, col=getColor("red", res[3]))
  text(295, 435, classes[2], cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=getColor("red", res[2]))
  rect(250, 305, 340, 365, col=getColor("green", res[4]))
  text(140, 400, classes[1], cex=1.2, srt=90)
  text(140, 335, classes[2], cex=1.2, srt=90)
  
  # add in the cm results
  text(195, 400, res[1], cex=1.6, font=2, col='black')
  text(195, 335, res[2], cex=1.6, font=2, col='black')
  text(295, 400, res[3], cex=1.6, font=2, col='black')
  text(295, 335, res[4], cex=1.6, font=2, col='black')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}

#####
#Preprocessing and understanding the dataset
#####
#Read in matlab file 
matdata2 <- readMat("OHRdatafe_2.mat")
matdata3 <- readMat("OHRdatafe_3.mat")
matdata4 <- readMat("OHRdatafe_4.mat")
matdata5 <- readMat("OHRdatafe_5.mat")
matdata6 <- readMat("OHRdatafe_6.mat")
matdata7 <- readMat("OHRdatafe_7.mat")
matdata8 <- readMat("OHRdatafe_8.mat")
#Convert to dataframes
data2 <- data.frame(matdata2$DATA)
data3 <- data.frame(matdata3$DATA)
data4 <- data.frame(matdata4$DATA)
data5 <- data.frame(matdata5$DATA)
data6 <- data.frame(matdata6$DATA)
data7 <- data.frame(matdata7$DATA)
data8 <- data.frame(matdata8$DATA)
#Combine dataframes
data <- rbind(data2, data3, data4,data5,data6, data7, data8)
#Rename columns
data <- data %>%
  rename(Gna = X1, Gto = X2, Pca = X3,Gkr= X4, Gks = X5,
         Gk1 = X6,Gncx = X7, Gpca = X8, Cycle=X9, APD90=X10,Class = X11)

#How many alternans are there? 
summary(as.factor(data$Class))

#Change to factors, rename levels and make EADs the positive class
data$Class[data$Class<0.5] <-10
data$Class <- as.factor(data$Class) 
levels(data$Class) <- c("EAD","EAD","EAD", "EAD", "Non_EAD")
levels(data$Class)
#####
#Look at EAD to Non_EAD proportion and structure of data
#####
table(data$Class)
prop.table(table(data$Class))*100
str(data)

#Plot class imbalance
dat <- data.frame("Class" = c("EAD", "Non EAD"),
                  "val" = c(2332, 209917),
                  "perc" = c(1.10, 98.90))

classimbplot <- ggplot(data = dat, aes(x=Class, y=val, fill=Class))+
  geom_bar(stat = 'identity')+
  geom_text(aes(label = paste0(round(perc), "%")),size=5, hjust=-0.1)+
  geom_text(aes(label = val,size=5, hjust=-0.1, vjust=-3))+
  xlab("EAD Class")+
  scale_y_continuous(limits=c(0, 250000))+
  ylab("Number of observations")+
  ggtitle("Number of observations in each class")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
  coord_flip() 

#####
#Splitting the data
#####
#Remove cycle from the data
classdata <- data[,-9]
#Set seed
set.seed(1235)
#Splitting the data
inTrain <- createDataPartition(y = classdata$Class, p = 0.75, list = FALSE)
#training set and its dimensions
training <- classdata[inTrain,]
dim(training)
#testing set and its dimensions
testing <- classdata[-inTrain,]
dim(testing)
######
#Exploring training set
#####
#Exploratory analysis will only be done on the training set rather than full dataset
#Summary of training set
summary(training$Class)
#summary for each AP class
training %>% split(.$Class) %>% map(summary)
#Summary of test set
summary(testing$Class)
summary(data)

#Scatterplots
#For visibility purposes, reoder factors so that EAD, which is the minority class, will
#be bigger and hence more visible. Do not do these on the training data, create 
#Duplicate only used for the next six scatterplots
training2 <- training
training2$Class <- factor(training2$Class, levels = c("Non_EAD", "EAD"))
levels(training2$Class) #Level order changed
#Gkr vs Pca
ggplot(training2, aes(x=Gkr, y=Pca)) +
geom_point(aes(colour=Class, size=Class,alpha=Class))+
ggtitle("Pca vs Gkr") +
scale_x_continuous(limits=c(0,.35))+
scale_y_continuous(limits=c(0.4,1.9))

#Gkr vs Gks
ggplot(training2, aes(x=Gkr, y=Gks)) +
geom_point(aes(colour=Class, size=Class,alpha=Class))+
scale_x_continuous(limits=c(0,.35))+
ggtitle("Gks vs Gkr")

#Pca vs Gncx
ggplot(training2, aes(x=Pca, y=Gncx)) +
geom_point(aes(colour=Class, size=Class,alpha=Class))+
scale_x_continuous(limits=c(0.4,1.9))+
ggtitle("Gncx vs Pca")

#Gkr vs Gncx
ggplot(training2, aes(x=Gkr, y=Gncx)) +
geom_point(aes(colour=Class, size=Class,alpha=Class))+
scale_x_continuous(limits=c(0,0.35))+
ggtitle("Gncx vs Gkr")


#Scatterplot matrix
matrixplot <- ggpairs(data = training, ggplot2::aes(colour=Class),
                      upper = list(continuous = wrap("cor", cex=2.6)),
                      lower = list(continuous = wrap("points", cex=0.001), 
                                   combo = wrap("box", cex=0.3), 
                                   discrete = wrap("facetbar", cex=0.1)), 
                      diag = list(continuous = wrap("densityDiag", alpha = 0.5, cex=0.1, limits=c(0,2.5)))) +
  ggtitle("Scatter plot matrix by EAD class")


#####
#Remove APD90 before training models for easier coding (will not be used as predictor variable)
#####
training <- training[,-9]
testing <- testing[,-9]

#####
###############################################Logistic Regression 
#####
##Creating sampled datasets

summary(training$Class)
#1897 AP and 9 EAD classes
# First try oversampling, increases the EAD records until matches AP  1749  157438 
# Here N= 157438*2
set.seed(1235)
over_training <- ovun.sample(Class ~ ., data = training, method="over", p=0.5)$data
#We can see that now the sizes are equal for both classes
summary(over_training$Class)
prop.table(table(over_training$Class))

#Undersampling will decrease AP to match EAD so data lost, N = 1749*2
set.seed(1235)
under_training <- ovun.sample(Class ~ ., data = training, method="under", N=3498)$data
#We can see that now the sizes are equal for both classes
summary(under_training$Class)
prop.table(table(under_training$Class))

# ROSE Sampling, generates artificial data instead of duplicate data.
rose_training <- ROSE(Class ~ ., data = training,  seed=1235)$data
summary(rose_training$Class)
prop.table(table(rose_training$Class))

# SMOTE(Synthetic Minority Over-sampling Technique) Sampling
# formula : how the dependent variable act based on other independent variables
# perc.over : increases the size of EAD class
# perc.under : decreases the size of non EAD class
set.seed(1235)
smote_training <-DMwR::SMOTE(Class ~ ., training, perc.over = 100, perc.under=200)
summary(smote_training$Class)
prop.table(table(smote_training$Class))

#When sampling through the ROSE package, it changed the factor levels
#imbalanced and SMOTE levels 
levels(training$Class)[1]
levels(smote_training$Class)[1]
#vs under, over and ROSE
levels(over_training$Class)[1]
levels(under_training$Class)[1]
levels(rose_training$Class)[1]

#Relevel under, over and ROSE to match original
over_training$Class <- factor(over_training$Class, levels = c("EAD", "Non_EAD"))
under_training$Class <- factor(under_training$Class, levels = c("EAD", "Non_EAD"))
rose_training$Class <- factor(rose_training$Class, levels = c("EAD", "Non_EAD"))

#####
#Logistic regression models
#####
#Logreg on imbalanced data
glmmodel <- glm(Class ~ ., family = binomial(link = "logit"), data= training)

#Model summary
summary(glmmodel)
#In R, EAD is coded in R as level one and Non_EAD is coded in R as level 2

#Oversampling
overglmmodel <- glm(Class ~ ., family=binomial(link = "logit"), data=over_training, control = list(maxit=50))
summary(overglmmodel)
anova(overglmmodel,test="Chisq")

#Undersampling
underglmmodel <- glm(Class ~ ., family=binomial, data=under_training, control = list(maxit=50))
summary(underglmmodel)
anova(underglmmodel,test="Chisq")

#Rose
roseglmmodel <- glm(Class ~ ., family=binomial, data=rose_training, control = list(maxit=50))
summary(roseglmmodel)
anova(roseglmmodel,test="Chisq")

#Smote
smoteglmmodel <- glm(Class ~ ., family=binomial, data=smote_training, control = list(maxit=50))
summary(smoteglmmodel)
anova(smoteglmmodel,test="Chisq")

#Predictions on training data 
glmpred=predict(glmmodel,newdata=training, type = "response")
over_glmpred=predict(overglmmodel,newdata=over_training, type = "response")
under_glmpred=predict(underglmmodel,newdata=under_training, type = "response")
rose_glmpred=predict(roseglmmodel,newdata=rose_training, type = "response")
smote_glmpred=predict(smoteglmmodel,newdata=smote_training, type = "response")


######Finding optimal thresholds for each model
#Glm(imbalanced)
glmpredictions <- prediction(glmpred, training$Class)
glmsens <- data.frame(x=unlist(performance(glmpredictions, "sens")@x.values), 
                      y=unlist(performance(glmpredictions, "sens")@y.values))
glmspec <- data.frame(x=unlist(performance(glmpredictions, "spec")@x.values), 
                      y=unlist(performance(glmpredictions, "spec")@y.values))
glmthresplot <- glmsens %>% ggplot(aes(x,y)) + 
  geom_line() + 
  geom_line(data=glmspec, aes(x,y,col="red")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Specificity")) +
  labs(x='Threshold', y="Sensitivity") +
  theme(axis.title.y.right = element_text(colour = "red"), legend.position="none") +
  scale_x_continuous(breaks=seq(0, 1, by = 0.05))+
  ggtitle("Glm thresholds - imbalanced data")+
  theme(axis.text.x=element_text(angle=45,hjust=1))
#Optimal threshold is 0.97 for best Sensitivity/Specificity



#Oversampling
overglmpredictions <- prediction(over_glmpred, over_training$Class)
overglmsens <- data.frame(x=unlist(performance(overglmpredictions, "sens")@x.values), 
                          y=unlist(performance(overglmpredictions, "sens")@y.values))
overglmspec <- data.frame(x=unlist(performance(overglmpredictions, "spec")@x.values), 
                          y=unlist(performance(overglmpredictions, "spec")@y.values))
overthresplot <- overglmsens %>% ggplot(aes(x,y)) + 
  geom_line() + 
  geom_line(data=overglmspec, aes(x,y,col="red")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Specificity")) +
  labs(x=' Threshold', y="Sensitivity") +
  theme(axis.title.y.right = element_text(colour = "red"), legend.position="none") +
  scale_x_continuous(breaks=seq(0, 1, by = 0.05))+
  ggtitle("Glm thresholds - oversampled data")+
  theme(axis.text.x=element_text(angle=45,hjust=1))
#Optimal threshold is 0.25 for best Sensitivity/Specificity

#Undersampling
underglmpredictions <- prediction(under_glmpred, under_training$Class)
underglmsens <- data.frame(x=unlist(performance(underglmpredictions, "sens")@x.values), 
                           y=unlist(performance(underglmpredictions, "sens")@y.values))
underglmspec <- data.frame(x=unlist(performance(underglmpredictions, "spec")@x.values), 
                           y=unlist(performance(underglmpredictions, "spec")@y.values))
underthresplot <- underglmsens %>% ggplot(aes(x,y)) + 
  geom_line() + 
  geom_line(data=underglmspec, aes(x,y,col="red")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Specificity")) +
  labs(x=' Threshold', y="Sensitivity") +
  theme(axis.title.y.right = element_text(colour = "red"), legend.position="none") +
  scale_x_continuous(breaks=seq(0, 1, by = 0.05))+
  ggtitle("Glm thresholds - undersampled data")+
  theme(axis.text.x=element_text(angle=45,hjust=1))
#Optimal threshold is 0.25 for best Sensitivity/Specificity

#ROSE Sampling
roseglmpredictions <- prediction(rose_glmpred, rose_training$Class)
roseglmsens <- data.frame(x=unlist(performance(roseglmpredictions, "sens")@x.values), 
                          y=unlist(performance(roseglmpredictions, "sens")@y.values))
roseglmspec <- data.frame(x=unlist(performance(roseglmpredictions, "spec")@x.values), 
                          y=unlist(performance(roseglmpredictions, "spec")@y.values))
rosethresplot <- roseglmsens %>% ggplot(aes(x,y)) + 
  geom_line() + 
  geom_line(data=roseglmspec, aes(x,y,col="red")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Specificity")) +
  labs(x=' Threshold', y="Sensitivity") +
  theme(axis.title.y.right = element_text(colour = "red"), legend.position="none") +
  scale_x_continuous(breaks=seq(0, 1, by = 0.05))+
  ggtitle("Glm thresholds - ROSE sampled data")+
  theme(axis.text.x=element_text(angle=45,hjust=1))
#Optimal threshold is 0.25 for best Sensitivity/Specificity

#Glm(SMOTE)
smoteglmpredictions <- prediction(smote_glmpred, smote_training$Class)
smoteglmsens <- data.frame(x=unlist(performance(smoteglmpredictions, "sens")@x.values), 
                           y=unlist(performance(smoteglmpredictions, "sens")@y.values))
smoteglmspec <- data.frame(x=unlist(performance(smoteglmpredictions, "spec")@x.values), 
                           y=unlist(performance(smoteglmpredictions, "spec")@y.values))
smoteglmthresplot <- smoteglmsens %>% ggplot(aes(x,y)) + 
  geom_line() + 
  geom_line(data=smoteglmspec, aes(x,y,col="red")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Specificity")) +
  labs(x='Threshold', y="Sensitivity") +
  theme(axis.title.y.right = element_text(colour = "red"), legend.position="none") +
  scale_x_continuous(breaks=seq(0, 1, by = 0.05))+
  ggtitle("Glm thresholds - SMOTE sampled data")+
  theme(axis.text.x=element_text(angle=45,hjust=1))
#Optimal threshold is 0.3 for best Sensitivity/Specificity
#Grid arrange of threshold plots
gridExtra::grid.arrange(glmthresplot,overthresplot,underthresplot,rosethresplot,smoteglmthresplot, 
                        ncol=2, heights=c(2,2,2))

#Predict on test data and produce CM

#Using model trained on imbalanced data
aa=predict(glmmodel, newdata = testing, type="response")
ag=ifelse(aa < 0.97,"EAD","Non_EAD")
cmtest_glm<-confusionMatrix(as.factor(ag), testing$Class, mode="everything") 

#Using model trained on over sampled data
overtestdataframe <- data.frame(actual = testing$Class,
                                predict(overglmmodel,newdata=testing, type = "response"))
overtestdataframe$predict <- as.factor(ifelse(overtestdataframe[,2] < 0.5, "EAD", "Non_EAD"))
cmtest_glmover <- confusionMatrix(overtestdataframe$predict,overtestdataframe$actual, mode="everything")

#Using model trained on undersampled data
undertestdataframe <- data.frame(actual = testing$Class,
                                 predict(underglmmodel,newdata=testing, type = "response"))
undertestdataframe$predict <- as.factor(ifelse(undertestdataframe[,2] < 0.5, "EAD", "Non_EAD"))
cmtest_glmunder <- confusionMatrix(undertestdataframe$predict,undertestdataframe$actual, mode="everything")

#Using model trained on ROSE data
rosetestdataframe <- data.frame(actual = testing$Class,
                                predict(roseglmmodel,newdata=testing, type = "response"))
rosetestdataframe$predict <- as.factor(ifelse(rosetestdataframe[,2] < 0.25, "EAD", "Non_EAD"))
cmtest_glmrose <- confusionMatrix(rosetestdataframe$predict,rosetestdataframe$actual, mode="everything")

#Using model trained on SMOTE data
af=predict(smoteglmmodel, newdata = testing, type="response")
al=ifelse(af < 0.5,"EAD","Non_EAD")
cmtest_glmsmote <- confusionMatrix(as.factor(al), testing$Class, mode="everything")

#####
#Produce CM test plots
#####
draw_confusion_matrix(cmtest_glm) 
draw_confusion_matrix(cmtest_glmover) 
draw_confusion_matrix(cmtest_glmunder) 
draw_confusion_matrix(cmtest_glmrose) 
draw_confusion_matrix(cmtest_glmsmote) 

#####
#Create plots for confusion matrices - testing
#####
#Create list of models
models <- list(glm = glmmodel,
               glmover = overglmmodel,
               glmunder = underglmmodel,
               glmrose = roseglmmodel,
               glmsmote = smoteglmmodel)


comparison <- data.frame(model = names(models))


for (name in names(models)) {
  model <- get(paste0("cmtest_", name))
  comparison[comparison$model == name, "Sensitivity"] <- model$byClass[["Sensitivity"]]
  comparison[comparison$model == name, "Precision"] <- model$byClass[["Precision"]]
  comparison[comparison$model == name, "Neg Pred Value"] <- model$byClass[["Neg Pred Value"]]
  comparison[comparison$model == name, "F1"] <- model$byClass[["F1"]]
  comparison[comparison$model == name, "Balanced Accuracy"] <- model$byClass[["Balanced Accuracy"]]
  comparison[comparison$model == name, "Specificity"] <- model$byClass[["Specificity"]]
}


comparison %>%
  gather(x, y, Sensitivity: Specificity) %>%
  ggplot(aes(x = x, y = y, color = model)) +
  labs(x="Measure of Performance", y="Value") + 
  ggtitle("GLM performance on test set - sampling comparison") + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 3) + 
  scale_y_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1), limits=c(0,1.05))

######
##############################Decision Trees
#####
#Decision trees using caret package
#Different balancing techniques are compared and the cp is used to post-prune trees
#####
#Setting up control function  and weights
#####
library(MLeval)
library(MLmetrics)
#AUPRC function
library(pROC)

ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 5,
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE)

##Class weights for Class calculated
#number of observations in each class
training %>% count(Class)
#weights argument 
cwt <- ifelse(training$Class == "EAD",
              (1/table(training$Class)[1]) * 0.5,
              (1/table(training$Class)[2]) * 0.5)

#####
##Building model on imbalanced and balanced data
#####
set.seed(1235)
dectreesoriginalfit <- train(Class ~ ., data = training,
                             method = "rpart",
                             metric = "ROC",
                             tuneLength = 30,
                             trControl = ctrl)
summary(dectreesoriginalfit)
plot(dectreesoriginalfit, main = "Model AUC vs Complexity Parameter -Decision Trees without sampling")
#Building weighted Class model
set.seed(1235)
dectreesweightedfit <- train(Class ~ ., data = training,
                             method = "rpart",
                             weights = cwt,
                             metric = "ROC",
                             tuneLength = 30,
                             trControl = ctrl)

#Building  oversampling Class model
#Set seed
set.seed(1235)
#Add sampling method
ctrl$sampling = "up"
#Build model
dectreesoverfit <- train(Class ~ ., data = training,
                         method = "rpart",
                         metric = "ROC",
                         tuneLength = 30,
                         trControl = ctrl)


#Building  undersampling Class model
#Set seed
set.seed(1235)
#Change sampling method
ctrl$sampling <- "down"
#Build model
dectreesunderfit <- train(Class ~ ., data = training,
                          method = "rpart",
                          metric = "ROC",
                          tuneLength = 30,
                          trControl = ctrl)



#Building  rose sampling Class model
#Set seed
set.seed(1235)
#Change sampling method

ctrl$sampling <- "rose"
#Build model
dectreesrosefit <- train(Class ~ ., data = training,
                         method = "rpart",
                         metric = "ROC",
                         tuneLength = 30,
                         trControl = ctrl)


#Building  smote sampling Class model
#Set seed
set.seed(1235)
#Change sampling method
ctrl$sampling <- "smote"
#Build model
dectreessmotefit <- train(Class ~ ., data = training,
                          method = "rpart",
                          metric = "ROC",
                          tuneLength = 30,
                          trControl = ctrl)




#####
#Confusion matrices for test sets
#####
#Original 
dectreesoriginal.predict <- predict(dectreesoriginalfit, testing)
cmtest_dectrees <- confusionMatrix(dectreesoriginal.predict, testing$Class, mode="everything")
#Weights
dectreesweighted.predict <- predict(dectreesweightedfit, testing)
cmtest_dectrees_weight <- confusionMatrix(dectreesweighted.predict, testing$Class, mode="everything")
#oversampling
dectreesover.predict <- predict(dectreesoverfit, testing)
cmtest_dectrees_over <- confusionMatrix(dectreesover.predict, testing$Class, mode="everything")
#undersampling
dectreesunder.predict <- predict(dectreesunderfit, testing)
cmtest_dectrees_under <- confusionMatrix(dectreesunder.predict, testing$Class, mode="everything")
#Rose sampling
dectreesrose.predict <- predict(dectreesrosefit, testing)
cmtest_dectrees_rose <- confusionMatrix(dectreesrose.predict, testing$Class, mode="everything")
#Smote sampling
dectreessmote.predict <- predict(dectreessmotefit, testing)
cmtest_dectrees_smote <- confusionMatrix(dectreessmote.predict, testing$Class, mode="everything")


#####
#Produce CM test plots for dec trees
#####
draw_confusion_matrix(cmtest_dectrees) 
draw_confusion_matrix(cmtest_dectrees_weight)
draw_confusion_matrix(cmtest_dectrees_over) 
draw_confusion_matrix(cmtest_dectrees_under) 
draw_confusion_matrix(cmtest_dectrees_rose) 
draw_confusion_matrix(cmtest_dectrees_smote) 

#####
#Plot of confusion matrices for comparison - test set
#####
#Create list of models
models <- list(dectrees = dectreesoriginalfit,
               dectrees_weight = dectreesweightedfit,
               dectrees_over = dectreesoverfit,
               dectrees_under = dectreesunderfit,
               dectrees_rose = dectreesrosefit,
               dectrees_smote = dectreessmotefit)

#Comparing overall model performance
comparison <- data.frame(model = names(models))

for (name in names(models)) {
  model <- get(paste0("cmtest_", name))
  comparison[comparison$model == name, "Sensitivity"] <- model$byClass[["Sensitivity"]]
  comparison[comparison$model == name, "Pos Pred Value"] <- model$byClass[["Pos Pred Value"]]
  comparison[comparison$model == name, "Neg Pred Value"] <- model$byClass[["Neg Pred Value"]]
  comparison[comparison$model == name, "F1"] <- model$byClass[["F1"]]
  comparison[comparison$model == name, "Balanced Accuracy"] <- model$byClass[["Balanced Accuracy"]]
  comparison[comparison$model == name, "Specificity"] <- model$byClass[["Specificity"]]
}

#Create CM comparison plot
comparison %>%
  gather(x, y, Sensitivity: Specificity) %>%
  ggplot(aes(x = x, y = y, color = model)) +
  labs(x="Measure of Performance", y="Value") + 
  ggtitle("Decision Trees performance on test set - sampling comparison")+
  geom_jitter(width = 0.2, alpha = 0.5, size = 3)+
  scale_y_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1), limits=c(0,1))

#Not the chosen model so decision tree not generated

#####
############################Random Forests
#####
#Random forests using caret package
#Different balancing techniques are compared and 
#the cp is used to post-prune trees
#####
#Setting up control function  and weights
#####
ctrl <- trainControl(method = "cv",
                     number = 10,
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE)


##Class weights for Class calculated
#number of observations in each class
training %>% count(Class)
#weights argument 
cwt <- ifelse(training$Class == "EAD",
              (1/table(training$Class)[1]) * 0.5,
              (1/table(training$Class)[2]) * 0.5)

#####
##Building model on imbalanced and balanced data
#####
set.seed(1235)
rforestoriginalfit <- train(Class ~ ., data = training,
                            method = "rf",
                            metric = "ROC",
                            tuneLength = 7,
                            trControl = ctrl)
rfplot <- plot(rforestoriginalfit, main="Model AUC vs mtry - Random Forests without sampling",
               xlab="mtry")
#Building weighted Class model
set.seed(1235)
rforestweightedfit <- train(Class ~ ., data = training,
                            method = "rf",
                            weights = cwt,
                            metric = "ROC",
                            tuneLength = 7,
                            trControl = ctrl)

#Building  upsampling Class model
#Set seed
set.seed(1235)
#Add sampling method
ctrl$sampling = "up"
#Build model
rforestupfit <- train(Class ~ ., data = training,
                      method = "rf",
                      metric = "ROC",
                      trControl = ctrl)


#Building  downsampling Class model
#Set seed
set.seed(1235)
#Change sampling method
ctrl$sampling <- "down"
#Build model
rforestdownfit <- train(Class ~ ., data = training,
                        method = "rf",
                        metric = "ROC",
                        tuneLength = 7,
                        trControl = ctrl)



#Building  rose sampling Class model
#Set seed
set.seed(1235)
#Change sampling method
ctrl$sampling <- "rose"
#Build model
rforestrosefit <- train(Class ~ ., data = training,
                        method = "rf",
                        metric = "ROC",
                        tuneLength = 7,
                        trControl = ctrl)


#Building  smote sampling Class model
#Set seed
set.seed(1235)
#Change sampling method
ctrl$sampling <- "smote"
#Build model
rforestsmotefit <- train(Class ~ ., data = training,
                         method = "rf",
                         tuneLength = 7,
                         metric = "ROC",
                         trControl = ctrl)

#####
#Confusion matrices for test sets
#####
#Original 
rforestoriginal.predict <- predict(rforestoriginalfit, testing)
cmtest_rforest <- confusionMatrix(rforestoriginal.predict, testing$Class)
#Weights
rforestweighted.predict <- predict(rforestweightedfit, testing)
cmtest_rforest_weights <- confusionMatrix(rforestweighted.predict, testing$Class)
#Upsampling
rforestup.predict <- predict(rforestupfit, testing)
cmtest_rforest_over <- confusionMatrix(rforestup.predict, testing$Class)
#Downsampling
rforestdown.predict <- predict(rforestdownfit, testing)
cmtest_rforest_under <- confusionMatrix(rforestdown.predict, testing$Class)
#Rose sampling
rforestrose.predict <- predict(rforestrosefit, testing)
cmtest_rforest_rose <- confusionMatrix(rforestrose.predict, testing$Class)
#Smote sampling
rforestsmote.predict <- predict(rforestsmotefit, testing)
cmtest_rforest_smote <- confusionMatrix(rforestsmote.predict, testing$Class)

#####
#Produce CM test plots for dec trees
#####
draw_confusion_matrix(cmtest_rforest) 
draw_confusion_matrix(cmtest_rforest_weights)
draw_confusion_matrix(cmtest_rforest_over) 
draw_confusion_matrix(cmtest_rforest_under) 
draw_confusion_matrix(cmtest_rforest_rose) 
draw_confusion_matrix(cmtest_rforest_smote) 


#####
#Plot of confusion matrices for comparison - test set
#####
#Create list of models
models <- list(original = rforestoriginalfit,
               weights = rforestweightedfit,
               up = rforestupfit,
               down = rforestdownfit,
               rose = rforestrosefit,
               smote = rforestsmotefit)

#Comparing overall model performance
comparison <- data.frame(model = names(models))

for (name in names(models)) {
  model <- get(paste0("cmtest_", name))
  comparison[comparison$model == name, "Sensitivity"] <- model$byClass[["Sensitivity"]]
  comparison[comparison$model == name, "Precision"] <- model$byClass[["Precision"]]
  comparison[comparison$model == name, "Neg Pred Value"] <- model$byClass[["Neg Pred Value"]]
  comparison[comparison$model == name, "Balanced Accuracy"] <- model$byClass[["Balanced Accuracy"]]
  comparison[comparison$model == name, "F1"] <- model$byClass[["F1"]]
  comparison[comparison$model == name, "Specificity"] <- model$byClass[["Specificity"]]
}

#Create CM comparison plot
comparison %>%
  gather(x, y, Sensitivity: Specificity) %>%
  ggplot(aes(x = x, y = y, color = model)) +
  labs(x="Measure of Performance", y="Value") + 
  ggtitle("Random Forest performance on test set - sampling comparison")+
  geom_jitter(width = 0.2, alpha = 0.5, size = 3)+
  scale_y_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1), limits=c(0,1.05))












#####
##########################################XGBoost, the best model
#####
#Load packages
#####
#load packages and calculate scale_pos_weight
library(tidyverse)
library(funModeling)
library(xgboost)
library(caret)
library(data.table)
library(SHAPforxgboost)
library(ggpubr)
source("shap.R")

#Class balance
training %>% count(Class)
#Weights for scale_pos_weight =  number of neg samples/number of pos samples = 1749/157438
#0.01110914

#####
#Convert to dummy variables, convert Class to integer classes
#####
dmytr = dummyVars(" ~ .", data = training[,-9], fullRank=T)
xgbdata = predict(dmytr, newdata = training[,-9])
target_var=ifelse(as.character(training$Class)=="Non_EAD",1,0)

#Convert Class to an integer class beginnning with 0 for test
dmytr1 = dummyVars(" ~ .", data = testing[,-9], fullRank=T)
xgbdatatest = predict(dmytr1, newdata = testing[,-9])
target_vartest=ifelse(as.character(testing$Class)=="Non_EAD",1,0)


#####
#XGBoost train model without weights and without sampling 
#####
params1 <- list(booster = "gbtree", objective = "binary:logistic")
set.seed(1235)
#XGB built in CV
cv.res = xgb.cv(params = params1, data=xgbdata, label = target_var ,nfold=10, nrounds = 500,eval_metric = "auc",
                showsd = T, stratified = T,
                print_every_n = 10,early_stopping_rounds =5)
#Run model with 134 rounds, best rounds according to cv.res
set.seed(1235)
xgb1 <- xgboost(data = xgbdata,nrounds=134, 
                objective="binary:logistic", label=target_var)
#Predictions and confusion matrix
testpredxgb <- predict(xgb1,xgbdatatest)
testpredxgb <- as.numeric(testpredxgb>0.5)
cmtest_xgboost <- confusionMatrix(factor(testpredxgb),factor(target_vartest), 
                                  mode="everything")
#####
#XGBoost train model with weights and without sampling 
#####
#Train model with CV
params <- list(booster = "gbtree", objective = "binary:logistic",
               scale_pos_weight = 0.01110914)
#Compute new CV
cv.res = xgb.cv(params = params, data=xgbdata, label = target_var ,nfold=10, nrounds = 500,eval_metric = "auc",
                showsd = T, stratified = T,
                print_every_n = 10,early_stopping_rounds =5)
#Best iteration at 159 rounds
set.seed(1235)
xgb2 <- xgboost(data = xgbdata, params = params,nrounds=159, eval_metric="auc", label=target_var)
#Predictions and confusion matrix for test
testpredxgb2 <- predict(xgb2,xgbdatatest)
testpredxgb2 <- as.numeric(testpredxgb2>0.5)
cmtest_xgboost_weights <- confusionMatrix(factor(testpredxgb2),factor(target_vartest), mode="everything")



#####
#xgboost on oversampled data
#####
##Convert Class to integer classes
#Convert to dummy variables
over_dmytr = dummyVars(" ~ .", data = over_training[,-9], fullRank=T)
over_xgbdata = predict(dmytr, newdata = over_training[,-9])
over_target_var=ifelse(as.character(over_training$Class)=="Non_EAD",1,0)

#Train model with CV
#Weights removed for sampled datasets 
sampledparams <- list(booster = "gbtree", objective = "binary:logistic")
set.seed(1235)
over.cv.res = xgb.cv(params = sampledparams, data=over_xgbdata, label = over_target_var ,nfold=10, nrounds = 500,eval_metric = "auc",
                     showsd = T, stratified = T,
                     print_every_n = 10,early_stopping_rounds =5)
#Best iteration at 44 rounds
set.seed(1235)
over_xgb2 <- xgboost(data = over_xgbdata, params = sampledparams,nrounds=44, eval_metric="auc", label=over_target_var)
#Predictions and confusion matrix for test
over_testpredxgb2 <- predict(over_xgb2,xgbdatatest)
over_testpredxgb2 <- as.numeric(over_testpredxgb2>0.5)
cmtest_xgboost_over <- confusionMatrix(factor(over_testpredxgb2),factor(target_vartest), mode="everything")


#####
#xgboost on undersampled data
#####
##Convert Class to integer classes
#Convert to dummy variables
under_dmytr = dummyVars(" ~ .", data = under_training[,-9], fullRank=T)
under_xgbdata = predict(dmytr, newdata = under_training[,-9])
under_target_var=ifelse(as.character(under_training$Class)=="Non_EAD",1,0)

#Train model with CV
#Weights removed for sampled datasets 
sampledparams <- list(booster = "gbtree", objective = "binary:logistic")
set.seed(1235)
under.cv.res = xgb.cv(params = sampledparams, data=under_xgbdata, label = under_target_var ,nfold=10, nrounds = 500,eval_metric = "auc",
                      showsd = T, stratified = T,
                      print_every_n = 10,early_stopping_rounds =5)
#Best iteration at 24  rounds
set.seed(1235)
under_xgb2 <- xgboost(data = under_xgbdata, params = sampledparams,nrounds=24, eval_metric="auc", label=under_target_var)
#Predictions and confusion matrix for test
under_testpredxgb2 <- predict(under_xgb2,xgbdatatest)
under_testpredxgb2 <- as.numeric(under_testpredxgb2>0.5)
cmtest_xgboost_under <- confusionMatrix(factor(under_testpredxgb2),factor(target_vartest), mode="everything")




#####
#xgboost on rose sampled data
#####
##Convert Class to integer classes
#Convert to dummy variables
rose_dmytr = dummyVars(" ~ .", data = rose_training[,-9], fullRank=T)
rose_xgbdata = predict(dmytr, newdata = rose_training[,-9])
rose_target_var=ifelse(as.character(rose_training$Class)=="Non_EAD",1,0)

#Train model with CV
#Weights removed for sampled datasets 
sampledparams <- list(booster = "gbtree", objective = "binary:logistic")
set.seed(1235)
rose.cv.res = xgb.cv(params = sampledparams, data=rose_xgbdata, label = rose_target_var ,nfold=10, nrounds = 500,eval_metric = "auc",
                     showsd = T, stratified = T,
                     print_every_n = 10,early_stopping_rounds =5)
#Best iteration at 72 rounds
set.seed(1235)
rose_xgb2 <- xgboost(data = rose_xgbdata, params = sampledparams,nrounds=72, eval_metric="auc", label=rose_target_var)
#Predictions and confusion matrix for test
rose_testpredxgb2 <- predict(rose_xgb2,xgbdatatest)
rose_testpredxgb2 <- as.numeric(rose_testpredxgb2>0.5)
cmtest_xgboost_rose <- confusionMatrix(factor(rose_testpredxgb2),factor(target_vartest), mode="everything")




#####
#xgboost on smote sampled data
#####
##Convert Class to integer classes
#Convert to dummy variables
smote_dmytr = dummyVars(" ~ .", data = smote_training[,-9], fullRank=T)
smote_xgbdata = predict(dmytr, newdata = smote_training[,-9])
smote_target_var=ifelse(as.character(smote_training$Class)=="Non_EAD",1,0)

#Train model with CV
#Weights removed for sampled datasets 
sampledparams <- list(booster = "gbtree", objective = "binary:logistic")
set.seed(1235)
smote.cv.res = xgb.cv(params = sampledparams, data=smote_xgbdata, label = smote_target_var ,nfold=10, nrounds = 500,eval_metric = "auc",
                      showsd = T, stratified = T,
                      print_every_n = 10,early_stopping_rounds =5)
#Best iteration at 35 rounds
set.seed(1235)
smote_xgb2 <- xgboost(data = smote_xgbdata, params = sampledparams,nrounds=35, eval_metric="auc", label=smote_target_var)
#Predictions and confusion matrix for test
smote_testpredxgb2 <- predict(smote_xgb2,xgbdatatest)
smote_testpredxgb2 <- as.numeric(smote_testpredxgb2>0.5)
cmtest_xgboost_smote <- confusionMatrix(factor(smote_testpredxgb2),factor(target_vartest), mode="everything")



######
#Confusion matrix on test plots
#####
draw_confusion_matrix(cmtest_xgboost) 
draw_confusion_matrix(cmtest_xgboost_2) 
draw_confusion_matrix(cmtest_xgboost_over) 
draw_confusion_matrix(cmtest_xgboost_under) 
draw_confusion_matrix(cmtest_xgboost_rose) 
draw_confusion_matrix(cmtest_xgboost_smote) 

#####
#Plot of confusion matrices for comparison - test set
#####
#Create list of models

models <- list(xgboost = predxgb,
               xgboost_weights = predxgb2,
               xgboost_over = over_predxgb2,
               xgboost_under = under_predxgb2,
               xgboost_rose = rose_predxgb2,
               xgboost_smote = smote_predxgb2)

#Comparing overall model performance
comparison <- data.frame(model = names(models))

for (name in names(models)) {
  model <- get(paste0("cmtest_", name))
  comparison[comparison$model == name, "Sensitivity"] <- model$byClass[["Sensitivity"]]
  comparison[comparison$model == name, "Precision"] <- model$byClass[["Precision"]]
  comparison[comparison$model == name, "Neg Pred Value"] <- model$byClass[["Neg Pred Value"]]
  comparison[comparison$model == name, "F1"] <- model$byClass[["F1"]]
  comparison[comparison$model == name, "Balanced Accuracy"] <- model$byClass[["Balanced Accuracy"]]
  comparison[comparison$model == name, "Specificity"] <- model$byClass[["Specificity"]]
}

#Create CM comparison plot
comparison %>%
  gather(x, y, Sensitivity: Specificity) %>%
  ggplot(aes(x = x, y = y, color = model)) +
  labs(x="Measure of Performance", y="Value") + 
  ggtitle("XGboost performance on test set - sampling comparison")+
  geom_jitter(width = 0.2, alpha = 0.5, size = 3)+
  scale_y_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1), limits=c(0,1.05))
#####
##############Comparison of best models
#####
#Create list of models
models <- list(glm = glmmodel,
               dectrees = dectreesoriginalfit,
               rforest = rforestoriginalfit,
               xgboost = xgb1)

####Comparison on training
comparison <- data.frame(model = names(models))

for (name in names(models)) {
  model <- get(paste0("cmtest_", name))
  comparison[comparison$model == name, "Sensitivity"] <- model$byClass[["Sensitivity"]]
  comparison[comparison$model == name, "Pos Pred Value"] <- model$byClass[["Pos Pred Value"]]
  comparison[comparison$model == name, "Neg Pred Value"] <- model$byClass[["Neg Pred Value"]]
  comparison[comparison$model == name, "F1"] <- model$byClass[["F1"]]
  comparison[comparison$model == name, "Balanced Accuracy"] <- model$byClass[["Balanced Accuracy"]]
  comparison[comparison$model == name, "Specificity"] <- model$byClass[["Specificity"]]
}



comparison %>%
  gather(x, y, Sensitivity: Specificity) %>%
  ggplot(aes(x = x, y = y, color = model)) +
  labs(x="Measure of Performance", y="Value") + 
  ggtitle("Model performance on testing set - comparison of best models") + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 3)

#####
#SHAP Values for original xgboost
#####
## Prepare shap data
shap_values <- shap.values(xgb_model = xgb1,X_train = xgbdata)
shap_long <- SHAPforxgboost::shap.prep(shap_contrib = shap_values$shap_score, X_train =xgbdata)

shap_result <- shap.score.rank(xgb_model = xgb1, 
                               X_train = xgbdata,
                               shap_approx = F)

## Plot shap summary plot
shapsumplot <- plot.shap.summary(data_long = shap_long) 

##Plot shap dependence 
plot.shap <- xgb.plot.shap(data = xgbdata, 
                           model = xgb1, 
                           features = names(shap_result$mean_shap_score)[1:8], 
                           n_col = 4, 
                           plot_loess = T)



