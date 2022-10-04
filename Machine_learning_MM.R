#load necessary packages
fun_packages <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      #Si le package ne peut pas être chargé, on l'installe
      install.packages( i , dependencies = TRUE )
      #Charger le package après l'installation
      require( i , character.only = TRUE )
    }
    if (!requireNamespace("BiocManager", quietly = TRUE))
    {install.packages("BiocManager")}
  }
}

fun_packages(c("dplyr","randomForest","org.Hs.eg.db","clusterProfiler","caret"))
#BiocManager::install("randomForest")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("clusterProfiler")

message("The previous pre-processing step is done with jupyter, please output the file, our machine learning model is based directly on this result.")


Preprocess=function(data_mmrf){
  print("load the data that has been pre-processing")
df_preprocess=data.frame(data_mmrf)
  Correct_Colnames <- function(df) {
    
    delete.columns <- grep("(^X$)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
    
    if (length(delete.columns) > 0) {
      
      row.names(df) <- as.character(df[, grep("^X$", colnames(df))])
      #other data types might apply than character or 
      #introduction of a new separate column might be suitable
      
      df <- df[,-delete.columns]
      
      colnames(df) <- gsub("^X", "",  colnames(df))
      #X might be replaced by different characters, instead of being deleted
    }
    
    return(df)
  }
  
  df_correct=Correct_Colnames(df_preprocess)
  df_correct <- mutate_all(df_correct, function(x) as.numeric(as.character(x)))
  df_correct$HR_FLAG=factor(df_correct$HR_FLAG)
  df_withHR=subset(df_correct,select=HR_FLAG)
  df_withoutHR=df_correct[,-which(names(df_correct)%in%c("HR_FLAG"))]
  
  print("transform the entrez id to gene symbol")
  TransSymbol<-function(df_correct){
  trans_df=t(df_correct)
  rank_names <- bitr(rownames(trans_df), 
                     fromType="ENTREZID", 
                     toType=c("SYMBOL"),
                     OrgDb="org.Hs.eg.db")
  expression_filter=data.frame(trans_df[rank_names$ENTREZID,])
  expression_filter$SYMBOL=rank_names$SYMBOL
  merge_data=merge(rank_names,expression_filter,by="SYMBOL",all=F)
  rownames(merge_data)=make.names(merge_data$SYMBOL,unique=T)
  data_pret=merge_data[ , -which(colnames(merge_data) %in% c("SYMBOL","ENTREZID"))]
  df_pret=data.frame(t(data_pret))
  return(df_pret)
  }
  df_correct=TransSymbol(df_withoutHR)
  df_right=cbind(df_correct,df_withHR)
  return(df_right)}


# Preparation the train and test set for RandomForest
Train_verify_model=function(df_right){
  set.seed(123)
  print("separate the train dataset and test dataset: 70% train, 30% test for training")
  select_train <- sample(dim(df_right)[1], dim(df_right)[1]*0.7)
  dfMMRF_train <- df_right[select_train, ]
  dfMMRF_test <- df_right[-select_train, ]
  #dfMMRF_test
  
  dfMMRF_train.forest <- randomForest(x = dfMMRF_train[, colnames(dfMMRF_train) != "HR_FLAG"],
                                      y = dfMMRF_train$HR_FLAG)
  dfMMRF_train.forest
  
  plot(margin(dfMMRF_train.forest, dfMMRF_train$HR_FLAG), main = 'the plot of good prediction')
  
  #use train set to see the exacte of prediction
  train_predict <- predict(dfMMRF_train.forest, dfMMRF_train[, colnames(dfMMRF_train) != "HR_FLAG"])
  compare_train <- table(train_predict, dfMMRF_train$HR_FLAG)
  #compare_train
  #sum(diag(compare_train)/sum(compare_train))
  
  #train_predict   0   1
  #             0 289   0
  #             1   0 224
  
  print("use the test set to predict the high_HR patients and compare the result with the actual result")
  test_predict <- predict(dfMMRF_train.forest, dfMMRF_test[, colnames(dfMMRF_train) != "HR_FLAG"])
  compare_test <- table(dfMMRF_test$HR_FLAG, test_predict, dnn = c('Actual', 'Predicted'))
  print("The comparaison of prediction result and actual result for test set:")
  print(compare_test)
  #      Predicted
  #Actual   0   1
  #      0 104  20
  #      1  61  36
    return(list(dfMMRF_train.forest,dfMMRF_train))}




########################## Cross Validation #######################

Cross_validation=function(dfMMRF_train){
print("the cross validation for our model, attention!!! this procesure may take long times and might need to run in terminal")

# R program to implement
# repeated K-fold cross-validation

# setting seed to generate a
# reproducible random sampling
set.seed(1985)


control <- trainControl(method="repeatedcv", number=10, repeats=3)

seed <- 7
metric <- "Accuracy"

set.seed(seed)

mtry <- sqrt(ncol(dfMMRF_train))
tunegrid <- expand.grid(.mtry=mtry)
rf_group <- train(HR_FLAG~., data=dfMMRF_train, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)

print(rf_group)

#Random Forest 
#513 samples
#24048 predictors
#2 classes: '0', '1' 

#No pre-processing
#Resampling: Cross-Validated (10 fold, repeated 3 times) 
#Summary of sample sizes: 462, 461, 461, 461, 462, 462, ... 
#Resampling results:
  
#  Accuracy   Kappa    
#0.6258301  0.2310248
}




######################## classifier with Top 30 genes test ###################
Train_test=function(importance_mmrf){
print("classifier with Top 30 genes")

##
#chose top30 genes by sorting with importances
  mmrf_select <- rownames(importance_mmrf)[1:30]

  #train set and test set
  mmrf_train_top30 <- dfMMRF_train[ ,c(mmrf_select, 'HR_FLAG')]
  mmrf_test_top30 <- dfMMRF_test[ ,c(mmrf_select, 'HR_FLAG')]
  
  # randomForest
  set.seed(123)
  mmrf_train.forest_30 <- randomForest(HR_FLAG ~ ., data = mmrf_train_top30, importance = TRUE)
  mmrf_train.forest_30
  #Call:
    #randomForest(formula = HR_FLAG ~ ., data = mmrf_train_top30,      importance = TRUE) 
  #Type of random forest: classification
  #Number of trees: 500
  #No. of variables tried at each split: 5
  
  #OOB estimate of  error rate: 43.66%
  #Confusion matrix:
   # 0  1    class.error
  #0 215 74   0.2560554
  #1 150 74   0.6696429
  
  print("prediction with train set itself")
  #train set
  train_predict <- predict(mmrf_train.forest_30, mmrf_train_top30)
  compare_train <- table(train_predict, mmrf_train_top30$HR_FLAG)
  compare_train
  
  #train_predict   0   1
  #             0 289   0
  #             1   0 224
  
  print("prediction with test set")
  #test set evaluation
  test_predict <- predict(mmrf_train.forest_30, mmrf_test_top30)
  compare_test <- table(mmrf_test_top30$HR_FLAG, test_predict, dnn = c('Actual', 'Predicted'))
  compare_test
  Predicted
  #Actual  0  1
  #      0 88 36
  #      1 67 30
}
#####################################################################

#load data
data_mmrf <- read.csv2(file = 'df_preprocessing_and_merged.csv',header=T,sep=",")
df_right=Preprocess(data_mmrf)
#Model construction
Model_ele=Train_verify_model(df_right)
dfMMRF_train.forest=Model_ele[[1]]
dfMMRF_train=Model_ele[[2]]

##############  prediction of data ######################

#new_data=read.csv2(file = 'df_preprocessing_and_merged.csv',header=T,sep=",")
#pret_newdata=Preprocess(new_data)
#test_predict <- predict(dfMMRF_train.forest, new_data[, colnames(new_data) != "HR_FLAG"])

########################################################
#see the importances of each variable
summary(dfMMRF_train.forest)
importance_mmrf <- dfMMRF_train.forest$importance
importance_mmrf <- data.frame(importance(dfMMRF_train.forest))
#head(importance_mmrf)

print("The top20 gene that matters for our machine learning model ")
#sort with the value, the most important items are presented at the first
importance_mmrf_sort <-apply(importance_mmrf, 2, sort, decreasing=T)
importance_mmrf_sort[1:20,]

#SNORD89  AURKAIP1  RNVU1.18   RNA5S10   RNA5S16    RNU1.1    RNA5S7     ELANE 
#0.3366621 0.2741121 0.2673002 0.2242182 0.2070937 0.2050225 0.2023584 0.1943041 
#LINC02108   RNA5S17   RNVU1.7    RNA5S2   ANKRD11     NCOA3      CYLD      RNY1 
#0.1860030 0.1672151 0.1641209 0.1631687 0.1621456 0.1575083 0.1573002 0.1572664 
#RNA5S11     RFESD    RNA5S3   RNA5S14 
#0.1537682 0.1518332 0.1517751 0.1488466 

#export the table
#write.table(importance_mmrf_sort, 'importance_mmrf.txt', sep = '\t', col.names = NA, quote = FALSE)

#varImpPlot(dfMMRF_train.forest, n.var = min(30, nrow(dfMMRF_train.forest$importance)), main = 'Top 30 - variable importance')

############### cross validation and classifier #########################
#cross validation
#Cross_validation(dfMMRF_train)
# classifier with top 30 genes
#Train_test(importance_mmrf_sort)
