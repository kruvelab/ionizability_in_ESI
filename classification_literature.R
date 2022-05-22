library(magrittr)
library(caTools)
library(rpart)
library(randomForest)
library(dplyr)
library(caret)
library(e1071)
library(class)
library(readr)


dataset = read_delim('Fingerprints_literature_decoded.csv',
                     delim = ",",
                     col_names = TRUE,
                     trim_ws = TRUE)


#dataset = dataset %>%
  group_by(PubChemCID) %>%
  mutate(Class = max(Class)) %>%
  ungroup() %>%
  unique()


datasetlit2=dataset%>%
  dplyr::select(-c(V882, V883, V884, V885, V886, V887, V888, Datasetname, InChIKey, PubChemCID, X1, SDFlineStart, SDFID, SDFlineEnd, Fingerprint))

datasetlit2=datasetlit2%>%
  unique()

#------------------------------------
datasetlit2 = datasetlit2 %>%
  mutate(Class = as.factor(Class))


set.seed(123)
split <- sample.split(datasetlit2$'Class', SplitRatio = 0.75)
training_set <- subset(datasetlit2, split == TRUE) 
test_set <- subset(datasetlit2, split == FALSE)


gmodels::CrossTable(training_set$Class)
gmodels::CrossTable(test_set$Class)
gmodels::CrossTable(datasetlit2$Class)

training_set=training_set%>%
  dplyr::select(-SMILES)


#cross-validation
set.seed(123)
fitControl<-trainControl(
  method = "repeatedcv",
  number = 2,
  repeats = 5)

#knn
set.seed(123)
classifier_knn <- train (`Class`~.,
                     data=training_set%>%
                       select(-Name),
                     method="knn",
                     trControl=fitControl)

summary(classifier_knn)
classifier_knn

training_set_knn<-training_set %>%
  mutate(knn_pred=predict(classifier_knn, newdata = training_set)) %>%
  dplyr::select('Class', knn_pred, everything())

test_set_knn<-test_set %>%
  mutate(knn_pred=predict(classifier_knn, newdata = test_set, type = "raw")) %>%
  dplyr::select('Class', knn_pred, everything())


accuracy_knn = test_set_knn %>%
  summarise(accuracy_knn = mean(Class==knn_pred))

#accuracy 0.7865787

saveRDS(classifier_knn,
        file = "classifier_knn_literature.RDS")

#----------------------------------------------------
set.seed(123)
classifier_lit_xgblinear <- train (`Class`~.,
                               data=training_set%>%
                                 select(-c(Name)),
                               method="xgbLinear",
                               trControl=fitControl)

summary(classifier_lit_xgblinear)

#model = readRDS('classifier_xgblinear_literature.RDS')


training_set_xgblinear<-training_set %>%
  mutate(xgblinear_pred=predict(classifier_lit_xgblinear, newdata = training_set)) %>%
  dplyr::select('Class', xgblinear_pred, everything())

test_set_xgblinear<-test_set %>%
  mutate(xgblinear_pred=predict(classifier_lit_xgblinear, newdata = test_set, type = "raw")) %>%
  dplyr::select('Class', xgblinear_pred, everything())


accuracy_xgblinear = test_set_xgblinear %>%
  summarise(accuracy_xgblinear = mean(Class==xgblinear_pred))

#accuracy  0.7953795

test_set_filtered=test_set_xgblinear %>%
  filter(Class != xgblinear_pred)

data_filtered=test_set_filtered%>%
  select(Name, Class, xgblinear_pred)


write_delim(data_filtered,
            "falses_literature.csv",
            delim = ",")


gmodels::CrossTable(test_set_xgblinear$Class, test_set_xgblinear$xgblinear_pred)



# importance <- varImp(classifier_lit_xgblinear)
# importance
# or
# caret::varImp(model)
# 
# varimp <- caret::varImp(model)[[1]] %>%
#   rownames_to_column()
# 
# varimp <- varimp[10:1,]
# 
# barplot(varimp$Overall, main="Variable importance", horiz=TRUE,
#         names.arg=varimp$rowname,
#         las = 1)


saveRDS(classifier_lit_xgblinear,
        file = "classifier_xgblinear_literature.RDS")



----------------------------------------------------
#xgbTree
set.seed(123)
classifier_xgbTree <- train (`Class`~.,
                             data=training_set%>%
                               dplyr::select(-c(Name)),
                             method="xgbTree",
                             trControl=fitControl)

classifier_xgbTree

training_set_xgbtree<-training_set %>%
  mutate(xgbtree_pred=predict(classifier_xgbTree, newdata = training_set)) %>%
  dplyr::select('Class', xgbtree_pred, everything())


test_set_xgbtree<-test_set %>%
  mutate(xgbtree_pred=predict(classifier_xgbTree, newdata = test_set, type = "raw")) %>%
  dplyr::select('Class', xgbtree_pred, everything())


accuracy_xgbtree = test_set_xgbtree %>%
  summarise(accuracy_xgbtree = mean(Class==xgbtree_pred))

#accuracy 0.8162816


saveRDS(classifier_xgbTree,
        file = "classifier_xgbTree_literature.RDS")
        


gmodels::CrossTable(test_set_xgbtree$Class, test_set_xgbtree$xgbtree_pred)
gmodels::CrossTable(training_set_xgbtree$Class, training_set_xgbtree$xgbtree_pred)


        
--------------------------------------
#dtrees
set.seed(123)
classifier_dtrees <- train (`Class`~.,
                            data=training_set%>%
                              dplyr::select(-c(Name)),
                            method="rpart",
                            trControl=fitControl)

summary(classifier_dtrees)
classifier_dtrees


training_set_dtrees<-training_set %>%
  mutate(dtrees_pred=predict(classifier_dtrees, newdata = training_set)) %>%
  dplyr::select('Class', dtrees_pred, everything())

test_set_dtrees<-test_set %>%
  mutate(dtrees_pred=predict(classifier_dtrees, newdata = test_set, type = "raw")) %>%
  dplyr::select('Class', dtrees_pred, everything())


accuracy_dtrees = test_set_dtrees %>%
  summarise(accuracy_dtrees = mean(Class==dtrees_pred))

#accuracy 0.7491749


saveRDS(classifier_dtrees,
        file = "classifier_dtrees_literature.RDS")
--------------------------------
#xgbDART
set.seed(123)
classifier_xgbdart <- train (`Class`~.,
                             data=training_set%>%
                               select(-Name),
                             method="xgbDART",
                             trControl=fitControl)



-----------------------------
#svmLinear (not the newest)
set.seed(123)
classifier_svmlinear <- train (`Class`~.,
                               data=training_set%>%
                                 dplyr::select(-Name),
                               method="svmLinear",
                               trControl=fitControl)

classifier_svmlinear


training_set_svmlinear<-training_set %>%
  mutate(svmlinear_pred=predict(classifier_svmlinear, newdata = training_set)) %>%
  dplyr::select('Class', svmlinear_pred, everything())

test_set_svmlinear<-test_set %>%
  mutate(svmlinear_pred=predict(classifier_svmlinear, newdata = test_set, type = "raw")) %>%
  dplyr::select('Class', svmlinear_pred, everything())


accuracy_svmlinear = test_set_svmlinear %>%
  summarise(accuracy_svmlinear = mean(Class==svmlinear_pred))

#accuracy 0.7918051

----------------------------------
#svmPoly (not the newest)
set.seed(123)
classifier_svmpoly <- train (`Class`~.,
                             data=training_set%>%
                               dplyr::select(-Name),
                             method="svmPoly",
                             trControl=fitControl)

summary(classifier_svmpoly)
classifier_svmpoly

training_set_svmpoly<-training_set %>%
  mutate(svmpoly_pred=predict(classifier_svmpoly, newdata = training_set)) %>%
  dplyr::select('Class', svmpoly_pred, everything())

test_set_svmpoly<-test_set %>%
  mutate(svmpoly_pred=predict(classifier_svmpoly, newdata = test_set, type = "raw")) %>%
  dplyr::select('Class', svmpoly_pred, everything())


accuracy_svmpoly = test_set_svmpoly %>%
  summarise(accuracy_svmpoly = mean(Class==svmpoly_pred))

#accuracy 0.8161683


#--------------------------------
#balancing

gmodels::CrossTable(datasetlit2$Class)

dataset_up = sample_n(datasetlit2 %>%
                        filter(Class == 0),
                      size = 1800,
                      replace = TRUE)


dataset_down = sample_n(datasetlit2 %>%
                          filter(Class == 1),
                        size = 2000,
                        replace = FALSE)

dataset_new=dataset_up%>%
  bind_rows(dataset_down)


write_delim(dataset_new,
            "balanced_data_literature.csv",
            delim = ",")

#dataset_SMILES = read_delim('Fingerprints_literature_decoded.csv',
#                     delim = ",",
 #                    col_names = TRUE,
  #                   trim_ws = TRUE) %>%
  #select(Name, SMILES)

#dataset_new2 = dataset_new %>%
#  left_join(dataset_SMILES)


dataset_new = dataset_new %>%
  mutate(SMILES = as.factor(SMILES))


data_for_splitting = dataset_new %>%
  dplyr::select(SMILES) %>%
  unique()

split <- sample.split(data_for_splitting$SMILES, SplitRatio = 0.75)


data_for_splitting = data_for_splitting %>%
  mutate(split = split)

dataset_new = dataset_new %>%
  left_join(data_for_splitting)

training_set <- dataset_new %>%
  filter(split == TRUE) %>%
  dplyr::select(-split, -SMILES)

test_set <- dataset_new %>%
  filter(split == FALSE) %>%
  dplyr::select(-split, -SMILES)



gmodels::CrossTable(dataset_new$Class)
gmodels::CrossTable(training_set$Class)
gmodels::CrossTable(test_set$Class)


#--------
training_set = training_set %>%
  mutate(Class = as.factor(Class))

test_set = test_set %>%
  mutate(Class = as.factor(Class))
#---------------

set.seed(123)
fitControl<-trainControl(
  method = "repeatedcv",
  number = 2,
  repeats = 5)


classifier_knn <- train (`Class`~.,
                         data=training_set%>%
                           dplyr::select(-Name),
                         method="knn",
                         trControl=fitControl)

summary(classifier_knn)
classifier_knn

training_set_knn<-training_set %>%
  mutate(knn_pred=predict(classifier_knn, newdata = training_set)) %>%
  dplyr::select('Class', knn_pred, everything())

test_set_knn<-test_set %>%
  mutate(knn_pred=predict(classifier_knn, newdata = test_set, type = "raw")) %>%
  dplyr::select('Class', knn_pred, everything())


accuracy_knn = test_set_knn %>%
  summarise(accuracy_knn = mean(Class==knn_pred))
#0.7526205

accuracy_knn2 = training_set_knn %>%
  summarise((accuracy_knn2 = mean(Class==knn_pred)))
#0.7993675



#------------------

set.seed(123)
classifier_xgblinear <- train (`Class`~.,
                                   data=training_set%>%
                                     select(-c(Name)),
                                   method="xgbLinear",
                                   trControl=fitControl)

summary(classifier_xgblinear)
classifier_xgblinear


training_set_xgblinear<-training_set %>%
  mutate(xgblinear_pred=predict(classifier_xgblinear, newdata = training_set)) %>%
  dplyr::select('Class', xgblinear_pred, everything())

test_set_xgblinear<-test_set %>%
  mutate(xgblinear_pred=predict(classifier_xgblinear, newdata = test_set, type = "raw")) %>%
  dplyr::select('Class', xgblinear_pred, everything())


accuracy_xgblinear = test_set_xgblinear %>%
  summarise(accuracy_xgblinear = mean(Class==xgblinear_pred))
#0.78826

accuracy_xgblinear2 = training_set_xgblinear %>%
  summarise(accuracy_xgblinear2 = mean(Class==xgblinear_pred))
#0.902319



----------------------------------------------------
#xgbTree

classifier_xgbTree <- train (`Class`~.,
                             data=training_set%>%
                               dplyr::select(-c(Name)),
                             method="xgbTree",
                             trControl=fitControl)

classifier_xgbTree

training_set_xgbtree<-training_set %>%
  mutate(xgbtree_pred=predict(classifier_xgbTree, newdata = training_set)) %>%
  dplyr::select('Class', xgbtree_pred, everything())


test_set_xgbtree<-test_set %>%
  mutate(xgbtree_pred=predict(classifier_xgbTree, newdata = test_set, type = "raw")) %>%
  dplyr::select('Class', xgbtree_pred, everything())


accuracy_xgbtree = test_set_xgbtree %>%
  summarise(accuracy_xgbtree = mean(Class==xgbtree_pred))

#accuracy 0.7756813

accuracy_xgbtree2 = training_set_xgbtree %>%
  summarise(accuracy_xgbtree2 = mean(Class==xgbtree_pred))

#0.8826423


#----------------------------
classifier_dtrees <- train (`Class`~.,
                            data=training_set%>%
                              dplyr::select(-c(Name)),
                            method="rpart",
                            trControl=fitControl)

summary(classifier_dtrees)
classifier_dtrees


training_set_dtrees<-training_set %>%
  mutate(dtrees_pred=predict(classifier_dtrees, newdata = training_set)) %>%
  dplyr::select('Class', dtrees_pred, everything())

test_set_dtrees<-test_set %>%
  mutate(dtrees_pred=predict(classifier_dtrees, newdata = test_set, type = "raw")) %>%
  dplyr::select('Class', dtrees_pred, everything())


accuracy_dtrees = test_set_dtrees %>%
  summarise(accuracy_dtrees = mean(Class==dtrees_pred))

#accuracy 0.7327044

accuracy_dtrees2 = training_set_dtrees %>%
  summarise(accuracy_dtrees2 = mean(Class==dtrees_pred))

#0.7125791


