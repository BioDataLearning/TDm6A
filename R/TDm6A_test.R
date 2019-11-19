#' @title <TDm6A_test>
#'
#' @param x_test a 3D array
#' @param y_test a 2D array
#' @param cellType a string to be "A549", or "CD8T", or "HEK293"
#' @return NULL
#' @description This function test TDm6A models with test case
#' @export

TDm6A_test <- function(x_test,y_test,cellType){
  ### load the saved best model
  cat("\n\nLoading the TDm6A models...\n")

  # load the model
  cat("\n\nLoading the TDm6A model 1...\n")
  modelName <- paste("TDm6A_",cellType,"_1.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_1 = model %>% predict_proba(x_test);
  cat("\n\nLoading the TDm6A model 2...\n")
  modelName <- paste("TDm6A_",cellType,"_2.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_2 = model %>% predict_proba(x_test);
  cat("\n\nLoading the TDm6A model 3...\n")
  modelName <- paste("TDm6A_",cellType,"_3.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_3 = model %>% predict_proba(x_test);
  cat("\n\nLoading the TDm6A model 4...\n")
  modelName <- paste("TDm6A_",cellType,"_4.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_4 = model %>% predict_proba(x_test);
  cat("\n\nLoading the TDm6A model 5...\n")
  modelName <- paste("TDm6A_",cellType,"_5.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_5 = model %>% predict_proba(x_test);
  cat("\n\nLoading the TDm6A model 6...\n")
  modelName <- paste("TDm6A_",cellType,"_6.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_6 = model %>% predict_proba(x_test);
  cat("\n\nLoading the TDm6A model 7...\n")
  modelName <- paste("TDm6A_",cellType,"_7.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_7 = model %>% predict_proba(x_test);
  cat("\n\nLoading the TDm6A model 8...\n")
  modelName <- paste("TDm6A_",cellType,"_8.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_8 = model %>% predict_proba(x_test);
  cat("\n\nLoading the TDm6A model 9...\n")
  modelName <- paste("TDm6A_",cellType,"_9.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_9 = model %>% predict_proba(x_test);
  cat("\n\nLoading the TDm6A model 10...\n")
  modelName <- paste("TDm6A_",cellType,"_10.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_10 = model %>% predict_proba(x_test);
  cat("\n\nCalculating average...\n")
  classes = (classes_1 + classes_2 + classes_3 + classes_4 + classes_5 + classes_6 + classes_7 + classes_8 + classes_9 + classes_10)/10


  cat("\n\nTDm6A models loaded :D\n")

  ### generate predictions
  cat("\n\nPredicting test set...\n")

  ###
  pred <- prediction(classes[,2],y_test[,2])
  ### accuracy
  acc = performance(pred, "acc")
  accuracy = acc@y.values[[1]][max(which(acc@x.values[[1]] >= 0.5))]
  cat("\n\nAccuracy is:\n")
  print(accuracy)
  ### sensitivity
  sens = performance(pred, "sens")
  sensitivity = sens@y.values[[1]][max(which(sens@x.values[[1]] >= 0.5))]
  cat("\n\nSensitivity is:\n")
  print(sensitivity)
  ### Spectivity
  spec = performance(pred, "spec")
  specificity = spec@y.values[[1]][max(which(spec@x.values[[1]] >= 0.5))]
  cat("\n\nSpecificity is:\n")
  print(specificity)
  ### MCC
  mcc = performance(pred, "mat")
  MCC = mcc@y.values[[1]][max(which(mcc@x.values[[1]] >= 0.5))]
  cat("\n\nMCC is:\n")
  print(MCC)

  # draw ROCR curve
  perf <- performance(pred, "tpr","fpr");
  plot(perf, col = "blue", lty =3, lwd=3, main = "ROC Curve ");
  abline(a = 0, b = 1, lty = 3)
  # print the auc vaules
  auc <- performance(pred, "auc")
  cat("\n\nROC_AUC value is:\n")
  print(auc@y.values)
}

