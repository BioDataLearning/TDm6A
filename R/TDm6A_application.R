#' @title <TDm6A_application>
#'
#' @param faFile_Input full path of the input fasta file storing the ID and sequence information
#' @param plot_Output full path of the output .png file to show the predicted m6A sites on the given sequence
#' @param csvFile_Output full path of the output .csv file to store the predicted results
#' @param cellType a string to be "A549", or "CD8T", or "HEK293"
#' @return NULL
#' @description This function is used to apply TDm6A models to predict m6A sites on a given nucleotide sequence
#' @export


TDm6A_application <- function(faFile_Input,plot_Output,csvFile_Output, cellType){
  ### refine the dataset:
  fa = read.fasta(faFile_Input, as.string = TRUE, forceDNAtolower = FALSE,whole.header = TRUE)
  fa = as.data.frame(unlist(fa))
  fa$ID = rownames(fa)
  colnames(fa) = c("Sequence","ID")
  fa$Sequence = as.character(fa$Sequence)
  fa$Sequence = toupper(fa$Sequence)
  fa$Sequence = sub("\n","",fa$Sequence)
  #### find the all the DRACH sites in lncRNAs:
  print("Running step 1: Finding all adenosine sites conforming to the DRACH motif...\n")
  DRACH_all = gregexpr("[AGT][AG]AC[ACT]", fa$Sequence)
  ### create dataframe to save each DRACH site:
  myData = as.data.frame(NA)
  colnames(myData) = "DRACH_site"
  ### put into different rows:
  DRACH_all = unlist(DRACH_all)

  i = 1; # df row
  j = 1; # DRACH site row
  while (j <= length(DRACH_all)) {
    myData[i,1] = DRACH_all[j]
    i = i + 1;
    j = j + 1;
  }

  #### find the Adenosine site
  myData$A_site = NA;
  myData$A_site = myData$DRACH_site +2;

  #### extract the 500t*2 +1 sequences
  print("Running step 2: Input feature preparation...\n")
  myData$flanking_500nt = NA;

  i = 1;
  while ( i <= nrow(myData))
  {
    myData$flanking_500nt[i] = substr(fa$Sequence, myData$A_site[i]-500, myData$A_site[i] + 500);
    if (myData$A_site[i] <= 500)
    {
      ga = rep("X", 501-myData$A_site[i]);
      gap = paste(ga, collapse = "");
      myData$flanking_500nt[i] = paste(gap, myData$flanking_500nt[i], sep = "");
    }
    if (nchar(fa$Sequence) - myData$A_site[i] < 500)
    {
      ga = rep("X", 500 + myData$A_site[i] - nchar(fa$Sequence));
      gap = paste(ga, collapse = "");
      myData$flanking_500nt[i] = paste(myData$flanking_500nt[i], gap, sep = "");
    }
    i = i + 1;
  }

  if(mean(nchar(myData$flanking_500nt)) == 1001)
  {
    print("1001nt sequences checked.\n")
  }else{
    print("Error for extracting flanking sequences.\n")
    stop("Program stopped.\n")
  }

  ##### one-hot encoding:
  myData$encoded_seq = NA

  i = 1;
  while ( i <= nrow(myData))
  {
    code_seq_A = gsub("A", "1000", myData$flanking_500nt[i]);
    code_seq_C = gsub("C", "0100", code_seq_A);
    code_seq_G = gsub("G", "0010", code_seq_C);
    code_seq_T = gsub("[TU]", "0001", code_seq_G);
    code_seq_X = gsub("X", "0000", code_seq_T);
    myData$encoded_seq[i] = as.character(code_seq_X);
    i = i + 1;
  }

  if(mean(nchar(myData$encoded_seq)) == 4004)
  {
    print("one-hot-encoding checked.\n")
  }else{
    print("Error for one-hot-encoding.\n")
    stop("Program stopped.\n")
  }

  # convert to 3D array x_train:
  x_array = array(data = NA, dim = c(4,1001, nrow(myData)));
  print(dim(x_array));

  s = 1;
  r = 1;
  c = 1;
  j = 1;
  while (s <= nrow(myData) )
  {
    c = 1;
    while ( c <= 1001)
    {
      r = 1;
      while (r <= 4)
      {
        x_array[r,c,s] = as.integer(substr(myData$encoded_seq[s],j,j))
        r = r + 1;
        j = j + 1;
      }
      c = c + 1;
    }
    s = s+1;
    j = 1;
  }

  if(any(is.na(x_array)))
  {
    print("Error for matrix preparation.\n")
    stop("Program stopped.\n")
  }else{
    print("Input data checked.\n")
  }

  ### change the index order to the ones keras package can use:
  x_test = aperm(x_array, c(3,2,1));
  print(dim(x_test));

  ############################################################################## model prediction ###################################################
  print("TDm6A predicting...\n")
  ## load models and predict:
  modelName <- paste("TDm6A_",cellType,"_1.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_1 = model %>% predict_proba(x_test);

  modelName <- paste("TDm6A_",cellType,"_2.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_2 = model %>% predict_proba(x_test);

  modelName <- paste("TDm6A_",cellType,"_3.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_3 = model %>% predict_proba(x_test);

  modelName <- paste("TDm6A_",cellType,"_4.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_4 = model %>% predict_proba(x_test);

  modelName <- paste("TDm6A_",cellType,"_5.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_5 = model %>% predict_proba(x_test);

  modelName <- paste("TDm6A_",cellType,"_6.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_6 = model %>% predict_proba(x_test);

  modelName <- paste("TDm6A_",cellType,"_7.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_7 = model %>% predict_proba(x_test);

  modelName <- paste("TDm6A_",cellType,"_8.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_8 = model %>% predict_proba(x_test);

  modelName <- paste("TDm6A_",cellType,"_9.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_9 = model %>% predict_proba(x_test);

  modelName <- paste("TDm6A_",cellType,"_10.hdf5", sep = "")
  path <- system.file(modelName, package ="TDm6A")
  model <- load_model_hdf5(filepath = path)
  classes_10 = model %>% predict_proba(x_test);

  classes = (classes_1 + classes_2 + classes_3 + classes_4 + classes_5 + classes_6 + classes_7 + classes_8 + classes_9 + classes_10)/10

  print("TDm6A prediction finished.\n")
  ### read the position of lncRNA DRACH sites:
  classes = as.data.frame(classes)
  colnames(classes) = c("Pred_Neg_Prob","Pred_Posi_Prob")
  myData = cbind.data.frame(myData,classes)

  ############################################### Plot:
  print("Plotting...\n")
  #### label the TP
  myData$predict_class = "P"

  i = 1;
  while ( i <= nrow(myData))
  {
    if ( myData$Pred_Posi_Prob[i]<0.5)
    {
      myData$predict_class[i] = "N"
    }
    i = i + 1
  }

  colour <- c("black","red")
  col.list <- rep(0,length(myData$predict_class))
  col.list[myData$predict_class=="N"] <-1
  col.list[myData$predict_class=="P"]<-2

  #### save plot
  png(plot_Output,width = 1200,height = 800)
  par(mar=c(5,5,5,5))
  plot(myData$A_site,myData$Pred_Posi_Prob,
       ylim = c(0,1),xlim = c(1,nchar(fa$Sequence)),
       col =c(colour[col.list]), lwd=1,
       xlab = "Nucleotide position on sequence",
       ylab="Probability to be m6A site",
       main = fa$ID, type = "h",
       cex=1.5,
       cex.axis=1.5,
       cex.lab=1.5)
  legend("topright",legend = c("Positive m6A sites","Negative m6A sites"),
         col=c("red","black"),
         lty = 1:1,
         cex = 1.5)
  dev.off()
  #### save predict result:
  myData=myData[,c(2,5,6,7)]
  write.csv(myData,csvFile_Output,row.names = FALSE)

  print("\nDone\n")
}





