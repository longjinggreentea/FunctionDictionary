###############################################
#  This is Wenjun's function collections:
###############################################


###
#
###
#
#
#...................................................................
find_inf <- function(x) which(is.infinite(x)==TRUE)
#...................................................................



###
#
###
#
# Combine expression data and clinical data; 
# Remove patient id with expression data as inf; 
# Data for specifical gene is ready for OS analysis.
#...................................................................
my.stage.dependency.getdata <- function(mygene){
  
  find_inf <- function(x) which(is.infinite(x)==TRUE)
  
  myexpr.stage.df <- data.frame(t(log.scal.est.trim.update[mygene, comID.clin_scal.est.stage]), 
                                (LUAD.clin.os.working$stage)[which(rownames(LUAD.clin.os.working) %in% comID.clin_scal.est.stage)])
  
  colnames(myexpr.stage.df)<- c(mygene, "stage")
  
  
  # scal.est.trim.update[mytarget.gene,]
  rm.target.id <- rownames(myexpr.stage.df)[find_inf(myexpr.stage.df[,1])] # 126
  keep.target.id <- (rownames(myexpr.stage.df))[- (find_inf(myexpr.stage.df[,1]))] # 373
  
  sub.myexpr.stage.df <- myexpr.stage.df[keep.target.id,]
  
  if (length(rm.target.id) == 0 ){
    sub.myexpr.stage.df <- myexpr.stage.df
  }
  
  print(paste("The number of pateints to be removed is  ", c(length(rm.target.id)))) 
  print(paste("The number of pateints kept for analysis", nrow(sub.myexpr.stage.df)))
  # mygene.ID.sort <- list(rm.target.id, keep.target.id)
  # return(mygene.ID.sort)
  return(sub.myexpr.stage.df)
}
#...................................................................

###
#
###
#
#
#...................................................................
my.stage.dependency.boxplot<- function(mydata){ 
  # X11()
  my.plot.mygene<- ggboxplot(mydata, x=colnames(mydata)[2], y=mydata[1],
                             color = colnames(mydata)[2], palette = c(1:4), 
                             ylab = colnames(mydata)[1], xlab="",
                             add = "jitter", # font.label = list(size = 28, fact = "italic")
                             title = "",
                             size = 0.1, notch = FALSE
  ) 
  ggpar(my.plot.mygene, legend = "top", font.x = c(18, "bold"), font.y = c(18, "bold"))
}


my.stage.dependency.ANOVAassay<- function(mydata, mygene){   
  
  # ANOVA analysis
  #........... gene name needs to be updated manually.........
  res.aov.mygene <- aov(`HLA-DQB2` ~ stage, mydata); print(summary(res.aov.mygene))
  
  # X11(); par(mfrow=c(1,2)); 
  plot(res.aov.mygene, 1, main=mygene); plot(res.aov.mygene, 2, main=mygene)
  
  #........... gene name needs to be updated manually.........
  homo.var.test.mygene <- leveneTest(`HLA-DQB2` ~ stage, mydata); print(homo.var.test.mygene)
  normal.test.mygene <- shapiro.test(x=residuals(object = res.aov.mygene)); print(normal.test.mygene)   
  
  #........... gene name needs to be updated manually.........
  kruskal.res.mygene<- kruskal.test(`HLA-DQB2` ~ stage, data = mydata); print(kruskal.res.mygene)
  
}
#...................................................................



###
#
###
#
#
#...................................................................
match.check.rownames <- function(genData1, genData2){
  identical.res <- identical(rownames(genData1), rownames(genData2))
  if(identical.res == TRUE) {
    print("Rownames are in the same order")
  } else {
    print("Rownames don't match.")
    comm.sample.IDs <- intersect(rownames(genData1), rownames(Data2))
    return(comm.sample.IDs)
  }
}
#...................................................................


###
#
###
#
#
#........................................................................
## suppose there are missing data for gene expression
my.jointData.function <- function(geneData1, geneData2){
  # Double check that the rownames are not the same in lenghth 
  if (nrow(geneData1)!=nrow(geneData2)) {
    print("These two data sets are not the same in row number")                
  }
  # identify the common rownames (or sample ID)
  comm.Sample.ID <- intersect(rownames(geneData1), rownames(geneData2))
  print(paste("There are ", length(comm.Sample.ID), "common sample IDs between these two data set."))
  comm.df <- data.frame(geneData1[comm.Sample.ID, 1], geneData2[comm.Sample.ID, c(1:2)])
  rownames(comm.df) <- comm.Sample.ID
  colnames(comm.df) <- c(colnames(geneData1)[1], colnames(geneData2))
  return(comm.df)
}
#.........................................................................


###
#
###
#
#
#...................................................................
get.OS.Analy.ready <- function(targetDat, clinDat){
  # Note: colnames(targetDat)="geneName" "stage"; colnames(clinDat)="vital.status" "OS.time" "stage"
  # check patients from targetDat are all included in the clinDat
  question.id.len <- length(which(!(rownames(targetDat) %in% rownames(clinDat))))
  if (question.id.len > 0) {
    print(paste("There are", question.id.len, "patients without corresponding clinical data for OS analysis"))
  }
  # create dataset for OS.analysis
  OS.dat.targ <- data.frame(targetDat, clinDat[rownames(targetDat), c("vital.status", "OS.time")])
  OS.dat.targ$event <- ifelse(OS.dat.targ$vital.status=="dead", 1, 0)
  OS.dat.targ$time <- OS.dat.targ$OS.time/30.5
  colnames(OS.dat.targ) <-c(colnames(targetDat)[1:2], "vital.status", "OS.time", "event", "time.month")
  
  return(OS.dat.targ)
}
#...................................................................















