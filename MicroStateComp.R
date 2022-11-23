library('R.matlab')
library(stringr)

path = '/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis/UpperBeta/Microstate_Stats/'
#allGrps = list.dirs(path = path, full.names = TRUE, recursive = FALSE)
allGrps = c()
for (grp_ in c('Group-01', 'Group-02', 'Group-03', 'Group-04', 'Group-05', 'Group-06', 'Group-07', 'Group-08')){
  allGrps = c(allGrps, file.path(path, grp_))
}

#### Microstate Occurence
output <- matrix(ncol=6, nrow=24)
totalCnt = 1
for (path_ in allGrps){
  print(path_)
  Occ = readMat(file.path(path_, 'Occurence.mat'))
  maxArr = c(mean(Occ$Occ[,1]), mean(Occ$Occ[,2]), mean(Occ$Occ[,3]), mean(Occ$Occ[,4]))
  print(maxArr)
  maxWhich = which.max(maxArr)
  allIndxs = c(1:4)
  remainIdxs = allIndxs[-c(maxWhich)]
  for (i in remainIdxs){
    tRes = t.test(Occ$Occ[,maxWhich], Occ$Occ[,i], paired = TRUE, alternative='two.sided')
    print(c(tRes$statistic, tRes$parameter, tRes$p.value, tRes$estimate))
    output[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
    output[totalCnt, 2] = paste(as.character(maxWhich), '-', as.character(i))
    output[totalCnt, 3] = round(tRes$statistic, 2)
    output[totalCnt, 4] = as.integer(tRes$parameter)
    output[totalCnt, 5] = round(tRes$p.value*24, 5)
    output[totalCnt, 6] = round(tRes$estimate, 2)
    totalCnt = totalCnt + 1
  }  
}
output = data.frame(output)
colnames(output) <- c('Group','MapDiff','Stats', 'Parameter', 'p.Value', 'meanDiff')
write.csv(output, file.path(path, 'MapOccurence.csv'))

############## Comparison with Maximum Value 

#### Microstate Occurence
output <- matrix(ncol=6, nrow=24)
totalCnt = 1
for (path_ in allGrps){
  print(path_)
  Occ = readMat(file.path(path_, 'Occurence.mat'))
  maxArr = c(mean(Occ$Occ[,1]), mean(Occ$Occ[,2]), mean(Occ$Occ[,3]), mean(Occ$Occ[,4]))
  print(maxArr)
  maxWhich = which.max(maxArr)
  allIndxs = c(1:4)
  remainIdxs = allIndxs[-c(maxWhich)]
  for (i in remainIdxs){
    tRes = t.test(Occ$Occ[,maxWhich], Occ$Occ[,i], paired = TRUE, alternative='two.sided')
    print(c(tRes$statistic, tRes$parameter, tRes$p.value, tRes$estimate))
    output[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
    output[totalCnt, 2] = paste(as.character(maxWhich), '-', as.character(i))
    output[totalCnt, 3] = round(tRes$statistic, 2)
    output[totalCnt, 4] = as.integer(tRes$parameter)
    output[totalCnt, 5] = round(tRes$p.value*24, 5)
    output[totalCnt, 6] = round(tRes$estimate, 2)
    totalCnt = totalCnt + 1
  }  
}
output = data.frame(output)
colnames(output) <- c('Group','MapDiff','Stats', 'Parameter', 'p.Value', 'meanDiff')
write.csv(output, file.path(path, 'MapOccurence.csv'))

#### Global Field Potential
output <- matrix(ncol=6, nrow=24)
totalCnt = 1
for (path_ in allGrps){
  print(path_)
  Gfp = readMat(file.path(path_, 'GFP.mat'))
  maxArr = c(mean(Gfp$Gfp[,1]), mean(Gfp$Gfp[,2]), mean(Gfp$Gfp[,3]), mean(Gfp$Gfp[,4]))
  print(maxArr)
  maxWhich = which.max(maxArr)
  allIndxs = c(1:4)
  remainIdxs = allIndxs[-c(maxWhich)]
  for (i in remainIdxs){
    tRes = t.test(Gfp$Gfp[,maxWhich], Gfp$Gfp[,i], paired = TRUE, alternative='two.sided')
    print(c(tRes$statistic, tRes$parameter, tRes$p.value, tRes$estimate))
    output[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
    output[totalCnt, 2] = paste(as.character(maxWhich), '-', as.character(i))
    output[totalCnt, 3] = round(tRes$statistic, 2)
    output[totalCnt, 4] = as.integer(tRes$parameter)
    output[totalCnt, 5] = round(tRes$p.value*24, 5)
    output[totalCnt, 6] = round(tRes$estimate, 2)
    totalCnt = totalCnt + 1
  }  
}
output = data.frame(output)
colnames(output) <- c('Group','MapDiff','Stats', 'Parameter', 'p.Value', 'meanDiff')
write.csv(output, file.path(path, 'GFP.csv'))

#### Microstate Correlation
output <- matrix(ncol=6, nrow=24)
totalCnt = 1
for (path_ in allGrps){
  print(path_)
  Msp = readMat(file.path(path_, 'MspatCorr.mat'))
  maxArr = c(mean(Msp$Msp[,1]), mean(Msp$Msp[,2]), mean(Msp$Msp[,3]), mean(Msp$Msp[,4]))
  print(maxArr)
  maxWhich = which.max(maxArr)
  allIndxs = c(1:4)
  remainIdxs = allIndxs[-c(maxWhich)]
  for (i in remainIdxs){
    tRes = t.test(Msp$Msp[,maxWhich], Msp$Msp[,i], paired = TRUE, alternative='two.sided')
    print(c(tRes$statistic, tRes$parameter, tRes$p.value, tRes$estimate))
    output[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
    output[totalCnt, 2] = paste(as.character(maxWhich), '-', as.character(i))
    output[totalCnt, 3] = round(tRes$statistic, 2)
    output[totalCnt, 4] = as.integer(tRes$parameter)
    output[totalCnt, 5] = round(tRes$p.value*24, 5)
    output[totalCnt, 6] = round(tRes$estimate, 2)
    totalCnt = totalCnt + 1
  }  
}
output = data.frame(output)
colnames(output) <- c('Group','MapDiff','Stats', 'Parameter', 'p.Value', 'meanDiff')
write.csv(output, file.path(path, 'MspatCorr.csv'))

#### Global Explained Variance
output <- matrix(ncol=6, nrow=24)
totalCnt = 1
for (path_ in allGrps){
  print(path_)
  Gev = readMat(file.path(path_, 'GEV.mat'))
  maxArr = c(mean(Gev$Gev[,1]), mean(Gev$Gev[,2]), mean(Gev$Gev[,3]), mean(Gev$Gev[,4]))
  print(maxArr)
  maxWhich = which.max(maxArr)
  allIndxs = c(1:4)
  remainIdxs = allIndxs[-c(maxWhich)]
  for (i in remainIdxs){
    tRes = t.test(Gev$Gev[,maxWhich], Gev$Gev[,i], paired = TRUE, alternative='two.sided')
    print(c(tRes$statistic, tRes$parameter, tRes$p.value, tRes$estimate))
    output[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
    output[totalCnt, 2] = paste(as.character(maxWhich), '-', as.character(i))
    output[totalCnt, 3] = round(tRes$statistic, 2)
    output[totalCnt, 4] = as.integer(tRes$parameter)
    output[totalCnt, 5] = round(tRes$p.value*24, 5)
    output[totalCnt, 6] = round(tRes$estimate, 2)
    totalCnt = totalCnt + 1
  }  
}
output = data.frame(output)
colnames(output) <- c('Group','MapDiff','Stats', 'Parameter', 'p.Value', 'meanDiff')
write.csv(output, file.path(path, 'GEV.csv'))

#### Microstate Durations
output <- matrix(ncol=6, nrow=24)
totalCnt = 1
for (path_ in allGrps){
  print(path_)
  Dur = readMat(file.path(path_, 'Duration.mat'))
  maxArr = c(mean(Dur$Dur[,1]), mean(Dur$Dur[,2]), mean(Dur$Dur[,3]), mean(Dur$Dur[,4]))
  print(maxArr)
  maxWhich = which.max(maxArr)
  allIndxs = c(1:4)
  remainIdxs = allIndxs[-c(maxWhich)]
  for (i in remainIdxs){
    tRes = t.test(Dur$Dur[,maxWhich], Dur$Dur[,i], paired = TRUE, alternative='two.sided')
    print(c(tRes$statistic, tRes$parameter, tRes$p.value, tRes$estimate))
    output[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
    output[totalCnt, 2] = paste(as.character(maxWhich), '-', as.character(i))
    output[totalCnt, 3] = round(tRes$statistic, 2)
    output[totalCnt, 4] = as.integer(tRes$parameter)
    output[totalCnt, 5] = round(tRes$p.value*24, 5)
    output[totalCnt, 6] = round(tRes$estimate, 2)
    totalCnt = totalCnt + 1
  }  
}
output = data.frame(output)
colnames(output) <- c('Group','MapDiff','Stats', 'Parameter', 'p.Value', 'meanDiff')
write.csv(output, file.path(path, 'Duration.csv'))

#### Microstate Coverage
output <- matrix(ncol=6, nrow=24)
totalCnt = 1
for (path_ in allGrps){
  print(path_)
  Cov = readMat(file.path(path_, 'Coverage.mat'))
  maxArr = c(mean(Cov$Cov[,1]), mean(Cov$Cov[,2]), mean(Cov$Cov[,3]), mean(Cov$Cov[,4]))
  print(maxArr)
  maxWhich = which.max(maxArr)
  allIndxs = c(1:4)
  remainIdxs = allIndxs[-c(maxWhich)]
  for (i in remainIdxs){
    tRes = t.test(Cov$Cov[,maxWhich], Cov$Cov[,i], paired = TRUE, alternative='two.sided')
    print(c(tRes$statistic, tRes$parameter, tRes$p.value, tRes$estimate))
    output[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
    output[totalCnt, 2] = paste(as.character(maxWhich), '-', as.character(i))
    output[totalCnt, 3] = round(tRes$statistic, 2)
    output[totalCnt, 4] = as.integer(tRes$parameter)
    output[totalCnt, 5] = round(tRes$p.value*24, 5)
    output[totalCnt, 6] = round(tRes$estimate, 2)
    totalCnt = totalCnt + 1
  }  
}
output = data.frame(output)
colnames(output) <- c('Group','MapDiff','Stats', 'Parameter', 'p.Value', 'meanDiff')
write.csv(output, file.path(path, 'Coverage.csv'))

##########################################

###### Transition Matrix (Significance Test) ######
library(ggplot2)
basePath = file.path(path,  'Baseline', 'TransitionMatrices')
BaseTPFiles = list.files(path = basePath, full.names = TRUE, recursive = FALSE)
allPValues = c()
allComrs = c()

BaseSubOrder = c()
for (bas in BaseTPFiles){
  strBase = strsplit(bas, '/')[[1]]
  sub_ = strsplit(strBase[length(strBase)], '.mat')[[1]]
  BaseSubOrder = c(BaseSubOrder, sub_)
}

for (path_ in allGrps){
  newPath = file.path(path_, 'TransitionMatrices')
  TPFiles = list.files(path = newPath, full.names = TRUE, recursive = FALSE)
  EmtSubOrder = c()
  for (bas in TPFiles){

  }
  
  TpArrayEmt = array(rep(1, 4*4*length(TPFiles)), dim=c(4, 4, length(TPFiles)))
  TpArrayBase = array(rep(1, 4*4*length(TPFiles)), dim=c(4, 4, length(TPFiles)))
  fileCount = 1
  
  for (fil_ in TPFiles){
    TPData = readMat(fil_)
    TpArrayEmt[, , fileCount] = TPData$observed.TM
    ### Getting the Subject Name
    strBase = strsplit(fil_, '/')[[1]]
    sub_ = strsplit(strBase[length(strBase)], '_')[[1]][1]
    ###
    baseData = readMat(file.path(basePath, paste(sub_, '.mat', sep = '')))$observed.TM
    TpArrayBase[, , fileCount] = baseData
    fileCount = fileCount + 1
  }
  ### Calculating Significance Now
  for (row_ in c(1:4)){
    for (col_ in c(1:4)){
      if (row_!=col_){
        tRes = t.test(TpArrayEmt[row_, col_, ], TpArrayBase[row_, col_, ], paired = TRUE, alternative='two.sided') 
        allPValues = c(allPValues, tRes$p.value)
        allComrs = c(allComrs, paste(strsplit(path_, '//')[[1]][2], as.character(row_), as.character(col_), sep = '_'))
      }
    }
  }
}

  PValFrame = data.frame(allComrs, allPValues)
  write.csv(PValFrame, file.path(path, 'TransMatPVals.csv'))
  sortIndx = order(allPValues)
  dividend = 97
  corrPVals = c()
  
  for (idx_ in c(1:length(sortIndx))){
    print((dividend-idx_))
    corrPVals = c(corrPVals, allPValues[sortIndx[idx_]]/(dividend-idx_))
  }
  corrPVals<0.01
  signiTransMat = allComrs[sortIndx][corrPVals<0.01]
  
########### Transition Matrix Plotting: Difference with Baseline ############
  library(ggplot2)
  baseTM = readMat(file.path(path, 'Baseline', 'TransitionMatrix.mat'))
  baseTM = baseTM$TM
  for (path_ in allGrps){
    print(path_)
    TM = readMat(file.path(path_, 'TransitionMatrix.mat'))
    TransMat = TM$TM#/(max(TM$TM))
    diff_TM = TransMat-baseTM
    diff_TM = data.frame(diff_TM)
    TransMat = data.frame(TransMat)
    X = c('Map1', 'Map2', 'Map3', 'Map4')
    Y = c('Map1', 'Map2', 'Map3', 'Map4')
    ####################
    colnames(TransMat) = X
    rownames(TransMat) = Y
    heatData = expand.grid(Y=X,X=Y)
    heatData$val = round(as.vector(t(TransMat)), 2)
    ### Difference between Transition Matrices
    colnames(diff_TM) = X
    rownames(diff_TM) = Y
    heatDataDiff = expand.grid(Y=X, X=Y)
    #heatDataDiff$val = round(as.vector(t(diff_TM)), 2)
    heatData$diff_TM = round(as.vector(t(diff_TM)), 2)
    heatData$Sign = '0'
    heatData$Sign[heatData$diff_TM<0.0] = '-1'
    heatData$Sign[heatData$diff_TM>0.0] = '1'
    
    title_ = strsplit(path_, "//")[[1]][2]
    ### Significance Annotation
    for (id_ in c(1:16)){
      pattern = paste(title_, as.integer((id_-1)/4)+1, ((id_-1)%%4)+1, sep = '_')
      Presence = str_detect(signiTransMat, pattern)
      if (is.element(TRUE, Presence)){
        heatData[id_, 'Significant'] = '*'
      }else{
        heatData[id_, 'val'] = ''
        heatData[id_, 'Sign'] = '0'
        heatDataDiff[id_, 'val'] = ''
      }
    }
    
    plot_ = ggplot(heatData, aes(X, Y, fill=Sign), show.legend=FALSE) + geom_tile() + ggtitle(title_) + xlab("") + ylab("") + 
            geom_text(aes(label = val), size=10) + theme(text = element_text(size = 22)) #+
            #scale_fill_distiller(palette = "Spectral") + theme(legend.position = "none")
    ggsave(file.path(strsplit(path_, "//")[[1]][1], paste(title_, "-Baseline.png", sep = '')))
  }
  
  
  ########### Transition Matrix Plotting: Difference with Neutral ############
  library(ggplot2)
  baseTM = readMat(file.path(path, 'Neutral', 'TransitionMatrix.mat'))
  baseTM = baseTM$TM
  
  for (path_ in allGrps){
    print(path_)
    TM = readMat(file.path(path_, 'TransitionMatrix.mat'))
    TransMat = TM$TM#/(max(TM$TM))
    diff_TM = TransMat-baseTM
    diff_TM = data.frame(diff_TM)
    TransMat = data.frame(TransMat)
    X = c('Map1', 'Map2', 'Map3', 'Map4')
    Y = c('Map1', 'Map2', 'Map3', 'Map4')
    ####################
    colnames(TransMat) = X
    rownames(TransMat) = Y
    heatData = expand.grid(Y=X,X=Y)
    heatData$val = round(as.vector(t(TransMat)), 2)
    ### Difference between Transition Matrices
    colnames(diff_TM) = X
    rownames(diff_TM) = Y
    heatDataDiff = expand.grid(Y=X, X=Y)
    #heatDataDiff$val = round(as.vector(t(diff_TM)), 2)
    heatData$diff_TM = round(as.vector(t(diff_TM)), 2)
    heatData$Sign = '0'
    heatData$Sign[heatData$diff_TM<0.0] = '-1'
    heatData$Sign[heatData$diff_TM>0.0] = '1'
    
    title_ = strsplit(path_, "//")[[1]][2]
    ### Significance Annotation
    for (id_ in c(1:16)){
      pattern = paste(title_, as.integer((id_-1)/4)+1, ((id_-1)%%4)+1, sep = '_')
      Presence = str_detect(signiTransMat, pattern)
      if (is.element(TRUE, Presence)){
        heatData[id_, 'Significant'] = '*'
      }else{
        heatData[id_, 'val'] = ''
        heatData[id_, 'Sign'] = '0'
        heatDataDiff[id_, 'val'] = ''
      }
    }
    
    plot_ = ggplot(heatData, aes(X, Y, fill=Sign), show.legend=FALSE) + geom_tile() + ggtitle(title_) + xlab("") + ylab("") + 
      geom_text(aes(label = val), size=10) + theme(text = element_text(size = 22)) #+
    #scale_fill_distiller(palette = "Spectral") + theme(legend.position = "none")
    ggsave(file.path(strsplit(path_, "//")[[1]][1], paste(title_, "-Neutral.png", sep = '')))
  }
  
    #plot_ = ggplot(heatDataDiff, aes(X, Y, fill=val), show.legend=FALSE) + geom_tile() + ggtitle(title_) + xlab("") + ylab("") + 
    #  geom_tile(aes(fill = val)) + geom_text(aes(label = val), size=10)+
    #  theme(text = element_text(size = 22))+theme(legend.position = "none")
    #ggsave(file.path(strsplit(path_, "//")[[1]][1], paste(title_, "-Baseline.png", sep = '')))
    
      #plot_ = ggplot(heatData, aes(X, Y, fill=val)) + geom_tile() + ggtitle(title_) + xlab("") + ylab("") + 
      #geom_tile(aes(fill = val)) + geom_text(aes(label = val), size=10)+
      #theme(text = element_text(size = 22))
      
      #theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20))
      #geom_point(aes(size = val)) + scale_fill_distiller(palette = "Spectral") + guides(size = "legend") +
      #theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20), legend.text = element_text(size = 20),
      #      legend.title = element_text(size = 20))
  #}
  
######### Distance Between Transition Probabilities ###########
  distance <- matrix(0, ncol=8, nrow=8)
  p_DistSign <- matrix(0, ncol=8, nrow=8)

  firstCnt = 1
  for (path_ in allGrps){
    secondCnt = 1
   # print(path_)
    TM = readMat(file.path(path_, 'TransitionMatrix.mat'))
    TransMat_1 = TM$TM#/(max(TM$TM))
    
    for (path_ in allGrps){
      #print(path_)
      TM = readMat(file.path(path_, 'TransitionMatrix.mat'))
      TransMat_2 = TM$TM#/(max(TM$TM))
      poolA = c()
      poolB = c()
      
      for (i in c(1:4)){
        for (j in c(1:4)){
          if (i != j){
            poolA = c(poolA, TransMat_1[i, j])
            poolB = c(poolB, TransMat_2[i, j])
            distance[firstCnt, secondCnt] = distance[firstCnt, secondCnt] + (abs(TransMat_1[i, j] - TransMat_2[i, j])) 
          }
        }
      }
      print(poolA)
      print(poolB)
      print('-------------------')
      #signRes = wilcox.test(poolA, poolB, conf.int = TRUE)
      signRes = t.test(poolA, poolB, conf.int = TRUE)
      p_DistSign[firstCnt, secondCnt] = signRes$p.value
      secondCnt = secondCnt + 1
    }
    firstCnt = firstCnt + 1
  }
  
distance = round(distance, 3)
write.csv(distance, file.path(strsplit(path_, "//")[[1]][1], "Distance_TransitionProb.csv"))
# How is it different with respect to baseline.
# axis.title = element_text(size = 20), 
# scale_fill_continuous(limits = c(0,1), breaks = c(0, 0.5, 0.7, 0.8, 0.9, 1.0),
#guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE))