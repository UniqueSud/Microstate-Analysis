library('R.matlab')
library('stringr')

path = '/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis/UpperBeta/Microstate_Stats/'
allGrps = list.dirs(path = path, full.names = TRUE, recursive = FALSE)

baselinePath = file.path(path, 'Baseline')
fileBase = readMat(file.path(baselinePath, 'fileNames.mat'))
baseFiles = c()
for (i in c(1:length(fileBase$fileNames))){
  baseFiles = c(baseFiles, fileBase$fileNames[[i]][1][[1]][1])
}

neutralPath = file.path(path, 'Neutral')
fileNeutral = readMat(file.path(neutralPath, 'fileNames.mat'))
NeutFiles = c()
for (i in c(1:length(fileNeutral$fileNames))){
  NeutFiles = c(NeutFiles, strsplit(fileNeutral$fileNames[[i]][1][[1]][1], '_')[[1]][1])
}

# ----------------- Microstate Occurrence
maxIter = 500

for (mapIdx in c(1:4)){
  StatsFrame <- matrix(data = 0, ncol = 8, nrow = 8)
  row_ = 1
  col_ = 1
  for (path_ in allGrps){
    if (str_detect(strsplit(path_, split = '//')[[1]][2], 'Group')){
      print(path_)
      Occ = readMat(file.path(path_, 'Occurence.mat'))
      #firstGrp = as.vector(Occ$Occ)
      firstGrp = Occ$Occ[, mapIdx]
      for (path_ in allGrps){
        if (str_detect(strsplit(path_, split = '//')[[1]][2], 'Group')){
          print(strsplit(path_, split = '//')[[1]][2])
          Occ = readMat(file.path(path_, 'Occurence.mat'))
          #print(colMeans(Occ$Occ))
          #secondGrp = as.vector(Occ$Occ) 
          secondGrp = Occ$Occ[, mapIdx]
          minVal = min(c(length(firstGrp), length(secondGrp)))
          idxx = which(c(length(firstGrp), length(secondGrp))==minVal)
          if (length(idxx) < 2){
            countSign = 0
            for (j in c(1:maxIter)){ ## Because we are sampling minimum number of samples from one of the observation pool to match with the second pool.
              if (idxx == 1){
                secondGrp = sample(secondGrp, size = length(firstGrp))
              }else{
                firstGrp = sample(firstGrp, size = length(secondGrp))
              }
              t_test = t.test(firstGrp, secondGrp, paired = FALSE)
              if (t_test$p.value<0.05){
                countSign = countSign + 1
                diff_ = round(t_test$estimate[1]-t_test$estimate[2], 2)
                if (StatsFrame[row_, col_] < diff_){
                  StatsFrame[row_, col_] = diff_ 
                }
              } 
            }
            if (countSign < (maxIter/2)){
              StatsFrame[row_, col_] = 0
            }
          }else{
            t_test = t.test(firstGrp, secondGrp, paired = FALSE)
            if (t_test$p.value<0.05){
              #print(t_test$statistic)
              print(row_)
              print(col_)
              StatsFrame[row_, col_] = round(t_test$estimate[1]-t_test$estimate[2], 2)
            } 
          }
          col_ = col_ + 1
        }
      }
      row_ = row_ + 1
      col_ = 1
    }
  }
  rownames(StatsFrame) = c('G-1', 'G-2', 'G-3', 'G-4', 'G-5', 'G-6', 'G-7', 'G-8')
  colnames(StatsFrame) = c('G-1', 'G-2', 'G-3', 'G-4', 'G-5', 'G-6', 'G-7', 'G-8')
  write.csv(StatsFrame, file.path(path, paste('GroupWiseOccurrenceComparisonForMap-', as.character(mapIdx), '.csv', sep = '')))
}

# ----------------- Microstate Duration
maxIter = 500

for (mapIdx in c(1:4)){
  StatsFrame <- matrix(data = 0, ncol = 8, nrow = 8)
  row_ = 1
  col_ = 1
  for (path_ in allGrps){
    if (str_detect(strsplit(path_, split = '//')[[1]][2], 'Group')){
      print(path_)
      Dur = readMat(file.path(path_, 'Duration.mat'))
      #firstGrp = as.vector(Occ$Occ)
      firstGrp = Dur$Dur[, mapIdx]
      for (path_ in allGrps){
        if (str_detect(strsplit(path_, split = '//')[[1]][2], 'Group')){
          print(strsplit(path_, split = '//')[[1]][2])
          Dur = readMat(file.path(path_, 'Duration.mat'))
          #print(colMeans(Occ$Occ))
          #secondGrp = as.vector(Occ$Occ) 
          secondGrp = Dur$Dur[, mapIdx]
          minVal = min(c(length(firstGrp), length(secondGrp)))
          idxx = which(c(length(firstGrp), length(secondGrp))==minVal)
          if (length(idxx) < 2){
            countSign = 0
            for (j in c(1:maxIter)){ ## Because we are sampling minimum number of samples from one of the observation pool to match with the second pool.
              if (idxx == 1){
                secondGrp = sample(secondGrp, size = length(firstGrp))
              }else{
                firstGrp = sample(firstGrp, size = length(secondGrp))
              }
              t_test = t.test(firstGrp, secondGrp, paired = FALSE)
              if (t_test$p.value<0.05){
                countSign = countSign + 1
                diff_ = round(t_test$estimate[1]-t_test$estimate[2], 2)
                if (StatsFrame[row_, col_] < diff_){
                  StatsFrame[row_, col_] = diff_ 
                }
              } 
            }
            if (countSign < (maxIter/2)){
              StatsFrame[row_, col_] = 0
            }
          }else{
            t_test = t.test(firstGrp, secondGrp, paired = FALSE)
            if (t_test$p.value<0.05){
              #print(t_test$statistic)
              print(row_)
              print(col_)
              StatsFrame[row_, col_] = round(t_test$estimate[1]-t_test$estimate[2], 2)
            } 
          }
          col_ = col_ + 1
        }
      }
      row_ = row_ + 1
      col_ = 1
    }
  }
  rownames(StatsFrame) = c('G-1', 'G-2', 'G-3', 'G-4', 'G-5', 'G-6', 'G-7', 'G-8')
  colnames(StatsFrame) = c('G-1', 'G-2', 'G-3', 'G-4', 'G-5', 'G-6', 'G-7', 'G-8')
  write.csv(StatsFrame, file.path(path, paste('GroupWiseDurationComparisonForMap-', as.character(mapIdx), '.csv', sep = '')))  
}


############## ============================= Comparison with Neutral and Baseline only =========================================##################
#### Microstate Occurence
outputBase <- matrix(ncol=6, nrow=32)
outputNeut <- matrix(ncol=6, nrow=32)
totalCnt = 1
for (path_ in allGrps){
  if (str_detect(strsplit(path_, split = '//')[[1]][2], 'Group')){
    #print(path_)
    Occ = readMat(file.path(path_, 'Occurence.mat'))
    fileGrp = readMat(file.path(path_, 'fileNames.mat'))
    GrpFiles = c()
    for (i in c(1:length(fileGrp$fileNames))){
      GrpFiles = c(GrpFiles, strsplit(fileGrp$fileNames[[i]][1][[1]][1], split = '_')[[1]][1])
    }
    
    OccBase = readMat(file.path(path, 'Baseline', 'Occurence.mat'))
    BaseMapMat = matrix(, nrow = length(GrpFiles), ncol = 4, byrow = TRUE)
    for (fl_ in baseFiles){
      #print(fl_)
      idxs = which(fl_ == GrpFiles)
      idxsBase = which(fl_ == baseFiles)
      for (j in idxs){
        BaseMapMat[j, ] = OccBase$Occ[idxsBase, ]  
      }
    }
    
    OccNeutral = readMat(file.path(path, 'Neutral', 'Occurence.mat'))
    NeutMapMat = matrix(, nrow = length(GrpFiles), ncol = 4, byrow = TRUE)
    for (fl_ in NeutFiles){
      idxs = which(fl_ == GrpFiles)
      idxsNeut = which(fl_ == NeutFiles)
      for (j in idxs){
        NeutMapMat[j, ] = OccNeutral$Occ[idxsNeut, ]  
      }
    }
    
    for (i in c(1:4)){
      tRes = t.test(Occ$Occ[,i], BaseMapMat[,i], paired = TRUE, alternative='two.sided')
      outputBase[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
      outputBase[totalCnt, 2] = paste('Map-', as.character(i), sep = '')
      outputBase[totalCnt, 3] = round(tRes$statistic, 2)
      outputBase[totalCnt, 4] = as.integer(tRes$parameter)
      outputBase[totalCnt, 5] = round(tRes$p.value*24, 5)
      outputBase[totalCnt, 6] = round(tRes$estimate, 2)
      
      tRes = t.test(Occ$Occ[,i], NeutMapMat[,i], paired = TRUE, alternative='two.sided')
      outputNeut[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
      outputNeut[totalCnt, 2] = paste('Map-', as.character(i), sep = '')
      outputNeut[totalCnt, 3] = round(tRes$statistic, 2)
      outputNeut[totalCnt, 4] = as.integer(tRes$parameter)
      outputNeut[totalCnt, 5] = round(tRes$p.value*24, 5)
      outputNeut[totalCnt, 6] = round(tRes$estimate, 2)
      totalCnt = totalCnt + 1
    }  
  }
}

outputBase = data.frame(outputBase)
colnames(outputBase) <- c('Group', 'Map', 'Stats', 'Parameter', 'p.Value', 'meanDiff')
outputBase = outputBase[order(outputBase$p.Value), ] 
outputBase$adjPval = as.numeric(as.vector(outputBase$p.Value))*seq(32,1,by=-1)
write.csv(outputBase, file.path(path, 'Group-Base_MapOccurence.csv'))

outputNeut = data.frame(outputNeut)
colnames(outputNeut) <- c('Group', 'Map','Stats', 'Parameter', 'p.Value', 'meanDiff')
outputNeut = outputNeut[order(outputNeut$p.Value), ] 
outputNeut$adjPval = as.numeric(as.vector(outputNeut$p.Value))*seq(32,1,by=-1)
write.csv(outputNeut, file.path(path, 'Group-Neutral_MapOccurence.csv'))


#### Microstate Duration
outputBase <- matrix(ncol=6, nrow=32)
outputNeut <- matrix(ncol=6, nrow=32)

totalCnt = 1

for (path_ in allGrps){
  if (str_detect(strsplit(path_, split = '//')[[1]][2], 'Group')){
    Dur = readMat(file.path(path_, 'Duration.mat'))
    print(Dur$Dur)
  }
}
  
for (path_ in allGrps){
  if (str_detect(strsplit(path_, split = '//')[[1]][2], 'Group')){
    #print(path_)
    Dur = readMat(file.path(path_, 'Duration.mat'))
    fileGrp = readMat(file.path(path_, 'fileNames.mat'))
    GrpFiles = c()
    for (i in c(1:length(fileGrp$fileNames))){
      GrpFiles = c(GrpFiles, strsplit(fileGrp$fileNames[[i]][1][[1]][1], split = '_')[[1]][1])
    }
    
    DurBase = readMat(file.path(path, 'Baseline', 'Duration.mat'))
    BaseMapMat = matrix(, nrow = length(GrpFiles), ncol = 4, byrow = TRUE)
    for (fl_ in baseFiles){
      #print(fl_)
      idxs = which(fl_ == GrpFiles)
      idxsBase = which(fl_ == baseFiles)
      for (j in idxs){
        BaseMapMat[j, ] = DurBase$Dur[idxsBase, ]  
      }
    }
    
    DurNeutral = readMat(file.path(path, 'Neutral', 'Duration.mat'))
    NeutMapMat = matrix(, nrow = length(GrpFiles), ncol = 4, byrow = TRUE)
    for (fl_ in NeutFiles){
      idxs = which(fl_ == GrpFiles)
      idxsNeut = which(fl_ == NeutFiles)
      for (j in idxs){
        NeutMapMat[j, ] = DurNeutral$Dur[idxsNeut, ]  
      }
    }
    
    for (i in c(1:4)){
      tRes = t.test(Dur$Dur[,i], BaseMapMat[,i], paired = TRUE, alternative='two.sided')
      outputBase[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
      outputBase[totalCnt, 2] = paste('Map-', as.character(i), sep = '')
      outputBase[totalCnt, 3] = round(tRes$statistic, 2)
      outputBase[totalCnt, 4] = as.integer(tRes$parameter)
      outputBase[totalCnt, 5] = round(tRes$p.value*24, 5)
      outputBase[totalCnt, 6] = round(tRes$estimate, 2)
      
      tRes = t.test(Dur$Dur[,i], NeutMapMat[,i], paired = TRUE, alternative='two.sided')
      outputNeut[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
      outputNeut[totalCnt, 2] = paste('Map-', as.character(i), sep = '')
      outputNeut[totalCnt, 3] = round(tRes$statistic, 2)
      outputNeut[totalCnt, 4] = as.integer(tRes$parameter)
      outputNeut[totalCnt, 5] = round(tRes$p.value*24, 5)
      outputNeut[totalCnt, 6] = round(tRes$estimate, 2)
      totalCnt = totalCnt + 1
    }  
  }
}

outputBase = data.frame(outputBase)
colnames(outputBase) <- c('Group', 'Map', 'Stats', 'Parameter', 'p.Value', 'meanDiff')
outputBase = outputBase[order(outputBase$p.Value), ] 
outputBase$adjPval = as.numeric(as.vector(outputBase$p.Value))*seq(32,1,by=-1)
write.csv(outputBase, file.path(path, 'Group-Base_MapDuration.csv'))

outputNeut = data.frame(outputNeut)
colnames(outputNeut) <- c('Group', 'Map','Stats', 'Parameter', 'p.Value', 'meanDiff')
outputNeut = outputNeut[order(outputNeut$p.Value), ] 
outputNeut$adjPval = as.numeric(as.vector(outputNeut$p.Value))*seq(32,1,by=-1)
write.csv(outputNeut, file.path(path, 'Group-Neutral_MapDuration.csv'))


#### Microstate Coverage
outputBase <- matrix(ncol=6, nrow=32)
outputNeut <- matrix(ncol=6, nrow=32)

totalCnt = 1
for (path_ in allGrps){
  if (str_detect(strsplit(path_, split = '//')[[1]][2], 'Group')){
    #print(path_)
    Cov = readMat(file.path(path_, 'Coverage.mat'))
    fileGrp = readMat(file.path(path_, 'fileNames.mat'))
    GrpFiles = c()
    for (i in c(1:length(fileGrp$fileNames))){
      GrpFiles = c(GrpFiles, strsplit(fileGrp$fileNames[[i]][1][[1]][1], split = '_')[[1]][1])
    }
    
    CovBase = readMat(file.path(path, 'Baseline', 'Coverage.mat'))
    BaseMapMat = matrix(, nrow = length(GrpFiles), ncol = 4, byrow = TRUE)
    for (fl_ in baseFiles){
      #print(fl_)
      idxs = which(fl_ == GrpFiles)
      idxsBase = which(fl_ == baseFiles)
      for (j in idxs){
        BaseMapMat[j, ] = CovBase$Cov[idxsBase, ]  
      }
    }
    
    CovNeutral = readMat(file.path(path, 'Neutral', 'Coverage.mat'))
    NeutMapMat = matrix(, nrow = length(GrpFiles), ncol = 4, byrow = TRUE)
    for (fl_ in NeutFiles){
      idxs = which(fl_ == GrpFiles)
      idxsNeut = which(fl_ == NeutFiles)
      for (j in idxs){
        NeutMapMat[j, ] = CovNeutral$Cov[idxsNeut, ]  
      }
    }
    
    for (i in c(1:4)){
      tRes = t.test(Cov$Cov[,i], BaseMapMat[,i], paired = TRUE, alternative='two.sided')
      outputBase[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
      outputBase[totalCnt, 2] = paste('Map-', as.character(i), sep = '')
      outputBase[totalCnt, 3] = round(tRes$statistic, 2)
      outputBase[totalCnt, 4] = as.integer(tRes$parameter)
      outputBase[totalCnt, 5] = round(tRes$p.value*24, 5)
      outputBase[totalCnt, 6] = round(tRes$estimate, 2)
      
      tRes = t.test(Cov$Cov[,i], NeutMapMat[,i], paired = TRUE, alternative='two.sided')
      outputNeut[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
      outputNeut[totalCnt, 2] = paste('Map-', as.character(i), sep = '')
      outputNeut[totalCnt, 3] = round(tRes$statistic, 2)
      outputNeut[totalCnt, 4] = as.integer(tRes$parameter)
      outputNeut[totalCnt, 5] = round(tRes$p.value*24, 5)
      outputNeut[totalCnt, 6] = round(tRes$estimate, 2)
      totalCnt = totalCnt + 1
    }  
  }
}

outputBase = data.frame(outputBase)
colnames(outputBase) <- c('Group', 'Map', 'Stats', 'Parameter', 'p.Value', 'meanDiff')
outputBase = outputBase[order(outputBase$p.Value), ] 
outputBase$adjPval = as.numeric(as.vector(outputBase$p.Value))*seq(32,1,by=-1)
write.csv(outputBase, file.path(path, 'Group-Base_MapCoverage.csv'))

outputNeut = data.frame(outputNeut)
colnames(outputNeut) <- c('Group', 'Map','Stats', 'Parameter', 'p.Value', 'meanDiff')
outputNeut = outputNeut[order(outputNeut$p.Value), ] 
outputNeut$adjPval = as.numeric(as.vector(outputNeut$p.Value))*seq(32,1,by=-1)
write.csv(outputNeut, file.path(path, 'Group-Neutral_MapCoverage.csv'))


####################################################################################################################


#### Global Field Potential
output <- matrix(ncol=5, nrow=32)
totalCnt = 1
for (path_ in allGrps){
  if (str_detect(strsplit(path_, split = '//')[[1]][2], 'Group')){
    #print(path_)
    Gfp = readMat(file.path(path_, 'GFP.mat'))
    
    ### Knowing the sequence of Subjects in GFP files
    fileGrp = readMat(file.path(path_, 'fileNames.mat'))
    GrpFiles = c()
    for (i in c(1:length(fileGrp$fileNames))){
      GrpFiles = c(GrpFiles, strsplit(fileGrp$fileNames[[i]][1][[1]][1], split = '_')[[1]][1])
    }
    
    GfpBase = readMat(file.path(path, 'Baseline', 'GFP.mat'))
    BaseMapMat = matrix(, nrow = length(GrpFiles), ncol = 4, byrow = TRUE)
    for (fl_ in baseFiles){
      #print(fl_)
      idxs = which(fl_ == GrpFiles)
      idxsBase = which(fl_ == baseFiles)
      for (j in idxs){
        BaseMapMat[j, ] = GfpBase$Gfp[idxsBase, ]  
      }
    }
    for (i in c(1, 4)){
      tRes = t.test(Gfp$Gfp[,i], BaseMapMat[,i], paired = TRUE, alternative='two.sided')
      #print(c(tRes$statistic, tRes$parameter, tRes$p.value, tRes$estimate))
      output[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
      output[totalCnt, 2] = paste(as.character(maxWhich), '-', as.character(i))
      output[totalCnt, 3] = round(tRes$statistic, 2)
      output[totalCnt, 4] = as.integer(tRes$parameter)
      output[totalCnt, 5] = round(tRes$p.value*24, 5)
      output[totalCnt, 6] = round(tRes$estimate, 2)
      totalCnt = totalCnt + 1
    }  
  }
}
output = data.frame(output)
colnames(output) <- c('Group','Stats', 'Parameter', 'p.Value', 'meanDiff')
write.csv(output, file.path(path, 'Group-Base_GFP.csv'))

#### Microstate Correlation
output <- matrix(ncol=5, nrow=32)
totalCnt = 1
for (path_ in allGrps){
  if (str_detect(strsplit(path_, split = '//')[[1]][2], 'Group')){
    #print(path_)
    Msp = readMat(file.path(path_, 'MspatCorr.mat'))
    
    ### Knowing the sequence of Subjects in GFP files
    fileGrp = readMat(file.path(path_, 'fileNames.mat'))
    GrpFiles = c()
    for (i in c(1:length(fileGrp$fileNames))){
      GrpFiles = c(GrpFiles, strsplit(fileGrp$fileNames[[i]][1][[1]][1], split = '_')[[1]][1])
    }
    
    MspBase = readMat(file.path(path, 'Baseline', 'MspatCorr.mat'))
    BaseMapMat = matrix(, nrow = length(GrpFiles), ncol = 4, byrow = TRUE)
    for (fl_ in baseFiles){
      #print(fl_)
      idxs = which(fl_ == GrpFiles)
      idxsBase = which(fl_ == baseFiles)
      for (j in idxs){
        BaseMapMat[j, ] = MspBase$Msp[idxsBase, ]  
      }
    }
    for (i in c(1:4)){
      tRes = t.test(Msp$Msp[,i], BaseMapMat[,i], paired = TRUE, alternative='two.sided')
      #print(c(tRes$statistic, tRes$parameter, tRes$p.value, tRes$estimate))
      output[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
      output[totalCnt, 2] = round(tRes$statistic, 2)
      output[totalCnt, 3] = as.integer(tRes$parameter)
      output[totalCnt, 4] = round(tRes$p.value*24, 5)
      output[totalCnt, 5] = round(tRes$estimate, 2)
      totalCnt = totalCnt + 1
    }  
  }
}
output = data.frame(output)
colnames(output) <- c('Group','Stats', 'Parameter', 'p.Value', 'meanDiff')
write.csv(output, file.path(path, 'Group-Base_MspatCorr.csv'))

#### Global Explained Variance
output <- matrix(ncol=5, nrow=32)
totalCnt = 1
for (path_ in allGrps){
  if (str_detect(strsplit(path_, split = '//')[[1]][2], 'Group')){
    #print(path_)
    Gev = readMat(file.path(path_, 'GEV.mat'))
    
    ### Knowing the sequence of Subjects in GFP files
    fileGrp = readMat(file.path(path_, 'fileNames.mat'))
    GrpFiles = c()
    for (i in c(1:length(fileGrp$fileNames))){
      GrpFiles = c(GrpFiles, strsplit(fileGrp$fileNames[[i]][1][[1]][1], split = '_')[[1]][1])
    }
    
    GevBase = readMat(file.path(path, 'Baseline', 'GEV.mat'))
    BaseMapMat = matrix(, nrow = length(GrpFiles), ncol = 4, byrow = TRUE)
    for (fl_ in baseFiles){
      #print(fl_)
      idxs = which(fl_ == GrpFiles)
      idxsBase = which(fl_ == baseFiles)
      for (j in idxs){
        BaseMapMat[j, ] = GevBase$Gev[idxsBase, ]  
      }
    }
    for (i in c(1:4)){
      tRes = t.test(Gev$Gev[,i], BaseMapMat[,i], paired = TRUE, alternative='two.sided')
      #print(c(tRes$statistic, tRes$parameter, tRes$p.value, tRes$estimate))
      output[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
      output[totalCnt, 2] = round(tRes$statistic, 2)
      output[totalCnt, 3] = as.integer(tRes$parameter)
      output[totalCnt, 4] = round(tRes$p.value*24, 5)
      output[totalCnt, 5] = round(tRes$estimate, 2)
      totalCnt = totalCnt + 1
    }  
  }
}
output = data.frame(output)
colnames(output) <- c('Group','Stats', 'Parameter', 'p.Value', 'meanDiff')
write.csv(output, file.path(path, 'Group-Base_GEV.csv'))

#### Microstate Durations
output <- matrix(ncol=5, nrow=32)
totalCnt = 1
for (path_ in allGrps){
  if (str_detect(strsplit(path_, split = '//')[[1]][2], 'Group')){
    #print(path_)
    Dur = readMat(file.path(path_, 'Duration.mat'))
    
    ### Knowing the sequence of Subjects in GFP files
    fileGrp = readMat(file.path(path_, 'fileNames.mat'))
    GrpFiles = c()
    for (i in c(1:length(fileGrp$fileNames))){
      GrpFiles = c(GrpFiles, strsplit(fileGrp$fileNames[[i]][1][[1]][1], split = '_')[[1]][1])
    }
    
    DurBase = readMat(file.path(path, 'Baseline', 'Duration.mat'))
    BaseMapMat = matrix(, nrow = length(GrpFiles), ncol = 4, byrow = TRUE)
    for (fl_ in baseFiles){
      #print(fl_)
      idxs = which(fl_ == GrpFiles)
      idxsBase = which(fl_ == baseFiles)
      for (j in idxs){
        BaseMapMat[j, ] = DurBase$Dur[idxsBase, ]  
      }
    }
    for (i in c(1:4)){
    tRes = t.test(Dur$Dur[,i], BaseMapMat[,i], paired = TRUE, alternative='two.sided')
    #print(c(tRes$statistic, tRes$parameter, tRes$p.value, tRes$estimate))
    output[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
    output[totalCnt, 2] = round(tRes$statistic, 2)
    output[totalCnt, 3] = as.integer(tRes$parameter)
    output[totalCnt, 4] = round(tRes$p.value*24, 5)
    output[totalCnt, 5] = round(tRes$estimate, 2)
    totalCnt = totalCnt + 1
    }  
  }
}
output = data.frame(output)
colnames(output) <- c('Group','Stats', 'Parameter', 'p.Value', 'meanDiff')
write.csv(output, file.path(path, 'Group-Base_Duration.csv'))

#### Microstate Coverage
output <- matrix(ncol=5, nrow=32)
totalCnt = 1
for (path_ in allGrps){
  if (str_detect(strsplit(path_, split = '//')[[1]][2], 'Group')){
    #print(path_)
    Cov = readMat(file.path(path_, 'Coverage.mat'))
    
    ### Knowing the sequence of Subjects in GFP files
    fileGrp = readMat(file.path(path_, 'fileNames.mat'))
    GrpFiles = c()
    for (i in c(1:length(fileGrp$fileNames))){
      GrpFiles = c(GrpFiles, strsplit(fileGrp$fileNames[[i]][1][[1]][1], split = '_')[[1]][1])
    }
    
    CovBase = readMat(file.path(path, 'Baseline', 'Coverage.mat'))
    BaseMapMat = matrix(, nrow = length(GrpFiles), ncol = 4, byrow = TRUE)
    for (fl_ in baseFiles){
      #print(fl_)
      idxs = which(fl_ == GrpFiles)
      idxsBase = which(fl_ == baseFiles)
      for (j in idxs){
        BaseMapMat[j, ] = CovBase$Cov[idxsBase, ]  
      }
    }
    for (i in c(1:4)){
      tRes = t.test(Cov$Cov[,i], BaseMapMat[,i], paired = TRUE, alternative='two.sided')
      #print(c(tRes$statistic, tRes$parameter, tRes$p.value, tRes$estimate))
      output[totalCnt, 1] = strsplit(path_, '//')[[1]][2]
      output[totalCnt, 2] = round(tRes$statistic, 2)
      output[totalCnt, 3] = as.integer(tRes$parameter)
      output[totalCnt, 4] = round(tRes$p.value*24, 5)
      output[totalCnt, 5] = round(tRes$estimate, 2)
      totalCnt = totalCnt + 1
    }  
  }
}
output = data.frame(output)
colnames(output) <- c('Group','Stats', 'Parameter', 'p.Value', 'meanDiff')
write.csv(output, file.path(path, 'Group-Base_Coverage.csv'))

########### Transition Matrix
library(ggplot2)
for (path_ in allGrps){
  #print(path_)
  TM = readMat(file.path(path_, 'TransitionMatrix.mat'))
  TransMat = TM$TM/(max(TM$TM))
  TransMat = data.frame(TransMat)
  X = c('Map1', 'Map2', 'Map3', 'Map4')
  Y = c('Map1', 'Map2', 'Map3', 'Map4')
  colnames(TransMat) = X
  rownames(TransMat) = Y
  heatData = expand.grid(Y=X,X=Y)
  heatData$val = as.vector(t(TransMat))
  title_ = strsplit(path_, "//")[[1]][2]
  plot_ = ggplot(heatData, aes(X, Y, fill=val)) + geom_tile() + ggtitle(title_) + xlab("") + ylab("") + 
    geom_point(aes(size = val)) + scale_fill_distiller(palette = "Spectral") + guides(size = "legend") +
    theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20), legend.text = element_text(size = 16))
  
  ggsave(file.path(strsplit(path_, "//")[[1]][1], paste(title_, ".png", sep = '')))
}
# axis.title = element_text(size = 20), 
# scale_fill_continuous(limits = c(0,1), breaks = c(0, 0.5, 0.7, 0.8, 0.9, 1.0),
#guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE))
