library(reticulate)
np <- import("numpy")

statsPath = '/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis/UpperBeta/Microstate_Stats/MSStats'

for (mapId in c(0:3)){
  for (grp_ in c('Group-01', 'Group-02', 'Group-03', 'Group-04', 'Group-05', 'Group-06', 'Group-07', 'Group-08')) {
    p_Val = np$load(file.path(statsPath, paste('MicroMaps_PValue_', grp_, '-Base_map-', as.character(mapId), '.npy', sep = '')))
    print(p_Val)
    tstats = qt(1-p_Val, 50)
    if (p_Val <= 0.001){
      origStats = np$load(file.path(statsPath, paste('MicroMaps_orig_tstatistics_', grp_, '-Base_map-', as.character(mapId), '.npy', sep = '')))
      signiVoxels = which(abs(origStats)>tstats[1])
      file.path(statsPath, paste('MicroMaps_Sign_Voxels_', grp_, '-Base_map-', as.character(mapId), '.npy', sep = ''))
    }
  }
}
