import scipy.io as io
import numpy as np
import glob
import os
import nibabel as nib
import pandas as pd
#from nilearn.plotting import plot_stat_map
import pdb

Dir_MS = "/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis/UpperBeta/MSsequence"
Dir_Source = "/data_hpc/home/ust/sudhakar/Processed_Emotions/python_files/EGI/brain/Raw_Data/sLoreta_20.0_30.0/emotion_wise/Upper_Beta"
Tar_MSSrc = "/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis/UpperBeta/Source_MSsequence"

def creatingImages(meanGrpMap, grpName):

    sourcePath = '/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis'
    targetPath = '/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis/UpperBeta/sourceImages'

    mni_template_file = os.path.join(sourcePath,'mni_icbm152_gm_tal_nlin_asym_09c.nii.gz')
    mni_template = nib.load(mni_template_file)
    data = mni_template.get_data()
    aft = mni_template.affine
    iaft = np.linalg.pinv(aft)

    csv_file = pd.read_csv(os.path.join(sourcePath, 'MNI-BAs-6239-voxels.csv'), header=None)
    all_indexes = csv_file.loc[:, [0,1,2]].values

    for mapId, map_ in enumerate(meanGrpMap):
        for idxCnt, xyz in enumerate(all_indexes):
            xyz = xyz.tolist()
            xyz.extend([1])
            voxel_image_index = iaft.dot(xyz)

            for k_ind in range(5):
                for j_ind in range(5):
                    for i_ind in range(5):
                        data[int(voxel_image_index[0])-2+i_ind,   # Sort out the problem of np.round. Use no.round instead of int.
                        int(voxel_image_index[1])-2+j_ind,
                        int(voxel_image_index[2])-2+k_ind] = map_[idxCnt]

                        '''t2_val = max_val+(float(csv_file['stats'][all_indexes[0][indexes]])/float(40))*(200-max_val)
                        data_t2[int(voxel_image_index_t2[0])-2+i_ind,
                        int(voxel_image_index_t2[1])-2+j_ind,
                        int(voxel_image_index_t2[2])-2+k_ind] = t2_val'''

        #pdb.set_trace()
        #mni_temp = nib.dataobj_images.DataobjImage(data, affine=mni_template.affine)#header = mni_template.header(), file_map=mni_template.file_map)
        #mni_temp = nib.spatialimages.SpatialImage(data, affine=mni_template.affine)
        niftiImg  = nib.Nifti1Image(data, mni_template.affine)
        file_name = os.path.join(targetPath, grpName+'_map-'+str(mapId)+'.nii.gz')
        nib.save(niftiImg, file_name)

def main():
    for MSG in glob.glob(os.path.join(Dir_MS, 'Group-05*')): #Baseline
        grpName = MSG.split('/')[-1]
        cnt = 1
        sourceGrp = os.path.join(Dir_Source, grpName)
        totalSamples = 0

        for file_ in glob.glob(os.path.join(MSG, '*.mat')):
            print(file_)
            allSplit = file_.split('/')[-1].split('_')

            if 'durationArr.mat' not in allSplit:
                if grpName == 'Baseline':
                    sourceFlName = grpName+'_'+allSplit[5]+'-slor.txt'
                else:
                    sourceFlName = allSplit[11].split('.mat')[0]+'_'+allSplit[4]+'Trial-'+allSplit[8].split('-')[1]+'Click-'+allSplit[10]+'-slor.txt'

                if allSplit[5] == 'mit017':
                    pdb.set_trace()

                if not os.path.exists(os.path.join(Tar_MSSrc, MSG.split('/')[-1], 'Source_'+sourceFlName.split('.txt')[0]+'.npy')):
                    
                    try:
                        slorData = np.loadtxt(os.path.join(sourceGrp,sourceFlName))
                        MsSeqData = io.loadmat(file_)['MSsequence'][0]
                    except:
                        print(os.path.join(sourceGrp,sourceFlName))
                        pdb.set_trace()
                        continue

                    ### Considering Microstates one by one
                    for i in np.arange(1, 5):
                        indexes = np.where(MsSeqData == i)[0]
                        if len(indexes):
                            averageMap = np.average(slorData[indexes, :], 0)
                            map_ = np.reshape(np.append(averageMap, i), (1, 6240))
                            if i == 1:
                                sourceArr = map_
                            else:
                                sourceArr = np.concatenate((sourceArr, map_), axis=0)

                    np.save(os.path.join(Tar_MSSrc, MSG.split('/')[-1], 'Source_'+sourceFlName.split('.txt')[0]+'.npy'), sourceArr)
                else:
                    totalSamples = totalSamples + 1
                    sourceArr = np.load(os.path.join(Tar_MSSrc, MSG.split('/')[-1], 'Source_'+sourceFlName.split('.txt')[0]+'.npy'))
                    sourceArr = sourceArr.round(3)
                    if cnt == 1:
                        meanGrpMap = sourceArr
                        cnt = 0
                    else:
                        try:
                            meanGrpMap = meanGrpMap + sourceArr
                        except:
                            newMapArr = np.zeros((4, 6240))
                            for indCnt, mapId in enumerate(sourceArr[:, -1]):
                                mapId = int(mapId)
                                if indCnt == mapId-1:
                                    newMapArr[mapId-1, :] = sourceArr[indCnt, :]
                                else:
                                    newMapArr[mapId-2, -1] = mapId-1
                                    newMapArr[mapId-1, :] = sourceArr[indCnt, :]

                            meanGrpMap = meanGrpMap + newMapArr

        meanGrpMap = meanGrpMap/totalSamples
        meanGrpMap = meanGrpMap.round(3)[:, :-1]
        creatingImages(meanGrpMap, grpName)

def statsImages(meanGrpMap, indexes, grpName):

    print(grpName)
    sourcePath = '/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis'
    targetPath = '/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis/UpperBeta/sourceImages'

    mni_template_file = os.path.join(sourcePath,'mni_icbm152_gm_tal_nlin_asym_09c.nii.gz')
    mni_template = nib.load(mni_template_file)
    data = mni_template.get_data()
    aft = mni_template.affine
    iaft = np.linalg.pinv(aft)

    csv_file = pd.read_csv(os.path.join(sourcePath, 'MNI-BAs-6239-voxels.csv'), header=None)
    all_indexes = csv_file.loc[:, [0,1,2]].values[indexes]

    for mapId, map_ in enumerate(meanGrpMap):
        for idxCnt, xyz in enumerate(all_indexes):
            xyz = xyz.tolist()
            xyz.extend([1])
            voxel_image_index = iaft.dot(xyz)

            for k_ind in range(5):
                for j_ind in range(5):
                    for i_ind in range(5):
                        data[int(voxel_image_index[0])-2+i_ind,   # Sort out the problem of np.round. Use no.round instead of int.
                        int(voxel_image_index[1])-2+j_ind,
                        int(voxel_image_index[2])-2+k_ind] = map_[idxCnt]

        niftiImg  = nib.Nifti1Image(data, mni_template.affine)
        file_name = os.path.join(targetPath, grpName+'_Statistical_Map-'+str(mapId)+'.nii.gz')
        nib.save(niftiImg, file_name)

def statsSigni():
    tarStats = '/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis/UpperBeta/MSStats'
    noPermut = 1000
    defaulter = 0

    for MSG in glob.glob(os.path.join(Dir_MS, 'Group*')): #Baseline
        grpName = MSG.split('/')[-1]
        totalSamples = 0

        for mapId in np.arange(4):
            if not os.path.isfile(os.path.join(tarStats, 'MicroMaps_MaxStatsPermutationDistributionHistogram_%s-Base_map-%s.npy' %(grpName, str(mapId)))):
                for file_ in glob.glob(os.path.join(MSG, '*.mat')):
                    print(file_)
                    allSplit = file_.split('/')[-1].split('_')
                    subName = file_.split('/')[-1].split('_')[4]

                    if 'durationArr.mat' not in allSplit:
                        if grpName == 'Baseline':
                            sourceFlName = grpName+'_'+allSplit[5]+'-slor.txt'
                        else:
                            sourceFlName = allSplit[11].split('.mat')[0]+'_'+allSplit[4]+'Trial-'+allSplit[8].split('-')[1]+'Click-'+allSplit[10]+'-slor.txt'

                        try:
                            sourceArr = np.load(os.path.join(Tar_MSSrc, MSG.split('/')[-1], 'Source_'+sourceFlName.split('.txt')[0]+'.npy'))
                            BaselnArr = np.load(glob.glob(os.path.join(Tar_MSSrc, 'Baseline', '*%s*.npy' %subName))[0])
                        except:
                            defaulter = defaulter + 1
                            continue

                        try:
                            TDiffData = sourceArr[mapId, :-1] - BaselnArr[mapId, :-1]
                            TDiffData = np.reshape(TDiffData, (1, TDiffData.shape[0]))
                        except:
                            #pdb.set_trace()
                            defaulter = defaulter + 1
                            continue

                        if totalSamples == 0:
                            sampleSub = TDiffData
                        else:
                            sampleSub = np.concatenate((sampleSub, TDiffData), axis=0)

                        totalSamples = totalSamples + 1

                MeanCalc = np.mean(sampleSub, axis=0)
                stdError = np.std(sampleSub, axis=0)/np.sqrt(totalSamples)
                tstats = MeanCalc/stdError
                np.save(os.path.join(tarStats, 'MicroMaps_orig_tstatistics_%s-Base_map-%s.npy' %(grpName, str(mapId))), tstats)

                ### Max Statistics ######
                actMaxStat = np.max(np.abs(tstats))
                np.save(os.path.join(tarStats, 'MicroMaps_ActualMaxStatistics_%s-Base_map-%s.npy' %(grpName, str(mapId))), actMaxStat)

                perMaxStat = []
                size_dt = sampleSub.shape 

                for shuffleCount in np.arange(noPermut):

                    #MeanCalcPerm = np.zeros((size_dt[1],size_dt[2]), np.float32)
                    noSampToFlip = np.random.randint(int(size_dt[0]/4),int(size_dt[0]/2))
                    print(shuffleCount, noSampToFlip)
                    samplToFlip = np.random.randint(0,size_dt[0],noSampToFlip)
                    forProcessing = sampleSub.copy()

                    flag = 0
                    for i in np.arange(size_dt[0]):
                        if i in samplToFlip:
                            if flag == 0:
                                print(i)
                                print(forProcessing[i, :])
                                flag = 1
                            forProcessing[i, :] = -forProcessing[i, :]

                    MeanCalcPerm = np.mean(forProcessing, axis=0)
                    stdErrorPerm = np.std(forProcessing, axis=0)/np.sqrt(forProcessing.shape[0])
                    tstatsPerm = MeanCalcPerm/stdErrorPerm
                    perMaxStat.append(np.max(np.abs(tstatsPerm)))
                    #perMaxStat.append(np.max(MeanCalcPerm))

                histVal = np.histogram(perMaxStat, bins=100)
                np.save(os.path.join(tarStats, 'MicroMaps_MaxStatsPermutationDistributionHistogram_%s-Base_map-%s.npy' %(grpName, str(mapId))), histVal)
                np.save(os.path.join(tarStats, 'MicroMaps_MaxStatsPermutationAArrays_%s-Base_map-%s.npy' %(grpName, str(mapId))), perMaxStat)

            else:
                histVal = np.load(os.path.join(tarStats, 'MicroMaps_MaxStatsPermutationDistributionHistogram_%s-Base_map-%s.npy' %(grpName, str(mapId))), allow_pickle=True)
                perMaxStat = np.load(os.path.join(tarStats, 'MicroMaps_MaxStatsPermutationAArrays_%s-Base_map-%s.npy' %(grpName, str(mapId))), allow_pickle=True)
                actMaxStat = np.load(os.path.join(tarStats, 'MicroMaps_ActualMaxStatistics_%s-Base_map-%s.npy' %(grpName, str(mapId))))

            idxs = np.where(histVal[1]>actMaxStat)
            #np.save(os.path.join(targetDirEmt, 'phaseRandomConnAnalEEG_ActualMaxStatistics_%s_seg-%s.npy' %(bandName[bandId], str(segId))), actMaxStat)

            if (len(idxs[0])):
                p_val = sum(histVal[0][idxs[0][0]:])/noPermut
            else:
                p_val = 1/noPermut

            np.save(os.path.join(tarStats, 'MicroMaps_PValue_%s-Base_map-%s.npy' %(grpName, str(mapId))), p_val)

            if p_val <= 0.001:
                tstats = np.load(os.path.join(tarStats, 'MicroMaps_orig_tstatistics_%s-Base_map-%s.npy' %(grpName, str(mapId))))
                indexes = np.where(tstats > 3.5)[0]
                statsImages(tstats, indexes, grpName)
    pdb.set_trace()

def plotting():
    sourcePath = '/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis/UpperBeta/sourceImages'
    target_folder = os.path.join(sourcePath, 'Images')

    anatomy_file = os.path.join('/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis','mni_icbm152_gm_tal_nlin_asym_09c.nii.gz')
    plot_stat_map(stat_map_img=anatomy_file, bg_img=anatomy_file,
                        output_file=os.path.join(target_folder, 'anatomyFile.png'),
                        display_mode='z', colorbar=True, figure=None,
                        axes=None, title=None, threshold=1e-06, annotate=True, draw_cross=True,
                        black_bg='auto', symmetric_cbar='auto', dim='auto', vmax=None)

    for files_ in glob.glob(os.path.join(sourcePath, '*nii*')):
        brain_volume = nib.load(files_)
        anatomy_file = os.path.join('/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis','mni_icbm152_gm_tal_nlin_asym_09c.nii.gz')

        plot_stat_map(stat_map_img=brain_volume, bg_img=anatomy_file,
                            output_file=os.path.join(target_folder, files_.split('/')[-1].split('.nii')[0]+'.png'),
                            display_mode='z', colorbar=True, figure=None,
                            axes=None, title=None, threshold=1e-06, annotate=True, draw_cross=True,
                            black_bg='auto', symmetric_cbar='auto', dim='auto', vmax=None)
