import os
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import glob
import scipy.io
from scipy.io import loadmat
import pandas as pd
import seaborn as sns
import scipy
import pdb

#frqBands = ['AllBands', 'Delta', 'Theta', 'Alpha', 'LowerBeta', 'UpperBeta', 'Gamma']
frqBands = ['UpperBeta']
#Groups = ['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08']

#sourceDir = '/media/silp/7c7fe885-7b62-45e2-a4d9-e64f37f48bef/ForMSAAnalysis'
sourceDir = '/media/sudhakar/4080C90A80C9077E/698GBDrive/ForMSAAnalysis'
############# Classification Based on Cumulants
#subjectsForTesting = ['mit072', 'mit082', 'mit099', 'mit114']
subjectsForTesting = ['mit114']
groupDir = {}
noG = 2 
#groupdict_ = {'Baseline':0, 'Group-01':1, 'Group-02':2, 'Group-03':3, 'Group-04':4, 'Group-05':5, 'Group-06':6, 'Group-07':7, 'Group-08':8}
#groupdict_ = {'Baseline':0, 'Group-01':1, 'Group-02':1, 'Group-03':2, 'Group-04':3, 'Group-05':2, 'Group-06':3, 'Group-07':4, 'Group-08':4}
groupdict_ = {'Group-01':1, 'Group-02':2, 'Group-03':3, 'Group-04':4, 'Group-05':5, 'Group-06':6, 'Group-07':7, 'Group-08':8}
#groupdict_ = {'Group-01':1, 'Group-02':1, 'Group-03':2, 'Group-04':2, 'Group-05':2, 'Group-06':2, 'Group-07':3, 'Group-08':3}
#groupdict_ = {'Group-01':1, 'Group-02':2, 'Group-03':2, 'Group-04':2, 'Group-05':3, 'Group-06':2, 'Group-07':3, 'Group-08':4}

allGroups = ['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08']
RatioDf = pd.DataFrame(0, index=allGroups, columns=allGroups)
df_accKNN = pd.DataFrame(0, index=allGroups, columns=allGroups)
df_accSVMSig = pd.DataFrame(0, index=allGroups, columns=allGroups)
df_accSVMRbf = pd.DataFrame(0, index=allGroups, columns=allGroups)

for id1, gp1 in enumerate(allGroups):
    for gp2 in allGroups[id1+1:]:
        print('===================== New Pair ======================')
        print(gp1, gp2)
        for ii in np.arange(2,3):
            if os.path.isfile(os.path.join(sourceDir, 'ScalingExponent', 'ScalingExp_db-%s.png' %str(ii))):
                flag_ = 0
                subDir = {}
                for frq_ in frqBands:
                    grpSort = np.sort(glob.glob(os.path.join(sourceDir, frq_, '*')))
                    for grp_ in grpSort:
                        #if '.png' in grp_.split('/')[-1] or 'Baseline' in grp_ or 'Group-03' in grp_ or 'Group-04' in grp_ or 'Group-05' in grp_ or 'Group-06' in grp_:
                        #if '.png' in grp_.split('/')[-1] or 'Baseline' in grp_ or 'Group-03' in grp_ or 'Group-04' in grp_ or 'Group-06' in grp_:
                        #if '.png' in grp_.split('/')[-1] or 'Baseline' in grp_ or 'Group-03' in grp_ or 'Group-04' in grp_ :
                        if '.png' in grp_.split('/')[-1]:# or 'Baseline' in grp_ or 'Group-03' in grp_ or 'Group-04' in grp_ or 'Group-05' in grp_  or 'Group-06' in grp_:# or 'Group-07' in grp_:
                            continue

                        '''if gp2=='Group-03':
                            pdb.set_trace()'''

                        if not((gp1 in grp_) or (gp2 in grp_)):
                            continue

                        cumulants = scipy.io.loadmat(os.path.join(grp_, 'db%s_cumulants.mat' %str(ii)))                        
                        keys_ = [i for i in cumulants.keys()]
                        cumulant = cumulants[keys_[-1]]
                        Meancumulant = np.median(cumulant, 0)
                        Meancumulant = [np.round(i,2) for i in Meancumulant]
                        Eq_cumulants = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_cumulants.mat' %str(ii)))
                        keys_ = [i for i in Eq_cumulants.keys()]
                        Eq_cumulant = Eq_cumulants[keys_[-1]]
                        Eq_Meancumulant = np.median(Eq_cumulant, 0)
                        Eq_Meancumulant = [np.round(i,2) for i in Eq_Meancumulant]
                        Sh_cumulants = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_cumulants.mat' %str(ii)))
                        keys_ = [i for i in Sh_cumulants.keys()]
                        Sh_cumulant = Sh_cumulants[keys_[-1]]
                        Sh_Meancumulant = np.median(Sh_cumulant, 0)
                        Sh_Meancumulant = [np.round(i,2) for i in Sh_Meancumulant]
                        grpName = grp_.split('/')[-1].split('_')[0]
                        tar_ = groupdict_[grpName]
                        with open(os.path.join(grp_, 'subjectsIncluded.txt')) as fl_:
                            subjectsInc = fl_.readlines()
                        for _id, sub_ in enumerate(subjectsInc):
                            if _id > len(cumulant)/3: ## Taking only one kind of microstate label
                                break

                            subIdx_ = np.where(np.array(['mit' in sub_[i]+sub_[i+1]+sub_[i+2] for i in np.arange(len(sub_)-3)])==True)[0][0]
                            upper = subIdx_+6

                            featTar = np.append(cumulant[_id,:],tar_)
                            #if sub_[subIdx_:upper] not in subjectsForTesting:
                            toAdd = np.reshape(featTar, (1,4))
                            if sub_[subIdx_:upper] not in subDir:
                                subDir[sub_[subIdx_:upper]] = toAdd
                            else:
                                subDir[sub_[subIdx_:upper]] = np.concatenate((subDir[sub_[subIdx_:upper]], toAdd), axis=0)
                            #else:

                            if sub_[subIdx_:upper] not in groupDir:
                                groupDir[sub_[subIdx_:upper]] = [grpName]
                            else:
                                groupDir[sub_[subIdx_:upper]] = np.append(groupDir[sub_[subIdx_:upper]], grpName)
                            groupDir[sub_[subIdx_:upper]] = np.unique(groupDir[sub_[subIdx_:upper]])

                del cumulant
                del Meancumulant
                del Eq_cumulant
                del Eq_Meancumulant
                del Sh_cumulant
                del Sh_Meancumulant

                trainFlag = 0
                testFlag = 0
                for sub_ in subDir.keys():
                    if sub_ not in subjectsForTesting:
                        if trainFlag == 0:
                            Xtrain = subDir[sub_][:, 0:2]
                            Ytrain = subDir[sub_][:, 3]
                            trainFlag = 1
                        else:
                            Xtrain = np.concatenate((Xtrain, subDir[sub_][:, 0:2]), axis=0)
                            Ytrain = np.append(Ytrain, subDir[sub_][:, 3])
                    else:
                        if testFlag == 0:
                            Xtest = subDir[sub_][:, 0:2]
                            Ytest = subDir[sub_][:, 3]
                            testFlag = 1
                        else:
                            Xtest = np.concatenate((Xtest, subDir[sub_][:, 0:2]), axis=0)
                            Ytest = np.append(Ytest, subDir[sub_][:, 3])

                del subDir
                #del groupDir
                del featTar
                '''for g_ in np.arange(1, noG+1):
                    print('Train-%s' %str(len(np.where(Ytrain==g_)[0])))
                    print('Test-%s' %str(len(np.where(Ytest==g_)[0])))'''
                print('Train-%s' %str(len(np.where(Ytrain==groupdict_[gp1])[0])))
                print('Test-%s' %str(len(np.where(Ytest==groupdict_[gp1])[0])))
                print('Train-%s' %str(len(np.where(Ytrain==groupdict_[gp2])[0])))
                print('Test-%s' %str(len(np.where(Ytest==groupdict_[gp2])[0])))
                ratio = len(np.where(Ytest==groupdict_[gp1])[0])/len(np.where(Ytest==groupdict_[gp2])[0])
                RatioDf.loc[gp1, gp2] = ratio

                np.where(Ytrain == 2)[0]
                from sklearn.neighbors import KNeighborsClassifier
                nbrs = KNeighborsClassifier(n_neighbors=8, algorithm='ball_tree')
                nbrs.fit(Xtrain, Ytrain)
                predLabels = nbrs.predict(Xtest)
                diff_= abs(predLabels-Ytest)
                false_ = len(np.where(diff_!=0)[0])
                true_ = 100-((false_*100)/len(diff_))
                print(true_, false_)
                df_accKNN.loc[gp1, gp2] = true_
                ### SVM
                from sklearn import svm
                clf = svm.SVC(kernel='sigmoid', max_iter=-1, tol=1e-5)
                clf.fit(Xtrain, Ytrain)
                pred_labels = clf.predict(Xtest)
                diff_= abs(pred_labels-Ytest)
                false_ = len(np.where(diff_!=0)[0])
                true_ = 100-((false_*100)/len(diff_))
                print(true_, false_)
                df_accSVMSig.loc[gp1, gp2] = true_
                '''clf = svm.SVC(kernel='poly', degree=5.0, max_iter=-1, tol=1e-5)
                clf.fit(Xtrain, Ytrain)
                pred_labels = clf.predict(Xtest)
                diff_= abs(pred_labels-Ytest)
                false_ = len(np.where(diff_!=0)[0])
                true_ = 100-((false_*100)/len(diff_))
                print(true_, false_)'''

                clf = svm.SVC(kernel='rbf', max_iter=-1, tol=1e-5)
                clf.fit(Xtrain, Ytrain)
                pred_labels = clf.predict(Xtest)
                diff_= abs(pred_labels-Ytest)
                false_ = len(np.where(diff_!=0)[0])
                true_ = 100-((false_*100)/len(diff_))
                print(true_, false_)
                df_accSVMRbf.loc[gp1, gp2] = true_

if not os.path.exists(os.path.join(sourceDir, 'Singularity')):
    os.makedirs(os.path.join(sourceDir, 'Singularity'))

for ii in np.arange(2,6):
    #if not os.path.isfile(os.path.join(sourceDir, 'Singularity', 'Singularity_db-%s.png' %str(ii))):
    if not os.path.isfile(os.path.join(sourceDir, 'Singularity', 'Singularity_db-%s.png' %str(ii))):
        Or_hPeak = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Or_Width = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Or_Asy = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])

        Eq_hPeak = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Eq_Width = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Eq_Asy = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Sh_hPeak = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Sh_Width = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Sh_Asy = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])

        fig = plt.figure(figsize=(20, 10.5))
        nrow = 7
        ncol = 9
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), hspace=0.3, left=0.04, right=0.99, bottom=0.05, top=0.93)
        row_ = 0
        col_ = 0
        for frq_ in frqBands:
            grpSort = np.sort(glob.glob(os.path.join(sourceDir, frq_, '*')))            
            for grpIdx, grp_ in enumerate(grpSort[1:]):
                if '.png' in grp_.split('/')[-1]:
                    continue
                
                dhFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityDim.mat' %str(ii)))
                Eq_dhFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_singularityDim.mat' %str(ii)))
                Sh_dhFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_singularityDim.mat' %str(ii)))
                keys_ = [i for i in dhFiles.keys()]
                dhData = dhFiles[keys_[-1]] ## Normal Files
                '''Taking median here is right because we want to get a representative statistics 
                for the distribution with respect to each moment. But, in the case where we want 
                to do the statistical test of value of h where the spectrum got its peak value 
                there the peakIndex for each row (data line) should be calculated. As done in for creating PeakValDist_db-X.png'''

                dhMean = np.median(dhData,0) 
                Eq_dhData = Eq_dhFiles[keys_[-1]] ## Equalized Time Files
                Eq_dhMean = np.median(Eq_dhData,0)## This is right because we want to get a representative statistics for the distribution with respect to each moment.
                Sh_dhData = Sh_dhFiles[keys_[-1]] ## Shuffled Label File
                Sh_dhMean = np.median(Sh_dhData,0)## This is right because we want to get a representative statistics for the distribution with respect to each moment.
                peakIndexMedian = np.argmax(dhMean)
                Eq_peakIndexMedian = np.argmax(Eq_dhMean)
                Sh_peakIndexMedian = np.argmax(Sh_dhMean)
                '''There was inherent assumption here that 1 is the 
                largest value and what we get is the larget value if 
                compare with 0.999 which is clearly not the case 
                because 1 is not the maximum value.
                modeIndex = np.where(dhMean>0.999)[0][0]. Hence the implementation of peakIndex'''

                hFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityExponent.mat' %str(ii)))
                Eq_hFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_singularityExponent.mat' %str(ii)))
                Sh_hFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_singularityExponent.mat' %str(ii)))
                keys_ = [i for i in hFiles.keys()]
                hData = hFiles[keys_[-1]] ## Normal Files
                hMean = np.median(hData,0)
                hModeVal = np.round(hMean[peakIndexMedian],2)
                hmin = np.min(hMean)
                hmax = np.max(hMean)
                hpeak = hModeVal
                DeltaL = hpeak-hmin
                DeltaR = hmax-hpeak
                A = np.round((DeltaL-DeltaR)/(DeltaL+DeltaR),2)

                print(np.round(scipy.stats.ttest_rel(hData[0:28, 0], hData[28:56, 0]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 1], hData[28:56, 1]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 2], hData[28:56, 2]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 3], hData[28:56, 3]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 4], hData[28:56, 4]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 5], hData[28:56, 5]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 6], hData[28:56, 6]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 7], hData[28:56, 7]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 8], hData[28:56, 8]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 9], hData[28:56, 9]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 10], hData[28:56, 10]), 2))

                print(np.round(scipy.stats.ttest_rel(hData[28:56, 0], hData[56:84, 0]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[28:56, 1], hData[56:84, 1]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[28:56, 2], hData[56:84, 2]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[28:56, 3], hData[56:84, 3]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[28:56, 4], hData[56:84, 4]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[28:56, 5], hData[56:84, 5]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[28:56, 6], hData[56:84, 6]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[28:56, 7], hData[56:84, 7]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[28:56, 8], hData[56:84, 8]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[28:56, 9], hData[56:84, 9]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[28:56, 10], hData[56:84, 10]), 2))

                print(np.round(scipy.stats.ttest_rel(hData[0:28, 0], hData[56:84, 0]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 1], hData[56:84, 1]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 2], hData[56:84, 2]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 3], hData[56:84, 3]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 4], hData[56:84, 4]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 5], hData[56:84, 5]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 6], hData[56:84, 6]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 7], hData[56:84, 7]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 8], hData[56:84, 8]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 9], hData[56:84, 9]), 2))
                print(np.round(scipy.stats.ttest_rel(hData[0:28, 10], hData[56:84, 10]), 2))
                
                width_ = np.round(np.quantile(hMean, 0.75)-np.quantile(hMean, 0.25),2)

                Eq_hData = Eq_hFiles[keys_[-1]] ## Equalized Time Files
                Eq_hMean = np.median(Eq_hData,0)
                Eq_hModeVal = np.round(Eq_hMean[Eq_peakIndexMedian],2)
                hmin = np.min(Eq_hMean)
                hmax = np.max(Eq_hMean)
                hpeak = Eq_hModeVal
                DeltaL = hpeak-hmin
                DeltaR = hmax-hpeak
                Eq_A = np.round((DeltaL-DeltaR)/(DeltaL+DeltaR),2)
                Eq_width_ = np.round(np.quantile(Eq_hMean, 0.75)-np.quantile(Eq_hMean, 0.25),2)

                Sh_hData = Sh_hFiles[keys_[-1]] ## Shuffled Label File
                Sh_hMean = np.median(Sh_hData,0)
                Sh_hModeVal = np.round(Sh_hMean[Sh_peakIndexMedian],2)
                hmin = np.min(Sh_hMean)
                hmax = np.max(Sh_hMean)
                hpeak = Sh_hModeVal
                DeltaL = hpeak-hmin
                DeltaR = hmax-hpeak

                Sh_A = np.round((DeltaL-DeltaR)/(DeltaL+DeltaR),2)
                Sh_width_ = np.round(np.quantile(Sh_hMean, 0.75)-np.quantile(Sh_hMean, 0.25),2)

                cumulants = scipy.io.loadmat(os.path.join(grp_, 'db%s_cumulants.mat' %str(ii)))
                keys_ = [i for i in cumulants.keys()]
                cumulant = cumulants[keys_[-1]]
                Meancumulant = np.median(cumulant, 0)

                Or_hPeak.loc[frq_, grp_.split('/')[-1].split('_')[0]] = hModeVal
                Or_Width.loc[frq_, grp_.split('/')[-1].split('_')[0]] = width_
                Or_Asy.loc[frq_, grp_.split('/')[-1].split('_')[0]] = A
                Eq_hPeak.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Eq_hModeVal
                Eq_Width.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Eq_width_
                Eq_Asy.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Eq_A
                Sh_hPeak.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Sh_hModeVal
                Sh_Width.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Sh_width_
                Sh_Asy.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Sh_A

                ax = plt.subplot(gs[row_,col_])
                '''ax.plot(hMean, dhMean, 'r-', label='(%s,%s,%s)' %(str(hModeVal),str(width_),str(A)))
                ax.plot(Eq_hMean, Eq_dhMean, 'g--', label='(%s,%s,%s)' %(str(Eq_hModeVal),str(Eq_width_),str(Eq_A)))
                ax.plot(Sh_hMean, Sh_dhMean, 'b--', label='(%s,%s,%s)' %(str(Sh_hModeVal),str(Sh_width_),str(Sh_A)))'''

                from scipy.interpolate import make_interp_spline, BSpline
                # making the plot smooth here.
                #xnew = np.linspace(T.min(), T.max(), 300)  
                
                startFrom = np.where(dhMean<0)[0][-1]+1
                sh_startFrom = np.where(Sh_dhMean<0)[0][-1]+1

                ax.plot(hMean[startFrom:], dhMean[startFrom:], 'r-', label='(%s,%s)' %(str(hModeVal),str(width_)), linewidth=2)
                ax.plot(Eq_hMean, Eq_dhMean, 'g--', label='(%s,%s)' %(str(Eq_hModeVal),str(Eq_width_)), linewidth=2)
                ax.plot(Sh_hMean[sh_startFrom:], Sh_dhMean[sh_startFrom:], 'b.--', label='(%s,%s)' %(str(Sh_hModeVal),str(Sh_width_)), linewidth=2)

                #ax.legend(loc=1, prop={'size': 8})
                ax.legend(prop={'size': 10})
                #ax.text(np.min(hMean), np.min(dhMean), 'M-%s,%s,%s\nW-%s,%s,%s' %(str(hModeVal),str(Eq_hModeVal),str(Sh_hModeVal), str(width_),str(Eq_width_),str(Sh_width_)))            
                #ax.scatter(hMean, dhMean, s=5)
                if col_ == 0:
                    ax.set_ylabel(frq_)
                if row_ == 0:
                    ax.set_title(grp_.split('/')[-1].split('_')[0])
                col_ = col_ + 1
            row_ = row_ + 1
            col_ = 0

        Or_hPeak.to_csv('Original_hPeak-db%s.csv' %str(ii))
        Or_Width.to_csv('Original_width-db%s.csv' %str(ii))
        Or_Asy.to_csv('Original_Asymmetry-db%s.csv' %str(ii))
        Eq_hPeak.to_csv('Equalized_hPeak-db%s.csv' %str(ii))
        Eq_Width.to_csv('Equalized_width-db%s.csv' %str(ii))
        Eq_Asy.to_csv('Equalized_Asy-db%s.csv' %str(ii))
        Sh_hPeak.to_csv('Shuffled_hPeak-db%s.csv' %str(ii))
        Sh_Width.to_csv('Shuffled_width-db%s.csv' %str(ii))
        Sh_Asy.to_csv('Shuffled_Asy-db%s.csv' %str(ii))
        plt.suptitle('db-%s' %str(ii))
        pdb.set_trace()
        plt.savefig(os.path.join(sourceDir, 'Singularity', 'Singularity_db-%s.png' %str(ii)))
        plt.savefig(os.path.join(sourceDir, 'Singularity', 'Singularity_db-%s.pdf' %str(ii)))
        plt.close()
        plt.clf()

    if os.path.isfile(os.path.join(sourceDir, 'Singularity', 'UpperBeta_Singularity_db-%s.png' %str(ii))):
        frq_ = 'UpperBeta'
        plt.figure(figsize=(20, 10.5))
        #plt.rcParams['font.size'] = 30
        nrow = 3
        ncol = 3
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), wspace=0.2, hspace=0.4, left=0.04, right=0.99, bottom=0.05, top=0.93)
        row_ = 0
        col_ = 0
        grpSort = np.sort(glob.glob(os.path.join(sourceDir, frq_, '*')))
        for grpIdx, grp_ in enumerate(grpSort[1:]):
            if '.png' in grp_.split('/')[-1]:
                continue

            dhFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityDim.mat' %str(ii)))
            Eq_dhFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_singularityDim.mat' %str(ii)))
            Sh_dhFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_singularityDim.mat' %str(ii)))
            keys_ = [i for i in dhFiles.keys()]
            dhData = dhFiles[keys_[-1]] ## Normal Files
            '''Taking median here is right because we want to get a representative statistics 
            for the distribution with respect to each moment. But, in the case where we want 
            to do the statistical test of value of h where the spectrum got its peak value 
            there the peakIndex for each row (data line) should be calculated. As done in for creating PeakValDist_db-X.png'''
            dhMean = np.median(dhData,0) 
            Eq_dhData = Eq_dhFiles[keys_[-1]] ## Equalized Time Files
            Eq_dhMean = np.median(Eq_dhData,0)## This is right because we want to get a representative statistics for the distribution with respect to each moment.
            Sh_dhData = Sh_dhFiles[keys_[-1]] ## Shuffled Label File
            Sh_dhMean = np.median(Sh_dhData,0)## This is right because we want to get a representative statistics for the distribution with respect to each moment.
            peakIndexMedian = np.argmax(dhMean)
            Eq_peakIndexMedian = np.argmax(Eq_dhMean)
            Sh_peakIndexMedian = np.argmax(Sh_dhMean)
            '''There was inherent assumption here that 1 is the 
            largest value and what we get is the larget value if 
            compare with 0.999 which is clearly not the case 
            because 1 is not the maximum value.
            modeIndex = np.where(dhMean>0.999)[0][0]. Hence the implementation of peakIndex'''

            hFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityExponent.mat' %str(ii)))
            Eq_hFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_singularityExponent.mat' %str(ii)))
            Sh_hFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_singularityExponent.mat' %str(ii)))
            keys_ = [i for i in hFiles.keys()]
            hData = hFiles[keys_[-1]] ## Normal Files
            hMean = np.median(hData,0)
            hModeVal = np.round(hMean[peakIndexMedian],2)
            hmin = np.min(hMean)
            hmax = np.max(hMean)
            hpeak = hModeVal
            DeltaL = hpeak-hmin
            DeltaR = hmax-hpeak
            A = np.round((DeltaL-DeltaR)/(DeltaL+DeltaR),2)
            startFrom = np.where(dhMean<0)[0][-1]+1            
            #width_ = np.round(np.quantile(hMean, 0.75)-np.quantile(hMean, 0.25),2)
            width_ = np.round(hMean[startFrom]-hMean[-1], 2)

            Eq_hData = Eq_hFiles[keys_[-1]] ## Equalized Time Files
            Eq_hMean = np.median(Eq_hData,0)
            Eq_hModeVal = np.round(Eq_hMean[Eq_peakIndexMedian],2)
            hmin = np.min(Eq_hMean)
            hmax = np.max(Eq_hMean)
            hpeak = Eq_hModeVal
            DeltaL = hpeak-hmin
            DeltaR = hmax-hpeak
            Eq_A = np.round((DeltaL-DeltaR)/(DeltaL+DeltaR),2)
            Eq_width_ = np.round(np.quantile(Eq_hMean, 0.75)-np.quantile(Eq_hMean, 0.25),2)

            Sh_hData = Sh_hFiles[keys_[-1]] ## Shuffled Label File
            Sh_hMean = np.median(Sh_hData,0)
            Sh_hModeVal = np.round(Sh_hMean[Sh_peakIndexMedian],2)
            hmin = np.min(Sh_hMean)
            hmax = np.max(Sh_hMean)
            hpeak = Sh_hModeVal
            DeltaL = hpeak-hmin
            DeltaR = hmax-hpeak

            Sh_A = np.round((DeltaL-DeltaR)/(DeltaL+DeltaR),2)
            sh_startFrom = np.where(Sh_dhMean<0)[0][-1]+1
            #Sh_width_ = np.round(np.quantile(Sh_hMean, 0.75)-np.quantile(Sh_hMean, 0.25),2)
            Sh_width_ = np.round(Sh_hMean[sh_startFrom]-Sh_hMean[-1], 2)

            cumulants = scipy.io.loadmat(os.path.join(grp_, 'db%s_cumulants.mat' %str(ii)))
            keys_ = [i for i in cumulants.keys()]
            cumulant = cumulants[keys_[-1]]
            Meancumulant = np.median(cumulant, 0)

            '''Or_hPeak.loc[frq_, grp_.split('/')[-1].split('_')[0]] = hModeVal
            Or_Width.loc[frq_, grp_.split('/')[-1].split('_')[0]] = width_
            Or_Asy.loc[frq_, grp_.split('/')[-1].split('_')[0]] = A
            Eq_hPeak.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Eq_hModeVal
            Eq_Width.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Eq_width_
            Eq_Asy.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Eq_A
            Sh_hPeak.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Sh_hModeVal
            Sh_Width.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Sh_width_
            Sh_Asy.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Sh_A'''

            ax = plt.subplot(gs[row_,col_])
            ax.plot(hMean[startFrom:], dhMean[startFrom:], 'r-', label='(%s,%s)' %(str(hModeVal),str(width_)), linewidth=4)
            ax.plot(Eq_hMean, Eq_dhMean, 'g--', label='(%s,%s)' %(str(Eq_hModeVal),str(Eq_width_)), linewidth=4)
            ax.plot(Sh_hMean[sh_startFrom:], Sh_dhMean[sh_startFrom:], 'b.--', label='(%s,%s)' %(str(Sh_hModeVal),str(Sh_width_)), linewidth=4)
            ax.legend(fontsize=13, fancybox=False, framealpha=0.1, ncol=3, title='h_peak, width', columnspacing=0.5, markerscale=0.5, title_fontsize=15)
            #ax.legend(prop={'size': 25})
            #ax.text(np.min(hMean), np.min(dhMean), 'M-%s,%s,%s\nW-%s,%s,%s' %(str(hModeVal),str(Eq_hModeVal),str(Sh_hModeVal), str(width_),str(Eq_width_),str(Sh_width_)))            
            #ax.scatter(hMean, dhMean, s=5)
            col_ = col_ + 1

            ax.set_title(grp_.split('/')[-1].split('_')[0], fontsize=25)
            if (grpIdx + 1) %3 == 0:
                col_ = 0
                row_ = row_ + 1
            ax.tick_params(axis='both', which='both', labelsize=25)

        #plt.suptitle('db-%s' %str(ii))                 
        plt.savefig(os.path.join(sourceDir, 'Singularity', '%s_Singularity_db-%s.png' %(frq_, str(ii))), bbox_inches='tight')
        plt.savefig(os.path.join(sourceDir, 'Singularity', '%s_Singularity_db-%s.pdf' %(frq_, str(ii))), bbox_inches='tight')
        plt.close()
        plt.clf()

###################### Spectrum Width Comparison ###############

#for ii in np.arange(2,6):
ii = 2
for grp in ['Baseline', 'Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08']:        
    grp_ = glob.glob(os.path.join(sourceDir, 'UpperBeta', grp+'*'))[0]
    if '.png' in grp_.split('/')[-1]:
        continue

    hFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityExponent.mat' %str(ii)))
    Eq_hFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_singularityExponent.mat' %str(ii)))
    Sh_hFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_singularityExponent.mat' %str(ii)))

    tmp_ = scipy.io.loadmat(os.path.join(grp_, 'filesProcessed.mat'))
    UB_hSubFiles = [j[0][0] for j in tmp_[[i for i in tmp_.keys()][-1]]]
    tmp_ = scipy.io.loadmat(os.path.join(grp_, 'EQ_filesProcessed.mat'))
    UB_Eq_hSubFiles = [j[0][0] for j in tmp_[[i for i in tmp_.keys()][-1]]    ]
    tmp_ = scipy.io.loadmat(os.path.join(grp_, 'SH_filesProcessed.mat'))
    UB_Sh_hSubFiles = [j[0][0] for j in tmp_[[i for i in tmp_.keys()][-1]]    ]

    keys_ = [i for i in hFiles.keys()]

    hData = hFiles[keys_[-1]] ## Normal Files        
    UB_widthArr = np.array([np.quantile(i, 0.75) for i in hData])-np.array([np.quantile(i, 0.25) for i in hData])

    Eq_hData = Eq_hFiles[keys_[-1]] ## Equalized Time Files        
    UB_Eq_widthArr = np.array([np.quantile(i, 0.75) for i in Eq_hData])-np.array([np.quantile(i, 0.25) for i in Eq_hData])

    Sh_hData = Sh_hFiles[keys_[-1]] ## Equalized Time Files        
    UB_Sh_widthArr = np.array([np.quantile(i, 0.75) for i in Sh_hData])-np.array([np.quantile(i, 0.25) for i in Sh_hData])
    
    for frq_ in frqBands:

        if os.path.isfile(os.path.join(sourceDir, 'WidthDifference', grp+'_upperBeta-'+frq_+'-Orig.csv')):
            continue

        if frq_ == 'AllBands':
            continue

        grp_ = glob.glob(os.path.join(sourceDir, frq_, grp+'*'))[0]
        if '.png' in grp_.split('/')[-1]:
            continue

        hFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityExponent.mat' %str(ii)))
        Eq_hFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_singularityExponent.mat' %str(ii)))
        Sh_hFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_singularityExponent.mat' %str(ii)))

        tmp_ = scipy.io.loadmat(os.path.join(grp_, 'filesProcessed.mat'))
        hSubFiles = [j[0][0] for j in tmp_[[i for i in tmp_.keys()][-1]]]
        tmp_ = scipy.io.loadmat(os.path.join(grp_, 'EQ_filesProcessed.mat'))
        Eq_hSubFiles = [j[0][0] for j in tmp_[[i for i in tmp_.keys()][-1]]    ]
        tmp_ = scipy.io.loadmat(os.path.join(grp_, 'SH_filesProcessed.mat'))
        Sh_hSubFiles = [j[0][0] for j in tmp_[[i for i in tmp_.keys()][-1]]    ]

        if frq_ == 'Theta' and grp == 'Baseline':            
            hSubFiles = [j.split('corrected')[0]+'corrected_7Sec_'+'_'.join(j.split('_')[4:]) for j in hSubFiles] 
            Eq_hSubFiles = [j.split('corrected')[0]+'corrected_7Sec_'+'_'.join(j.split('_')[5:]) for j in Eq_hSubFiles] 
            Sh_hSubFiles = [j.split('corrected')[0]+'corrected_7Sec_'+'_'.join(j.split('_')[5:]) for j in Sh_hSubFiles] 

        keys_ = [i for i in hFiles.keys()]        
        hData = hFiles[keys_[-1]] ## Normal Files        
        widthArr = np.array([np.quantile(i, 0.75) for i in hData])-np.array([np.quantile(i, 0.25) for i in hData])

        Eq_hData = Eq_hFiles[keys_[-1]] ## Equalized Time Files        
        Eq_widthArr = np.array([np.quantile(i, 0.75) for i in Eq_hData])-np.array([np.quantile(i, 0.25) for i in Eq_hData])

        Sh_hData = Sh_hFiles[keys_[-1]] ## Equalized Time Files        
        Sh_widthArr = np.array([np.quantile(i, 0.75) for i in Sh_hData])-np.array([np.quantile(i, 0.25) for i in Sh_hData])
     
        widthFrame = pd.DataFrame(-1, index=UB_hSubFiles, columns=['UBOrig', frq_+'Orig'])
        for sub_ in UB_hSubFiles:
            if (sub_ in hSubFiles):                
                print(sub_)
                idxs = np.where(sub_ == np.array(UB_hSubFiles))[0]
                idxs2 = np.where(sub_ == np.array(hSubFiles))[0]                                
                widthFrame.loc[sub_, 'UBOrig'] = UB_widthArr[idxs]
                widthFrame.loc[sub_, frq_+'Orig'] = widthArr[idxs2]                
        widthFrame.drop(widthFrame.index.values[np.where(np.sum(widthFrame.values, axis=1)==-2)[0]], inplace=True)
        widthFrame.to_csv(os.path.join(sourceDir, 'WidthDifference', grp+'_upperBeta-'+frq_+'-Orig.csv'))
        #scipy.stats.ttest_rel(widthFrame['UBOrig'], widthFrame[frq_+'Orig'], alternative='less')

        Eq_widthFrame = pd.DataFrame(-1, index=UB_Eq_hSubFiles, columns=['UBEqual', frq_+'Equal'])
        for sub_ in UB_Eq_hSubFiles:
            if (sub_ in Eq_hSubFiles):
                print(sub_)
                idxs = np.where(sub_ == np.array(UB_Eq_hSubFiles))[0]
                idxs2 = np.where(sub_ == np.array(Eq_hSubFiles))[0]
                Eq_widthFrame.loc[sub_, 'UBEqual'] = UB_Eq_widthArr[idxs]
                Eq_widthFrame.loc[sub_, frq_+'Equal'] = Eq_widthArr[idxs2]
        Eq_widthFrame.drop(Eq_widthFrame.index.values[np.where(np.sum(Eq_widthFrame.values, axis=1)==-2)[0]], inplace=True)
        Eq_widthFrame.to_csv(os.path.join(sourceDir, 'WidthDifference', grp+'_upperBeta-'+frq_+'-Equal.csv'))
        #scipy.stats.ttest_rel(Eq_widthFrame['UBEqual'], Eq_widthFrame[frq_+'Equal'], alternative='greater')
        
        Sh_widthFrame = pd.DataFrame(-1, index=UB_Sh_hSubFiles, columns=['UBShuff', frq_+'Shuff'])
        for sub_ in UB_Sh_hSubFiles:
            if (sub_ in Sh_hSubFiles):
                print(sub_)
                idxs = np.where(sub_ == np.array(UB_Sh_hSubFiles))[0]
                idxs2 = np.where(sub_ == np.array(Sh_hSubFiles))[0]
                Sh_widthFrame.loc[sub_, 'UBShuff'] = UB_Sh_widthArr[idxs]
                Sh_widthFrame.loc[sub_, frq_+'Shuff'] = Sh_widthArr[idxs2]
        Sh_widthFrame.drop(Sh_widthFrame.index.values[np.where(np.sum(Sh_widthFrame.values, axis=1)==-2)[0]], inplace=True)
        Sh_widthFrame.to_csv(os.path.join(sourceDir, 'WidthDifference', grp+'_upperBeta-'+frq_+'-Shuff.csv'))
        #scipy.stats.ttest_rel(Sh_widthFrame['UBShuff'], Sh_widthFrame[frq_+'Shuff'], alternative='greater')
        
################## BarPlot for Statistics

if not os.path.exists(os.path.join(sourceDir, 'hpeak-width-asymmetry')):
    os.makedirs(os.path.join(sourceDir, 'hpeak-width-asymmetry'))

for ii in np.arange(2,6):
    print('hpeak-width-asymmetry-db-%s' %str(ii))
    #if os.path.isfile(os.path.join(sourceDir, 'hpeak-width-asymmetry', 'hpeak-width-asymmetry-db-%s.png' %str(ii))):
    if not os.path.isfile(os.path.join(sourceDir, 'hpeak-width-asymmetry', 'hpeak-width-asymmetry-db-%s.png' %str(ii))):    
        fig = plt.figure(figsize=(20, 10.5))
        nrow = 7
        ncol = 5
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), wspace=0.15, hspace=0.14, left=0.04, right=0.995, bottom=0.075, top=0.94)
        row_ = 0
        col_ = 0

        Or_wilcoxS = pd.read_csv('Or-Eq-WilcoxStats-db%s.csv' %str(ii), index_col=0)
        Sh_wilcoxS = pd.read_csv('Or-Sh-WilcoxStats-db%s.csv' %str(ii), index_col=0)
        Or_wilcoxP = pd.read_csv('Or-Eq-WilcoxPval-db%s.csv' %str(ii), index_col=0)
        Sh_wilcoxP = pd.read_csv('Or-Sh-WilcoxPval-db%s.csv' %str(ii), index_col=0)
        Or_hPeak= pd.read_csv('Original_hPeak-db%s.csv' %str(ii), index_col=0)
        Or_Width= pd.read_csv('Original_width-db%s.csv' %str(ii), index_col=0)
        Or_Asy= pd.read_csv('Original_Asymmetry-db%s.csv' %str(ii), index_col=0)
        Eq_hPeak= pd.read_csv('Equalized_hPeak-db%s.csv' %str(ii), index_col=0)
        Eq_Width= pd.read_csv('Equalized_width-db%s.csv' %str(ii), index_col=0)
        Eq_Asy= pd.read_csv('Equalized_Asy-db%s.csv' %str(ii), index_col=0)
        Sh_hPeak= pd.read_csv('Shuffled_hPeak-db%s.csv' %str(ii), index_col=0)
        Sh_Width= pd.read_csv('Shuffled_width-db%s.csv' %str(ii), index_col=0)
        Sh_Asy= pd.read_csv('Shuffled_Asy-db%s.csv' %str(ii), index_col=0)

        for frq_ in frqBands:
            #for statArr in list([Or_hPeak, Eq_hPeak, Sh_hPeak, Or_Width, Eq_Width, Sh_Width, Or_Asy, Eq_Asy, Sh_Asy]):
            ax = plt.subplot(gs[row_,0])
            newDf = pd.DataFrame([], columns=['Original', 'Equalized', 'Shuffled'])
            newDf['Original'] = Or_hPeak.loc[frq_]
            newDf['Equalized'] = Eq_hPeak.loc[frq_]
            newDf['Shuffled'] = Sh_hPeak.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('h(max(spectrum))')
            if row_ != nrow-1:
                ax.set_xticklabels('')
            ax.legend('')
            ax.set_ylabel(frq_)

            ax = plt.subplot(gs[row_,1])
            newDf = pd.DataFrame([], columns=['Original', 'Equalized', 'Shuffled'])
            newDf['Original'] = Or_Width.loc[frq_]
            newDf['Equalized'] = Eq_Width.loc[frq_]
            newDf['Shuffled'] = Sh_Width.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Spectrum Width')
            if row_ != nrow-1:
                ax.set_xticklabels('')
            ax.legend('')

            ax = plt.subplot(gs[row_,2])
            newDf = pd.DataFrame([], columns=['Original', 'Equalized', 'Shuffled'])
            newDf['Original'] = Or_Asy.loc[frq_]
            newDf['Equalized'] = Eq_Asy.loc[frq_]
            newDf['Shuffled'] = Sh_Asy.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Spectrum Asymmetry')
            if row_ != nrow-1:
                ax.set_xticklabels('')
                ax.legend('')
            else:
                ax.legend(loc=2, fontsize=8)

            ax = plt.subplot(gs[row_,3])
            newDf = pd.DataFrame([], columns=['Orig-Shuf', 'Orig>Equa'])
            newDf['Orig>Equa'] = Or_wilcoxP.loc[frq_]
            newDf['Orig-Shuf'] = Sh_wilcoxP.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Wilcoxon P-val of h_peak')
            if row_ != nrow-1:
                ax.set_xticklabels('')
            ax.legend('')

            ax = plt.subplot(gs[row_,4])
            newDf = pd.DataFrame([], columns=['Orig-Shuf', 'Orig>Equa'])
            newDf['Orig>Equa'] = Or_wilcoxS.loc[frq_]
            newDf['Orig-Shuf'] = Sh_wilcoxS.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Wilcoxon Stats of h_peak')
            if row_ != nrow-1:
                ax.set_xticklabels('')
                ax.legend('')
            else:
                ax.legend(loc=2, fontsize=8)
            ax.set_ylim([0, 10000])

            row_ = row_ + 1

        plt.suptitle('db-%s' %str(ii))
        plt.savefig(os.path.join(sourceDir, 'hpeak-width-asymmetry', 'hpeak-width-asymmetry-db-%s.png' %str(ii)))
        plt.savefig(os.path.join(sourceDir, 'hpeak-width-asymmetry', 'hpeak-width-asymmetry-db-%s.pdf' %str(ii)))
        plt.close()
        plt.clf()        
    
    if not os.path.isfile(os.path.join(sourceDir, 'hpeak-width-asymmetry', 'UpperBeta_hPeakStatsTest-db-%s.png' %str(ii))):
        frq_ = 'UpperBeta'
        newDf = pd.DataFrame([], columns=['Orig-Shuf', 'Orig>Equa'])
        newDf['Orig>Equa'] = Or_wilcoxP.loc[frq_]
        newDf['Orig-Shuf'] = Sh_wilcoxP.loc[frq_]
        ax = newDf.plot.bar(fontsize=30, figsize=(20,10))
        #ax.set_title('Wilcoxon P-val of h_peak', fontsize=30)
        ax.legend(ncol=3, fontsize=25, fancybox=True, framealpha=0.6, title='Alternative H',loc='upper left', columnspacing=0.5, markerscale=0.5, title_fontsize=20)
        ax.tick_params(axis='both', which='both', labelsize=35)
        plt.savefig(os.path.join(sourceDir, 'hpeak-width-asymmetry', '%s_hPeakStatsTest-db-%s.png' %(frq_, str(ii))), bbox_inches='tight')
        plt.savefig(os.path.join(sourceDir, 'hpeak-width-asymmetry', '%s_hPeakStatsTest-db-%s.pdf' %(frq_, str(ii))), bbox_inches='tight')
        plt.close()
        plt.clf()

################## BarPlot for Local Singularity
'''GroupsAre = ['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08']
for grpIdx, grps_ in enumerate([GroupsAre[3]]):
    fig = plt.figure(figsize=(20, 10.5))
    nrow = 7
    ncol = 4
    gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), wspace = 0.05, hspace=0.1, left=0.04, right=0.99, bottom=0.1, top=0.92)
    row_ = 0
    col_ = 0
    for frq_ in frqBands:
        grp_ = glob.glob(os.path.join(sourceDir, frq_, grps_+"*"))[0]
        categ = np.arange(0, 3+(3/25), 3/25)
        categ = [round(i,3) for i in categ]
        for ii in np.arange(2,6):
            newArr_ = pd.DataFrame(0, index=categ, columns=['Original'])
            LSData = loadmat(os.path.join(grp_, 'db%s_LocalSingularity.mat' %str(ii)))['alphaLocalSingularity%s' %str(ii)][0]
            print('%s-%s Number of Local Signularity greater than 3 for db-%s is %s' %(frq_, grp_, str(ii), str(sum(LSData>3))))
            for d_ in LSData:
                if d_ <= 3:
                    dfIndex = np.argmin(abs(categ-d_))
                    newArr_.loc[categ[dfIndex], 'Original'] = newArr_.loc[categ[dfIndex], 'Original'] + 1

            Eq_LSData = loadmat(os.path.join(grp_, 'EQ_db%s_LocalSingularity.mat' %str(ii)))['alphaLocalSingularity%s' %str(ii)][0]
            print('%s-%s-Eq Number of Local Signularity greater than 3 for db-%s is %s' %(frq_, grp_, str(ii), str(sum(Eq_LSData>3))))
            newArr_['Equalized'] = 0
            for d_ in Eq_LSData:
                if d_ <= 3:
                    dfIndex = np.argmin(abs(categ-d_))
                    newArr_.loc[categ[dfIndex], 'Equalized'] = newArr_.loc[categ[dfIndex], 'Equalized'] + 1

            Sh_LSData = loadmat(os.path.join(grp_, 'SH_db%s_LocalSingularity.mat' %str(ii)))['alphaLocalSingularity%s' %str(ii)][0]
            print('%s-%s-Sh Number of Local Signularity greater than 3 for db-%s is %s' %(frq_, grp_, str(ii), str(sum(Sh_LSData>3))))
            newArr_['Shuffled'] = 0
            for d_ in Sh_LSData:
                if d_ <= 3:
                    dfIndex = np.argmin(abs(categ-d_))
                    newArr_.loc[categ[dfIndex], 'Shuffled'] = newArr_.loc[categ[dfIndex], 'Shuffled'] + 1

            ax = plt.subplot(gs[row_,col_])
            newArr_.plot.bar(ax=ax)

            if row_ == 0:
                ax.set_title('db-%s' %str(ii))
            if row_ != nrow-1:
                ax.set_xticklabels('')
            if (row_ == nrow-1) and (ii== 5):
                ax.legend(loc=2, fontsize=8)
            else: 
                ax.legend('')

            if col_ == 0:
                ax.set_ylabel(frq_)
            else:
                ax.set_yticklabels('')

            ax.set_ylim([0, 5000])
            col_ = col_ + 1

        row_ = row_ + 1
        col_ = 0

    pdb.set_trace()
    plt.suptitle(grps_)
    plt.savefig('%s-LocalSingularity.png' %str(grps_))
    plt.close()
    plt.clf()'''

################## BarPlot for Local Singularity: Only Original
'''if not os.path.exists(os.path.join(sourceDir, 'OriginalLocalSingularity')):
    os.makedirs(os.path.join(sourceDir, 'OriginalLocalSingularity'))

GroupsAre = ['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08']
for grpIdx, grps_ in enumerate(GroupsAre):
    if not os.path.isfile(os.path.join(sourceDir, 'OriginalLocalSingularity', '%s-Original_LocalSingularity.png' %str(grps_))):    
        fig = plt.figure(figsize=(20, 10.5))
        nrow = 7
        ncol = 4
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), wspace = 0.05, hspace=0.1, left=0.04, right=0.99, bottom=0.1, top=0.92)
        row_ = 0
        col_ = 0
        for frq_ in frqBands:
            grp_ = glob.glob(os.path.join(sourceDir, frq_, grps_+"*"))[0]
            categ = np.arange(0, 3+(3/25), 3/25)
            categ = [round(i,3) for i in categ]
            for ii in np.arange(2,6):
                newArr_ = pd.DataFrame(0, index=categ, columns=['Original'])
                LSData = loadmat(os.path.join(grp_, 'db%s_LocalSingularity.mat' %str(ii)))['alphaLocalSingularity%s' %str(ii)][0]
                print('%s-%s Number of Local Signularity greater than 3 for db-%s is %s' %(frq_, grp_, str(ii), str(sum(LSData>3))))
                for d_ in LSData:
                    if d_ <= 3:
                        dfIndex = np.argmin(abs(categ-d_))
                        newArr_.loc[categ[dfIndex], 'Original'] = newArr_.loc[categ[dfIndex], 'Original'] + 1

                ax = plt.subplot(gs[row_,col_])
                newArr_.plot.bar(ax=ax)

                if row_ == 0:
                    ax.set_title('db-%s' %str(ii))
                if row_ != nrow-1:
                    ax.set_xticklabels('')
                if (row_ == nrow-1) and (ii== 5):
                    ax.legend(loc=2, fontsize=8)
                else: 
                    ax.legend('')

                if col_ == 0:
                    ax.set_ylabel(frq_)
                #else:
                #    ax.set_yticklabels('')

                #ax.set_ylim([0, 5000])
                col_ = col_ + 1

            row_ = row_ + 1
            col_ = 0

        plt.suptitle(grps_)
        plt.savefig(os.path.join(sourceDir, 'OriginalLocalSingularity', '%s-Original_LocalSingularity.png' %str(grps_)))
        plt.savefig(os.path.join(sourceDir, 'OriginalLocalSingularity', '%s-Original_LocalSingularity.pdf' %str(grps_)))
        plt.close()
        plt.clf()

################## Width Distribution Plot
if not os.path.exists(os.path.join(sourceDir, 'WidthDist')):
    os.makedirs(os.path.join(sourceDir, 'WidthDist'))

for ii in np.arange(2,6):
    if not os.path.isfile(os.path.join(sourceDir, 'WidthDist', 'WidthDist_db-%s.png' %str(ii))):
        fig = plt.figure(figsize=(20, 10.5))
        nrow = 7
        ncol = 9
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), hspace=0.3, left=0.04, right=0.99, bottom=0.05, top=0.93)
        row_ = 0
        col_ = 0
        for frq_ in frqBands:
            grpSort = np.sort(glob.glob(os.path.join(sourceDir, frq_, '*')))
            for grp_ in grpSort:
                if '.png' in grp_.split('/')[-1]:
                    continue

                hFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityExponent.mat' %str(ii)))
                Eq_hFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_singularityExponent.mat' %str(ii)))
                Sh_hFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_singularityExponent.mat' %str(ii)))
                keys_ = [i for i in hFiles.keys()]
                hData = hFiles[keys_[-1]] ## Normal Files
                hMean = np.median(hData,0)
                hModeVal = np.round(5,2)
                width_ = np.round(np.quantile(hData, 0.75, axis=1)-np.quantile(hData, 0.25, axis=1),2)
                df_1 = pd.DataFrame(width_, columns=['Width'])
                df_1['Name'] = 'Original'
                Eq_hData = Eq_hFiles[keys_[-1]] ## Equalized Time Files
                Eq_hMean = np.median(Eq_hData,0)
                Eq_hModeVal = np.round(5,2)
                Eq_width_ = np.round(np.quantile(Eq_hData, 0.75, axis=1)-np.quantile(Eq_hData, 0.25, axis=1),2)
                df_2 = pd.DataFrame(Eq_width_, columns=['Width'])
                df_2['Name'] = 'Equalized'
                Sh_hData = Sh_hFiles[keys_[-1]] ## Shuffled Label File
                Sh_hMean = np.median(Sh_hData,0)
                Sh_hModeVal = np.round(5,2)
                Sh_width_ = np.round(np.quantile(Sh_hData, 0.75, axis=1)-np.quantile(Sh_hData, 0.25, axis=1),2)
                df_3 = pd.DataFrame(Sh_width_, columns=['Width'])
                df_3['Name'] = 'Shuffled'
                df_ = pd.concat((df_1, df_2), ignore_index=True)
                df_ = pd.concat((df_, df_3), ignore_index=True)

                ax = plt.subplot(gs[row_,col_])
                for name in ['Original', 'Equalized', 'Shuffled']:
                    # Subset to the airline
                    subset = df_[df_['Name'] == name]
                    
                    # Draw the density plot
                    sns.distplot(subset['Width'], hist = False, kde = True,
                                 kde_kws = {'linewidth': 3},
                                 label = name, ax=ax)
                # Plot formatting
                ax.legend(loc=1, prop={'size': 8})
                #ax.text(np.min(hMean), np.min(dhMean), 'M-%s,%s,%s\nW-%s,%s,%s' %(str(hModeVal),str(Eq_hModeVal),str(Sh_hModeVal), str(width_),str(Eq_width_),str(Sh_width_)))            
                #ax.scatter(hMean, dhMean, s=5)
                if col_ == 0:
                    ax.set_ylabel(frq_)
                else:
                    ax.set_ylabel([])
                if row_ == 0:
                    ax.set_title(grp_.split('/')[-1].split('_')[0])
                col_ = col_ + 1
            row_ = row_ + 1
            col_ = 0
        plt.suptitle('db-%s' %str(ii))
        plt.savefig(os.path.join(sourceDir, 'WidthDist', 'WidthDist_db-%s.png' %str(ii)))
        plt.savefig(os.path.join(sourceDir, 'WidthDist', 'WidthDist_db-%s.pdf' %str(ii)))
        plt.close()
        plt.clf()

################## H Peak Value Distribution Plot
if not os.path.exists(os.path.join(sourceDir, 'PeakValDist')):
    os.makedirs(os.path.join(sourceDir, 'PeakValDist'))

for ii in np.arange(2,6):
    if not os.path.isfile(os.path.join(sourceDir, 'PeakValDist', 'PeakValDist_db-%s.png' %str(ii))):
        Or_wilcoxS = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Eq_wilcoxS = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Sh_wilcoxS = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Or_wilcoxP = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Eq_wilcoxP = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Sh_wilcoxP = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])

        Or_Stats1 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Or_Stats2 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Or_Stats3 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Or_Pval1 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Or_Pval2 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Or_Pval3 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])

        Sh_Stats1 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Sh_Stats2 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Sh_Stats3 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Sh_Pval1 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Sh_Pval2 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Sh_Pval3 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])

        Eq_Stats1 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Eq_Stats2 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Eq_Stats3 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Eq_Pval1 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Eq_Pval2 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        Eq_Pval3 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])

        OrEq_Stats1 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        OrEq_Pval1 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        OrSh_Stats1 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        OrSh_Pval1 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        OrEq_Stats2 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        OrEq_Pval2 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        OrSh_Stats2 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        OrSh_Pval2 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        OrEq_Stats3 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        OrEq_Pval3 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        OrSh_Stats3 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        OrSh_Pval3 = pd.DataFrame(0, index=frqBands, columns=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])

        fig = plt.figure(figsize=(20, 10.5))
        nrow = 7
        ncol = 9
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), hspace=0.3, left=0.04, right=0.99, bottom=0.05, top=0.93)
        row_ = 0
        col_ = 0
        for frq_ in frqBands:
            grpSort = np.sort(glob.glob(os.path.join(sourceDir, frq_, '*')))
            for grp_ in grpSort:
                if '.png' in grp_.split('/')[-1]:
                    continue

                dhFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityDim.mat' %str(ii)))
                keys_ = [i for i in dhFiles.keys()]
                dhData = dhFiles[keys_[-1]] ## Normal Files
                peakIndex = [np.argmax(i) for i in dhData]
                peakVal = [np.max(i) for i in dhData]

                Eq_dhFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_singularityDim.mat' %str(ii))) ## Equalized Time Files
                Eq_dhData = Eq_dhFiles[keys_[-1]]
                Eq_peakIndex = [np.argmax(i) for i in Eq_dhData]
                Eq_peakVal = [np.max(i) for i in Eq_dhData]

                Sh_dhFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_singularityDim.mat' %str(ii))) ## Shuffled Label File
                Sh_dhData = Sh_dhFiles[keys_[-1]]
                Sh_peakIndex = [np.argmax(i) for i in Sh_dhData]
                Sh_peakVal = [np.max(i) for i in Sh_dhData]

                hFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityExponent.mat' %str(ii)))
                keys_ = [i for i in hFiles.keys()]
                hData = hFiles[keys_[-1]] ## Normal Files
                hModeVal = [hData[i, j] for i,j in enumerate(peakIndex)]
                ####### Just to cross check if the peakIndexes extracted from singularity spectrum are really representing the max peak value or not
                dhPeakVal = [dhData[i, j] for i,j in enumerate(peakIndex)]
                if sum(np.array(peakVal)-np.array(dhPeakVal)):
                    raise ValueError('There is some Problem')
                df_1 = pd.DataFrame(hModeVal, columns=['Mode'])
                df_1['Name'] = 'Original'

                Eq_hFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_singularityExponent.mat' %str(ii)))
                Eq_hData = Eq_hFiles[keys_[-1]] ## Equalized Time Files
                Eq_hModeVal = [Eq_hData[i, j] for i,j in enumerate(Eq_peakIndex)]
                ####### Just to cross check if the peakIndexes extracted from singularity spectrum are really representing the max peak value or not
                Eq_dhPeakVal = [Eq_dhData[i, j] for i,j in enumerate(Eq_peakIndex)]
                if sum(np.array(Eq_peakVal)-np.array(Eq_dhPeakVal)):
                    raise ValueError('There is some Problem')
                df_2 = pd.DataFrame(Eq_hModeVal, columns=['Mode'])
                df_2['Name'] = 'Equalized'

                Sh_hFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_singularityExponent.mat' %str(ii)))
                Sh_hData = Sh_hFiles[keys_[-1]] ## Shuffled Label File
                Sh_hModeVal = [Sh_hData[i, j] for i,j in enumerate(Sh_peakIndex)]
                ####### Just to cross check if the peakIndexes extracted from singularity spectrum are really representing the max peak value or not
                Sh_dhPeakVal = [Sh_dhData[i, j] for i,j in enumerate(Sh_peakIndex)]
                if sum(np.array(Sh_peakVal)-np.array(Sh_dhPeakVal)):
                    raise ValueError('There is some Problem')
                df_3 = pd.DataFrame(Sh_hModeVal, columns=['Mode'])
                df_3['Name'] = 'Shuffled'

                df_ = pd.concat((df_1, df_2), ignore_index=True)
                df_ = pd.concat((df_, df_3), ignore_index=True)

                stats, pval = scipy.stats.wilcoxon(hModeVal, Eq_hModeVal, alternative='greater')
                statsSH, pvalSH = scipy.stats.wilcoxon(hModeVal, Sh_hModeVal, alternative='two-sided')

                Or_wilcoxS.loc[frq_, grp_.split('/')[-1].split('_')[0]] = stats
                Sh_wilcoxS.loc[frq_, grp_.split('/')[-1].split('_')[0]] = statsSH
                Or_wilcoxP.loc[frq_, grp_.split('/')[-1].split('_')[0]] = pval
                Sh_wilcoxP.loc[frq_, grp_.split('/')[-1].split('_')[0]] = pvalSH

                stats = int(stats)
                pval = np.round(pval,4)
                statsSH = int(statsSH)
                pvalSH = np.round(pvalSH,4)
                ax = plt.subplot(gs[row_,col_])
                for name in ['Original', 'Equalized', 'Shuffled']:
                    # Subset to the airline
                    subset = df_[df_['Name'] == name]
                    
                    # Draw the density plot
                    sns.distplot(subset['Mode'], hist = False, kde = True,
                                 kde_kws = {'linewidth': 3},
                                 label = name, ax=ax)

                ax.text(0, ax.get_ylim()[1]-0.65,'(%s,%s\n%s,%s)' %(str(stats),str(pval),str(statsSH),str(pvalSH)), fontsize=8)
                ax.set_xlim([0, 1.5])
                ax.legend(loc=4, prop={'size': 8})
                if col_ == 0:
                    ax.set_ylabel(frq_)
                else:
                    ax.set_ylabel('')
                if row_ == 0:
                    ax.set_title(grp_.split('/')[-1].split('_')[0])
                col_ = col_ + 1

                cumulants = scipy.io.loadmat(os.path.join(grp_, 'db%s_cumulants.mat' %str(ii)))
                keys_ = [i for i in cumulants.keys()]
                cumulant = cumulants[keys_[-1]]
                stats1, pval1 = scipy.stats.wilcoxon(cumulant[:,0], np.ones(len(cumulant))*0.6, alternative='greater')
                stats2, pval2 = scipy.stats.wilcoxon(cumulant[:,1], np.ones(len(cumulant))*(-0.1), alternative='less')
                stats3, pval3 = scipy.stats.wilcoxon(cumulant[:,2], np.ones(len(cumulant))*0.1, alternative='greater')
                #stats2, pval2 = scipy.stats.wilcoxon(cumulant[:,1], np.ones(len(cumulant))*0, alternative='less')
                #stats3, pval3 = scipy.stats.wilcoxon(cumulant[:,2], np.ones(len(cumulant))*0, alternative='greater')
                Or_Pval1.loc[frq_, grp_.split('/')[-1].split('_')[0]] = pval1
                Or_Pval2.loc[frq_, grp_.split('/')[-1].split('_')[0]] = pval2
                Or_Pval3.loc[frq_, grp_.split('/')[-1].split('_')[0]] = pval3

                Eq_cumulants = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_cumulants.mat' %str(ii)))
                keys_ = [i for i in Eq_cumulants.keys()]
                Eq_cumulant = Eq_cumulants[keys_[-1]]
                Eq_stats1, Eq_pval1 = scipy.stats.wilcoxon(Eq_cumulant[:,0], np.ones(len(Eq_cumulant))*0.6, alternative='greater')
                Eq_stats2, Eq_pval2 = scipy.stats.wilcoxon(Eq_cumulant[:,1], np.ones(len(Eq_cumulant))*(-0.1), alternative='less')
                Eq_stats3, Eq_pval3 = scipy.stats.wilcoxon(Eq_cumulant[:,2], np.ones(len(Eq_cumulant))*0.1, alternative='greater')
                #Eq_stats2, Eq_pval2 = scipy.stats.wilcoxon(Eq_cumulant[:,1], np.ones(len(Eq_cumulant))*0, alternative='less')
                #Eq_stats3, Eq_pval3 = scipy.stats.wilcoxon(Eq_cumulant[:,2], np.ones(len(Eq_cumulant))*0, alternative='greater')
                Eq_Pval1.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Eq_pval1
                Eq_Pval2.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Eq_pval2
                Eq_Pval3.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Eq_pval3

                Sh_cumulants = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_cumulants.mat' %str(ii)))
                keys_ = [i for i in Sh_cumulants.keys()]
                Sh_cumulant = Sh_cumulants[keys_[-1]]
                Sh_stats1, Sh_pval1 = scipy.stats.wilcoxon(Sh_cumulant[:,0], np.ones(len(Sh_cumulant))*0.6, alternative='greater')
                Sh_stats2, Sh_pval2 = scipy.stats.wilcoxon(Sh_cumulant[:,1], np.ones(len(Sh_cumulant))*(-0.1), alternative='less')
                Sh_stats3, Sh_pval3 = scipy.stats.wilcoxon(Sh_cumulant[:,2], np.ones(len(Sh_cumulant))*0.1, alternative='greater')
                #Sh_stats2, Sh_pval2 = scipy.stats.wilcoxon(Sh_cumulant[:,1], np.ones(len(Sh_cumulant))*0, alternative='less')
                #Sh_stats3, Sh_pval3 = scipy.stats.wilcoxon(Sh_cumulant[:,2], np.ones(len(Sh_cumulant))*0, alternative='greater')
                Sh_Pval1.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Sh_pval1
                Sh_Pval2.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Sh_pval2
                Sh_Pval3.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Sh_pval3

                OrEq_stats1, OrEq_pval1 = scipy.stats.wilcoxon(cumulant[:,0], Eq_cumulant[:,0], alternative='greater')
                OrSh_stats1, OrSh_pval1 = scipy.stats.wilcoxon(cumulant[:,0], Sh_cumulant[:,0], alternative='two-sided')
                OrEq_stats2, OrEq_pval2 = scipy.stats.wilcoxon(cumulant[:,1], Eq_cumulant[:,1], alternative='less')
                OrSh_stats2, OrSh_pval2 = scipy.stats.wilcoxon(cumulant[:,1], Sh_cumulant[:,1], alternative='two-sided')
                OrEq_stats3, OrEq_pval3 = scipy.stats.wilcoxon(cumulant[:,2], Eq_cumulant[:,2], alternative='greater')
                OrSh_stats3, OrSh_pval3 = scipy.stats.wilcoxon(cumulant[:,2], Sh_cumulant[:,2], alternative='two-sided')

                OrEq_Stats1.loc[frq_, grp_.split('/')[-1].split('_')[0]] = OrEq_stats1
                OrEq_Pval1.loc[frq_, grp_.split('/')[-1].split('_')[0]] = OrEq_pval1
                OrSh_Stats1.loc[frq_, grp_.split('/')[-1].split('_')[0]] = OrSh_stats1
                OrSh_Pval1.loc[frq_, grp_.split('/')[-1].split('_')[0]] = OrSh_pval1
                OrEq_Stats2.loc[frq_, grp_.split('/')[-1].split('_')[0]] = OrEq_stats2
                OrEq_Pval2.loc[frq_, grp_.split('/')[-1].split('_')[0]] = OrEq_pval2
                OrSh_Stats2.loc[frq_, grp_.split('/')[-1].split('_')[0]] = OrSh_stats2
                OrSh_Pval2.loc[frq_, grp_.split('/')[-1].split('_')[0]] = OrSh_pval2
                OrEq_Stats3.loc[frq_, grp_.split('/')[-1].split('_')[0]] = OrEq_stats3
                OrEq_Pval3.loc[frq_, grp_.split('/')[-1].split('_')[0]] = OrEq_pval3
                OrSh_Stats3.loc[frq_, grp_.split('/')[-1].split('_')[0]] = OrSh_stats3
                OrSh_Pval3.loc[frq_, grp_.split('/')[-1].split('_')[0]] = OrSh_pval3

            row_ = row_ + 1
            col_ = 0

        Or_wilcoxS.to_csv('Or-Eq-WilcoxStats-db%s.csv' %str(ii))
        Sh_wilcoxS.to_csv('Or-Sh-WilcoxStats-db%s.csv' %str(ii))
        Or_wilcoxP.to_csv('Or-Eq-WilcoxPval-db%s.csv' %str(ii))
        Sh_wilcoxP.to_csv('Or-Sh-WilcoxPval-db%s.csv' %str(ii))
        Or_Pval1.to_csv('Original-Wilcox-Cumulant-1-db%s.csv' %str(ii))
        Or_Pval2.to_csv('Original-Wilcox-Cumulant-2-db%s.csv' %str(ii))
        Or_Pval3.to_csv('Original-Wilcox-Cumulant-3-db%s.csv' %str(ii))
        Eq_Pval1.to_csv('Equalized-Wilcox-Cumulant-1-db%s.csv' %str(ii))
        Eq_Pval2.to_csv('Equalized-Wilcox-Cumulant-2-db%s.csv' %str(ii))
        Eq_Pval3.to_csv('Equalized-Wilcox-Cumulant-3-db%s.csv' %str(ii))
        Sh_Pval1.to_csv('Shuffled-Wilcox-Cumulant-1-db%s.csv' %str(ii))
        Sh_Pval2.to_csv('Shuffled-Wilcox-Cumulant-2-db%s.csv' %str(ii))
        Sh_Pval3.to_csv('Shuffled-Wilcox-Cumulant-3-db%s.csv' %str(ii))

        OrEq_Stats1.to_csv('stats_OrEq-Wilcox-Cumulant-1-db%s.csv' %str(ii))
        OrEq_Pval1.to_csv('pval_OrEq-Wilcox-Cumulant-1-db%s.csv' %str(ii))
        OrSh_Stats1.to_csv('stats_OrSh-Wilcox-Cumulant-1-db%s.csv' %str(ii))
        OrSh_Pval1.to_csv('pval_OrSh-Wilcox-Cumulant-1-db%s.csv' %str(ii))
        OrEq_Stats2.to_csv('stats_OrEq-Wilcox-Cumulant-2-db%s.csv' %str(ii))
        OrEq_Pval2.to_csv('pval_OrEq-Wilcox-Cumulant-2-db%s.csv' %str(ii))
        OrSh_Stats2.to_csv('stats_OrSh-Wilcox-Cumulant-2-db%s.csv' %str(ii))
        OrSh_Pval2.to_csv('pval_OrSh-Wilcox-Cumulant-2-db%s.csv' %str(ii))
        OrEq_Stats3.to_csv('stats_OrEq-Wilcox-Cumulant-3-db%s.csv' %str(ii))
        OrEq_Pval3.to_csv('pval_OrEq-Wilcox-Cumulant-3-db%s.csv' %str(ii))
        OrSh_Stats3.to_csv('stats_OrSh-Wilcox-Cumulant-3-db%s.csv' %str(ii))
        OrSh_Pval3.to_csv('pval_OrSh-Wilcox-Cumulant-3-db%s.csv' %str(ii))

        plt.suptitle('db-%s' %str(ii))
        plt.savefig(os.path.join(sourceDir, 'PeakValDist', 'PeakValDist_db-%s.png' %str(ii)))
        plt.savefig(os.path.join(sourceDir, 'PeakValDist', 'PeakValDist_db-%s.pdf' %str(ii)))
        plt.close()
        plt.clf()

################## How emotion multifractal parameters are different from the baseline multifractal parameters.
if not os.path.exists(os.path.join(sourceDir, 'Singularity')):
    os.makedirs(os.path.join(sourceDir, 'Singularity'))
for ii in np.arange(2,6):
    if not os.path.isfile(os.path.join(sourceDir, 'Singularity', 'SingularitySpectrumShapeComparisonBaseGroupComp_db-%s.png' %str(ii))):
        BaseGComp_wilcoxGS = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        BaseGComp_wilcoxLS = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        BaseGComp_wilcoxGP = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        BaseGComp_wilcoxLP = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        WBaseGComp_wilcoxGS = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        WBaseGComp_wilcoxLS = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        WBaseGComp_wilcoxGP = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        WBaseGComp_wilcoxLP = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])

        fig = plt.figure(figsize=(20, 10.5))
        nrow = 7
        ncol = 2
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), hspace=0.3, left=0.04, right=0.99, bottom=0.1, top=0.93)
        row_ = 0
        col_ = 0
        for frq_ in frqBands:
            ## Calculating peak value of h for baseline
            baseline = glob.glob(os.path.join(sourceDir, frq_, '*Baseline*'))[0]
            with open(os.path.join(baseline, 'subjectsIncluded.txt')) as f_:
                filesIncludedBase = f_.readlines()
            f_.close()
            nobaseFileOneBifur = int(len(filesIncludedBase)/3)

            Base_dhFiles = scipy.io.loadmat(os.path.join(baseline, 'db%s_singularityDim.mat' %str(ii)))
            keys_ = [i for i in Base_dhFiles.keys()]
            Base_dhData = Base_dhFiles[keys_[-1]] ## Normal Files
            peakIndex = [np.argmax(i) for i in Base_dhData]
            Base_hFiles = scipy.io.loadmat(os.path.join(baseline, 'db%s_singularityExponent.mat' %str(ii)))
            keys_ = [i for i in Base_hFiles.keys()]
            Base_hData = Base_hFiles[keys_[-1]] ## Normal Files
            Base_width = np.quantile(Base_hData, 0.75, axis=1)-np.quantile(Base_hData, 0.25, axis=1)
            Base_hModeVal = [Base_hData[i, j] for i,j in enumerate(peakIndex)]
            del peakIndex
            del keys_

            plotF = pd.DataFrame([], index=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
            WplotF = pd.DataFrame([], index=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
            ## Calculating peak value of h for each group and comparing it with every emotion group.
            grpSort = np.sort(glob.glob(os.path.join(sourceDir, frq_, '*Group*')))
            for grp_ in grpSort:
                if '.png' in grp_.split('/')[-1]:
                    continue
                with open(os.path.join(grp_, 'subjectsIncluded.txt')) as f_:
                    filesIncludedGroup = f_.readlines()
                f_.close()
                noGroupFileOneBifur = int(len(filesIncludedGroup)/3)

                dhFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityDim.mat' %str(ii)))
                keys_ = [i for i in dhFiles.keys()]
                dhData = dhFiles[keys_[-1]] ## Normal Files
                peakIndex = [np.argmax(i) for i in dhData]
                hFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityExponent.mat' %str(ii)))
                keys_ = [i for i in hFiles.keys()]
                hData = hFiles[keys_[-1]] ## Normal Files
                Group_width = np.quantile(hData, 0.75, axis=1)-np.quantile(hData, 0.25, axis=1)
                hModeVal = [hData[i, j] for i,j in enumerate(peakIndex)]

                for idx_, val_ in enumerate(hModeVal): ## Enumerating over all the values of hModeVal
                    subName = filesIncludedGroup[idx_].split('_')[1] ## Extrating the subjectName so that index for same subject baseline is extracted in the coming codes.
                    bifurIndex = idx_/noGroupFileOneBifur # reaching to the bifurcation pool. Three kind of microstate bifurcations: [1&2:-1;3&4:+1],[1&3:-1;2&4:+1],[1&4:-1;2&3:+1]
                    lower, upper = [bifurIndex*nobaseFileOneBifur, (bifurIndex+1)*nobaseFileOneBifur]
                    BaselineFiles = filesIncludedBase[int(lower):int(upper)]
                    for flIdx, fl_ in enumerate(BaselineFiles):
                        if subName in fl_: # Getting the index of same subject baseline file from the appropriate bifurcation pool.
                            break;

                    CorrBaseFileIndex = int((bifurIndex*nobaseFileOneBifur)+flIdx); # Reaching to the same subject baseline file in the particular bifurcation pool
                    if subName in filesIncludedBase[CorrBaseFileIndex]:
                        print(subName, CorrBaseFileIndex, filesIncludedGroup[idx_], filesIncludedBase[CorrBaseFileIndex])
                        temp_ = np.reshape(np.array([hModeVal[idx_], Base_hModeVal[CorrBaseFileIndex]]), (1,2))
                        temp_2 = np.reshape(np.array([Group_width[idx_], Base_width[CorrBaseFileIndex]]), (1,2))
                        if idx_ == 0:
                            grpBase_ = temp_
                            grpBaseW = temp_2
                        else:
                            grpBase_ = np.concatenate((grpBase_, temp_), axis=0)
                            grpBaseW = np.concatenate((grpBaseW, temp_2), axis=0)

                df_ = pd.DataFrame(grpBase_, columns=['Group', 'Base'])
                dfW = pd.DataFrame(grpBaseW, columns=['Group', 'Base'])
                Gstats, Gpval = scipy.stats.wilcoxon(df_['Group'], df_['Base'], alternative='greater')
                Lstats, Lpval = scipy.stats.wilcoxon(df_['Group'], df_['Base'], alternative='less')
                Tstats, Tpval = scipy.stats.wilcoxon(df_['Group'], df_['Base'], alternative='two-sided')
                GstatsW, GpvalW = scipy.stats.wilcoxon(dfW['Group'], dfW['Base'], alternative='greater')
                LstatsW, LpvalW = scipy.stats.wilcoxon(dfW['Group'], dfW['Base'], alternative='less')
                TstatsW, TpvalW = scipy.stats.wilcoxon(dfW['Group'], dfW['Base'], alternative='two-sided')

                BaseGComp_wilcoxGS.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Gstats
                BaseGComp_wilcoxLS.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Lstats
                BaseGComp_wilcoxGP.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Gpval
                BaseGComp_wilcoxLP.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Lpval

                WBaseGComp_wilcoxGS.loc[frq_, grp_.split('/')[-1].split('_')[0]] = GstatsW
                WBaseGComp_wilcoxLS.loc[frq_, grp_.split('/')[-1].split('_')[0]] = LstatsW
                WBaseGComp_wilcoxGP.loc[frq_, grp_.split('/')[-1].split('_')[0]] = GpvalW
                WBaseGComp_wilcoxLP.loc[frq_, grp_.split('/')[-1].split('_')[0]] = LpvalW

                plotF.loc[grp_.split('/')[-1].split('_')[0], 'Group>Base'] = Gpval
                plotF.loc[grp_.split('/')[-1].split('_')[0], 'Group<Base'] = Lpval
                plotF.loc[grp_.split('/')[-1].split('_')[0], 'Group-Base'] = Tpval
                WplotF.loc[grp_.split('/')[-1].split('_')[0], 'Group>Base'] = GpvalW
                WplotF.loc[grp_.split('/')[-1].split('_')[0], 'Group<Base'] = LpvalW
                WplotF.loc[grp_.split('/')[-1].split('_')[0], 'Group-Base'] = TpvalW
                del df_
                del dfW

            ax = plt.subplot(gs[row_,0])
            plotF.plot.bar(ax=ax)
            ax.legend(loc=4, prop={'size': 8})
            if col_ == 0:
                ax.set_ylabel(frq_)
            else:
                ax.set_ylabel('')
            if row_ == 0:
                ax.set_title('Singularity Exponent Comparison')
            if row_ != 6:
                ax.set_xticklabels([])
            ax.set_ylim([0, 0.1])

            ax = plt.subplot(gs[row_,1])
            WplotF.plot.bar(ax=ax)
            ax.legend(loc=4, prop={'size': 8})
            ax.set_ylabel('')
            if row_ == 0:
                ax.set_title('Width Comparison')
            if row_ != 6:
                ax.set_xticklabels([])
            ax.set_ylim([0, 0.1])
            row_ = row_ + 1

        plt.suptitle('db-%s' %str(ii))
        plt.savefig(os.path.join(sourceDir, 'Singularity', 'SingularitySpectrumShapeComparisonBaseGroupComp_db-%s.png' %str(ii)))
        plt.savefig(os.path.join(sourceDir, 'Singularity', 'SingularitySpectrumShapeComparisonBaseGroupComp_db-%s.pdf' %str(ii)))
        plt.close()
        plt.clf()

################## Bifurcation Wise: How emotion multifractal parameters are different from the baseline multifractal parameters. 
if not os.path.exists(os.path.join(sourceDir, 'Singularity')):
    os.makedirs(os.path.join(sourceDir, 'Singularity'))

for ii in np.arange(2,6):
    if not os.path.isfile(os.path.join(sourceDir, 'Singularity', 'BifurcationWise_SingularitySpectrumShapeComparisonBaseGroupComp_db-%s.png' %str(ii))):
        BaseGComp_wilcoxGS = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        BaseGComp_wilcoxLS = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        BaseGComp_wilcoxGP = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        BaseGComp_wilcoxLP = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        WBaseGComp_wilcoxGS = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        WBaseGComp_wilcoxLS = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        WBaseGComp_wilcoxGP = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        WBaseGComp_wilcoxLP = pd.DataFrame(0, index=frqBands, columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])

        fig = plt.figure(figsize=(20, 10.5))
        nrow = 7
        ncol = 6
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), hspace=0.3, left=0.04, right=0.99, bottom=0.1, top=0.93)
        row_ = 0
        col_ = 0
        for frq_ in frqBands:
            ## Calculating peak value of h for baseline
            baseline = glob.glob(os.path.join(sourceDir, frq_, '*Baseline*'))[0]
            with open(os.path.join(baseline, 'subjectsIncluded.txt')) as f_:
                filesIncludedBase = f_.readlines()
            f_.close()
            nobaseFileOneBifur = int(len(filesIncludedBase)/3)

            for bfIdx in np.arange(3):
                lowerBase, upperBase = [int(bfIdx*nobaseFileOneBifur), int((bfIdx+1)*nobaseFileOneBifur)]
                Base_dhFiles = scipy.io.loadmat(os.path.join(baseline, 'db%s_singularityDim.mat' %str(ii)))
                keys_ = [i for i in Base_dhFiles.keys()]
                Base_dhData = Base_dhFiles[keys_[-1]][int(lowerBase):int(upperBase), :] ## Normal Files
                peakIndex = [np.argmax(i) for i in Base_dhData]
                Base_hFiles = scipy.io.loadmat(os.path.join(baseline, 'db%s_singularityExponent.mat' %str(ii)))
                keys_ = [i for i in Base_hFiles.keys()]
                Base_hData = Base_hFiles[keys_[-1]][int(lowerBase):int(upperBase), :] ## Normal Files
                Base_width = np.quantile(Base_hData, 0.75, axis=1)-np.quantile(Base_hData, 0.25, axis=1)
                Base_hModeVal = [Base_hData[i, j] for i,j in enumerate(peakIndex)]
                del peakIndex
                del keys_

                plotF = pd.DataFrame([], index=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
                WplotF = pd.DataFrame([], index=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
                ## Calculating peak value of h for each group and comparing it with every emotion group.
                grpSort = np.sort(glob.glob(os.path.join(sourceDir, frq_, '*Group*')))
                for grp_ in grpSort:
                    if '.png' in grp_.split('/')[-1]:
                        continue
                    with open(os.path.join(grp_, 'subjectsIncluded.txt')) as f_:
                        filesIncludedGroup = f_.readlines()
                    f_.close()
                    noGroupFileOneBifur = int(len(filesIncludedGroup)/3)
                    lowerGrp, upperGrp = [int(bfIdx*noGroupFileOneBifur), int((bfIdx+1)*noGroupFileOneBifur)]
                    dhFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityDim.mat' %str(ii)))
                    keys_ = [i for i in dhFiles.keys()]
                    dhData = dhFiles[keys_[-1]][int(lowerGrp):int(upperGrp), :] ## Normal Files
                    peakIndex = [np.argmax(i) for i in dhData]
                    hFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityExponent.mat' %str(ii)))
                    keys_ = [i for i in hFiles.keys()]
                    hData = hFiles[keys_[-1]][int(lowerGrp):int(upperGrp), :] ## Normal Files
                    Group_width = np.quantile(hData, 0.75, axis=1)-np.quantile(hData, 0.25, axis=1)
                    hModeVal = [hData[i, j] for i,j in enumerate(peakIndex)]

                     # Selecting the appropriate pool both for baseline and groups
                    BaselineFiles = filesIncludedBase[int(lowerBase):int(upperBase)]
                    GroupFiles = filesIncludedGroup[int(lowerGrp):int(upperGrp)]
                    for idx_, val_ in enumerate(hModeVal):
                        subName = GroupFiles[idx_].split('_')[1]
                        for flIdx, fl_ in enumerate(BaselineFiles):
                            if subName in fl_:
                                break

                        if subName in BaselineFiles[flIdx]:
                            print(subName, flIdx, GroupFiles[idx_], BaselineFiles[flIdx])
                            temp_ = np.reshape(np.array([hModeVal[idx_], Base_hModeVal[flIdx]]), (1,2))
                            temp_2 = np.reshape(np.array([Group_width[idx_], Base_width[flIdx]]), (1,2))
                            if idx_ == 0:
                                grpBase_ = temp_
                                grpBaseW = temp_2
                            else:
                                grpBase_ = np.concatenate((grpBase_, temp_), axis=0)
                                grpBaseW = np.concatenate((grpBaseW, temp_2), axis=0)

                    df_ = pd.DataFrame(grpBase_, columns=['Group', 'Base'])
                    dfW = pd.DataFrame(grpBaseW, columns=['Group', 'Base'])
                    Gstats, Gpval = scipy.stats.wilcoxon(df_['Group'], df_['Base'], alternative='greater')
                    Lstats, Lpval = scipy.stats.wilcoxon(df_['Group'], df_['Base'], alternative='less')
                    Tstats, Tpval = scipy.stats.wilcoxon(df_['Group'], df_['Base'], alternative='two-sided')
                    GstatsW, GpvalW = scipy.stats.wilcoxon(dfW['Group'], dfW['Base'], alternative='greater')
                    LstatsW, LpvalW = scipy.stats.wilcoxon(dfW['Group'], dfW['Base'], alternative='less')
                    TstatsW, TpvalW = scipy.stats.wilcoxon(dfW['Group'], dfW['Base'], alternative='two-sided')

                    BaseGComp_wilcoxGS.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Gstats
                    BaseGComp_wilcoxLS.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Lstats
                    BaseGComp_wilcoxGP.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Gpval
                    BaseGComp_wilcoxLP.loc[frq_, grp_.split('/')[-1].split('_')[0]] = Lpval

                    WBaseGComp_wilcoxGS.loc[frq_, grp_.split('/')[-1].split('_')[0]] = GstatsW
                    WBaseGComp_wilcoxLS.loc[frq_, grp_.split('/')[-1].split('_')[0]] = LstatsW
                    WBaseGComp_wilcoxGP.loc[frq_, grp_.split('/')[-1].split('_')[0]] = GpvalW
                    WBaseGComp_wilcoxLP.loc[frq_, grp_.split('/')[-1].split('_')[0]] = LpvalW

                    plotF.loc[grp_.split('/')[-1].split('_')[0], 'Group>Base'] = Gpval
                    plotF.loc[grp_.split('/')[-1].split('_')[0], 'Group<Base'] = Lpval
                    plotF.loc[grp_.split('/')[-1].split('_')[0], 'Group-Base'] = Tpval
                    WplotF.loc[grp_.split('/')[-1].split('_')[0], 'Group>Base'] = GpvalW
                    WplotF.loc[grp_.split('/')[-1].split('_')[0], 'Group<Base'] = LpvalW
                    WplotF.loc[grp_.split('/')[-1].split('_')[0], 'Group-Base'] = TpvalW
                    del df_
                    del dfW

                ax = plt.subplot(gs[row_,bfIdx])
                plotF.plot.bar(ax=ax)
                ax.legend([])
                if bfIdx == 0:
                    ax.set_ylabel(frq_)
                else:
                    ax.set_ylabel('')
                    ax.set_yticklabels([])

                if row_ == 0:
                    ax.set_title('Singularity Exponent-%s' %str(bfIdx+1))
                if row_ != 6:
                    ax.set_xticklabels([])
                ax.set_ylim([0, 0.1])

                ax = plt.subplot(gs[row_,bfIdx+3])
                WplotF.plot.bar(ax=ax)
                ax.set_ylabel('')
                ax.set_yticklabels([])

                if row_ == 0:
                    ax.set_title('Width-%s' %str(bfIdx+1))
                if row_ != 6:
                    ax.set_xticklabels([])
                ax.set_ylim([0, 0.1])
                ax.legend([])
            row_ = row_ + 1
        ax.legend(loc=4, prop={'size': 8})
        plt.suptitle('db-%s' %str(ii))
        plt.savefig(os.path.join(sourceDir, 'Singularity', 'BifurcationWise_SingularitySpectrumShapeComparisonBaseGroupComp_db-%s.png' %str(ii)))
        plt.savefig(os.path.join(sourceDir, 'Singularity', 'BifurcationWise_SingularitySpectrumShapeComparisonBaseGroupComp_db-%s.pdf' %str(ii)))
        plt.close()
        plt.clf()
'''
################## Bifurcation Wise:  Statistical Calculation
''' Two things has to be shown.
1. They are multifractal (statistical significance with equalized signals)
2. differences between them is trivial (show only the statistical representations eg.median peak for different bipartition scheme)
'''

if not os.path.exists(os.path.join(sourceDir, 'PeakValDist')):
    os.makedirs(os.path.join(sourceDir, 'PeakValDist'))

for ii in np.arange(2,6):
    if not os.path.isfile(os.path.join(sourceDir, 'PeakValDist', 'Bipart_MultifractalityCheck_db-%s.png' %str(ii))):
        Stats = pd.DataFrame(0, index=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'], columns=['h:Or>Eq','W:Or>Eq'])
        Pval = pd.DataFrame(0, index=['Baseline','Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'], columns=['h:Or>Eq','W:Or>Eq'])

        fig = plt.figure(figsize=(20, 10.5))
        #plt.locator_params(axis="y", tight='False',nbins=6)
        nrow = 7
        ncol = 3
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), wspace=0.1, hspace=0.3, left=0.04, right=0.99, bottom=0.1, top=0.93)
        row_ = 0
        col_ = 0
        for frq_ in frqBands:
            grpSort = np.sort(glob.glob(os.path.join(sourceDir, frq_, '*')))
            for bfIdx in np.arange(3):
                for grp_ in grpSort:
                    with open(os.path.join(grp_, 'subjectsIncluded.txt')) as f_:
                        filesIncludedGroup = f_.readlines()
                    f_.close()
                    noGroupFileOneBifur = int(len(filesIncludedGroup)/3)
                    lowerGrp, upperGrp = [int(bfIdx*noGroupFileOneBifur), int((bfIdx+1)*noGroupFileOneBifur)]
                    ##### Reading Singularity Dimension Files
                    dhFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityDim.mat' %str(ii)))
                    keys_ = [i for i in dhFiles.keys()]
                    dhData = dhFiles[keys_[-1]][int(lowerGrp):int(upperGrp), :] ## Normal Files
                    peakIndex = [np.argmax(i) for i in dhData]
                    dhMean = np.median(dhData,0) 
                    peakIndexMedian = np.argmax(dhMean)

                    Eq_dhFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_singularityDim.mat' %str(ii)))
                    Eq_keys_ = [i for i in Eq_dhFiles.keys()]
                    Eq_dhData = Eq_dhFiles[keys_[-1]][int(lowerGrp):int(upperGrp), :] ## Normal Files
                    Eq_peakIndex = [np.argmax(i) for i in Eq_dhData]
                    Eq_dhMean = np.median(Eq_dhData,0) 
                    Eq_peakIndexMedian = np.argmax(Eq_dhMean)

                    Sh_dhFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_singularityDim.mat' %str(ii)))
                    Sh_keys_ = [i for i in Sh_dhFiles.keys()]
                    Sh_dhData = Sh_dhFiles[keys_[-1]][int(lowerGrp):int(upperGrp), :] ## Normal Files
                    Sh_peakIndex = [np.argmax(i) for i in Sh_dhData]
                    Sh_dhMean = np.median(Sh_dhData,0) 
                    Sh_peakIndexMedian = np.argmax(Sh_dhMean)

                    ##### Reading Singularity Exponent Files Now
                    hFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityExponent.mat' %str(ii)))
                    keys_ = [i for i in hFiles.keys()]
                    hData = hFiles[keys_[-1]][int(lowerGrp):int(upperGrp), :] ## Normal Files
                    #### Stats on Arr
                    widthArr = np.quantile(hData, 0.75, axis=1)-np.quantile(hData, 0.25, axis=1)
                    hPeakArr = [hData[i, j] for i,j in enumerate(peakIndex)]
                    hminArr = np.min(hData, axis=1)
                    hmaxArr = np.max(hData, axis=1)
                    DeltaLArr = hPeakArr-hminArr
                    DeltaRArr = hmaxArr-hPeakArr
                    AsyArr = np.round((DeltaLArr-DeltaRArr)/(DeltaLArr+DeltaRArr),2)
                    #### Stats on Representative Values
                    hMean = np.median(hData,0)
                    hModeVal = np.round(hMean[peakIndexMedian],2)
                    hmin = np.min(hMean)
                    hmax = np.max(hMean)
                    hpeak = hModeVal
                    DeltaL = hpeak-hmin
                    DeltaR = hmax-hpeak
                    A = np.round((DeltaL-DeltaR)/(DeltaL+DeltaR),2)
                    width_ = np.round(np.quantile(hMean, 0.75)-np.quantile(hMean, 0.25),2)
    ################## Equalzed Time File
                    Eq_hFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_singularityExponent.mat' %str(ii)))
                    Eq_keys_ = [i for i in Eq_hFiles.keys()]
                    Eq_hData = Eq_hFiles[keys_[-1]][int(lowerGrp):int(upperGrp), :] ## Normal Files
                    #### Stats on Arr
                    Eq_widthArr = np.quantile(Eq_hData, 0.75, axis=1)-np.quantile(Eq_hData, 0.25, axis=1)
                    Eq_hPeakArr = [Eq_hData[i, j] for i,j in enumerate(Eq_peakIndex)]
                    Eq_hminArr = np.min(Eq_hData, axis=1)
                    Eq_hmaxArr = np.max(Eq_hData, axis=1)
                    Eq_DeltaLArr = Eq_hPeakArr-Eq_hminArr
                    Eq_DeltaRArr = Eq_hmaxArr-Eq_hPeakArr
                    Eq_AsyArr = np.round((Eq_DeltaLArr-Eq_DeltaRArr)/(Eq_DeltaLArr+Eq_DeltaRArr),2)
                    #### Stats on Representative Values
                    Eq_hMean = np.median(Eq_hData,0)
                    Eq_hModeVal = np.round(Eq_hMean[Eq_peakIndexMedian],2)
                    Eq_hmin = np.min(Eq_hMean)
                    Eq_hmax = np.max(Eq_hMean)
                    Eq_hpeak = Eq_hModeVal
                    Eq_DeltaL = Eq_hpeak-Eq_hmin
                    Eq_DeltaR = Eq_hmax-Eq_hpeak
                    Eq_A = np.round((Eq_DeltaL-Eq_DeltaR)/(Eq_DeltaL+Eq_DeltaR),2)
                    Eq_width_ = np.round(np.quantile(Eq_hMean, 0.75)-np.quantile(Eq_hMean, 0.25),2)
    ################## Shuffled Time File
                    Sh_hFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_singularityExponent.mat' %str(ii)))
                    Sh_keys_ = [i for i in Sh_hFiles.keys()]
                    Sh_hData = Sh_hFiles[keys_[-1]][int(lowerGrp):int(upperGrp), :] ## Normal Files
                    #### Stats on Arr
                    Sh_widthArr = np.quantile(Sh_hData, 0.75, axis=1)-np.quantile(Sh_hData, 0.25, axis=1)
                    Sh_hPeakArr = [Sh_hData[i, j] for i,j in enumerate(Sh_peakIndex)]
                    Sh_hminArr = np.min(Sh_hData, axis=1)
                    Sh_hmaxArr = np.max(Sh_hData, axis=1)
                    Sh_DeltaLArr = Sh_hPeakArr-Sh_hminArr
                    Sh_DeltaRArr = Sh_hmaxArr-Sh_hPeakArr
                    Sh_AsyArr = np.round((Sh_DeltaLArr-Sh_DeltaRArr)/(Sh_DeltaLArr+Sh_DeltaRArr),2)
                    #### Stats on Representative Values
                    Sh_hMean = np.median(Sh_hData,0)
                    Sh_hModeVal = np.round(Sh_hMean[Sh_peakIndexMedian],2)
                    Sh_hmin = np.min(Sh_hMean)
                    Sh_hmax = np.max(Sh_hMean)
                    Sh_hpeak = Sh_hModeVal
                    Sh_DeltaL = Sh_hpeak-Sh_hmin
                    Sh_DeltaR = Sh_hmax-Sh_hpeak
                    Sh_A = np.round((Sh_DeltaL-Sh_DeltaR)/(Sh_DeltaL+Sh_DeltaR),2)
                    Sh_width_ = np.round(np.quantile(Sh_hMean, 0.75)-np.quantile(Sh_hMean, 0.25),2)

                    stats, pval = scipy.stats.wilcoxon(hPeakArr, Eq_hPeakArr, alternative='greater')
                    Wstats, Wpval = scipy.stats.wilcoxon(widthArr, Eq_widthArr, alternative='greater')

                    Stats.loc[grp_.split('/')[-1].split('_')[0],'h:Or>Eq'] = stats
                    Pval.loc[grp_.split('/')[-1].split('_')[0],'h:Or>Eq'] = pval
                    Stats.loc[grp_.split('/')[-1].split('_')[0],'W:Or>Eq'] = Wstats
                    Pval.loc[grp_.split('/')[-1].split('_')[0],'W:Or>Eq'] = Wpval
                    #### No need to check shuffling scheme here
                    #Stats.loc[grp_.split('/')[-1].split('_')[0],'Or-Sh'] = statsSH
                    #Pval.loc[grp_.split('/')[-1].split('_')[0],'Or-Sh'] = pvalSH

                Pval.to_csv(os.path.join(sourceDir, 'PeakValDist', 'Frequency-%s_Bipart-%s_MultifractalityCheck_db-%s.csv' %(frq_, str(bfIdx), str(ii))))
                ax = plt.subplot(gs[row_,col_])
                Pval.plot.bar(ax=ax)
                ax.set_ylim([0, 0.15])
                #ax.set_yticks(np.arange(0,4))
                #ax.set_yticklabels(np.arange(0, 0.2, 0.05))
                if row_ == 0:
                    ax.set_title('BiPart-%s: Singularity Exponent' %str(bfIdx+1))
                if row_ != nrow-1:
                    ax.set_xticklabels('')
                if col_ == 0:
                    ax.set_ylabel(frq_)
                    #ax.set_ylim([0, 0.15])
                else:
                    ax.set_ylabel('')
                    #ax.set_yticklabels([])
                if (col_==(ncol-1) and row_==(nrow-1)):
                    print('show legend')
                else:
                    ax.legend('')

                col_ = col_ + 1

            row_ = row_ + 1
            col_ = 0

        pdb.set_trace()
        plt.suptitle('db-%s' %str(ii))
        plt.savefig(os.path.join(sourceDir, 'PeakValDist', 'Bipart_MultifractalityCheck_db-%s.png' %str(ii)))
        plt.savefig(os.path.join(sourceDir, 'PeakValDist', 'Bipart_MultifractalityCheck_db-%s.pdf' %str(ii)))
        plt.close()
        plt.clf()


    if not os.path.isfile(os.path.join(sourceDir, 'PeakValDist', 'Frequency-UpperBeta_Bipart_MultifractalityCheck_db-%s.png' %str(ii))):    
        frq_ = 'UpperBeta'
        fig = plt.figure(figsize=(20, 6))
        #plt.locator_params(axis="y", tight='False',nbins=6)
        nrow = 1
        ncol = 3
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), wspace=0.2, hspace=0.3, left=0.05, right=0.99, bottom=0.3, top=0.93)
        row_ = 0
        col_ = 0

        for bfIdx in np.arange(3):
            Pval = pd.read_csv(os.path.join(sourceDir, 'PeakValDist', 'Frequency-%s_Bipart-%s_MultifractalityCheck_db-%s.csv' %(frq_, str(bfIdx), str(ii))), index_col=0)
            ax = plt.subplot(gs[row_,bfIdx])
            Pval.plot.bar(ax=ax)
            ax.set_ylim([0, 0.06])
            #ax.set_yticks(np.arange(0,4))
            #ax.set_yticklabels(np.arange(0, 0.2, 0.05))
            if row_ == 0:
                ax.set_title('BiPart-%s: Singularity Exponent' %str(bfIdx+1), fontsize=25)

            ax.tick_params(axis='both', which='both', labelsize=25)
            ax.legend(ncol=2, fontsize=25, fancybox=True, framealpha=0.6, title='Alternative H',loc='upper left', columnspacing=0.5, markerscale=0.5, title_fontsize=20)            
        
        plt.savefig(os.path.join(sourceDir, 'PeakValDist', 'Frequency-%s_Bipart_MultifractalityCheck_db-%s.png' %(frq_, str(ii))))
        plt.savefig(os.path.join(sourceDir, 'PeakValDist', 'Frequency-%s_Bipart_MultifractalityCheck_db-%s.pdf' %(frq_, str(ii))))

'''################## Bifurcation Wise: BarPlot for Statistics
if not os.path.exists(os.path.join(sourceDir, 'hpeak-width-asymmetry')):
    os.makedirs(os.path.join(sourceDir, 'hpeak-width-asymmetry'))

for ii in np.arange(2,6):
    if not os.path.isfile(os.path.join(sourceDir, 'hpeak-width-asymmetry', 'Bipartition_hpeak-width-asymmetry-db-%s.png' %str(ii))):
        fig = plt.figure(figsize=(20, 10.5))
        nrow = 7
        ncol = 5
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), wspace=0.15, hspace=0.14, left=0.04, right=0.995, bottom=0.075, top=0.94)
        row_ = 0
        col_ = 0

        Or_wilcoxS = pd.read_csv('Or-Eq-WilcoxStats-db%s.csv' %str(ii), index_col=0)
        Sh_wilcoxS = pd.read_csv('Or-Sh-WilcoxStats-db%s.csv' %str(ii), index_col=0)
        Or_wilcoxP = pd.read_csv('Or-Eq-WilcoxPval-db%s.csv' %str(ii), index_col=0)
        Sh_wilcoxP = pd.read_csv('Or-Sh-WilcoxPval-db%s.csv' %str(ii), index_col=0)
        Or_hPeak= pd.read_csv('Original_hPeak-db%s.csv' %str(ii), index_col=0)
        Or_Width= pd.read_csv('Original_width-db%s.csv' %str(ii), index_col=0)
        Or_Asy= pd.read_csv('Original_Asymmetry-db%s.csv' %str(ii), index_col=0)
        Eq_hPeak= pd.read_csv('Equalized_hPeak-db%s.csv' %str(ii), index_col=0)
        Eq_Width= pd.read_csv('Equalized_width-db%s.csv' %str(ii), index_col=0)
        Eq_Asy= pd.read_csv('Equalized_Asy-db%s.csv' %str(ii), index_col=0)
        Sh_hPeak= pd.read_csv('Shuffled_hPeak-db%s.csv' %str(ii), index_col=0)
        Sh_Width= pd.read_csv('Shuffled_width-db%s.csv' %str(ii), index_col=0)
        Sh_Asy= pd.read_csv('Shuffled_Asy-db%s.csv' %str(ii), index_col=0)

        for frq_ in frqBands:
            #for statArr in list([Or_hPeak, Eq_hPeak, Sh_hPeak, Or_Width, Eq_Width, Sh_Width, Or_Asy, Eq_Asy, Sh_Asy]):
            ax = plt.subplot(gs[row_,0])
            newDf = pd.DataFrame([], columns=['Original', 'Equalized', 'Shuffled'])
            newDf['Original'] = Or_hPeak.loc[frq_]
            newDf['Equalized'] = Eq_hPeak.loc[frq_]
            newDf['Shuffled'] = Sh_hPeak.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('h(max(spectrum))')
            if row_ != nrow-1:
                ax.set_xticklabels('')
            ax.legend('')
            ax.set_ylabel(frq_)

            ax = plt.subplot(gs[row_,1])
            newDf = pd.DataFrame([], columns=['Original', 'Equalized', 'Shuffled'])
            newDf['Original'] = Or_Width.loc[frq_]
            newDf['Equalized'] = Eq_Width.loc[frq_]
            newDf['Shuffled'] = Sh_Width.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Spectrum Width')
            if row_ != nrow-1:
                ax.set_xticklabels('')
            ax.legend('')

            ax = plt.subplot(gs[row_,2])
            newDf = pd.DataFrame([], columns=['Original', 'Equalized', 'Shuffled'])
            newDf['Original'] = Or_Asy.loc[frq_]
            newDf['Equalized'] = Eq_Asy.loc[frq_]
            newDf['Shuffled'] = Sh_Asy.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Spectrum Asymmetry')
            if row_ != nrow-1:
                ax.set_xticklabels('')
                ax.legend('')
            else:
                ax.legend(loc=2, fontsize=8)

            ax = plt.subplot(gs[row_,3])
            newDf = pd.DataFrame([], columns=['Orig-Shuf', 'Orig>Equa'])
            newDf['Orig>Equa'] = Or_wilcoxP.loc[frq_]
            newDf['Orig-Shuf'] = Sh_wilcoxP.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Wilcoxon P-val of h_peak')
            if row_ != nrow-1:
                ax.set_xticklabels('')
            ax.legend('')

            ax = plt.subplot(gs[row_,4])
            newDf = pd.DataFrame([], columns=['Orig-Shuf', 'Orig>Equa'])
            newDf['Orig>Equa'] = Or_wilcoxS.loc[frq_]
            newDf['Orig-Shuf'] = Sh_wilcoxS.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Wilcoxon Stats of h_peak')
            if row_ != nrow-1:
                ax.set_xticklabels('')
                ax.legend('')
            else:
                ax.legend(loc=2, fontsize=8)
            ax.set_ylim([0, 10000])

            row_ = row_ + 1

        pdb.set_trace()
        plt.suptitle('db-%s' %str(ii))
        plt.savefig(os.path.join(sourceDir, 'hpeak-width-asymmetry', 'hpeak-width-asymmetry-db-%s.png' %str(ii)))
        plt.savefig(os.path.join(sourceDir, 'hpeak-width-asymmetry', 'hpeak-width-asymmetry-db-%s.pdf' %str(ii)))
        plt.close()
        plt.clf()

    if not os.path.isfile(os.path.join(sourceDir, 'hpeak-width-asymmetry', '%s_hPeakStatsTest-db-%s.pdf' %(frq_, str(ii)))):
        frq_ = 'UpperBeta'
        newDf = pd.DataFrame([], columns=['Orig-Shuf', 'Orig>Equa'])
        newDf['Orig>Equa'] = Or_wilcoxP.loc[frq_]
        newDf['Orig-Shuf'] = Sh_wilcoxP.loc[frq_]
        ax = newDf.plot.bar(fontsize=30, figsize=(20,10))
        ax.set_title('Wilcoxon P-val of h_mode', fontsize=30)
        ax.legend(ncol=3, fontsize=25, fancybox=True, framealpha=0.6, title='Alternative H',loc='upper left', columnspacing=0.5, markerscale=0.5, title_fontsize=20)
        pdb.set_trace()
        plt.savefig(os.path.join(sourceDir, 'hpeak-width-asymmetry', '%s_hPeakStatsTest-db-%s.png' %(frq_, str(ii))), bbox_inches='tight')
        plt.savefig(os.path.join(sourceDir, 'hpeak-width-asymmetry', '%s_hPeakStatsTest-db-%s.pdf' %(frq_, str(ii))), bbox_inches='tight')
        plt.close()
        plt.clf()        '''

################## BarPlot for Cumulants Statistics
if not os.path.exists(os.path.join(sourceDir, 'Cumulants')):
    os.makedirs(os.path.join(sourceDir, 'Cumulants'))

for ii in np.arange(2,6):
    if not os.path.isfile(os.path.join(sourceDir, 'Cumulants', 'Cumulants-db-%s.png' %str(ii))):
        fig = plt.figure(figsize=(20, 10.5))
        nrow = 7
        ncol = 3
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), wspace=0.15, hspace=0.14, left=0.04, right=0.995, bottom=0.075, top=0.94)
        row_ = 0
        col_ = 0

        Or_Pval1 = pd.read_csv('Original-Wilcox-Cumulant-1-db%s.csv' %str(ii), index_col=0)
        Or_Pval2 = pd.read_csv('Original-Wilcox-Cumulant-2-db%s.csv' %str(ii), index_col=0)
        Or_Pval3 = pd.read_csv('Original-Wilcox-Cumulant-3-db%s.csv' %str(ii), index_col=0)
        Eq_Pval1 = pd.read_csv('Equalized-Wilcox-Cumulant-1-db%s.csv' %str(ii), index_col=0)
        Eq_Pval2 = pd.read_csv('Equalized-Wilcox-Cumulant-2-db%s.csv' %str(ii), index_col=0)
        Eq_Pval3 = pd.read_csv('Equalized-Wilcox-Cumulant-3-db%s.csv' %str(ii), index_col=0)
        Sh_Pval1 = pd.read_csv('Shuffled-Wilcox-Cumulant-1-db%s.csv' %str(ii), index_col=0)
        Sh_Pval2 = pd.read_csv('Shuffled-Wilcox-Cumulant-2-db%s.csv' %str(ii), index_col=0)
        Sh_Pval3 = pd.read_csv('Shuffled-Wilcox-Cumulant-3-db%s.csv' %str(ii), index_col=0)

        for frq_ in frqBands:
            #for statArr in list([Or_hPeak, Eq_hPeak, Sh_hPeak, Or_Width, Eq_Width, Sh_Width, Or_Asy, Eq_Asy, Sh_Asy]):
            ax = plt.subplot(gs[row_,0])
            newDf = pd.DataFrame([], columns=['Original', 'Equalized', 'Shuffled'])
            newDf['Original'] = Or_Pval1.loc[frq_]
            newDf['Equalized'] = Eq_Pval1.loc[frq_]
            newDf['Shuffled'] = Sh_Pval1.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Cumulant-1')
            if row_ != nrow-1:
                ax.set_xticklabels('')
            ax.legend('')
            ax.set_ylim([0,0.5])
            plt.axhline(y=0.001, color='r', linestyle='--')
            plt.axhline(y=0.01, color='b', linestyle='--')
            plt.axhline(y=0.05, color='g', linestyle='--')
            ax.set_ylabel(frq_)

            ax = plt.subplot(gs[row_,1])
            newDf = pd.DataFrame([], columns=['Original', 'Equalized', 'Shuffled'])
            newDf['Original'] = Or_Pval2.loc[frq_]
            newDf['Equalized'] = Eq_Pval2.loc[frq_]
            newDf['Shuffled'] = Sh_Pval2.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Cumulant-2')
            if row_ != nrow-1:
                ax.set_xticklabels('')
            ax.legend('')
            ax.set_ylim([0,0.5])
            plt.axhline(y=0.001, color='r', linestyle='--')
            plt.axhline(y=0.01, color='b', linestyle='--')
            plt.axhline(y=0.05, color='g', linestyle='--')
  
            ax = plt.subplot(gs[row_,2])
            newDf = pd.DataFrame([], columns=['Original', 'Equalized', 'Shuffled'])
            newDf['Original'] = Or_Pval3.loc[frq_]
            newDf['Equalized'] = Eq_Pval3.loc[frq_]
            newDf['Shuffled'] = Sh_Pval3.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Cumulant-3')
            if row_ != nrow-1:
                ax.set_xticklabels('')
                ax.legend('')
            else:
                ax.legend(loc=2, fontsize=8)
            ax.set_ylim([0,0.5])
            plt.axhline(y=0.001, color='r', linestyle='--')
            plt.axhline(y=0.01, color='b', linestyle='--')
            plt.axhline(y=0.05, color='g', linestyle='--')

            row_ = row_ + 1

        plt.suptitle('db-%s' %str(ii))
        plt.savefig(os.path.join(sourceDir, 'Cumulants', 'Cumulants-db-%s.png' %str(ii)))
        plt.savefig(os.path.join(sourceDir, 'Cumulants', 'Cumulants-db-%s.pdf' %str(ii)))
        plt.close()
        plt.clf()

################## BarPlot for Cumulants Statistics

if not os.path.exists(os.path.join(sourceDir, 'Cumulants')):
    os.makedirs(os.path.join(sourceDir, 'Cumulants'))

### Creating a Table Here
for ii in np.arange(2,6):
    if not os.path.isfile(os.path.join(sourceDir, 'Cumulants', 'Comparison_Cumulants-db-%s.csv' %str(ii))):

        OrEq_Stats1= pd.read_csv('stats_OrEq-Wilcox-Cumulant-1-db%s.csv' %str(ii), index_col=0)
        OrEq_Pval1= pd.read_csv('pval_OrEq-Wilcox-Cumulant-1-db%s.csv' %str(ii), index_col=0)
        OrSh_Stats1= pd.read_csv('stats_OrSh-Wilcox-Cumulant-1-db%s.csv' %str(ii), index_col=0)
        OrSh_Pval1= pd.read_csv('pval_OrSh-Wilcox-Cumulant-1-db%s.csv' %str(ii), index_col=0)
        OrEq_Stats2= pd.read_csv('stats_OrEq-Wilcox-Cumulant-2-db%s.csv' %str(ii), index_col=0)
        OrEq_Pval2= pd.read_csv('pval_OrEq-Wilcox-Cumulant-2-db%s.csv' %str(ii), index_col=0)
        OrSh_Stats2= pd.read_csv('stats_OrSh-Wilcox-Cumulant-2-db%s.csv' %str(ii), index_col=0)
        OrSh_Pval2= pd.read_csv('pval_OrSh-Wilcox-Cumulant-2-db%s.csv' %str(ii), index_col=0)

        forTable = pd.DataFrame([], index=['1:Os>Es', '1:Os-Ss', '2:Os>Es', '2:Os-Ss'], columns=['Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08'])
        for clmns in OrEq_Stats1.columns.values:
################# For equalized data
            if OrEq_Pval1.loc['UpperBeta', clmns] < 0.001:
                forTable.loc['1:Os>Es', clmns] = str(int(OrEq_Stats1.loc['UpperBeta', clmns])) + '*'
            elif OrEq_Pval1.loc['UpperBeta', clmns] < 0.01:
                forTable.loc['1:Os>Es', clmns] = str(int(OrEq_Stats1.loc['UpperBeta', clmns])) + '#'
            else:
                forTable.loc['1:Os>Es', clmns] = int(OrEq_Stats1.loc['UpperBeta', clmns])

            if OrEq_Pval2.loc['UpperBeta', clmns] < 0.001:
                forTable.loc['2:Os>Es', clmns] = str(int(OrEq_Stats2.loc['UpperBeta', clmns])) + '*'
            elif OrEq_Pval2.loc['UpperBeta', clmns] < 0.01:
                forTable.loc['2:Os>Es', clmns] = str(int(OrEq_Stats2.loc['UpperBeta', clmns])) + '#'
            else:
                forTable.loc['2:Os>Es', clmns] = int(OrEq_Stats2.loc['UpperBeta', clmns])

################# For shuffled data
            if OrSh_Pval1.loc['UpperBeta', clmns] < 0.001:
                forTable.loc['1:Os-Ss', clmns] = str(int(OrSh_Stats1.loc['UpperBeta', clmns])) + '*'
            elif OrSh_Pval1.loc['UpperBeta', clmns] < 0.01:
                forTable.loc['1:Os-Ss', clmns] = str(int(OrSh_Stats1.loc['UpperBeta', clmns])) + '#'                
            else:
                forTable.loc['1:Os-Ss', clmns] = int(OrSh_Stats1.loc['UpperBeta', clmns])

            if OrSh_Pval2.loc['UpperBeta', clmns] < 0.001:
                forTable.loc['2:Os-Ss', clmns] = str(int(OrSh_Stats2.loc['UpperBeta', clmns])) + '*'
            elif OrSh_Pval2.loc['UpperBeta', clmns] < 0.01:
                forTable.loc['2:Os-Ss', clmns] = str(int(OrSh_Stats2.loc['UpperBeta', clmns])) + '#'  
            else:
                forTable.loc['2:Os-Ss', clmns] = int(OrSh_Stats2.loc['UpperBeta', clmns])

        forTable.to_csv(os.path.join(sourceDir, 'Cumulants', 'Comparison_Cumulants-db-%s.csv' %str(ii)))

for ii in np.arange(2,6):
    if not os.path.isfile(os.path.join(sourceDir, 'Cumulants', 'Comparison_Cumulants-db-%s.png' %str(ii))):
        fig = plt.figure(figsize=(20, 10.5))
        nrow = 7
        ncol = 5
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), wspace=0.15, hspace=0.14, left=0.04, right=0.995, bottom=0.075, top=0.94)
        row_ = 0
        col_ = 0

        OrEq_Stats1= pd.read_csv('stats_OrEq-Wilcox-Cumulant-1-db%s.csv' %str(ii), index_col=0)
        OrEq_Pval1= pd.read_csv('pval_OrEq-Wilcox-Cumulant-1-db%s.csv' %str(ii), index_col=0)
        OrSh_Stats1= pd.read_csv('stats_OrSh-Wilcox-Cumulant-1-db%s.csv' %str(ii), index_col=0)
        OrSh_Pval1= pd.read_csv('pval_OrSh-Wilcox-Cumulant-1-db%s.csv' %str(ii), index_col=0)
        OrEq_Stats2= pd.read_csv('stats_OrEq-Wilcox-Cumulant-2-db%s.csv' %str(ii), index_col=0)
        OrEq_Pval2= pd.read_csv('pval_OrEq-Wilcox-Cumulant-2-db%s.csv' %str(ii), index_col=0)
        OrSh_Stats2= pd.read_csv('stats_OrSh-Wilcox-Cumulant-2-db%s.csv' %str(ii), index_col=0)
        OrSh_Pval2= pd.read_csv('pval_OrSh-Wilcox-Cumulant-2-db%s.csv' %str(ii), index_col=0)
        OrEq_Stats3= pd.read_csv('stats_OrEq-Wilcox-Cumulant-3-db%s.csv' %str(ii), index_col=0)
        OrEq_Pval3= pd.read_csv('pval_OrEq-Wilcox-Cumulant-3-db%s.csv' %str(ii), index_col=0)
        OrSh_Stats3= pd.read_csv('stats_OrSh-Wilcox-Cumulant-3-db%s.csv' %str(ii), index_col=0)
        OrSh_Pval3= pd.read_csv('pval_OrSh-Wilcox-Cumulant-3-db%s.csv' %str(ii), index_col=0)

        for frq_ in frqBands:
            #for statArr in list([Or_hPeak, Eq_hPeak, Sh_hPeak, Or_Width, Eq_Width, Sh_Width, Or_Asy, Eq_Asy, Sh_Asy]):
            ax = plt.subplot(gs[row_,0])
            newDf = pd.DataFrame([], columns=['Orig>Equal', 'Orig-Shuf'])
            newDf['Orig>Equal'] = OrEq_Stats1.loc[frq_]
            newDf['Orig-Shuf'] = OrSh_Stats1.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Cumulant-1')
            if row_ != nrow-1:
                ax.set_xticklabels('')
            ax.legend('')
            ax.set_ylabel(frq_)

            ax = plt.subplot(gs[row_,1])
            newDf = pd.DataFrame([], columns=['Orig>Equal', 'Orig-Shuf'])
            newDf['Orig>Equal'] = OrEq_Stats2.loc[frq_]
            newDf['Orig-Shuf'] = OrSh_Stats2.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Stats: Cumulant-2')
            if row_ != nrow-1:
                ax.set_xticklabels('')
            ax.legend('')
  
            ax = plt.subplot(gs[row_,2])
            newDf = pd.DataFrame([], columns=['Orig>Equal', 'Orig-Shuf'])
            newDf['Orig>Equal'] = OrEq_Stats3.loc[frq_]
            newDf['Orig-Shuf'] = OrSh_Stats3.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Stats: Cumulant-3')
            if row_ != nrow-1:
                ax.set_xticklabels('')
            ax.legend('')

            ax = plt.subplot(gs[row_,3])
            newDf = pd.DataFrame([], columns=['Orig>Equal', 'Orig-Shuf'])
            newDf['Orig>Equal'] = OrEq_Pval1.loc[frq_]
            newDf['Orig-Shuf'] = OrSh_Pval1.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Pval: Cumulant-1')
            if row_ != nrow-1:
                ax.set_xticklabels('')
            ax.legend('')

            ax = plt.subplot(gs[row_,4])
            newDf = pd.DataFrame([], columns=['Orig>Equal', 'Orig-Shuf'])
            newDf['Orig>Equal'] = OrEq_Pval2.loc[frq_]
            newDf['Orig-Shuf'] = OrSh_Pval2.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Pval: Cumulant-2')
            if row_ != nrow-1:
                ax.set_xticklabels('')
            ax.legend('')

            '''ax = plt.subplot(gs[row_,5])
            newDf = pd.DataFrame([], columns=['Orig>Equal', 'Orig-Shuf'])
            newDf['Orig>Equal'] = OrEq_Pval3.loc[frq_]
            newDf['Orig-Shuf'] = OrSh_Pval3.loc[frq_]
            newDf.plot.bar(ax=ax)
            if row_ == 0:
                ax.set_title('Pval: Cumulant-3')
            if row_ != nrow-1:
                ax.set_xticklabels('')
                ax.legend('')
            else:
                ax.legend(loc=2, fontsize=8)
            row_ = row_ + 1'''

        plt.suptitle('db-%s' %str(ii))
        plt.savefig(os.path.join(sourceDir, 'Cumulants', 'Comparison_Cumulants-db-%s.png' %str(ii)))
        plt.savefig(os.path.join(sourceDir, 'Cumulants', 'Comparison_Cumulants-db-%s.pdf' %str(ii)))
        plt.close()
        plt.clf()


    if not os.path.isfile(os.path.join(sourceDir, 'Cumulants', 'UpperBeta_Comparison_Cumulants-db-%s.png' %str(ii))):
        frq_ = 'UpperBeta'
        plt.figure(figsize=(20, 10.5))
        #plt.rcParams['font.size'] = 30
        nrow = 1
        ncol = 2
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), wspace=0.15, hspace=0.14, left=0.04, right=0.995, bottom=0.2, top=0.94)
        row_ = 0
        col_ = 0

        #for statArr in list([Or_hPeak, Eq_hPeak, Sh_hPeak, Or_Width, Eq_Width, Sh_Width, Or_Asy, Eq_Asy, Sh_Asy]):
        ax = plt.subplot(gs[row_,0])
        newDf = pd.DataFrame([], columns=['Orig>Equal', 'Orig-Shuf'])
        newDf['Orig>Equal'] = OrEq_Pval1.loc[frq_]
        newDf['Orig-Shuf'] = OrSh_Pval1.loc[frq_]
        newDf.plot.bar(ax=ax, figsize=(20,10))
        ax.set_title('Pval: Cumulant-1', fontsize=25)
        ax.tick_params(axis='both', which='both', labelsize=25)
        ax.legend(ncol=2, fontsize=25, fancybox=True, framealpha=0.6, title='Alternative H',loc='upper left', columnspacing=0.5, markerscale=0.5, title_fontsize=20)            

        ax = plt.subplot(gs[row_,1])
        newDf = pd.DataFrame([], columns=['Orig>Equal', 'Orig-Shuf'])
        newDf['Orig>Equal'] = OrEq_Pval2.loc[frq_]
        newDf['Orig-Shuf'] = OrSh_Pval2.loc[frq_]
        newDf.plot.bar(ax=ax, figsize=(20,10))
        ax.set_title('Pval: Cumulant-2', fontsize=25)
        ax.tick_params(axis='both', which='both', labelsize=25)
        ax.legend(ncol=2, fontsize=25, fancybox=True, framealpha=0.6, title='Alternative H',loc='upper left', columnspacing=0.5, markerscale=0.5, title_fontsize=20)            

        plt.savefig(os.path.join(sourceDir, 'Cumulants', '%s_Comparison_Cumulants-db-%s.png' %(frq_, str(ii))))
        plt.savefig(os.path.join(sourceDir, 'Cumulants', '%s_Comparison_Cumulants-db-%s.pdf' %(frq_, str(ii))))
        plt.close()
        plt.clf()

############# Scaling Graph
if not os.path.exists(os.path.join(sourceDir, 'ScalingExponent')):
    os.makedirs(os.path.join(sourceDir, 'ScalingExponent'))

for ii in np.arange(2,6):
    if not os.path.isfile(os.path.join(sourceDir, 'ScalingExponent', 'ScalingExp_db-%s.png' %str(ii))):
        fig = plt.figure(figsize=(20, 10.5))
        nrow = 7
        ncol = 9
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), hspace=0.3, left=0.04, right=0.99, bottom=0.05, top=0.93)
        row_ = 0
        col_ = 0
        for frq_ in frqBands:
            grpSort = np.sort(glob.glob(os.path.join(sourceDir, frq_, '*')))
            for grp_ in grpSort:
                if '.png' in grp_.split('/')[-1]:
                    continue

                ScFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_scalingExpoent.mat' %str(ii)))
                Eq_ScFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_scalingExpoent.mat' %str(ii)))
                Sh_ScFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_scalingExpoent.mat' %str(ii)))                
                keys_ = [i for i in ScFiles.keys()]
                ScData = ScFiles[keys_[-1]] ## Normal Files
                ScMean = np.median(ScData,0)
                Eq_ScData = Eq_ScFiles[keys_[-1]] ## Equalized Time Files
                Eq_ScMean = np.median(Eq_ScData,0)
                Sh_ScData = Sh_ScFiles[keys_[-1]] ## Shuffled Label File
                Sh_ScMean = np.median(Sh_ScData,0)

                cumulants = scipy.io.loadmat(os.path.join(grp_, 'db%s_cumulants.mat' %str(ii)))
                keys_ = [i for i in cumulants.keys()]
                cumulant = cumulants[keys_[-1]]
                Meancumulant = np.median(cumulant, 0)
                Meancumulant = [np.round(i,2) for i in Meancumulant]
                Eq_cumulants = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_cumulants.mat' %str(ii)))
                keys_ = [i for i in Eq_cumulants.keys()]
                Eq_cumulant = Eq_cumulants[keys_[-1]]
                Eq_Meancumulant = np.median(Eq_cumulant, 0)
                Eq_Meancumulant = [np.round(i,2) for i in Eq_Meancumulant]
                Sh_cumulants = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_cumulants.mat' %str(ii)))
                keys_ = [i for i in Sh_cumulants.keys()]
                Sh_cumulant = Sh_cumulants[keys_[-1]]
                Sh_Meancumulant = np.median(Sh_cumulant, 0)
                Sh_Meancumulant = [np.round(i,2) for i in Sh_Meancumulant]

                ax = plt.subplot(gs[row_,col_])
                '''ax.plot(ScMean, 'r-', label='(%s,%s,%s)' %(str(Meancumulant[0]),str(Meancumulant[1]),str(Meancumulant[2])))
                ax.plot(Eq_ScMean, 'g--', label='(%s,%s,%s)' %(str(Eq_Meancumulant[0]),str(Eq_Meancumulant[1]),str(Eq_Meancumulant[2])))
                ax.plot(Sh_ScMean, 'b--', label='(%s,%s,%s)' %(str(Sh_Meancumulant[0]),str(Sh_Meancumulant[1]),str(Sh_Meancumulant[2])))'''
                
                ax.plot(ScMean, 'r-', label='(%s,%s)' %(str(Meancumulant[0]),str(Meancumulant[1])), linewidth=2)
                ax.plot(Eq_ScMean, 'g--', label='(%s,%s)' %(str(Eq_Meancumulant[0]),str(Eq_Meancumulant[1])), linewidth=2)
                ax.plot(Sh_ScMean, 'b.--', label='(%s,%s)' %(str(Sh_Meancumulant[0]),str(Sh_Meancumulant[1])), linewidth=2)
                #ax.set_xticks([1,2,3,4,5], ['-5','-5','0','5','10'])
                ax.set_xticklabels(['-5','-5','0','5','10'])                

                del ScData
                del ScMean
                del Eq_ScData
                del Eq_ScMean
                del Sh_ScData
                del Sh_ScMean
                del cumulant
                del Meancumulant
                del Eq_cumulant
                del Eq_Meancumulant
                del Sh_cumulant
                del Sh_Meancumulant
                ax.legend(loc=4, prop={'size': 8})
                if col_ == 0:
                    ax.set_ylabel(frq_)
                if row_ == 0:
                    ax.set_title(grp_.split('/')[-1].split('_')[0])
                col_ = col_ + 1
            row_ = row_ + 1
            col_ = 0

        plt.suptitle('db-%s' %str(ii))
        plt.savefig(os.path.join(sourceDir, 'ScalingExponent', 'ScalingExp_db-%s.png' %str(ii)))
        plt.savefig(os.path.join(sourceDir, 'ScalingExponent', 'ScalingExp_db-%s.pdf' %str(ii)))
        plt.close()
        plt.clf()

    if os.path.isfile(os.path.join(sourceDir, 'ScalingExponent', 'UpperBeta_ScalingExp_db-%s.png' %str(ii))):
        frq_ = 'UpperBeta'
        plt.figure(figsize=(20, 10.5))
        #plt.rcParams['font.size'] = 30
        nrow = 3
        ncol = 3
        gs = gridspec.GridSpec(nrow, ncol, height_ratios=np.ones(nrow), width_ratios=np.ones(ncol), wspace=0.13, hspace=0.3, left=0.06, right=0.99, bottom=0.05, top=0.96)
        row_ = 0
        col_ = 0
        frq_ = 'UpperBeta'
        grpSort = np.sort(glob.glob(os.path.join(sourceDir, frq_, '*')))
        for pltCnt, grp_ in enumerate(grpSort[1:]):
            if '.png' in grp_.split('/')[-1]:
                continue

            ScFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_scalingExpoent.mat' %str(ii)))
            Eq_ScFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_scalingExpoent.mat' %str(ii)))
            Sh_ScFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_scalingExpoent.mat' %str(ii)))                
            keys_ = [i for i in ScFiles.keys()]
            ScData = ScFiles[keys_[-1]] ## Normal Files
            ScMean = np.median(ScData,0)
            Eq_ScData = Eq_ScFiles[keys_[-1]] ## Equalized Time Files
            Eq_ScMean = np.median(Eq_ScData,0)
            Sh_ScData = Sh_ScFiles[keys_[-1]] ## Shuffled Label File
            Sh_ScMean = np.median(Sh_ScData,0)

            cumulants = scipy.io.loadmat(os.path.join(grp_, 'db%s_cumulants.mat' %str(ii)))
            keys_ = [i for i in cumulants.keys()]
            cumulant = cumulants[keys_[-1]]
            Meancumulant = np.median(cumulant, 0)
            Meancumulant = [np.round(i,2) for i in Meancumulant]
            Eq_cumulants = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_cumulants.mat' %str(ii)))
            keys_ = [i for i in Eq_cumulants.keys()]
            Eq_cumulant = Eq_cumulants[keys_[-1]]
            Eq_Meancumulant = np.median(Eq_cumulant, 0)
            Eq_Meancumulant = [np.round(i,2) for i in Eq_Meancumulant]
            Sh_cumulants = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_cumulants.mat' %str(ii)))
            keys_ = [i for i in Sh_cumulants.keys()]
            Sh_cumulant = Sh_cumulants[keys_[-1]]
            Sh_Meancumulant = np.median(Sh_cumulant, 0)
            Sh_Meancumulant = [np.round(i,2) for i in Sh_Meancumulant]

            print(row_, col_)
            ax = plt.subplot(gs[row_,col_])
            '''ax.plot(ScMean, 'r-', label='(%s,%s,%s)' %(str(Meancumulant[0]),str(Meancumulant[1]),str(Meancumulant[2])))
            ax.plot(Eq_ScMean, 'g--', label='(%s,%s,%s)' %(str(Eq_Meancumulant[0]),str(Eq_Meancumulant[1]),str(Eq_Meancumulant[2])))
            ax.plot(Sh_ScMean, 'b--', label='(%s,%s,%s)' %(str(Sh_Meancumulant[0]),str(Sh_Meancumulant[1]),str(Sh_Meancumulant[2])))'''

            ax.plot(ScMean, 'r-', label='(%s,%s)' %(str(Meancumulant[0]),str(Meancumulant[1])), linewidth=4)
            ax.plot(Eq_ScMean, 'g--', label='(%s,%s)' %(str(Eq_Meancumulant[0]),str(Eq_Meancumulant[1])), linewidth=4)
            ax.plot(Sh_ScMean, 'b.--', label='(%s,%s)' %(str(Sh_Meancumulant[0]),str(Sh_Meancumulant[1])), linewidth=4)
            
            ax.tick_params(axis='both', which='both', labelsize=22)
            ax.set_xticklabels(['-7','-5','-3','-1','1','3','5','7'])

            del ScData
            del ScMean
            del Eq_ScData
            del Eq_ScMean
            del Sh_ScData
            del Sh_ScMean
            del cumulant
            del Meancumulant
            del Eq_cumulant
            del Eq_Meancumulant
            del Sh_cumulant
            del Sh_Meancumulant
    
            if col_ == 0:
                ax.set_ylabel('Scaling Exponent', fontsize=25)
            ax.set_title(grp_.split('/')[-1].split('_')[0], fontsize=25)
            col_ = col_ + 1

            if (pltCnt+1) % 3 == 0:
                col_ = 0
                row_ = row_ + 1                        
            
            ax.legend(fontsize=20, fancybox=True, framealpha=0.6, title='Cumulants-1 & 2', columnspacing=0.5, markerscale=0.5, title_fontsize=20)            
            
        plt.savefig(os.path.join(sourceDir, 'ScalingExponent', '%s_ScalingExp_db-%s.png' %(frq_, str(ii))))
        plt.savefig(os.path.join(sourceDir, 'ScalingExponent', '%s_ScalingExp_db-%s.pdf' %(frq_, str(ii))))
        plt.close()
        plt.clf()

pdb.set_trace()
###################### Cumulants Comparison ###############
ii = 2
for grp in ['Baseline', 'Group-01','Group-02','Group-03','Group-04','Group-05','Group-06','Group-07','Group-08']:        
    grp_ = glob.glob(os.path.join(sourceDir, 'UpperBeta', grp+'*'))[0]
    if '.png' in grp_.split('/')[-1]:
        continue

    UB_hFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_cumulants.mat' %str(ii)))
    keys_ = [i for i in UB_hFiles.keys()]
    pdb.set_trace()
    UB_hFiles = UB_hFiles[keys_[-1]][:, 1]
    UB_Eq_hFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_cumulants.mat' %str(ii)))    
    UB_Eq_hFiles = UB_Eq_hFiles[keys_[-1]][:, 1]
    UB_Sh_hFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_cumulants.mat' %str(ii)))    
    UB_Sh_hFiles = UB_Sh_hFiles[keys_[-1]][:, 1]

    hFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_singularityExponent.mat' %str(ii)))
    Eq_hFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_singularityExponent.mat' %str(ii)))
    Sh_hFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_singularityExponent.mat' %str(ii)))

    tmp_ = scipy.io.loadmat(os.path.join(grp_, 'filesProcessed.mat'))
    UB_hSubFiles = [j[0][0] for j in tmp_[[i for i in tmp_.keys()][-1]]]
    tmp_ = scipy.io.loadmat(os.path.join(grp_, 'EQ_filesProcessed.mat'))
    UB_Eq_hSubFiles = [j[0][0] for j in tmp_[[i for i in tmp_.keys()][-1]]    ]
    tmp_ = scipy.io.loadmat(os.path.join(grp_, 'SH_filesProcessed.mat'))
    UB_Sh_hSubFiles = [j[0][0] for j in tmp_[[i for i in tmp_.keys()][-1]]    ]
    
    for frq_ in frqBands:

        if frq_ == 'AllBands':
            continue

        grp_ = glob.glob(os.path.join(sourceDir, frq_, grp+'*'))[0]
        if '.png' in grp_.split('/')[-1]:
            continue

        hFiles = scipy.io.loadmat(os.path.join(grp_, 'db%s_cumulants.mat' %str(ii)))
        keys_ = [i for i in hFiles.keys()]
        hFiles = hFiles[keys_[-1]][:, 1]
        Eq_hFiles = scipy.io.loadmat(os.path.join(grp_, 'EQ_db%s_cumulants.mat' %str(ii)))    
        Eq_hFiles = Eq_hFiles[keys_[-1]][:, 1]
        Sh_hFiles = scipy.io.loadmat(os.path.join(grp_, 'SH_db%s_cumulants.mat' %str(ii)))    
        Sh_hFiles = Sh_hFiles[keys_[-1]][:, 1]

        tmp_ = scipy.io.loadmat(os.path.join(grp_, 'filesProcessed.mat'))
        hSubFiles = [j[0][0] for j in tmp_[[i for i in tmp_.keys()][-1]]]
        tmp_ = scipy.io.loadmat(os.path.join(grp_, 'EQ_filesProcessed.mat'))
        Eq_hSubFiles = [j[0][0] for j in tmp_[[i for i in tmp_.keys()][-1]]    ]
        tmp_ = scipy.io.loadmat(os.path.join(grp_, 'SH_filesProcessed.mat'))
        Sh_hSubFiles = [j[0][0] for j in tmp_[[i for i in tmp_.keys()][-1]]    ]

        if frq_ == 'Theta' and grp == 'Baseline':            
            hSubFiles = [j.split('corrected')[0]+'corrected_7Sec_'+'_'.join(j.split('_')[4:]) for j in hSubFiles] 
            Eq_hSubFiles = [j.split('corrected')[0]+'corrected_7Sec_'+'_'.join(j.split('_')[5:]) for j in Eq_hSubFiles] 
            Sh_hSubFiles = [j.split('corrected')[0]+'corrected_7Sec_'+'_'.join(j.split('_')[5:]) for j in Sh_hSubFiles] 
     
        widthFrame = pd.DataFrame(-1, index=UB_hSubFiles, columns=['UBOrig', frq_+'Orig'])
        for sub_ in UB_hSubFiles:
            if (sub_ in hSubFiles):
                print(sub_)
                idxs = np.where(sub_ == np.array(UB_hSubFiles))[0]
                idxs2 = np.where(sub_ == np.array(hSubFiles))[0]                                
                widthFrame.loc[sub_, 'UBOrig'] = UB_hFiles[idxs]
                widthFrame.loc[sub_, frq_+'Orig'] = hFiles[idxs2]                
        widthFrame.drop(widthFrame.index.values[np.where(np.sum(widthFrame.values, axis=1)==-2)[0]], inplace=True)
        widthFrame.to_csv(os.path.join(sourceDir, 'cumulantComp', grp+'_upperBeta-'+frq_+'-Orig.csv'))
        #scipy.stats.ttest_rel(widthFrame['UBOrig'], widthFrame[frq_+'Orig'], alternative='less')

        Eq_widthFrame = pd.DataFrame(-1, index=UB_Eq_hSubFiles, columns=['UBEqual', frq_+'Equal'])
        for sub_ in UB_Eq_hSubFiles:
            if (sub_ in Eq_hSubFiles):
                print(sub_)
                idxs = np.where(sub_ == np.array(UB_Eq_hSubFiles))[0]
                idxs2 = np.where(sub_ == np.array(Eq_hSubFiles))[0]
                Eq_widthFrame.loc[sub_, 'UBEqual'] = UB_Eq_hFiles[idxs]
                Eq_widthFrame.loc[sub_, frq_+'Equal'] = Eq_hFiles[idxs2]
        Eq_widthFrame.drop(Eq_widthFrame.index.values[np.where(np.sum(Eq_widthFrame.values, axis=1)==-2)[0]], inplace=True)
        Eq_widthFrame.to_csv(os.path.join(sourceDir, 'cumulantComp', grp+'_upperBeta-'+frq_+'-Equal.csv'))
        #scipy.stats.ttest_rel(Eq_widthFrame['UBEqual'], Eq_widthFrame[frq_+'Equal'], alternative='greater')
        
        Sh_widthFrame = pd.DataFrame(-1, index=UB_Sh_hSubFiles, columns=['UBShuff', frq_+'Shuff'])
        for sub_ in UB_Sh_hSubFiles:
            if (sub_ in Sh_hSubFiles):
                print(sub_)
                idxs = np.where(sub_ == np.array(UB_Sh_hSubFiles))[0]
                idxs2 = np.where(sub_ == np.array(Sh_hSubFiles))[0]
                Sh_widthFrame.loc[sub_, 'UBShuff'] = UB_Sh_hFiles[idxs]
                Sh_widthFrame.loc[sub_, frq_+'Shuff'] = Sh_hFiles[idxs2]
        Sh_widthFrame.drop(Sh_widthFrame.index.values[np.where(np.sum(Sh_widthFrame.values, axis=1)==-2)[0]], inplace=True)
        Sh_widthFrame.to_csv(os.path.join(sourceDir, 'cumulantComp', grp+'_upperBeta-'+frq_+'-Shuff.csv'))
        #scipy.stats.ttest_rel(Sh_widthFrame['UBShuff'], Sh_widthFrame[frq_+'Shuff'], alternative='greater')