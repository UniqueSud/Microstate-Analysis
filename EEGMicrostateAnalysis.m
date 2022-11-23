addpath('/media/forAll/Processed_Emotions/ToTransfer_Zuddler/eeglab2019_2/')
addpath('/home/zuddler/Downloads/logfit/')
%addpath('/home/shrikant/ForMSAAnalysis/')
%% 3 Tutorial: EEG microstates analysis on spontaneous EEG Data
%
% This script executes the analysis steps of the Tutorial described in
% detail in section 3.3 to 3.8 of:htop
% Poulsen, A. T., Pedroni, A., Langer, N., & Hansen, L. K. (2018).
% Microstate EEGlab toolbox: An introductionary guide. bioRxiv.
%
% Authors:
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Zurich, Psychologisches Institut, Methoden der
% Plastizitaetsforschung.
close all
clc
clear all 
bandInform = 'Alpha';
rootDir = fullfile('/media/forAll/Processed_Emotions/ToTransfer_Zuddler/ForMSAAnalysis', bandInform);
cd(rootDir)
%groups = dir('Baseline10Sec*');
%groups = dir('Baseline*');
%groups = dir('Group-0*');
groups = dir('Neutral*');
gi = 1;
%while gi < length(groups)
for gi = 1 : length(groups)
    subName_ = groups(gi).name;
    noIter = 100;
    % start EEGLAB to load all dependent paths
    %% set the path to the directory with the EEG files
    % change this path to the folder where the EEG files are saved

    %samplesToRemove = [9, 10, 12, 13, 16, 17, 23, 24, 29, 32, 35, 36, 43, 44, 45]; %HappyAmusedJoyous more than 100 peaks
    %samplesToRemove = [12, 13, 16, 24, 29, 32, 36, 43, 45]; %Baseline more than 50 peaks
    %samplesToRemove = [12, 18, 37]; %neutral more than 50 peaks
    %samplesToRemove = [2, 7,17,19,37]; %TenseFrustratedDistress more than 50 peaks
    %samplesToRemove = [3, 6,12,13]; %Angry more than 50 peaks
    samplesToRemove = [4,  6, 22, 27, 28, 34, 62, 63]; %OnlyBaseline more than 50 peaks
    %samplesToRemove = [28, 41]; % For sad
    %samplesToRemove = [3]; % For AlarmedStartled
    %samplesToRemove = [6, 47]; % For HateDisgust
    %subName_ = 'Baseline';
    EEGdir = fullfile(rootDir,subName_);
    targetdir = fullfile(EEGdir, 'MSAAnalysis_WithMoreThan50Peaks');       

    if ~ isdir(targetdir)
        mkdir(targetdir)
    end
    targetPath = fullfile(targetdir, 'fitBack');    
    if ~ isdir(targetPath)
        mkdir(targetPath)
    end    
    cd(EEGdir)
    EEGFiles = dir(fullfile(EEGdir, 'corrected*.set'));
    GFPfiles = dir(fullfile(targetdir, 'GFPEEGMSA*.set'));

    eeglab
    %% 3.3 Data selection and aggregation
    %% 3.3.1 Loading datasets in EEGLAB
    noFiles = 0;
    %% It has to be done manually after looking at the peakIndx in pop_micro_selectdata module.
    %%
    noFiles = 0;

    if length(GFPfiles) ~= noIter
        for i=1:length(EEGFiles)
            if sum(i == samplesToRemove) == 1
                disp(strcat('Sample = ', num2str(i), ' Should not be considered................................................'))
                continue
            else
                noFiles = noFiles + 1;
                EEG = pop_loadset('filename',EEGFiles(i).name,'filepath',EEGdir);
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                eeglab redraw % updates EEGLAB datasets            
            end
        end

    end

    if noFiles == 0
        noFiles = length(EEGFiles)-length(samplesToRemove);
    end

    if noFiles < 30
        error('Number of samples are very less')
    end

    for grdIdx = 1 : noIter 
        if ~ isfile(fullfile(targetdir, strcat('GFPEEGMSA_', num2str(grdIdx), '.set')))            % retrieve a list of all EEG Files in EEGdir

            %% 3.3.2 Select data for microstate analysis
            %[EEG, ALLEEG] = pop_micro_selectdata( EEG, ALLEEG, 'datatype', 'spontaneous', 'avgref', 1, 'normalise', 1, 'MinPeakDist', 10, 'Npeaks', 5000, 'GFPthresh', 0, 'dataset_idx', 1:40 );

            [EEG, ALLEEG] = pop_micro_selectdata( EEG, ALLEEG, 'datatype', 'spontaneous',...
            'avgref', 1, ...
            'normalise', 0, ...
            'MinPeakDist', 10, ...
            'Npeaks', 10000, ...
            'GFPthresh', 0, ...
            'dataset_idx', 1:noFiles);
            % store data in a new EEG structure
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            eeglab redraw % updates EEGLAB datasets
            %% 3.4 Microstate segmentation
            % select the "GFPpeak" dataset and make it the active set
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, noFiles,'retrieve', noFiles+1,'study',0);
            eeglab redraw

            % Perform the microstate segmentation
            EEG = pop_micro_segment( EEG, 'algorithm', 'modkmeans', ...
            'sorting', 'Global explained variance', ...
            'Nmicrostates', 2:6, ... % 'Nmicrostates', 2:12
            'verbose', 1, ...
            'normalise', 1, ...
            'Nrepetitions', 50, ...
            'max_iterations', 1000, ...
            'threshold', 1e-08, ...
            'fitmeas', 'CV',...
            'optimised',1);
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            %% 3.5.1 Plot microstate prototype topographies
            %%figure;MicroPlotTopo( EEG, 'plot_range', [] );
            %% 3.5.2 Select active number of microstates
            %EEG = pop_micro_selectNmicro( EEG);
            %[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            %% Import microstate prototypes from other dataset to the datasets that should be
            %back-fitted
            % note that dataset number 5 is the GFPpeaks dataset with the microstate
            % prototypes

            %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',noFiles+1,'study',0);
            pop_saveset(EEG, 'filename',strcat('GFPEEGMSA_', num2str(grdIdx)),'filepath', targetdir);    
            ALLEEG = pop_delset(ALLEEG, noFiles+1);
            eeglab redraw

        end
    end
    
    ALLEEG = [];
    CURRENTSET = 0;
    CURRENTSTUDY = 0;    
    
    maxfl_ = dir(fullfile(targetdir, 'MAX_GFPEEGMSA*.set'));
    if ~ length(maxfl_)
        max_ = 0;
        for grdIdx = 1 : noIter
            fileName = fullfile(targetdir, strcat('GFPEEGMSA_', num2str(grdIdx), '.set'));
            EEG = pop_loadset('filename',fileName);

            if sum(EEG.microstate.Res.GEV) > max_
                [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);            
                max_ = sum(EEG.microstate.Res.GEV); %% max_ is very less. Check Here.
                maxFile = strcat('GFPEEGMSA_', num2str(grdIdx), '.set');
            end    

            if grdIdx < noIter
                close all
                clear EEG
            end        
            disp('================End Here===============')
        end
        pop_saveset(ALLEEG(1,1), 'filename', strcat('MAX_', maxFile),'filepath', targetdir);    
    else
        maxFile = maxfl_.name;
        EEG = pop_loadset('filename', maxfl_.name,'filepath', targetdir); 
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);
    end
    eeglab redraw;
    if ~strmatch(maxFile, ALLEEG(1,1).filename)
        error('file is not matching')
    end
    cd(EEGdir)
    EEGFiles = dir(fullfile(EEGdir, '*.set'));
    %% 3.3 Data selection and aggregation
    %% 3.3.1 Loading datasets in EEGLAB

    %for micInd_ = 2 : 9
    for micInd_ = 2 : 5
        %% Actually We have calculated the Microstates starting from 2 clusters. So, if I am saying micInd_=2, actually I am referring to Microstates with 3 clusters.
        fitbackSetPath = fullfile(targetPath, 'OnlySetFiles', strcat('With_',num2str(micInd_+1),'_MS'));
        if ~ isdir(fitbackSetPath)
            mkdir(fitbackSetPath)
        end   
        if ~ strmatch(ALLEEG(1,1).filename, maxFile, 'exact')
            error('Not with the best Microstate file')
        end

        ALLEEG(1,1).microstate.prototypes = ALLEEG(1,1).microstate.Res.A_all{micInd_};        
        if ~isfile(strcat('Topoplot_', num2str(micInd_+1), '.mat'))
            tp_ = ALLEEG(1,1).microstate.prototypes;
            save(strcat('Topoplot_', num2str(micInd_+1)), 'tp_')

            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',1,'study',0);    
            eeglab redraw;
            if ~ strmatch(ALLEEG(1,1).filename, maxFile, 'exact')
                error('Not with the best Microstate file')
            end        
            s = get(0, 'ScreenSize'); 
            figure('Position', [0 0 s(3) s(4)]); 
            MicroPlotTopo( EEG, 'plot_range', micInd_+1 );    
            saveas(gca, strcat('Topoplot-', num2str(micInd_+1), subName_), 'jpeg')
        end    
        close all;

        subEEGGEVtotal = [];
        subEEGGEV = [];

        if ~ isfile(fullfile(targetdir, strcat('subjectWiseGEVTotal_for_Microstates_', num2str(micInd_+1), '.mat')))
            labelisnotthere = 0;
            for i=1:length(EEGFiles)           

                EEG = pop_loadset('filename',EEGFiles(i).name,'filepath',EEGdir);    
                if ~ strmatch(EEG.filename, EEGFiles(i).name, 'exact')
                    error('Not with the best Microstate file')
                end            
                if ~ strmatch(ALLEEG(1,1).filename, maxFile, 'exact')
                    error('Not with the best Microstate file')
                end                    
                %eeglab redraw
                EEG = pop_micro_import_proto( EEG, ALLEEG, 1);
                % 3.6 Back-fit microstates on EEG
                EEG = pop_micro_fit( EEG, 'polarity', 0 );
                % 3.7 Temporally smooth microstates labels
                EEG = pop_micro_smooth( EEG, 'label_type', 'backfit', 'smooth_type', 'reject segments', 'minTime', 30, 'polarity', 0 );        
                EEG = pop_micro_stats( EEG, 'label_type', 'backfit', 'polarity', 0 );  
                subEEGGEVtotal = [subEEGGEVtotal EEG.microstate.stats.GEVtotal];
                subEEGGEV = [subEEGGEV; EEG.microstate.stats.GEV];

                MSsequence = EEG.microstate.fit.labels;

                unq_ = unique(MSsequence);
                newInd = micInd_;
                if unq_(1)~=1
                    newInd = newInd-1;
                end
    % Why +1: Actually We have calculated the Microstates starting from 2 clusters. So, if I am saying micInd_=2, actually I am referring to Microstates with 3 clusters.                        
                newInd = newInd-((micInd_+1)-unq_(end));
                %if unq_(end)~=(micInd_+1)
                %    newInd = newInd-1;
                %end            
                if newInd ~= micInd_
                    labelisnotthere = labelisnotthere + 1;
                end
                if sum(unq_(2:end)-unq_(1:end-1)) ~= newInd
                    error('Some Problem In Calculation')
                end

                if ~ strmatch(EEG.filename, EEGFiles(i).name, 'exact')
                    error('Not with the best Microstate file')
                end            
                if ~ strmatch(ALLEEG(1,1).filename, maxFile, 'exact')
                    error('Not with the best Microstate file')
                end

                pop_saveset(EEG, 'filename',strcat('Microstates_', num2str(micInd_+1), '_EEGMSA_', EEGFiles(i).name),'filepath', fitbackSetPath);  
                %[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
                %[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                %eeglab redraw % updates EEGLAB datasets
                close all

            end
            if labelisnotthere > (length(EEGFiles)/2)
                error('Some Problem In Calculation')
            end
            save(fullfile(targetdir, strcat('subjectWiseGEVTotal_for_Microstates_', num2str(micInd_+1), '.mat')), 'subEEGGEVtotal')
            save(fullfile(targetdir, strcat('subjectWiseGEV_for_Microstates_', num2str(micInd_+1), '.mat')), 'subEEGGEV')
        end
    end

    figure('units','normalized','outerposition',[0 0 1 1]);
    %for micInd_ = 2 : 9
    for micInd_ = 2 : 5
        load(fullfile(targetdir, strcat('subjectWiseGEVTotal_for_Microstates_', num2str(micInd_+1), '.mat')));
        [h,p,ci,stats] = ttest(subEEGGEVtotal, 0.5, 'alpha', 0.001);
        if h == 0
            error('The GEV Value is not significant')
        end
        disp(strcat('==========For Microstates = ', num2str(micInd_+1), '============='))
        disp(mean(subEEGGEVtotal))
        disp(std(subEEGGEVtotal))
        disp(mean(subEEGGEVtotal)/std(subEEGGEVtotal))

        subplot(3, 3, micInd_-1)
        hist(subEEGGEVtotal)
        hold on
        xline(ci(1), 'Color', 'r')
        hold on
        try
            xline(mean(subEEGGEVtotal), 'Color', 'g')
        catch 
            disp('stop here')
        end
        hold on
        xline(ci(2), 'Color', 'r')
        title(strcat('Prototypes-', num2str(micInd_+1)))
        xlabel('Global Explained Variable')
        xticks([0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9])
        ylabel('Number of Subjects')
        %load(strcat('subjectWiseGEV_for_Microstates_', num2str(micInd_+1), '.mat'));
    end
    saveas(gca, strcat('Histogram ', subName_), 'jpeg')
    close all
    if ~ strmatch(ALLEEG(1,1).filename, maxFile, 'exact')
        error('Not with the best Microstate file')
    end

    % Fractal Analysis
    %{
    1. Positive Negative Assignments
    2. Embedded Random Walk
    3. Wavelet and Fractal Analysis 
    %}
%{
    EQFlag = 0;
    MSsequenceCalc(targetPath, EEGFiles, EEG, subName_, EQFlag)
    EQFlag = 1;
    MSsequenceCalc(targetPath, EEGFiles, EEG, subName_, EQFlag)
    EQFlag = 2;
    MSsequenceCalc(targetPath, EEGFiles, EEG, subName_, EQFlag)
    %}

    EQFlag = 0;
    meanTThrs = 25;
    meanTThrs = waveletPlotting(targetPath, ALLEEG, EEGFiles, EEG, subName_, EQFlag, meanTThrs)
    EQFlag = 1;
    meanTThrs = waveletPlotting(targetPath, ALLEEG, EEGFiles, EEG, subName_, EQFlag, meanTThrs)
    EQFlag = 2;
    meanTThrs = waveletPlotting(targetPath, ALLEEG, EEGFiles, EEG, subName_, EQFlag, meanTThrs)    
    %clear all

    clearvars -except gi bandInform
    close all
    
    rootDir = fullfile('/media/forAll/Processed_Emotions/ToTransfer_Zuddler/ForMSAAnalysis', bandInform);
    cd(rootDir)
    groups = dir('Group*');        
end