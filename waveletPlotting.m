function [meanTThrs] = waveletPlotting(targetPath, ALLEEG, EEGFiles, EEG, subName_, EQFlag, meanTThrs)
    %meanTThrs = 25;
    %meanTThrs = 35;
    if ~strmatch(ALLEEG(1,1).filename, 'MAX_')
        error('First file is not Max Record')
    end
    try
        allMicrostates = ALLEEG.microstate.Res.A_all{3,1};
    catch
        disp('')
    end
    
    if EQFlag == 1
        prefix = 'EQ_';
    elseif EQFlag == 2
        prefix = 'SH_';
    else
        prefix = '';
    end
    noLevels = 10;
    classAssGrp = [1,2,3,4;1,3,2,4;1,4,2,3];
    %for micInd_ = 4 : 4 % Don't let
    micInd_ = 4;
    tmpPath = fullfile(targetPath, 'OnlyMat', strcat('With_',num2str(micInd_),'_MS'));
       
    %%%%%%%%%%%%% This segment should go out and in if doing not group wise and groupwise respectively
    totalMT = 0;
    cnt_ = 1;
    %finalPath = fullfile(targetPath, strcat('embeddingGroup-', num2str(grp_)), strcat('FinalPlots_', num2str(meanTThrs)));

    %%%%%%%%%%%%% 
    firstPart = zeros(ALLEEG.nbchan, length(EEGFiles));
    secondPart = zeros(ALLEEG.nbchan, length(EEGFiles));

    %%%%%%%%%%%% Files To Consider
    subjectsNotToCons = {'mit017', 'mit096', 'mit097', 'mit104', 'mit108', 'mit111', 'mit122'}; %% These subjects are not considered in any analysis due to high frequency noise in power ratio
    subjectsNotNoisy = {};
    for i = 1 : length(EEGFiles)            
        tmpName = split(EEGFiles(i).name, '.set');    
        found = 0;
        for subIdx = 1 : length(subjectsNotToCons)
            if strfind(tmpName{1,1}, subjectsNotToCons{1, subIdx})
                found = 1;
                break
            end            
        end
        if found == 0
            subjectsNotNoisy(end+1) = {tmpName{1,1}};
        end
    end
    
    if EQFlag == 0
        while 1            
            disp(meanTThrs)
            partNotIncluded = 0;
            labelisnotthere = 0;
            for i = 1 : length(subjectsNotNoisy)            
                tmpName = split(subjectsNotNoisy(i), '.set');
                tmpName = tmpName{1};
                meanT = load(fullfile(tmpPath, strcat('MeanTime_MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')));
                meanT = meanT.meanT;            
                if meanT > meanTThrs                
                    partNotIncluded = partNotIncluded + 1;
                    partNotIncludedArr(partNotIncluded) = {strcat(tmpName, '-', num2str(meanT))};
                    continue
                end 
            end

            if partNotIncluded > (0.25*length(subjectsNotNoisy))
                meanTThrs = meanTThrs + 1;            
                clear partNotIncludedArr
            else
                break
            end    
        end
    end
    finalPath = fullfile(tmpPath, 'Results', strcat(subName_,'_FinalPlots_', num2str(meanTThrs)));
    if ~ isdir(finalPath)
        mkdir(finalPath)        
    end     

    alphaLocalSingularity2 = [];
    alphaLocalSingularity3 = [];
    alphaLocalSingularity4 = [];
    alphaLocalSingularity5 = [];
    
    subjectsIncluded = {};
    filesProcessed = [];
    for grp_ = 1 : size(classAssGrp,1)
        partNotIncludedArr = {};
        partNotIncluded = 0;
        labelisnotthere = 0;
        for i = 1 : length(subjectsNotNoisy)            
            tmpName = split(subjectsNotNoisy(i), '.set');
            tmpName = tmpName{1};
            meanT = load(fullfile(tmpPath, strcat('MeanTime_MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')));
            meanT = meanT.meanT;   
            if meanT > meanTThrs                
                partNotIncluded = partNotIncluded + 1;
                partNotIncludedArr(partNotIncluded) = {strcat(tmpName, '-', num2str(meanT))};
                continue
            end 
            subjectsIncluded(end+1, 1) = {tmpName};
            totalMT = totalMT + meanT;
            MSsequence = load(fullfile(tmpPath, strcat(prefix,'MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')));
            filesProcessed = [filesProcessed; {strcat(prefix,'MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')}];
            flds = fieldnames(MSsequence);
            MSsequence = getfield(MSsequence, flds{1,1});
            
            %% Here find out the average of microstate plots for firstPart 
            %segment, secondPart segment. And, then see which electrodes are consistently active across all the samples for one emotion. 
            % This will cross-validate our finding about coreAssociative
            % and CoreAffective Networks. Though, it should be done
            % frequency wise. 
            if grp_ == 1
                for seqId_ = 1 : length(MSsequence)
                    if seqId_ < 776
                        firstPart(:, i) = firstPart(:, i) + allMicrostates(:, MSsequence(seqId_));
                        firstCount = seqId_;
                    else                        
                        secondPart(:, i) = secondPart(:, i) + allMicrostates(:, MSsequence(seqId_));
                    end
                end 
                firstPart(:, i) = firstPart(:, i)/firstCount;
                secondPart(:, i) = secondPart(:, i)/(seqId_-firstCount);
            end
            %%
            unq_ = unique(MSsequence);
            newInd = micInd_;
            if unq_(1)~=1
                newInd = newInd-1;
            end
% Why +1: Actually We have calculated the Microstates starting from 2 clusters. So, if I am saying micInd_=2, actually I am referring to Microstates with 3 clusters.                        
            newInd = newInd-((micInd_)-unq_(end));
            %if unq_(end)~=(micInd_+1)
            %    newInd = newInd-1;
            %end            
            if newInd ~= micInd_
                labelisnotthere = labelisnotthere + 1;
            end
            if sum(unq_(2:end)-unq_(1:end-1)) ~= newInd-1
                error('Some Problem In Calculation')
            end
            try   
                embedding = load(fullfile(tmpPath, strcat(num2str(grp_), '_',prefix,'EmbeddingWithArraySum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')));
            catch
                disp('sudhakar')
            end
            embedding2 = load(fullfile(tmpPath, strcat(num2str(grp_), '_',prefix,'EmbeddingWithScalarSum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')));
            embedding = embedding.embedding;
            embedding2 = embedding2.embedding2;

            %{
            previous = MSsequence(1);
            countArr(1) = 1;
            statOrder(1) = MSsequence(1);
            for st_ = 2 : length(MSsequence)
                new = MSsequence(st_);
                if previous == new
                   countArr(end) = countArr(end)+1;                       
                else
                    statOrder(end+1) = new;
                    countArr(end+1) = 1;
                end
                previous = new;
            end
            meanTCalc = round(mean(countArr));
            if meanTCalc ~= meanT
                error('Mean Time is not equal')
            end
            figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(2,2,1);plot(MSsequence);grid on; title('MSsequence');
            subplot(2,2,2);plot(MSsequence);grid on; title('MSsequence');
            subplot(2,2,3);plot(embedding);grid on; title('Embedding');
            subplot(2,2,4);plot(embedding2);grid on; title('Embedding');  
            suptitle(strcat('Mean-Time is-', num2str(meanT)))
            tmpPath2 = fullfile(tmpPath, 'RecreatedImagesDuringFinalPlotting');            
            if ~ isdir(tmpPath2)
                mkdir(tmpPath2)
            end
            saveas(gca, fullfile(tmpPath2, strcat(num2str(grp_), '_MSseqEmbeArraySum_', num2str(micInd_), '_EEGMSA_', tmpName)), 'jpeg')                
            close all
            clear previous countArr statOrder meanTCalc
            %}
            
            [c, l] = wavedec(embedding2, noLevels, 'db2');
            approx = appcoef(c,l,'db2');
            [cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10] = detcoef(c,l,[1 : noLevels]);
            %[cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9] = detcoef(c,l,[1 : noLevels]);
            %[cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10, cd11, cd12] = detcoef(c,l,[1 : noLevels]);
            energy2(cnt_, 1) = log2((sum((abs(cd1).^2))/length(cd1)));
            energy2(cnt_, 2) = log2((sum((abs(cd2).^2))/length(cd2)));
            energy2(cnt_, 3) = log2((sum((abs(cd3).^2))/length(cd3)));
            energy2(cnt_, 4) = log2((sum((abs(cd4).^2))/length(cd4)));
            energy2(cnt_, 5) = log2((sum((abs(cd5).^2))/length(cd5)));
            energy2(cnt_, 6) = log2((sum((abs(cd6).^2))/length(cd6)));
            energy2(cnt_, 7) = log2((sum((abs(cd7).^2))/length(cd7)));
            energy2(cnt_, 8) = log2((sum((abs(cd8).^2))/length(cd8)));            
            energy2(cnt_, 9) = log2((sum((abs(cd9).^2))/length(cd9)));
            energy2(cnt_, 10) = log2((sum((abs(cd10).^2))/length(cd10)));
            %energy2(cnt_, 11) = log2((sum((abs(cd10).^2))/length(cd11)));
            %energy2(cnt_, 12) = log2((sum((abs(cd10).^2))/length(cd12)));
            reconsEmbedding_db2 = waverec(c(1:sum(l(1:10))), l(1:10), 'db2');
            clear c l approx cd1 cd2 cd3 cd4 cd5 cd6 cd7 cd8 cd9 cd10 cd11 cd12
            
            [c, l] = wavedec(embedding2, noLevels, 'db3');
            approx = appcoef(c,l,'db3');
            [cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10] = detcoef(c,l,[1 : noLevels]);
            %[cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9] = detcoef(c,l,[1 : noLevels]);
            %[cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10, cd11, cd12] = detcoef(c,l,[1 : noLevels]);
            energy3(cnt_, 1) = log2((sum((abs(cd1).^2))/length(cd1)));
            energy3(cnt_, 2) = log2((sum((abs(cd2).^2))/length(cd2)));
            energy3(cnt_, 3) = log2((sum((abs(cd3).^2))/length(cd3)));
            energy3(cnt_, 4) = log2((sum((abs(cd4).^2))/length(cd4)));
            energy3(cnt_, 5) = log2((sum((abs(cd5).^2))/length(cd5)));
            energy3(cnt_, 6) = log2((sum((abs(cd6).^2))/length(cd6)));
            energy3(cnt_, 7) = log2((sum((abs(cd7).^2))/length(cd7)));
            energy3(cnt_, 8) = log2((sum((abs(cd8).^2))/length(cd8)));            
            energy3(cnt_, 9) = log2((sum((abs(cd9).^2))/length(cd9)));
            energy3(cnt_, 10) = log2((sum((abs(cd10).^2))/length(cd10)));
            %energy3(cnt_, 11) = log2((sum((abs(cd10).^2))/length(cd11)));
            %energy3(cnt_, 12) = log2((sum((abs(cd10).^2))/length(cd12)));
            reconsEmbedding_db3 = waverec(c(1:sum(l(1:10))), l(1:10), 'db3');
            clear c l approx cd1 cd2 cd3 cd4 cd5 cd6 cd7 cd8 cd9 cd10 cd11 cd12
            
            [c, l] = wavedec(embedding2, noLevels, 'db4');
            approx = appcoef(c,l,'db4');
            [cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10] = detcoef(c,l,[1 : noLevels]);
            %[cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9] = detcoef(c,l,[1 : noLevels]);
            %[cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10, cd11, cd12] = detcoef(c,l,[1 : noLevels]);
            energy4(cnt_, 1) = log2((sum((abs(cd1).^2))/length(cd1)));
            energy4(cnt_, 2) = log2((sum((abs(cd2).^2))/length(cd2)));
            energy4(cnt_, 3) = log2((sum((abs(cd3).^2))/length(cd3)));
            energy4(cnt_, 4) = log2((sum((abs(cd4).^2))/length(cd4)));
            energy4(cnt_, 5) = log2((sum((abs(cd5).^2))/length(cd5)));
            energy4(cnt_, 6) = log2((sum((abs(cd6).^2))/length(cd6)));
            energy4(cnt_, 7) = log2((sum((abs(cd7).^2))/length(cd7)));
            energy4(cnt_, 8) = log2((sum((abs(cd8).^2))/length(cd8)));            
            energy4(cnt_, 9) = log2((sum((abs(cd9).^2))/length(cd9)));
            energy4(cnt_, 10) = log2((sum((abs(cd10).^2))/length(cd10)));
            %energy4(cnt_, 11) = log2((sum((abs(cd10).^2))/length(cd11)));
            %energy4(cnt_, 12) = log2((sum((abs(cd10).^2))/length(cd12)));
            reconsEmbedding_db4 = waverec(c(1:sum(l(1:10))), l(1:10), 'db4');
            clear c l approx cd1 cd2 cd3 cd4 cd5 cd6 cd7 cd8 cd9 cd10 cd11 cd12
            
            [c, l] = wavedec(embedding2, noLevels, 'db5');
            approx = appcoef(c,l,'db5');
            [cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10] = detcoef(c,l,[1 : noLevels]);
            %[cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9] = detcoef(c,l,[1 : noLevels]);
            %[cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10, cd11, cd12] = detcoef(c,l,[1 : noLevels]);
            energy5(cnt_, 1) = log2((sum((abs(cd1).^2))/length(cd1)));
            energy5(cnt_, 2) = log2((sum((abs(cd2).^2))/length(cd2)));
            energy5(cnt_, 3) = log2((sum((abs(cd3).^2))/length(cd3)));
            energy5(cnt_, 4) = log2((sum((abs(cd4).^2))/length(cd4)));
            energy5(cnt_, 5) = log2((sum((abs(cd5).^2))/length(cd5)));
            energy5(cnt_, 6) = log2((sum((abs(cd6).^2))/length(cd6)));
            energy5(cnt_, 7) = log2((sum((abs(cd7).^2))/length(cd7)));
            energy5(cnt_, 8) = log2((sum((abs(cd8).^2))/length(cd8)));            
            energy5(cnt_, 9) = log2((sum((abs(cd9).^2))/length(cd9)));
            energy5(cnt_, 10) = log2((sum((abs(cd10).^2))/length(cd10))); % 10th level. Match cd10 with the c(11:20) they should match.
            %energy5(cnt_, 11) = log2((sum((abs(cd10).^2))/length(cd12)));
            %energy5(cnt_, 12) = log2((sum((abs(cd10).^2))/length(cd12)));
            reconsEmbedding_db5 = waverec(c(1:sum(l(1:10))), l(1:10), 'db5');     
            clear c l approx cd1 cd2 cd3 cd4 cd5 cd6 cd7 cd8 cd9 cd10 cd11 cd12

            [dh2(cnt_,:),h2(cnt_,:),cp2(cnt_,:),tauq2(cnt_,:),leaders2,structfunc2,Localtmp_2, lenLocaltmp2, zkqq2] = dwtleader(reconsEmbedding_db2, 'db2');                           
            [dh3(cnt_,:),h3(cnt_,:),cp3(cnt_,:),tauq3(cnt_,:),leaders3,structfunc3,Localtmp_3, lenLocaltmp3, zkqq3] = dwtleader(reconsEmbedding_db3, 'db3');            
            [dh4(cnt_,:),h4(cnt_,:),cp4(cnt_,:),tauq4(cnt_,:),leaders4,structfunc4,Localtmp_4, lenLocaltmp4, zkqq4] = dwtleader(reconsEmbedding_db4, 'db4');
            [dh5(cnt_,:),h5(cnt_,:),cp5(cnt_,:),tauq5(cnt_,:),leaders5,structfunc5,Localtmp_5, lenLocaltmp5, zkqq5] = dwtleader(reconsEmbedding_db5, 'db5');            

            newL2 = [1, cumsum(lenLocaltmp2)];
            for zl = 1 : length(newL2)-1
                logVal(:,zl, cnt_) = log2(mean((zkqq2(:, newL2(zl):newL2(zl+1))),2));
            end    
    
            capacity = round(median(h2(cnt_,6)),2);
            if capacity > 1
                capacity = 1;
            end
            try
                FBMWithSameCapacity = wfbm(capacity, 1750, 8);
            catch
                FBMWithSameCapacity = wfbm(capacity+0.01, 1750, 8);
            end
            [dhFBMArr(cnt_,:),hFBMArr(cnt_,:),cpFBMArr(cnt_,:),tauqFBMArr(cnt_,:),leadersFBM, structfuncFBM,Localtmp_fbm, lenLocaltmpfbm, zkqqfbm] = dwtleader(FBMWithSameCapacity, 'db2');
            newLfbm = [1, cumsum(lenLocaltmpfbm)];
            for zl = 1 : length(newLfbm)-1
                logValfbm(:,zl, cnt_) = log2(mean((zkqqfbm(:, newLfbm(zl):newLfbm(zl+1))),2));
            end    
            
            alphaLocalSingularity2 = [alphaLocalSingularity2, Localtmp_2];
            alphaLocalSingularity3 = [alphaLocalSingularity3, Localtmp_3];
            alphaLocalSingularity4 = [alphaLocalSingularity4, Localtmp_4];
            alphaLocalSingularity5 = [alphaLocalSingularity5, Localtmp_5];
            clear Localtmp_2 Localtmp_3 Localtmp_4 Localtmp_5 
            if cnt_ > 1
                if length(lenLocaltmp2) < size(lenLocalSing2,2)
                    lenLocalSing2(cnt_, :) = [lenLocaltmp2, zeros(1,size(lenLocalSing2,2)-length(lenLocaltmp2))];
                else
                    lenLocalSing2(cnt_, :) = lenLocaltmp2;
                end
                if length(lenLocaltmp3) < size(lenLocalSing3,2)
                    lenLocalSing3(cnt_, :) = [lenLocaltmp3, zeros(1,size(lenLocalSing3,2)-length(lenLocaltmp3))];
                else
                    lenLocalSing3(cnt_, :) = lenLocaltmp3;
                end
                if length(lenLocaltmp4) < size(lenLocalSing4,2)
                    lenLocalSing4(cnt_, :) = [lenLocaltmp4, zeros(1,size(lenLocalSing4,2)-length(lenLocaltmp4))];
                else
                    lenLocalSing4(cnt_, :) = lenLocaltmp4;
                end
                if length(lenLocaltmp5) < size(lenLocalSing5,2)
                    lenLocalSing5(cnt_, :) = [lenLocaltmp5, zeros(1,size(lenLocalSing5,2)-length(lenLocaltmp5))];
                else
                    lenLocalSing5(cnt_, :) = lenLocaltmp5;
                end
            else
                lenLocalSing2(cnt_, :) = lenLocaltmp2;
                lenLocalSing3(cnt_, :) = lenLocaltmp3;
                lenLocalSing4(cnt_, :) = lenLocaltmp4;
                lenLocalSing5(cnt_, :) = lenLocaltmp5;
            end
            clear lenLocaltmp2 lenLocaltmp3 lenLocaltmp4 lenLocaltmp5
            
            min2(cnt_) = min(h2(cnt_, :));
            min3(cnt_) = min(h3(cnt_, :));
            min4(cnt_) = min(h4(cnt_, :));
            min5(cnt_) = min(h5(cnt_, :));

            max2(cnt_) = max(h2(cnt_, :));
            max3(cnt_) = max(h3(cnt_, :));
            max4(cnt_) = max(h4(cnt_, :));
            max5(cnt_) = max(h5(cnt_, :));
            
            width2(cnt_) = max(h2(cnt_, :))-min(h2(cnt_, :));
            width3(cnt_) = max(h3(cnt_, :))-min(h3(cnt_, :));
            width4(cnt_) = max(h4(cnt_, :))-min(h4(cnt_, :));
            width5(cnt_) = max(h5(cnt_, :))-min(h5(cnt_, :));
            
            %subplot(2,1,1);plot(embedding);
            %subplot(2,1,2);plot(embedding2);
            cnt_ = cnt_ + 1;
            close all
                
        end
    end

    %m_dh = median(dh2);
    %m_h = median(h2);
    %trueIdx = m_h < 1;
    %plot(m_h(trueIdx), m_dh(trueIdx))
    
    save(fullfile(finalPath, strcat(prefix, 'filesProcessed.mat')), 'filesProcessed');
    newlogVal = median(logVal,3);
    newlogVal = logVal(:,:,1);
    %plot([6:-1:1],newlogVal(1,:))
    %hold on ; plot([6:-1:1],newlogVal(2,:))
    %hold on ; plot([6:-1:1],newlogVal(3,:))
    %hold on ; plot([6:-1:1],newlogVal(4,:))
    %hold on ; plot([6:-1:1],newlogVal(5,:))
    %hold on ; plot([6:-1:1],newlogVal(6,:))
    %hold on ; plot([6:-1:1],newlogVal(7,:))
    %hold on ; plot([6:-1:1],newlogVal(8,:))
    %hold on ; plot([6:-1:1],newlogVal(9,:))
    %hold on ; plot([6:-1:1],newlogVal(10,:))
    %hold on ; plot([6:-1:1],newlogVal(11,:)) 
    %saveas(gca, fullfile(finalPath, strcat(prefix, 'structureFunctionSCale ', subName_)), 'png')
        
    newlogVal = median(logValfbm,3);
    plot([size(newlogVal,2):-1:1],newlogVal(1,:))
    hold on ; plot([size(newlogVal,2):-1:1],newlogVal(2,:))
    hold on ; plot([size(newlogVal,2):-1:1],newlogVal(3,:))
    hold on ; plot([size(newlogVal,2):-1:1],newlogVal(4,:))
    hold on ; plot([size(newlogVal,2):-1:1],newlogVal(5,:))
    hold on ; plot([size(newlogVal,2):-1:1],newlogVal(6,:))
    hold on ; plot([size(newlogVal,2):-1:1],newlogVal(7,:))
    hold on ; plot([size(newlogVal,2):-1:1],newlogVal(8,:))
    hold on ; plot([size(newlogVal,2):-1:1],newlogVal(9,:))
    hold on ; plot([size(newlogVal,2):-1:1],newlogVal(10,:))
    hold on ; plot([size(newlogVal,2):-1:1],newlogVal(11,:)) 
    saveas(gca, fullfile(finalPath, strcat(prefix, 'structureFunctionSCaleFBM ', subName_)), 'png')
    
    save(fullfile(finalPath, strcat(prefix, 'db2_LocalSingularity.mat')), 'alphaLocalSingularity2');
    save(fullfile(finalPath, strcat(prefix, 'db3_LocalSingularity.mat')), 'alphaLocalSingularity3');
    save(fullfile(finalPath, strcat(prefix, 'db4_LocalSingularity.mat')), 'alphaLocalSingularity4');
    save(fullfile(finalPath, strcat(prefix, 'db5_LocalSingularity.mat')), 'alphaLocalSingularity5');
        
    save(fullfile(finalPath, strcat(prefix, 'db2_LengthLocalSingularity.mat')), 'lenLocalSing2');
    save(fullfile(finalPath, strcat(prefix, 'db3_LengthLocalSingularity.mat')), 'lenLocalSing3');
    save(fullfile(finalPath, strcat(prefix, 'db4_LengthLocalSingularity.mat')), 'lenLocalSing4');
    save(fullfile(finalPath, strcat(prefix, 'db5_LengthLocalSingularity.mat')), 'lenLocalSing5');

    alphaLocalSingularity2(find(alphaLocalSingularity2>3))=[];
    alphaLocalSingularity3(find(alphaLocalSingularity3>3))=[];
    alphaLocalSingularity4(find(alphaLocalSingularity4>3))=[];
    alphaLocalSingularity5(find(alphaLocalSingularity5>3))=[];

    %%%%%%%%%%%%%%%% Spread of Local Singularities
    figure('units','normalized','outerposition',[0 0 1 1]);
    suptitle('PointWise Singularity To Check The Spread of Singularity');
    subplot(2,2,1)
    statParam = strcat('Mean-',num2str(mean(alphaLocalSingularity2)), '\nMedian-',num2str(median(alphaLocalSingularity2)));
    hist(alphaLocalSingularity2, 35);xticks([0:3/35:3]);xtickangle(90);
    title('db2', 'FontSize', 10);
    grid on
    xlabel('Counts')
    ylabel('Pointwise Singularity')
    %legend(statParam, 'Location','northeast', 'FontSize', 11)        
    %----------
    subplot(2,2,2)
    statParam = strcat('Mean-',num2str(mean(alphaLocalSingularity3)), '\nMedian-',num2str(median(alphaLocalSingularity3)));
    hist(alphaLocalSingularity3, 35);xticks([0:3/35:3]);xtickangle(90);    
    title('db3', 'FontSize', 10);
    grid on
    xlabel('Counts')
    ylabel('Pointwise Singularity')
    %legend(statParam, 'Location','northeast', 'FontSize', 11)   
    %----------
    subplot(2,2,3)
    statParam = strcat('Mean-',num2str(mean(alphaLocalSingularity4)), '\nMedian-',num2str(median(alphaLocalSingularity4)));
    hist(alphaLocalSingularity4, 35);xticks([0:3/35:3]);xtickangle(90);
    title('db4', 'FontSize', 10);
    grid on
    xlabel('Counts')
    ylabel('Pointwise Singularity')
    %legend(statParam, 'Location','northeast', 'FontSize', 11)   
    %----------
    subplot(2,2,4)
    statParam = strcat('Mean-',num2str(mean(alphaLocalSingularity5)), '\nMedian-',num2str(median(alphaLocalSingularity5)));
    hist(alphaLocalSingularity5, 35);xticks([0:3/35:3]);xtickangle(90);
    title('db5', 'FontSize', 10);
    grid on
    xlabel('Counts')
    ylabel('Pointwise Singularity')    
    %legend(statParam, 'Location','northeast', 'FontSize', 11)   
    saveas(gca, fullfile(finalPath, strcat(prefix, 'Local_scaling_singularity ', subName_)), 'png')
    close all
    
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'subjectsNotIncluded.txt')),'w');
    if length(partNotIncludedArr) ~= 0
        fprintf(fileID,'%s \r\n', partNotIncludedArr{1, :});
    else
        fprintf(fileID,'%s \r\n', '');
    end
    fclose(fileID);
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'subjectsIncluded.txt')),'w');
    fprintf(fileID,'%s \r\n', subjectsIncluded{:, 1});
    fclose(fileID);

    save(fullfile(finalPath, strcat(prefix, 'db2_singularityDimFBMArr.mat')), 'dhFBMArr');
    save(fullfile(finalPath, strcat(prefix, 'db2_singularityExponentFBMArr.mat')), 'hFBMArr');
    save(fullfile(finalPath, strcat(prefix, 'db2_scalingExpoentFBMArr.mat')), 'tauqFBMArr');
    save(fullfile(finalPath, strcat(prefix, 'db2_cumulantsFBMArr.mat')), 'cpFBMArr');    
    
    save(fullfile(finalPath, strcat(prefix, 'db2_singularityDim.mat')), 'dh2');
    save(fullfile(finalPath, strcat(prefix, 'db2_singularityExponent.mat')), 'h2');
    save(fullfile(finalPath, strcat(prefix, 'db2_scalingExpoent.mat')), 'tauq2');
    save(fullfile(finalPath, strcat(prefix, 'db2_cumulants.mat')), 'cp2');
    [pval1,hpt1,stats1] = signtest(cp2(:, 1), 0.6, 0.01);
    [pval2,hpt2,stats2] = signtest(cp2(:, 2), 0, 0.01);
    [pval3,hpt3,stats3] = signtest(cp2(:, 3), 0, 0.01);    
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'db2_PValues.txt')),'w');
    fprintf(fileID,'%f \r\n', [pval1, pval2, pval3]);
    fclose(fileID);
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'db2_HValues.txt')),'w');
    fprintf(fileID,'%f \r\n', [hpt1, hpt2, hpt3]);
    fclose(fileID);
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'db2_Stats.txt')),'w');
    fprintf(fileID,'%f \r\n', [stats1.zval, stats2.zval, stats3.zval]);
    fclose(fileID);
    
    save(fullfile(finalPath, strcat(prefix, 'db3_singularityDim.mat')), 'dh3');
    save(fullfile(finalPath, strcat(prefix, 'db3_singularityExponent.mat')), 'h3');
    save(fullfile(finalPath, strcat(prefix, 'db3_scalingExpoent.mat')), 'tauq3');
    save(fullfile(finalPath, strcat(prefix, 'db3_cumulants.mat')), 'cp3');
    [pval1,hpt1,stats1] = signtest(cp3(:, 1), 0.6, 0.01);
    [pval2,hpt2,stats2] = signtest(cp3(:, 2), 0, 0.01);
    [pval3,hpt3,stats3] = signtest(cp3(:, 3), 0, 0.01);    
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'db3_PValues.txt')),'w');
    fprintf(fileID,'%f \r\n', [pval1, pval2, pval3]);
    fclose(fileID);
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'db3_HValues.txt')),'w');
    fprintf(fileID,'%f \r\n', [hpt1, hpt2, hpt3]);
    fclose(fileID);
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'db3_Stats.txt')),'w');
    fprintf(fileID,'%f \r\n', [stats1.zval, stats2.zval, stats3.zval]);
    fclose(fileID);
    
    save(fullfile(finalPath, strcat(prefix, 'db4_singularityDim.mat')), 'dh4');
    save(fullfile(finalPath, strcat(prefix, 'db4_singularityExponent.mat')), 'h4');
    save(fullfile(finalPath, strcat(prefix, 'db4_scalingExpoent.mat')), 'tauq4');
    save(fullfile(finalPath, strcat(prefix, 'db4_cumulants.mat')), 'cp4');
    [pval1,hpt1,stats1] = signtest(cp4(:, 1), 0.6, 0.01);
    [pval2,hpt2,stats2] = signtest(cp4(:, 2), 0, 0.01);
    [pval3,hpt3,stats3] = signtest(cp4(:, 3), 0, 0.01);    
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'db4_PValues.txt')),'w');
    fprintf(fileID,'%f \r\n', [pval1, pval2, pval3]);
    fclose(fileID);
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'db4_HValues.txt')),'w');
    fprintf(fileID,'%f \r\n', [hpt1, hpt2, hpt3]);
    fclose(fileID);
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'db4_Stats.txt')),'w');
    fprintf(fileID,'%f \r\n', [stats1.zval, stats2.zval, stats3.zval]);
    fclose(fileID);
    
    save(fullfile(finalPath, strcat(prefix, 'db5_singularityDim.mat')), 'dh5');
    save(fullfile(finalPath, strcat(prefix, 'db5_singularityExponent.mat')), 'h5');
    save(fullfile(finalPath, strcat(prefix, 'db5_scalingExpoent.mat')), 'tauq5');
    save(fullfile(finalPath, strcat(prefix, 'db5_cumulants.mat')), 'cp5');        
    [pval1,hpt1,stats1] = signtest(cp5(:, 1), 0.6, 0.01);
    [pval2,hpt2,stats2] = signtest(cp5(:, 2), 0, 0.01);
    [pval3,hpt3,stats3] = signtest(cp5(:, 3), 0, 0.01);    
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'db5_PValues.txt')),'w');
    fprintf(fileID,'%f \r\n', [pval1, pval2, pval3]);
    fclose(fileID);
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'db5_HValues.txt')),'w');
    fprintf(fileID,'%f \r\n', [hpt1, hpt2, hpt3]);
    fclose(fileID);
    fileID = fopen(fullfile(finalPath, strcat(prefix, 'db5_Stats.txt')),'w');
    fprintf(fileID,'%f \r\n', [stats1.zval, stats2.zval, stats3.zval]);
    fclose(fileID);
    
    %return
    %% db2, db3, db4 and db5 plots.
    %db2    
    figure('units','normalized','outerposition',[0 0 1 1]);
    suptitle('db2');
    subplot(4,2,1)
    plot([1:11],mean(tauq2),'bo--')
    xticks([1 2 3 4 5 6 7 8 9 10 11])
    xticklabels([-5 -4 -3 -2 -1 0 1 2 3 4 5])
    title('Scaling Exponents', 'FontSize', 10);
    grid on
    xlabel('qth Moments')
    ylabel('\tau(q)')
    
    subplot(4,2,2)
    width = strcat('Mean Width-', num2str(mean(width2)));
    plot(mean(h2),mean(dh2),'bo--', 'DisplayName', width)    
    xlabel('H')
    ylabel('D(H)')
    title({'Singularity Spectrum'; ['cp-1:' num2str(mean(cp2(:,1))) strcat('(',num2str(std(cp2(:,1))),')')]}, 'FontSize', 10);
    legend('Location','southwest', 'FontSize', 11)
    grid on
    legend('Location','southwest', 'FontSize', 11)

    numberSubConsd = length(subjectsNotNoisy) - partNotIncluded;
    plotCnt = 3;
    for grp_ = 1 : size(classAssGrp,1)
        lower_ = ((grp_-1)*numberSubConsd)+1;
        upper_ = grp_*numberSubConsd;
        subplot(4,2,plotCnt)
        plotCnt = plotCnt + 1;
        plot([1:11],mean(tauq2(lower_:upper_, :)),'bo--')
        xticks([1 2 3 4 5 6 7 8 9 10 11])
        xticklabels([-5 -4 -3 -2 -1 0 1 2 3 4 5])
        %title('Scaling Exponents', 'FontSize', 10);
        grid on
        xlabel('qth Moments')
        ylabel(strcat('C-', num2str(grp_),'\tau(q)'))

        subplot(4,2,plotCnt)
        plotCnt = plotCnt + 1;
        width = strcat('Mean Width-', num2str(mean(width2(lower_:upper_))));
        plot(mean(h2(lower_:upper_, :)),mean(dh2(lower_:upper_, :)),'bo--', 'DisplayName', width)    
        xlabel('H')
        ylabel(strcat('C-', num2str(grp_),'D(H)'))
        title({['cp-1:' num2str(mean(cp2(lower_:upper_,1))) strcat('(',num2str(std(cp2(lower_:upper_,1))),')')]}, 'FontSize', 10);
        legend('Location','southwest', 'FontSize', 11)
        grid on
        legend('Location','southwest', 'FontSize', 11)        
    end
    saveas(gca, fullfile(finalPath, strcat(prefix, 'db2_scaling_singularity ', subName_)), 'jpeg')
    close all
    
    %db3
    figure('units','normalized','outerposition',[0 0 1 1]);
    suptitle('db3');
    subplot(4,2,1)
    plot([1:11],mean(tauq3),'bo--')
    xticks([1 2 3 4 5 6 7 8 9 10 11])
    xticklabels([-5 -4 -3 -2 -1 0 1 2 3 4 5])
    title('Scaling Exponents', 'FontSize', 10);
    grid on
    xlabel('qth Moments')
    ylabel('\tau(q)')
    
    subplot(4,2,2)
    width = strcat('Mean Width-', num2str(mean(width3)));
    plot(mean(h3),mean(dh3),'bo--', 'DisplayName', width)    
    xlabel('H')
    ylabel('D(H)')
    title({'Singularity Spectrum'; ['cp-1:' num2str(mean(cp3(:,1))) strcat('(',num2str(std(cp3(:,1))),')')]}, 'FontSize', 10);
    legend('Location','southwest', 'FontSize', 11)
    grid on
    legend('Location','southwest', 'FontSize', 11)
    
    numberSubConsd = length(subjectsNotNoisy) - partNotIncluded;
    plotCnt = 3;
    for grp_ = 1 : size(classAssGrp,1)
        lower_ = ((grp_-1)*numberSubConsd)+1;
        upper_ = grp_*numberSubConsd;
        subplot(4,2,plotCnt)
        plotCnt = plotCnt + 1;
        plot([1:11],mean(tauq3(lower_:upper_, :)),'bo--')
        xticks([1 2 3 4 5 6 7 8 9 10 11])
        xticklabels([-5 -4 -3 -2 -1 0 1 2 3 4 5])
        %title('Scaling Exponents', 'FontSize', 10);
        grid on
        xlabel('qth Moments')
        ylabel(strcat('C-', num2str(grp_),'\tau(q)'))

        subplot(4,2,plotCnt)
        plotCnt = plotCnt + 1;
        width = strcat('Mean Width-', num2str(mean(width3(lower_:upper_))));
        plot(mean(h3(lower_:upper_, :)),mean(dh3(lower_:upper_, :)),'bo--', 'DisplayName', width)    
        xlabel('H')
        ylabel(strcat('C-', num2str(grp_),'D(H)'))
        title({['cp-1:' num2str(mean(cp3(lower_:upper_,1))) strcat('(',num2str(std(cp3(lower_:upper_,1))),')')]}, 'FontSize', 10);
        legend('Location','southwest', 'FontSize', 11)
        grid on
        legend('Location','southwest', 'FontSize', 11)        
    end
    saveas(gca, fullfile(finalPath, strcat(prefix, 'db3_scaling_singularity ', subName_)), 'jpeg')
    close all
    % =============Commented For the time Being above one
    
    % db4
    figure('units','normalized','outerposition',[0 0 1 1]);
    suptitle('db4');
    subplot(4,2,1)
    plot([1:11],mean(tauq4),'bo--')
    xticks([1 2 3 4 5 6 7 8 9 10 11])
    xticklabels([-5 -4 -3 -2 -1 0 1 2 3 4 5])
    title('Scaling Exponents', 'FontSize', 10);
    grid on
    xlabel('qth Moments')
    ylabel('\tau(q)')
    
    subplot(4,2,2)
    width = strcat('Mean Width-', num2str(mean(width4)));
    plot(mean(h4),mean(dh4),'bo--', 'DisplayName', width)    
    xlabel('H')
    ylabel('D(H)')
    title({'Singularity Spectrum'; ['cp-1:' num2str(mean(cp4(:,1))) strcat('(',num2str(std(cp4(:,1))),')')]}, 'FontSize', 10);
    legend('Location','southwest', 'FontSize', 11)
    grid on
    legend('Location','southwest', 'FontSize', 11)
    
    numberSubConsd = length(subjectsNotNoisy) - partNotIncluded;
    plotCnt = 3;
    for grp_ = 1 : size(classAssGrp,1)
        lower_ = ((grp_-1)*numberSubConsd)+1;
        upper_ = grp_*numberSubConsd;
        subplot(4,2,plotCnt)
        plotCnt = plotCnt + 1;
        plot([1:11],mean(tauq4(lower_:upper_, :)),'bo--')
        xticks([1 2 3 4 5 6 7 8 9 10 11])
        xticklabels([-5 -4 -3 -2 -1 0 1 2 3 4 5])
        %title('Scaling Exponents', 'FontSize', 10);
        grid on
        xlabel('qth Moments')
        ylabel(strcat('C-', num2str(grp_),'\tau(q)'))

        subplot(4,2,plotCnt)
        plotCnt = plotCnt + 1;
        width = strcat('Mean Width-', num2str(mean(width4(lower_:upper_))));
        plot(mean(h4(lower_:upper_, :)),mean(dh4(lower_:upper_, :)),'bo--', 'DisplayName', width)    
        xlabel('H')
        ylabel(strcat('C-', num2str(grp_),'D(H)'))
        title({['cp-1:' num2str(mean(cp4(lower_:upper_,1))) strcat('(',num2str(std(cp4(lower_:upper_,1))),')')]}, 'FontSize', 10);
        legend('Location','southwest', 'FontSize', 11)
        grid on
        legend('Location','southwest', 'FontSize', 11)        
    end
    saveas(gca, fullfile(finalPath, strcat(prefix, 'db4_scaling_singularity ', subName_)), 'jpeg')
    
    % db5
    figure('units','normalized','outerposition',[0 0 1 1]);
    suptitle('db5');
    subplot(4,2,1)
    plot([1:11],mean(tauq5),'bo--')
    xticks([1 2 3 4 5 6 7 8 9 10 11])
    xticklabels([-5 -4 -3 -2 -1 0 1 2 3 4 5])
    title('Scaling Exponents', 'FontSize', 10);
    grid on
    xlabel('qth Moments')
    ylabel('\tau(q)')
    
    subplot(4,2,2)
    width = strcat('Mean Width-', num2str(mean(width5)));
    plot(mean(h5),mean(dh5),'bo--', 'DisplayName', width)    
    xlabel('H')
    ylabel('D(H)')
    title({'Singularity Spectrum'; ['cp-1:' num2str(mean(cp5(:,1))) strcat('(',num2str(std(cp5(:,1))),')')]}, 'FontSize', 10);
    legend('Location','southwest', 'FontSize', 11)
    grid on
    legend('Location','southwest', 'FontSize', 11)
    
    numberSubConsd = length(subjectsNotNoisy) - partNotIncluded;
    plotCnt = 3;
    for grp_ = 1 : size(classAssGrp,1)
        lower_ = ((grp_-1)*numberSubConsd)+1;
        upper_ = grp_*numberSubConsd;
        subplot(4,2,plotCnt)
        plotCnt = plotCnt + 1;
        plot([1:11],mean(tauq5(lower_:upper_, :)),'bo--')
        xticks([1 2 3 4 5 6 7 8 9 10 11])
        xticklabels([-5 -4 -3 -2 -1 0 1 2 3 4 5])
        %title('Scaling Exponents', 'FontSize', 10);
        grid on
        xlabel('qth Moments')
        ylabel(strcat('C-', num2str(grp_),'\tau(q)'))

        subplot(4,2,plotCnt)
        plotCnt = plotCnt + 1;
        width = strcat('Mean Width-', num2str(mean(width5(lower_:upper_))));
        plot(mean(h5(lower_:upper_, :)),mean(dh5(lower_:upper_, :)),'bo--', 'DisplayName', width)    
        xlabel('H')
        ylabel(strcat('C-', num2str(grp_),'D(H)'))
        title({['cp-1:' num2str(mean(cp5(lower_:upper_,1))) strcat('(',num2str(std(cp5(lower_:upper_,1))),')')]}, 'FontSize', 10);
        legend('Location','southwest', 'FontSize', 11)
        grid on
        legend('Location','southwest', 'FontSize', 11)        
    end
    saveas(gca, fullfile(finalPath, strcat(prefix, 'db5_scaling_singularity ', subName_)), 'jpeg')
    
    %return
    %{
    for chanIdx = 1 : ALLEEG.nbchan
        [pchanFirst(chanIdx), hChanFirst(chanIdx), statsChanFirst(chanIdx)] = signtest(firstPart(chanIdx, :), secondPart(chanIdx, :), 'tail', 'right','alpha',0.01);
        [pchanSecond(chanIdx), hChanSecond(chanIdx), statsChanSecond(chanIdx)] = signtest(secondPart(chanIdx, :), firstPart(chanIdx, :), 'tail', 'right','alpha',0.01);        
        if hChanFirst(chanIdx) == 0
            segWiseStatePlotFirst(chanIdx) = 0;
        else
            if (median(firstPart(chanIdx, :))>=median(secondPart(chanIdx, :)))
                segWiseStatePlotFirst(chanIdx) = median(firstPart(chanIdx, :));
            else
                if EQFlag == 0
                    error('Error in Calculation')
                end
            end
        end
        if hChanSecond(chanIdx) == 0
            segWiseStatePlotSecond(chanIdx) = 0;
        else
            if (median(firstPart(chanIdx, :))<=median(secondPart(chanIdx, :)))
                segWiseStatePlotSecond(chanIdx) = median(secondPart(chanIdx, :));
            else
                if EQFlag == 0
                    error('Error in Calculation')
                end
            end            
        end        
    end    
    save(fullfile(finalPath, strcat(prefix, 'Seg(1-4) ', subName_)), 'segWiseStatePlotFirst')
    save(fullfile(finalPath, strcat(prefix, 'Seg(5-9) ', subName_)), 'segWiseStatePlotSecond')
    %}
    
    %fid = fopen(FilePath,'w'); for rows = 1:size(Cell_in,1) fprintf(fid,'%s\n',Cell_in{rows,:}); end. fclose(fid)
    
    totalMT = totalMT/cnt_;
    close all
    %% 
    % =============Commented For the time Being
    
    % replace 1750 with len(embedding2)
    %FBMWithSameCapacity = wfbm(0.25, 1750, 8, 'db2');
    %FBM05 = wfbm(0.5, 1750, 8, 'db2');
    %FBM075 = wfbm(0.75, 1750, 8, 'db2');
    %FBM1 = wfbm(1.0, 1750, 8, 'db2');
    %WGN = wgn(1750,1,0);
    disp(median(h2(:,6)))
    %FBMWithSameCapacity = wfbm(1, 1750, 8);  % For Baseline condition in Alpha Frequency Band
    try
        FBMWithSameCapacity = wfbm(round(median(h2(:,6)),2), 1750, 8);
    catch
        FBMWithSameCapacity = wfbm(round(median(h2(:,6)),2)+0.01, 1750, 8);
    end
    
    FBM05 = wfbm(0.5, 1750, 8);
    FBM075 = wfbm(0.75, 1750, 8);
    FBM1 = wfbm(1.0, 1750, 8);
    WGN = wgn(1750,1,0);    
    %FBM15 = wfbm(1.5, 1750, 8, 'db4');
    %FBM20 = wfbm(2.0, 1750, 8, 'db4');
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    %suptitle(strcat('Mean Time is - ', num2str(round(totalMT, 2))));
    subplot(6,1,1)
    plot(FBMWithSameCapacity); title(strcat('FBM-', num2str(median(h2(:,6)))), 'FontSize', 14);
    subplot(6,1,2)
    plot(FBM05); title('FBM05', 'FontSize', 14);
    subplot(6,1,3)
    plot(FBM075); title('FBM075', 'FontSize', 14);
    subplot(6,1,4)
    plot(FBM1); title('FBM1', 'FontSize', 14);
    subplot(6,1,5)
    plot(WGN); title('WGN', 'FontSize', 14);
    subplot(6,1,6)
    plot(embedding2); title('Embedding', 'FontSize', 14);    
    saveas(gca, fullfile(finalPath, strcat(prefix, 'FBM WGN Plot ', subName_)), 'jpeg')
    %close all
    
    [dhFBMWithSameCapacity,hFBMWithSameCapacity,cpFBMWithSameCapacity,tauqFBMWithSameCapacity,leadersFBMWithSameCapacity, structfuncFBMWithSameCapacity] = dwtleader(FBMWithSameCapacity, 'db2');
    [dhfbm05,hfbm05,cpfbm05,tauqfbm05,leadersfbm05, structfuncfbm05] = dwtleader(FBM05, 'db2');
    [dhfbm075,hfbm075,cpfbm075,tauqfbm075,leadersfbm075, structfuncfbm075] = dwtleader(FBM075, 'db2');
    [dhfbm1,hfbm1,cpfbm1,tauqfbm1,leadersfbm1, structfuncfbm1] = dwtleader(FBM1, 'db2');
    [dhwgn,hwgn,cpwgn,tauqwgn,leaderswgn, structfuncwgn] = dwtleader(WGN, 'db2');
    save(fullfile(finalPath, strcat(prefix, 'db2_singularityDimFBM.mat')), 'dhFBMWithSameCapacity');
    save(fullfile(finalPath, strcat(prefix, 'db2_singularityExponentFBM.mat')), 'hFBMWithSameCapacity');
    save(fullfile(finalPath, strcat(prefix, 'db2_scalingExpoentFBM.mat')), 'tauqFBMWithSameCapacity');
    save(fullfile(finalPath, strcat(prefix, 'db2_cumulantsFBM.mat')), 'cpFBMWithSameCapacity');    
    %[dhfbm15,hfbm15,cpfbm15,tauqfbm15,leadersfbm15, structfuncfbm15] = dwtleader(FBM15, 'db4');
    %[dhfbm20,hfbm20,cpfbm20,tauqfbm20,leadersfbm20, structfuncfbm20] = dwtleader(FBM20, 'db4');

    %%
    % Scaling Exponents With Fractional Brownian Motion and White Gaussian Noise

    figure('units','normalized','outerposition',[0 0 1 1]);
    %suptitle(strcat('Mean Time is - ', num2str(round(totalMT, 2))));
    subplot(2,2,1)
    plot([1:11],tauqFBMWithSameCapacity,'bo--', 'DisplayName', strcat('FBM-', num2str(median(h2(:,6)))))
    hold on
    plot([1:11],tauqwgn, 'DisplayName', 'WGN')
    legend('Location','southeast', 'FontSize', 11);
    title(strcat('Scaling Exponents hurst=', num2str(median(h2(:,6)))), 'FontSize', 14);
    grid on
    xlabel('qth Moments')
    ylabel('\tau(q)')
            
    subplot(2,2,2)
    plot([1:11],tauqfbm05,'bo--', 'DisplayName', 'FBM05')
    hold on
    plot([1:11],tauqwgn, 'DisplayName', 'WGN')
    legend('Location','southeast', 'FontSize', 11);
    title('Scaling Exponents for fBM=0.5', 'FontSize', 14);
    grid on
    xlabel('qth Moments')
    ylabel('\tau(q)')

    subplot(2,2,3)
    plot([1:11],tauqfbm075,'bo--', 'DisplayName', 'FBM075')
    hold on
    plot([1:11],tauqwgn, 'DisplayName', 'WGN')
    legend('Location','southeast', 'FontSize', 11);
    title('Scaling Exponents for fBM=0.75', 'FontSize', 14);
    grid on
    xlabel('qth Moments')
    ylabel('\tau(q)')

    subplot(2,2,4)
    plot([1:11],tauqfbm1,'bo--', 'DisplayName', 'FBM1')
    hold on
    plot([1:11],tauqwgn, 'DisplayName', 'WGN')
    legend('Location','southeast', 'FontSize', 11);
    title('Scaling Exponents for fBM=1.0', 'FontSize', 14);
    grid on
    xlabel('qth Moments')
    ylabel('\tau(q)')
    
    saveas(gca, fullfile(finalPath, strcat(prefix, 'Scaling Exponent FBM ', subName_)), 'jpeg')

    %%
    % Singularity Spectrum With Fractional Brownian Motion and White Gaussian Noise

    figure('units','normalized','outerposition',[0 0 1 1]);
    %suptitle(strcat('Mean Time is - ', num2str(round(totalMT, 2))));
    subplot(2,2,1)
    plot((hFBMWithSameCapacity),(dhFBMWithSameCapacity),'bo--', 'DisplayName', strcat('FBM-',num2str(median(h2(:,6)))))
    hold on
    plot(hwgn, dhwgn, 'DisplayName', 'WGN')
    legend('Location','southeast', 'FontSize', 11);
    xlabel('h')
    ylabel('D(h): Singularity Spectrum')
    title({'db2 Singularity Spectrum'; ['First Cumulant ' num2str(cpFBMWithSameCapacity(1))]}, 'FontSize', 14);
    grid on

    subplot(2,2,2)
    plot((hfbm05),(dhfbm05),'bo--', 'DisplayName', 'FBM0.5')
    hold on
    plot(hwgn, dhwgn, 'DisplayName', 'WGN')
    legend('Location','southeast', 'FontSize', 11);
    xlabel('h')
    ylabel('D(h): Singularity Spectrum')
    title({'db2 Singularity Spectrum'; ['First Cumulant ' num2str(cpfbm05(1))]}, 'FontSize', 14);
    grid on

    subplot(2,2,3)
    plot((hfbm075),(dhfbm075),'bo--', 'DisplayName', 'FBM0.75')
    hold on
    plot(hwgn, dhwgn, 'DisplayName', 'WGN')
    legend('Location','southeast', 'FontSize', 11);
    xlabel('h')
    ylabel('D(h): Singularity Spectrum')
    title({'db2 Singularity Spectrum'; ['First Cumulant ' num2str(cpfbm075(1))]}, 'FontSize', 14);
    grid on

    subplot(2,2,4)
    plot((hfbm1),(hfbm1),'bo--', 'DisplayName', 'FBM1')
    hold on
    plot(hwgn, dhwgn, 'DisplayName', 'WGN')
    legend('Location','southeast', 'FontSize', 11);
    xlabel('h')
    ylabel('D(h): Singularity Spectrum')
    title({'db2 Singularity Spectrum'; ['First Cumulant ' num2str(cpfbm1(1))]}, 'FontSize', 14);
    grid on
    
    saveas(gca, fullfile(finalPath, strcat(prefix, 'Singularity Spectrum FBM ', subName_)), 'jpeg')
    close all
    %%
    % Energy Plot With Fractional Brownian Motion and White Gaussian Noise
    
    [c, l] = wavedec(WGN, noLevels, 'db4');
    [cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10] = detcoef(c,l,[1 : noLevels]);
    energywgn(1) = log2((sum((abs(cd1).^2))/length(cd1)));
    energywgn(2) = log2((sum((abs(cd2).^2))/length(cd2)));
    energywgn(3) = log2((sum((abs(cd3).^2))/length(cd3)));
    energywgn(4) = log2((sum((abs(cd4).^2))/length(cd4)));
    energywgn(5) = log2((sum((abs(cd5).^2))/length(cd5)));
    energywgn(6) = log2((sum((abs(cd6).^2))/length(cd6)));
    energywgn(7) = log2((sum((abs(cd7).^2))/length(cd7)));
    energywgn(8) = log2((sum((abs(cd8).^2))/length(cd8)));            
    energywgn(9) = log2((sum((abs(cd9).^2))/length(cd9)));
    energywgn(10) = log2((sum((abs(cd10).^2))/length(cd10)));

    [c, l] = wavedec(FBMWithSameCapacity, noLevels, 'db4');
    [cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10] = detcoef(c,l,[1 : noLevels]);
    FBMWithSameCapacity(1) = log2((sum((abs(cd1).^2))/length(cd1)));
    FBMWithSameCapacity(2) = log2((sum((abs(cd2).^2))/length(cd2)));
    FBMWithSameCapacity(3) = log2((sum((abs(cd3).^2))/length(cd3)));
    FBMWithSameCapacity(4) = log2((sum((abs(cd4).^2))/length(cd4)));
    FBMWithSameCapacity(5) = log2((sum((abs(cd5).^2))/length(cd5)));
    FBMWithSameCapacity(6) = log2((sum((abs(cd6).^2))/length(cd6)));
    FBMWithSameCapacity(7) = log2((sum((abs(cd7).^2))/length(cd7)));
    FBMWithSameCapacity(8) = log2((sum((abs(cd8).^2))/length(cd8)));            
    FBMWithSameCapacity(9) = log2((sum((abs(cd9).^2))/length(cd9)));
    FBMWithSameCapacity(10) = log2((sum((abs(cd10).^2))/length(cd10)));

    [c, l] = wavedec(FBM05, noLevels, 'db4');
    [cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10] = detcoef(c,l,[1 : noLevels]);
    fbm05(1) = log2((sum((abs(cd1).^2))/length(cd1)));
    fbm05(2) = log2((sum((abs(cd2).^2))/length(cd2)));
    fbm05(3) = log2((sum((abs(cd3).^2))/length(cd3)));
    fbm05(4) = log2((sum((abs(cd4).^2))/length(cd4)));
    fbm05(5) = log2((sum((abs(cd5).^2))/length(cd5)));
    fbm05(6) = log2((sum((abs(cd6).^2))/length(cd6)));
    fbm05(7) = log2((sum((abs(cd7).^2))/length(cd7)));
    fbm05(8) = log2((sum((abs(cd8).^2))/length(cd8)));            
    fbm05(9) = log2((sum((abs(cd9).^2))/length(cd9)));
    fbm05(10) = log2((sum((abs(cd10).^2))/length(cd10)));

    [c, l] = wavedec(FBM075, noLevels, 'db4');
    [cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10] = detcoef(c,l,[1 : noLevels]);
    fbm075(1) = log2((sum((abs(cd1).^2))/length(cd1)));
    fbm075(2) = log2((sum((abs(cd2).^2))/length(cd2)));
    fbm075(3) = log2((sum((abs(cd3).^2))/length(cd3)));
    fbm075(4) = log2((sum((abs(cd4).^2))/length(cd4)));
    fbm075(5) = log2((sum((abs(cd5).^2))/length(cd5)));
    fbm075(6) = log2((sum((abs(cd6).^2))/length(cd6)));
    fbm075(7) = log2((sum((abs(cd7).^2))/length(cd7)));
    fbm075(8) = log2((sum((abs(cd8).^2))/length(cd8)));            
    fbm075(9) = log2((sum((abs(cd9).^2))/length(cd9)));
    fbm075(10) = log2((sum((abs(cd10).^2))/length(cd10)));

    [c, l] = wavedec(FBM1, noLevels, 'db4');
    [cd1,cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10] = detcoef(c,l,[1 : noLevels]);
    fbm1(1) = log2((sum((abs(cd1).^2))/length(cd1)));
    fbm1(2) = log2((sum((abs(cd2).^2))/length(cd2)));
    fbm1(3) = log2((sum((abs(cd3).^2))/length(cd3)));
    fbm1(4) = log2((sum((abs(cd4).^2))/length(cd4)));
    fbm1(5) = log2((sum((abs(cd5).^2))/length(cd5)));
    fbm1(6) = log2((sum((abs(cd6).^2))/length(cd6)));
    fbm1(7) = log2((sum((abs(cd7).^2))/length(cd7)));
    fbm1(8) = log2((sum((abs(cd8).^2))/length(cd8)));            
    fbm1(9) = log2((sum((abs(cd9).^2))/length(cd9)));
    fbm1(10) = log2((sum((abs(cd10).^2))/length(cd10)));
       
    x1=linspace(1,noLevels, noLevels);
    [slopewgn, interceptwgn] = logfit(x1(1:noLevels), energywgn(1:noLevels), 'linear');
    yApproxNoise = slopewgn*x1+interceptwgn;

    %{
    figure('units','normalized','outerposition',[0 0 1 1]);
    suptitle(strcat('Mean Time is - ', num2str(round(totalMT, 2))));
    subplot(2,2,1)
    mnEnergy = FBMWithSameCapacity;
    [slope, intercept] = logfit(x1(1:noLevels), FBMWithSameCapacity(1:noLevels), 'linear');
    yApprox = slope*x1+intercept;
    plot([1:noLevels], FBMWithSameCapacity,'-o', 'DisplayName', 'Energy Spectrum'); 
    hold on
    plot([1:noLevels], yApprox,'-+', 'DisplayName',strcat('Fitted Slope:  ', num2str(slope)))
    hold on
    plot([1:noLevels], yApproxNoise,'-x', 'DisplayName',strcat('Noise Slope:  ', num2str(slopewgn)))
    xlabel('scale'); ylabel('log(energy)'); grid on;
    legend('Location','southeast', 'FontSize', 11);
    title('Hurst=0.25', 'FontSize', 14);

    subplot(2,2,2)
    mnEnergy = fbm05;
    [slope, intercept] = logfit(x1(1:noLevels), fbm05(1:noLevels), 'linear');
    yApprox = slope*x1+intercept;

    plot([1:noLevels], fbm05,'-o', 'DisplayName', 'Energy Spectrum'); 
    hold on
    plot([1:noLevels], yApprox,'-+', 'DisplayName',strcat('Fitted Slope:  ', num2str(slope)))
    hold on
    plot([1:noLevels], yApproxNoise,'-x', 'DisplayName',strcat('Noise Slope:  ', num2str(slopewgn)))
    xlabel('scale'); ylabel('log(energy)'); grid on;
    legend('Location','southeast', 'FontSize', 11);
    title('Hurst=0.5', 'FontSize', 14);

    subplot(2,2,3)
    mnEnergy = fbm075;
    [slope, intercept] = logfit(x1(1:noLevels), fbm075(1:noLevels), 'linear');
    yApprox = slope*x1+intercept;

    plot([1:noLevels], fbm075,'-o', 'DisplayName', 'Energy Spectrum'); 
    hold on
    plot([1:noLevels], yApprox,'-+', 'DisplayName',strcat('Fitted Slope:  ', num2str(slope)))
    hold on
    plot([1:noLevels], yApproxNoise,'-x', 'DisplayName',strcat('Noise Slope:  ', num2str(slopewgn)))
    xlabel('scale'); ylabel('log(energy)'); grid on;
    legend('Location','southeast', 'FontSize', 11);
    title('Hurst=0.75', 'FontSize', 14);

    subplot(2,2,4)
    mnEnergy = fbm1;
    [slope, intercept] = logfit(x1(1:noLevels), fbm1(1:noLevels), 'linear');
    yApprox = slope*x1+intercept;
    plot([1:noLevels], fbm1,'-o', 'DisplayName', 'Energy Spectrum'); 
    hold on
    plot([1:noLevels], yApprox,'-+', 'DisplayName',strcat('Fitted Slope:  ', num2str(slope)))
    hold on
    plot([1:noLevels], yApproxNoise,'-x', 'DisplayName',strcat('Noise Slope:  ', num2str(slopewgn)))
    xlabel('scale'); ylabel('log(energy)'); grid on;
    legend('Location','southeast', 'FontSize', 11);
    title('Hurst=1', 'FontSize', 14);
    
    saveas(gca, fullfile(finalPath, strcat(prefix, 'Simulation ', subName_)), 'jpeg')
    %}
    %%
    % Now My signal With Fractional Brownian Motion and White Gaussian Noise
    figure('units','normalized','outerposition',[0 0 1 1]);
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');  
    nexttile
    %suptitle(strcat('Mean Time is - ', num2str(round(totalMT, 2))));
    %subplot(2,2,1)
    x1=linspace(1,noLevels, noLevels);
    mnEnergy = mean(energy2);
    [slope, intercept] = logfit(x1(5:noLevels), mnEnergy(5:noLevels), 'linear');
    yApprox = slope*x1+intercept;
    plot([1:noLevels], median(energy2),'ro-', 'DisplayName', 'Energy Spectrum','LineWidth',3); 
    hold on
    plot([1:noLevels], yApprox,'b-+', 'DisplayName',strcat('Fitted Slope:  ', num2str(round(slope, 2))),'LineWidth',3); 
    hold on
    plot([1:noLevels], fbm075,'k-x', 'DisplayName', 'Energy Spectrum FBM','LineWidth',3);  
    hold on
    plot([1:noLevels], yApproxNoise,'m-s', 'DisplayName',strcat('WGN Noise Slope:  ', num2str(round(slopewgn,2))),'LineWidth',3); 
    xlabel('scale'); ylabel('log(energy)'); grid on;
    legend('Location','northwest', 'FontSize', 30);
    ax = gca;
    ax.FontSize = 30;   
    xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    xticklabels({'8ms','16','32','64','128','256','512','1s','2s','4s'})
    %title('db2', 'FontSize', 14);

    %{
    subplot(2,2,2)
    mnEnergy = mean(energy3);
    [slope, intercept] = logfit(x1(5:noLevels), mnEnergy(5:noLevels), 'linear');
    yApprox = slope*x1+intercept;
    plot([1:noLevels], mean(energy3),'ro-', 'DisplayName', 'Energy Spectrum'); 
    hold on
    plot([1:noLevels], yApprox,'b-+', 'DisplayName',strcat('Fitted Slope:  ', num2str(slope)))
    hold on
    plot([1:noLevels], fbm075,'g-x', 'DisplayName', 'Energy Spectrum FBM'); 
    hold on
    plot([1:noLevels], yApproxNoise,'c-s', 'DisplayName',strcat('WGN Noise Slope:  ', num2str(slopewgn)))
    xlabel('scale'); ylabel('log(energy)'); grid on;
    legend('Location','southeast', 'FontSize', 11);
    title('db3', 'FontSize', 14);

    subplot(2,2,3)
    mnEnergy = mean(energy4);
    [slope, intercept] = logfit(x1(5:noLevels), mnEnergy(5:noLevels), 'linear');
    yApprox = slope*x1+intercept;
    plot([1:noLevels], mean(energy4),'ro-', 'DisplayName', 'Energy Spectrum'); 
    hold on
    plot([1:noLevels], yApprox,'b-+', 'DisplayName',strcat('Fitted Slope:  ', num2str(slope)))
    hold on
    plot([1:noLevels], fbm075,'g-x', 'DisplayName', 'Energy Spectrum FBM'); 
    hold on
    plot([1:noLevels], yApproxNoise,'c-s', 'DisplayName',strcat('WGN Noise Slope:  ', num2str(slopewgn)))
    xlabel('scale'); ylabel('log(energy)'); grid on;
    legend('Location','southeast', 'FontSize', 11);
    title('db4', 'FontSize', 14);

    subplot(2,2,4)
    mnEnergy = mean(energy5);
    [slope, intercept] = logfit(x1(5:noLevels), mnEnergy(5:noLevels), 'linear');
    yApprox = slope*x1+intercept;
    plot([1:noLevels], mean(energy5),'ro-', 'DisplayName', 'Energy Spectrum'); 
    hold on
    plot([1:noLevels], yApprox,'b-+', 'DisplayName',strcat('Fitted Slope:  ', num2str(slope)))
    hold on
    plot([1:noLevels], fbm075,'g-x', 'DisplayName', 'Energy Spectrum FBM'); 
    %hold on
    %plot([1:noLevels], yApproxNoise,'c-s', 'DisplayName',strcat('WGN Noise Slope:  ', num2str(slopewgn)))
    %xlabel('scale'); ylabel('log(energy)'); grid on;
    %legend('Location','southeast', 'FontSize', 40);    
    %title('db5', 'FontSize', 14);
    %}
    
    saveas(gca, fullfile(finalPath, strcat(prefix, subName_, '_PowerLaw')), 'png')

    %close all

    %{
    ax1 = gca; % current axes
    ax1_pos = ax1.Position; % position of first axes
    %ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation','right', 'Color','none');
    hold on; plot([1 12], [0 mnEnergy(12)], '--');
    ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation','right', 'Color','none', 'YTick',[]);
    ax2.XTick = [0:1/6:1];
    xticklabels(ax2, {'','16ms','64ms','256ms','1s','4s','16s'})
    %}

    figure('units','normalized','outerposition',[0 0 1 1]);
    suptitle(strcat('Mean Time is - ', num2str(round(totalMT, 2))));
    subplot(2,2,1)
    plot([1:11],mean(tauq2),'bo--')
    title('Estimated Scaling Exponents db2', 'FontSize', 14);
    grid on
    xlabel('qth Moments')
    ylabel('\tau(q)')
            
    subplot(2,2,2)
    plot([1:11],mean(tauq3),'bo--')
    title('Estimated Scaling Exponents db3', 'FontSize', 14);
    grid on
    xlabel('qth Moments')
    ylabel('\tau(q)')

    subplot(2,2,3)
    plot([1:11],mean(tauq4),'bo--')
    title('Estimated Scaling Exponents db4', 'FontSize', 14);
    grid on
    xlabel('qth Moments')
    ylabel('\tau(q)')

    %subplot(2,2,4)
    %plot([1:11],mean(tauq5),'bo--')
    %title('Estimated Scaling Exponents db5', 'FontSize', 14);
    %grid on
    %xlabel('qth Moments')
    %ylabel('\tau(q)')    
    saveas(gca, fullfile(finalPath, strcat(prefix, 'Scaling Exponent ', subName_)), 'jpeg')
    close all

    figure('units','normalized','outerposition',[0 0 1 1]);
    suptitle(strcat('Mean Time is - ', num2str(round(totalMT, 2))));
    subplot(2,2,1)
    width = strcat('Mean Width-', num2str(mean(width2)));
    plot(mean(h2),mean(dh2),'bo--', 'DisplayName', width)    
    xlabel('h')
    ylabel('D(h): Singularity Spectrum')
    title({'db2 Singularity Spectrum'; ['First Cumulant ' num2str(cp2(1))]}, 'FontSize', 14);
    legend('Location','southwest', 'FontSize', 11)
    grid on

    subplot(2,2,2)
    width = strcat('Mean Width-', num2str(mean(width3)));
    plot(mean(h3),mean(dh3),'bo--', 'DisplayName', width)    
    xlabel('h')
    ylabel('D(h): Singularity Spectrum')
    title({'db3 Singularity Spectrum'; ['First Cumulant ' num2str(cp3(1))]}, 'FontSize', 14);
    legend('Location','southwest', 'FontSize', 11)
    grid on

    subplot(2,2,3)
    width = strcat('Mean Width-', num2str(mean(width4)));
    plot(mean(h4),mean(dh4),'bo--', 'DisplayName', width)    
    xlabel('h')
    ylabel('D(h): Singularity Spectrum')
    title({'db4 Singularity Spectrum'; ['First Cumulant ' num2str(cp4(1))]}, 'FontSize', 14);
    legend('Location','southwest', 'FontSize', 11)
    grid on

    %subplot(2,2,4)
    %plot(mean(h5),mean(dh5),'bo--')
    %xlabel('h')
    %ylabel('D(h): Singularity Spectrum')
    %title({'db5 Singularity Spectrum'; ['First Cumulant ' num2str(cp5(1))]}, 'FontSize', 14);
    %grid on    
    %xt=xticks(gca);
    %yt=yticks(gca);
    %width=max(xt)-min(xt);
    %dim = [.1 .1 .3 .3];
    %annotation('textbox', dim, strcat('width-',num2str(width)),'FitBoxToText','on')
    saveas(gca, fullfile(finalPath, strcat(prefix, 'Singularity Spectrum ', subName_)), 'jpeg')
    close all
    %}
    
    %%mean dh and std plots subject wise.
    disp(partNotIncludedArr)
    disp(partNotIncluded)
end