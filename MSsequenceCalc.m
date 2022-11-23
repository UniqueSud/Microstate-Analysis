function [] = MSsequenceCalc(targetPath, EEGFiles, EEG, subName_, EQFlag)
    if EQFlag == 0
        classAssGrp = [1,2,3,4;1,3,2,4;1,4,2,3];
        for micInd_ = 4 : 4
            fitbackSetPath = fullfile(targetPath, 'OnlySetFiles', strcat('With_',num2str(micInd_),'_MS'));
            if ~ isdir(fitbackSetPath)
                mkdir(fitbackSetPath)
            end              
            tmpPath = fullfile(targetPath, 'OnlyMat', strcat('With_',num2str(micInd_),'_MS'));
            if ~ isdir(tmpPath)
                mkdir(tmpPath)
            end   
            if micInd_==5
                disp('stop here')
            end
            for grp_ = 1 : size(classAssGrp,1)
                cpArr = [];
                labelisnotthere = 0;
                for i = 1 : length(EEGFiles)
                    tmpName = split(EEGFiles(i).name, '.set');
                    tmpName = tmpName{1};
                    if (~isfile(fullfile(tmpPath, strcat('MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')))) || (~isfile(fullfile(tmpPath, strcat(num2str(grp_), '_EmbeddingWithArraySum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')))) ||(~isfile(fullfile(tmpPath, strcat(num2str(grp_), '_EmbeddingWithScalarSum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat'))))
                        figure('units','normalized','outerposition',[0 0 1 1]);
                        EEG = pop_loadset('filename',strcat('Microstates_', num2str(micInd_), '_EEGMSA_', EEGFiles(i).name),'filepath', fitbackSetPath);  
                        MSsequence = EEG.microstate.fit.labels;
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
                        meanT = round(mean(countArr));
                        save(fullfile(tmpPath, strcat('MeanTime_MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')), 'meanT');
                        save(fullfile(tmpPath, strcat('MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')), 'MSsequence');
                
                        classAss = [];
                        classAss2 = 0;
                        embedding = [];
                        embedding2 = [];
                        for index_ = 1 : length(MSsequence)
                            if ~ any(classAssGrp(grp_, :)>=(MSsequence(index_)))
                                disp('********************************Problem MSsequenceIndex is greater.************')
                                break;                                
                            end
                            if (MSsequence(index_) == classAssGrp(grp_, 1)) || (MSsequence(index_) == classAssGrp(grp_, 2))
                                classAss2 = classAss2 + 1;
                                classAss = [classAss 1];
                            end

                            if (MSsequence(index_) == classAssGrp(grp_, 3)) || (MSsequence(index_) == classAssGrp(grp_, 4))
                                classAss2 = classAss2-1;
                                classAss = [classAss -1];
                            end        
                            embedding2 = [embedding2 classAss2];
                            embedding = [embedding sum(classAss)];
                        end
                        if index_ ~= length(MSsequence)
                            close all
                            continue                            
                        end
                        
                        EEG.microstate.fit.classAss = classAss;
                        EEG.microstate.fit.embedding2 = embedding2;
                        EEG.microstate.fit.embedding = embedding;

                        subplot(2,2,1);plot(MSsequence);grid on; title('MSsequence');
                        subplot(2,2,2);plot(MSsequence);grid on; title('MSsequence');
                        subplot(2,2,3);plot(embedding);grid on; title('Embedding');
                        subplot(2,2,4);plot(embedding2);grid on; title('Embedding');  
                        suptitle(strcat('Mean-Time is-', num2str(meanT)))
                        tmpPath2 = fullfile(tmpPath, 'Images');
                        if ~ isdir(tmpPath2)
                            mkdir(tmpPath2)
                        end
                        saveas(gca, fullfile(tmpPath2, strcat(num2str(grp_), '_MSseqEmbeArraySum_', num2str(micInd_), '_EEGMSA_', tmpName)), 'jpeg')                
                        close all
                        save(fullfile(tmpPath, strcat(num2str(grp_), '_EmbeddingWithArraySum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')), 'embedding');
                        save(fullfile(tmpPath, strcat(num2str(grp_), '_EmbeddingWithScalarSum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')), 'embedding2');
                        clear previous countArr statOrder MSsequence
                    end
                end
                if labelisnotthere > (length(EEGFiles)/2)
                    error('Some Problem In Calculation')
                end              
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%% Equalizing Time Here %%%%%%%%%%%%%%%%%%%%%
    elseif EQFlag == 1
        processedFile = [];
        meanTArray = [];       

        classAssGrp = [1,2,3,4;1,3,2,4;1,4,2,3];
        for micInd_ = 4 : 4
            fitbackSetPath = fullfile(targetPath, 'OnlySetFiles', strcat('With_',num2str(micInd_),'_MS'));
            if ~ isdir(fitbackSetPath)
                mkdir(fitbackSetPath)
            end            
            tmpPath2 = fullfile(targetPath, 'OnlyMat', strcat('With_',num2str(micInd_),'_MS'));
            if ~ isdir(tmpPath2)
                mkdir(tmpPath2)
            end   
            for grp_ = 1 : size(classAssGrp,1)
                if grp_ == 2
                    disp('Stop here')
                end
                cpArr = [];
                labelisnotthere = 0;
                for i = 1 : length(EEGFiles)
                    tmpName = split(EEGFiles(i).name, '.set');
                    tmpName = tmpName{1};                        
                    %if (~isfile(fullfile(tmpPath2, strcat('EQ_MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')))) || (~isfile(fullfile(tmpPath2, strcat(num2str(grp_), '_EQ_EmbeddingWithArraySum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')))) ||(~isfile(fullfile(tmpPath2, strcat(num2str(grp_), '_EQ_EmbeddingWithScalarSum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat'))))
                    if (isfile(fullfile(tmpPath2, strcat('EQ_MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')))) || (isfile(fullfile(tmpPath2, strcat(num2str(grp_), '_EQ_EmbeddingWithArraySum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')))) ||(isfile(fullfile(tmpPath2, strcat(num2str(grp_), '_EQ_EmbeddingWithScalarSum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat'))))
                        figure('units','normalized','outerposition',[0 0 1 1]);
                        EEG = pop_loadset('filename',strcat('Microstates_', num2str(micInd_), '_EEGMSA_', EEGFiles(i).name),'filepath', fitbackSetPath);  
                        MSsequence = EEG.microstate.fit.labels;

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

                        save(fullfile(tmpPath2, strcat('MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '_durationArr.mat')), 'countArr');
                        meanT = round(mean(countArr));
                        processedFile = [processedFile; {strcat('EQ_MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')}];
                        meanTArray = [meanTArray, meanT];

                        for st_ = 1 : length(statOrder)
                            lower = ((st_-1)*meanT) + 1;
                            upper = lower+meanT-1;
                            if st_ == length(statOrder)
                                upper = length(MSsequence);
                                EqMSsequence(1, lower:upper) = repmat(statOrder(st_), 1, (upper-lower+1));
                            else
                                EqMSsequence(1, lower:upper) = repmat(statOrder(st_), 1, (upper-lower+1));
                            end
                        end                        

                        Eqprevious = EqMSsequence(1);
                        EqcountArr(1) = 1;
                        EqstatOrder(1) = EqMSsequence(1);
                        for st_ = 2 : length(EqMSsequence)
                            Eqnew = EqMSsequence(st_);
                            if Eqprevious == Eqnew
                               EqcountArr(end) = EqcountArr(end)+1;                       
                            else
                                EqstatOrder(end+1) = Eqnew;
                                EqcountArr(end+1) = 1;
                            end
                            Eqprevious = Eqnew;
                        end

                        save(fullfile(tmpPath2, strcat('EQ_MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '_durationArr.mat')), 'EqcountArr');
                        %Uncoment save(fullfile(tmpPath2, strcat('EQ_MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')), 'EqMSsequence');

                        classAss = [];
                        classAss2 = 0;
                        embedding = [];
                        embedding2 = [];
                        for index_ = 1 : length(EqMSsequence)
                            disp(index_)
                            try
                            if ~ any(classAssGrp(grp_, :)>=(EqMSsequence(index_)))
                                disp('********************************Problem MSsequenceIndex is greater.************')
                                break;                                
                            end
                            catch
                                disp('Bheed')
                            end
                            if (EqMSsequence(index_) == classAssGrp(grp_, 1)) || (EqMSsequence(index_) == classAssGrp(grp_, 2))
                                classAss2 = classAss2 + 1;
                                classAss = [classAss 1];
                            end

                            if (EqMSsequence(index_) == classAssGrp(grp_, 3)) || (EqMSsequence(index_) == classAssGrp(grp_, 4))
                                classAss2 = classAss2-1;
                                classAss = [classAss -1];
                            end        
                            embedding2 = [embedding2 classAss2];
                            embedding = [embedding sum(classAss)];
                        end
                        if index_ ~= length(EqMSsequence)
                            close all
                            continue
                        end
                        
                        EEG.microstate.fit.classAssEQ = classAss;
                        EEG.microstate.fit.embedding2EQ = embedding2;
                        EEG.microstate.fit.embeddingEQ = embedding;
                        EEG.microstate.fit.MicrostatesMeanTime = meanT;

                        subplot(2,2,1);plot(MSsequence);grid on; title('MSsequence');
                        subplot(2,2,2);plot(EqMSsequence);grid on; title('MSsequence After Equalizing');
                        subplot(2,2,3);plot(embedding);grid on; title('Embedding with Equalization');
                        subplot(2,2,4);plot(embedding2);grid on; title('Embedding with Equalization2');
                        suptitle(strcat('Mean-Time is-', num2str(meanT)))
                        tmpPath3 = fullfile(tmpPath2, 'Images');
                        if ~ isdir(tmpPath3)
                            mkdir(tmpPath3)
                        end
                        tmpName = split(tmpName, '.set');
                        tmpName = tmpName{1};
                        %Uncoment saveas(gca, fullfile(tmpPath3, strcat(num2str(grp_), '_EQ_EmbeddingWithArraySum_', num2str(micInd_), '_EEGMSA_', tmpName)), 'jpeg')
                        close all
                        %Uncoment save(fullfile(tmpPath2, strcat(num2str(grp_), '_EQ_EmbeddingWithArraySum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')), 'embedding');
                        %Uncoment save(fullfile(tmpPath2, strcat(num2str(grp_), '_EQ_EmbeddingWithScalarSum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')), 'embedding2');
                        clear MSsequence previous countArr statOrder new classAss classAss2 embedding embedding2 EqMSsequence
                    end
                end
                if labelisnotthere > (length(EEGFiles)/2)
                    error('Some Problem In Calculation')
                end                              
            end
        end
        save(fullfile(tmpPath2, 'ProcessedFiles.mat'), 'processedFile')
        save(fullfile(tmpPath2, 'ArrayofMeanTime.mat'), 'meanTArray')
        %%%%%%%%%%%%%%%%%%%%%%%%%% Shuffling The Sequence Order %%%%%%%%%%%%%%%%%%%%%
    elseif EQFlag == 2
        classAssGrp = [1,2,3,4;1,3,2,4;1,4,2,3];
        for micInd_ = 4 : 4
            fitbackSetPath = fullfile(targetPath, 'OnlySetFiles', strcat('With_',num2str(micInd_),'_MS'));
            if ~ isdir(fitbackSetPath)
                mkdir(fitbackSetPath)
            end            
            tmpPath2 = fullfile(targetPath, 'OnlyMat', strcat('With_',num2str(micInd_),'_MS'));
            if ~ isdir(tmpPath2)
                mkdir(tmpPath2)
            end   
            for grp_ = 1 : size(classAssGrp,1)
                if grp_ == 2
                    disp('Stop here')
                end
                cpArr = [];
                labelisnotthere = 0;
                for i = 1 : length(EEGFiles)
                    tmpName = split(EEGFiles(i).name, '.set');
                    tmpName = tmpName{1};                        
                    if (~isfile(fullfile(tmpPath2, strcat('SH_MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')))) || (~isfile(fullfile(tmpPath2, strcat(num2str(grp_), '_SH_EmbeddingWithArraySum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')))) ||(~isfile(fullfile(tmpPath2, strcat(num2str(grp_), '_SH_EmbeddingWithScalarSum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat'))))
                        figure('units','normalized','outerposition',[0 0 1 1]);
                        EEG = pop_loadset('filename',strcat('Microstates_', num2str(micInd_), '_EEGMSA_', EEGFiles(i).name),'filepath', fitbackSetPath);  
                        MSsequence = EEG.microstate.fit.labels;
                        EqMSsequence = zeros(1, length(MSsequence));
                        
                        unq_ = unique(MSsequence);
                        unqShuffle = shuffle(unq_);       
                        for iii=1:length(unqShuffle)
                        EqMSsequence(find(MSsequence==unq_(iii)))=unqShuffle(iii);
                        end
                        %EqMSsequence(find(MSsequence==unq_(2)))=unqShuffle(2);
                        %EqMSsequence(find(MSsequence==unq_(3)))=unqShuffle(3);
                        %try
                        %EqMSsequence(find(MSsequence==unq_(4)))=unqShuffle(4);
                        %catch
                        %    disp('Why Stopped Here')
                        %end
     
                        %Uncoment save(fullfile(tmpPath2, strcat('SH_MSsequence_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')), 'EqMSsequence');

                        classAss = [];
                        classAss2 = 0;
                        embedding = [];
                        embedding2 = [];
                        for index_ = 1 : length(EqMSsequence)
                            disp(index_)
                            try
                            if ~ any(classAssGrp(grp_, :)>=(EqMSsequence(index_)))
                                disp('********************************Problem MSsSHuenceIndex is greater.************')
                                break;                                
                            end
                            catch
                                disp('Bheed')
                            end
                            if (EqMSsequence(index_) == classAssGrp(grp_, 1)) || (EqMSsequence(index_) == classAssGrp(grp_, 2))
                                classAss2 = classAss2 + 1;
                                classAss = [classAss 1];
                            end

                            if (EqMSsequence(index_) == classAssGrp(grp_, 3)) || (EqMSsequence(index_) == classAssGrp(grp_, 4))
                                classAss2 = classAss2-1;
                                classAss = [classAss -1];
                            end        
                            embedding2 = [embedding2 classAss2];
                            embedding = [embedding sum(classAss)];
                        end
                        if index_ ~= length(EqMSsequence)
                            close all
                            continue
                        end
                        
                        EEG.microstate.fit.classAssSH = classAss;
                        EEG.microstate.fit.embedding2SH = embedding2;
                        EEG.microstate.fit.embeddingSH = embedding;
                        %EEG.microstate.fit.MicrostatesMeanTime = meanT;

                        subplot(2,2,1);plot(MSsequence);grid on; title('MSsequence');
                        subplot(2,2,2);plot(EqMSsequence);grid on; title('MSsequence After Shuffling');
                        subplot(2,2,3);plot(embedding);grid on; title('Embedding with Shuffling');
                        subplot(2,2,4);plot(embedding2);grid on; title('Embedding with Shuffling2');
                        %suptitle(strcat('Mean-Time is-', num2str(meanT)))
                        tmpPath3 = fullfile(tmpPath2, 'Images');
                        if ~ isdir(tmpPath3)
                            mkdir(tmpPath3)
                        end
                        tmpName = split(tmpName, '.set');
                        tmpName = tmpName{1};
                        %Uncoment saveas(gca, fullfile(tmpPath3, strcat(num2str(grp_), '_SH_EmbeddingWithArraySum_', num2str(micInd_), '_EEGMSA_', tmpName)), 'jpeg')
                        close all
                        %Uncoment save(fullfile(tmpPath2, strcat(num2str(grp_), '_SH_EmbeddingWithArraySum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')), 'embedding');
                        %Uncoment save(fullfile(tmpPath2, strcat(num2str(grp_), '_SH_EmbeddingWithScalarSum_', num2str(micInd_), '_EEGMSA_', tmpName, '.mat')), 'embedding2');
                        clear MSsequence previous countArr statOrder new classAss classAss2 embedding embedding2 EqMSsequence
                    end
                end
                if labelisnotthere > (length(EEGFiles)/2)
                    error('Some Problem In Calculation')
                end                              
            end
        end        
    end
end
