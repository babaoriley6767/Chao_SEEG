% loop for CCEP analysis

clear
% path manipulation
[ScriptFolder,SubjectFolder] = PathManipulation('BaotianZ620');

% ScriptFolder = 'D:\EmotionalFaces';
cd(ScriptFolder)

% SubjectFolder = 'E:\';
cd(SubjectFolder)

load('CCEPChannels.mat')
cd('E:\CCEP')

CCEPDir = dir(pwd);
CCEPDir = CCEPDir(3:end);

for Subjects = 1:length(CCEPDir)
    cd(CCEPDir(Subjects).name)
    edfFiles = dir('*.edf');
    for i = 1:length(edfFiles)
        
        mkdir(edfFiles(i).name(1:end-4));
        movefile(edfFiles(i).name,edfFiles(i).name(1:end-4));
        cd(edfFiles(i).name(1:end-4));
        
        %% File IO
        % load the raw edf data and convert it to SPM format
        index = IndexInUse(Subjects).Channels;
        D = edf_TO_SPM_converter_GUI(edfFiles(i).name,index,'meeg_');
        
        % load and convert the DC channel
        loc_DC = IndexInUse(Subjects).DCChannel;
        DC = edf_TO_SPM_converter_GUI(edfFiles(i).name,loc_DC,'DC_');
        
        %% rename the electrode
        Channel_Labels_Raw = D.chanlabels';
        
        Channel_Labels_New = Deblank_Names(Channel_Labels_Raw);
        
        Pattern = '-Ref';
        Channel_Labels_New = Remove_Name_Pattern(Channel_Labels_New,Pattern);
        
        Channel_Labels_New = cellfun(@(x) x(3+1:end),Channel_Labels_New,'UniformOutput',false);
        
        D = struct(D);
        for i = 1:length(Channel_Labels_New)
            D.channels(i).label = Channel_Labels_New{i};
        end
        D = meeg(D);
        save(D);
        
        %% Downsampling
        % Downsample the data to 1000 if > 1000Hz
        if D.fsample > 1003
            clear S
            S.D = D;
            S.fsample_new = 1000;
            D = spm_eeg_downsample(S);
        end
        if DC.fsample > 1003
            clear S
            S.D = DC;
            S.fsample_new = 1000;
            DC = spm_eeg_downsample(S);
        end
        
        %% epoch
        timeStampNew = FindCCEPTriggers(DC);
        close;
        % define the trl
        for i = 1:length(timeStampNew)
            trial(i,1)=timeStampNew(i)*DC.fsample - 500;
            trial(i,2)=timeStampNew(i)*DC.fsample + 1000;
            trial(i,3)=-500;
        end
        
        % Epoching
        clear S
        S.D = D;
        S.bc = 0;
        S.trl = trial;
        S.prefix = 'e';
        D = spm_eeg_epochs(S);
        
        
        %% baseline correction
        clear S
        S.D = D;
        D = spm_eeg_bc(S);    % when epoching we have the negtive time , so dont need to set the timewindow
        
        
        %% filter
        clear S
        S.D = D;
        S.band = 'bandpass';
        S.freq = [1 300];
        D = spm_eeg_filter(S);
        
        %% Filter the signal 4 times, minimal preprocessing the raw resting data
        % 1st, 2nd, 3rd, 4th are bandstop around 50Hz, 100Hz, 150Hz, 200Hz,
        % 250Hz
        % respectively
        for i = 1:5
            clear S
            S.D              = D;
            S.band           = 'stop';
            S.freq           = [50*i-3 50*i+3];
            S.order          = 5;
            S.prefix         = 'f';
            D = spm_eeg_filter(S);
        end
        
        %% montage bipolar
        clear S
        S = D;
        D = SPM_bipolar_montage(S,'BipM_');
        
        %% plot
        mkdir('CCEP_Figs')
        cd('CCEP_Figs')
        parfor i = 1:length(D.chanlabels)
            ChannelInd = i;
            Data = squeeze(D(ChannelInd,:,:));
            figure
            plot(D.time,Data,'Color',[0.5 0.5 0.5]);
            grid on
            hold on
            plot(D.time,mean(Data,2),'LineWidth',2,'Color','r')
            %str1 = char(D.chanlabels(i))
            title(D.chanlabels(i));
            set(gca,'FontSize',14)
            %set(findobj(gca,'type','line'),'linew',4)
            set(gcf,'Position',[0 100 1920 600])
            print(D.chanlabels{i},'-dpng')
            close
        end
        cd ..
        cd ..
    end
    cd ..
end