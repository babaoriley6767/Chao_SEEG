

cd('/Users/Baotian/Desktop/Data/CCEP_data')
cd('/Users/Baotian/Desktop/ToolBoxes/CCEP_Pipeline/CCEP_Analysis_Pipeline')
spm('defaults', 'EEG');
%% Parameters settings
TimeWindow = [-500 1500];

%% File IO
% load the raw edf data and convert it to SPM format
D = edf_TO_SPM_converter_GUI([],[],'meeg_');

% load and convert the DC channel
DC = edf_TO_SPM_converter_GUI([],[],'DC_');

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

% define the trl
for i = 1:length(timeStampNew)
    trl(i,1)=timeStampNew(i)*DC.fsample - 500;
    trl(i,2)=timeStampNew(i)*DC.fsample + 1000;
    trl(i,3)=-500;
end

% Epoching
clear S
S.D = D;
S.bc = 0;
S.trl = trl;
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
% 1st, 2nd, 3rd, 4th are bandstop around 50Hz, 100Hz, 150Hz, 200Hz
% respectively
for i = 1:4
    clear S
    S.D              = D;
    S.band           = 'stop';
    S.freq           = [50*i-3 50*i+3];
    S.order          = 5;
    S.prefix         = 'f';
    D = spm_eeg_filter(S);
end

%% rename the electrode
Channel_Renaming_UI;  

                     % then delete the '-Ref' artificially
%% montage bipolar
S = D;
D = SPM_bipolar_montage(S,'BipM_');






%% these are for test
D = spm_eeg_load();
B = clone(D,'E:\mat\eeg\b_blank.mat');    %get the information from the entire edf file
S.D = D;
S.timewin = [19000 130000];
S.prefix = 'a1a2'
B  =spm_eeg_crop(S);

D=meeg(D);
save(D);




figure
plot(B(:,10000:190000));
S.D = B;
S.bc = 0;

S.trl =[1000 2000 -200;3000 4000 -200];
D = spm_eeg_epochs(S);
figure;
plot(D(1,:,1));

S.trl =[1000 2000;3000 4000];
D = spm_eeg_epochs(S);
figure;
plot(D(1,:,1));

plot(D(20,18000:111000))



D = spm_eeg_load();
for i = 1:144 
ChannelInd = i;
Data = squeeze(D(ChannelInd,:,:));
figure
plot(Data)
axis tight
end
idx = kmeans(Data',2)

figure
hold on
for i = 1:40
    if idx(i) == 1
        plot(Data(:,i),'r')
    elseif idx(i) == 2
        plot(Data(:,i),'g')
    end
end
