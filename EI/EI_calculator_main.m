%% calculate EI of different channels
tic;
clc;
clear;
EEG = m00_import('yena_onset.m00');
% calculate Nd(i) and ER(i) for each channel, i is the channel index
for i = 1:EEG.Channels;
    [Nd(i),ER(i,:)] = EI_calc(EEG,EEG.data(i,:),i);
end

% calculate N_0 and EI for each channel
N_0 = min(Nd);
H = 5;
for i = 1:EEG.Channels
    ER_temp = ER(i,:);
       if Nd(i) <= length(ER_temp)/100 - H % acturally it is H * fs / step_length * fs
           EI(i) = (1 / (Nd(i) - N_0 + 1)) * sum(ER_temp(round(Nd(i)*100):(round(Nd(i)*100) + H*100)));
       else
           EI(i) = (1 / (Nd(i) - N_0 + 1));
       end
end
% normalization EI index
for i = 1:EEG.Channels
    EI(i) = EI(i)/max(EI);
end
toc;
 
