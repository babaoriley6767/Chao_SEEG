function wave = noiseFiltData(globalVar, wave)
% Filtering 60 Hz line noise
% Dependencies:notch.m .eeglab toolbox
% Writen by Mohammad Dastjerdi, Parvizi Lab, Stanford
% Revision date SEP,2009

% filtering 60 Hz
if strcmp(globalVar.center, 'Stanford')
    [wave]= notch(wave, globalVar.iEEG_rate, 59, 61,1);
    [wave]= notch(wave, globalVar.iEEG_rate, 118,122,1); % Second harmonic of 60
    [wave]= notch(wave, globalVar.iEEG_rate, 178,182,1); % Third harmonic of 60
elseif strcmp(globalVar.center, 'China')
    [wave]= notch(wave, globalVar.iEEG_rate, 49, 51,1);
    [wave]= notch(wave, globalVar.iEEG_rate, 98,102,1); % Second harmonic of 60
    [wave]= notch(wave, globalVar.iEEG_rate, 148,152,1); % Third harmonic of 60
else
end


end