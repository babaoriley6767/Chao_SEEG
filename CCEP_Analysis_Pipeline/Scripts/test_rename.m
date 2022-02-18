spm('defaults','EEG');
D = spm_eeg_load();
S.D = D;
S.Channel_Labels_Raw = S.D.chanlabels';

S.Channel_Labels_New = Deblank_Names(S.Channel_Labels_Raw);

Pattern = '-Ref';
S.Channel_Labels_New = Remove_Name_Pattern(S.Channel_Labels_New,Pattern);

S.Channel_Labels_New = cellfun(@(x) x(3+1:end),S.Channel_Labels_New,'UniformOutput',false);

S.D = struct(S.D);
for i = 1:length(S.Channel_Labels_New)
    S.D.channels(i).label = S.Channel_Labels_New{i};
end
S.D = meeg(S.D);
save(S.D);




