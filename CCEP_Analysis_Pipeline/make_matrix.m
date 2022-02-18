CCEPdir = dir();
CCEPdir = CCEPdir(3:end);
for i = 1:length(CCEPdir)
    cd(CCEPdir(i).name);
    edfFiles = dir('*.mat');
    temp = 0;max_num= 0 ;
    for j = 1:length(edfFiles)
        if length(edfFiles(j).name)>temp
            temp = length(edfFiles(j).name);
            max_num = j;
        end
    end
    edfFile = edfFiles(max_num);
    D = spm_eeg_load(edfFile.name);
    Pattern = '-';
    chanlabels_new = D.chanlabels;
    chanlabels_new = Remove_Name_Pattern(chanlabels_new,Pattern);
    Evoked = find(ismember(chanlabels_new,CCEPdir(i).name));
    for j= 1:length(D.chanlabels)
        datas = squeeze(D(j,:,:));
        data = mean(datas,2);
        M(Evoked,j) = Calc_RMS(data,516,801);
    end
    cd ..
end
save('PTXXX_Matrix.mat','M');
save('PTXXX_chanlabels.mat','chanlabels_new');
