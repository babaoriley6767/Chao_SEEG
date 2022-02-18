
spm('defaults', 'EEG');
edfFiles = dir('*.edf');

for i = 1:length(edfFiles)
    mkdir(edfFiles(i).name(1:end-4));
    movefile(edfFiles(i).name,edfFiles(i).name(1:end-4));
    cd(edfFiles(i).name(1:end-4));
    
    %% File IO
    % load the raw edf data and convert it to SPM format
    index =[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,39,40,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,67,68,69,70,71,72,73,74,75,76,77,78,79,80];
    D = edf_TO_SPM_converter_GUI(edfFiles(i).name,index,'meeg_');

    % load and convert the DC channel
    loc_DC = [38];
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
    timeStampNew =[];
    timeStampNew = FindCCEPTriggers(DC);
    close;
    trl =[];
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


                     
    %% montage bipolar
    clear S
    S = D;
    D = SPM_bipolar_montage(S,'BipM_');
    close;

    for jj =1:3   %1:length(D.chanlabels)
    mkdir(D.chanlabels{jj});
    cd(D.chanlabels{jj});
    for ii =1:length(D.events)
        clear DD;clear DDD;clear tDD;clear P;clear PP;clear l_PP;clear l_PP_baseline;clear tPP;
        DD = D(jj,:,ii);
        tDD = remove_art(DD,500);
        t = -0.5:1/1000:1;
        freq = 1:300;
        P = morlet_transform(tDD, t, freq);
        PP = squeeze(P);
        o_PP{ii} = PP;
        %imagesc(PP);
        %imagesc(log10(PP));   axis xy;    set(gca,'CLim',[0 10]);    print(num2str(i),'-dpng');        close;
        %l_PP =log10(PP);
        l_PP = PP;
        l_PP_baseline = l_PP(:,200:300);
        tPP = (l_PP - mean(l_PP_baseline,2)) ./mean(l_PP_baseline,2);
        imagesc(tPP);        axis xy;        set(gca,'CLim',[0 10]); print(num2str(ii),'-dpng');  close;
        baselined_PP(:,:,ii) = tPP;
        %tPP =l_PP;
        gamma_n1(ii) = sum(sum(tPP(30:100,510:550)))/71/41;
        gamma_n2(ii) = sum(sum(tPP(30:100,600:700)))/71/101;
        ripple_n1(ii) = sum(sum(tPP(100:200,510:550)))/101/41;
        ripple_n2(ii) = sum(sum(tPP(100:200,600:700)))/101/101;
        
    end
    save('gamma_n1.mat','gamma_n1');
    save('gamma_n2.mat','gamma_n2');
    save('ripple_n1.mat','ripple_n1');
    save('ripple_n2.mat','ripple_n2');
    save('baselined_PP.mat','baselined_PP');
    cd ..;
    clear gamma_n1;clear gamma_n2;clear ripple_n1;clear ripple_n2;clear baselined_PP;
    end
    
   cd ..  
end











