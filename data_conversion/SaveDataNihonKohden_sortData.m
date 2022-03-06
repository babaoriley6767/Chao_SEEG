function SaveDataNihonKohden(sbj_name, project_name, block_name, dirs, refChan, epiChan, emptyChan)

%% load the data to define and eliminate bad channels

% Loop across blocks
for i = 1:length(block_name)
    bn = block_name{i};
    
    % Load globalVar
    fn = sprintf('%s/originalData/%s/global_%s_%s_%s.mat',dirs.data_root,sbj_name,project_name,sbj_name,bn);
    load(fn,'globalVar');
    
    data_dir = [dirs.original_data '/' sbj_name '/' bn]; % directory for saving data
    fname =  [dirs.original_data '/' sbj_name '/' bn '/' bn '.edf'];
    
    % exeption should be list in here
    % Sometimes the electrode slot of the previous patient is not closed, and the number
    % of electrodes of the previous patient is greater than the current one. In this case,
    % the electrode name may be duplicated.Check the raw electrode data and cut off the electrode data could fix the error.
    if strcmp(sbj_name, 'C18_33')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:105]);
    elseif strcmp(sbj_name,'2017_FT_SP_018_xiangxiang')&strcmp(project_name,'SEEG_test')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:69]);
    elseif strcmp(sbj_name,'group1_SJY')&strcmp(project_name,'SEEG_test')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:60]);
    elseif strcmp(sbj_name,'group1_WXL')&strcmp(project_name,'SEEG_test')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:64]);
    elseif strcmp(sbj_name,'group2_CBB')&strcmp(project_name,'SEEG_test')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:40,43:84]);
    elseif strcmp(sbj_name,'2015_FT_hankunpeng')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:140]);
    elseif strcmp(sbj_name,'2015_FT_huangzutu')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:43,46:69]);
    elseif strcmp(sbj_name,'2015_FT_songziwei')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:43,46:67,72:103]);
    elseif strcmp(sbj_name,'2015_FT_wangmingming')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:119]);
    elseif strcmp(sbj_name,'2015_FT_wangzheng')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:40,43:64,69:78]);%get rid of the right side
    elseif strcmp(sbj_name,'2015_FT_xumengjiao')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:62,69:96]);
    elseif strcmp(sbj_name,'2015_FT_yangbihong')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:98]);
    elseif strcmp(sbj_name,'2015_FT_yinfengting')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:40,43:64,69:80]);
    elseif strcmp(sbj_name,'2015_FT_yuanjinrui')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:40,43:58]);
    elseif strcmp(sbj_name,'2015_FT_zhaoxichun')
        [hdr,D] = edfread_China(fname,'targetsignals',[1:89]);
    else
        [hdr, D] = edfread_China(fname);
    end
    
    
    
    %% Add Exception for when channels don' have labels
    hdr.label =  hdr.label(~strcmp(hdr.label, 'POL')); %get rid of the pure 'POL' naming
    D = D(~strcmp(hdr.label, 'POL'),:); %get rid of the pure 'POL' naming
    
    fs = size(D,2)/(hdr.records * hdr.duration);
    % hdr.records = number of chuncks
    % hdr.duration = duration of each chunck
    % fs = total points/total time
    
    % Downsampling parameters
    target_fs = 1000; % defalt 1000Hz
    %     target_fs_comp = round(target_fs/5); % reduced fs for spectral data
    
    if fs <= target_fs
        ecog_ds = 1;
    else
        ecog_ds = round(fs/target_fs); % decimate factor
    end
    
    %     pdio_ds = 1; %downsample for photodiode signals %%%%
    %     pdio_oldinds = find(contains(hdr.label,
    %     'DC'));%????NihonKohden??????????
    %     pdio_newinds = [1:length(pdio_oldinds)];
    
    % Take the indices of the channels of interest
    ecog_oldinds = find(~contains(hdr.label, 'EKG') & ~contains(hdr.label, 'DC') ...
        & ~contains(hdr.label, 'REF') & ~contains(hdr.label, 'Annotations')& ~strcmp(hdr.label, 'POLE')...
        & ~contains(hdr.label, 'POLBP')& ~contains(hdr.label, 'POL$')& ~contains(hdr.label, 'ECG')& ~contains(hdr.label, 'EMG')& ~contains(hdr.label, 'POL0'));
    idx_ecog_oldinds = ~contains(hdr.label, 'EKG') & ~contains(hdr.label, 'DC') ...
        & ~contains(hdr.label, 'REF') & ~contains(hdr.label, 'Annotations')& ~strcmp(hdr.label, 'POLE')...
        & ~contains(hdr.label, 'POLBP')& ~contains(hdr.label, 'POL$')& ~contains(hdr.label, 'ECG')& ~contains(hdr.label, 'EMG')& ~contains(hdr.label, 'POL0');
    % deal with the data accordingly
    D = D(idx_ecog_oldinds,:);
    
    % Loop across channels and clean the channel name
    ecog_newinds = 1:length(ecog_oldinds);
    channame = cell(size(ecog_newinds));
    for ei = 1:length(ecog_oldinds)
        % Clean channel name 'POL'
        channame_tpm = hdr.label{ecog_oldinds(ei)};
        channame_tpm = strrep(channame_tpm,'POL','');
        % This is to correct some chan labels from China: 'EEG','-Ref','Ref'
        if contains(channame_tpm, 'EEG') && contains (channame_tpm, '-')
            channame_tpm = strrep(channame_tpm,'EEG','');
            channame{ei} = strrep(channame_tpm,'-Ref','');
        else
            channame{ei} = strrep(channame_tpm,'Ref','');
        end
    end
    
    % zero padding in the channame (channame can only be sorted correctly after zero padding)
    channame_zero = cell(size(channame));
    capExpr = '[a-zA-Z]';
    numExpr = '\d*';
    for ei = 1:length(channame)
        channame_cap = regexp(channame{ei},capExpr,'match');
        channame_num = regexp(channame{ei},numExpr,'match');
        channame_zero{ei} = strcat(channame_cap{:},sprintf('%02d', str2num(channame_num{:})));
    end
    
    [channame_zero_sort,izero,~] = unique(channame_zero);
    if length(channame_zero) > length(channame_zero_sort)
        error(['The name is duplicated, most likely that the previous electrode slot has not been closed!',...
            'Go to check on line 16'])
    end
    % deal with the data accordingly!!!
    D = D(izero,:);
    
    % un-zero padding in the channame
    channame_sort = cell(size(channame_zero_sort));
    for ei = 1:length(channame_sort)
        channame_cap = regexp(channame_zero_sort{ei},capExpr,'match');
        channame_num = regexp(channame_zero_sort{ei},numExpr,'match');
        channame_sort{ei} = strcat(channame_cap{:},sprintf('%d', str2num(channame_num{:})));
    end
    
    % save the single electrode data in the originalData folder and
    % dicimate
    for ei = 1:length(channame_sort)
        if ei<10
            chanlbl = ['0',num2str(ei)];
        else
            chanlbl = num2str(ei);
        end
        fp = sprintf('%s/iEEG%s_%s.mat',data_dir,bn,chanlbl);
        wave = squeeze(D(ei,:,1));
        if (ecog_ds > 1)
            wave = decimate(double(wave),ecog_ds);
        end
        channame_sort_tpm = channame_sort{ei};
        save(fp,'wave','fs','channame_sort_tpm')
        disp(['Saving chan ',chanlbl,' ',channame_sort{ei}]);
    end
    
    
    %     for pi = 1:length(pdio_oldinds)
    %         chanlbl = ['0',num2str(pdio_newinds(pi))];
    %         fp = sprintf('%s/Pdio%s_%s.mat',data_dir,bn,chanlbl);
    %         anlg = squeeze(D(pdio_oldinds(pi),:,1));
    %         if (pdio_ds > 1)
    %             anlg = decimate(double(anlg),pdio_ds);
    %         end
    %         save(fp,'anlg','fs')
    %         disp(['Saving pdio ',chanlbl])
    %     end
    
    %% Update global variable
    
    globalVar.iEEG_rate = fs/ecog_ds;
    %     globalVar.Pdio_rate = fs/pdio_ds;
    %     globalVar.fs_comp= target_fs_comp;
    globalVar.channame = channame_sort;
    globalVar.chanLength = length(wave);
    globalVar.nchan = length(globalVar.channame);
    globalVar.refChan = refChan;
    globalVar.epiChan = epiChan;
    globalVar.emptyChan = emptyChan;
    
    save(fn,'globalVar');
    disp('globalVar updated')
    
end