
function [subjVar_volume_BR, subjVar_volume_BR_created]  = CreateSubjVar_volume_BR(sbj_name, comp_root, server_root, code_root,center);


[fs_iEEG, data_format] = GetFSdataFormat(center);
dirs = InitializeDirs('SEEG_test', sbj_name, comp_root, server_root, code_root);
subDir = dirs.brainvisa;


% Load a given globalVar
if strcmp(data_format, 'edf')
    gv_dir = dir(fullfile([dirs.data_root filesep 'originalData/' sbj_name]));
    gv_inds = arrayfun(@(x) contains(x.name, 'global') && ~contains(x.name, '._'), gv_dir);
    fn_tmp = gv_dir(gv_inds);
    load([fn_tmp(1).folder filesep fn_tmp(1).name])
else
end


%% Creat a new folder for the Bipolor Brainvisa
folder_Brainvisa = [server_root,filesep,sbj_name,filesep,'Brainvisa'];
folder_Brainvisa_BR = [server_root,filesep,sbj_name,filesep,'Brainvisa_BR'];
if ~exist(folder_Brainvisa_BR)
    mkdir(folder_Brainvisa_BR);
    all_files = dir(fullfile(folder_Brainvisa));
    if isempty(all_files)
        warning('You don''t have the orignal Brainvisa files')
    else
        for i = 1:length(all_files)
            if contains(all_files(i).name, 'txt')
                copyfile([folder_Brainvisa,filesep,all_files(i).name], folder_Brainvisa_BR);
                display(sprintf('Copied txt file %s to %s', all_files(i).name, folder_Brainvisa_BR));
            else
            end
        end
        for i = 1:length(all_files)
            BrainMask(i) = contains(all_files(i).name, 'BrainMask.nii');
        end
        for i = 1:length(all_files)
            nii(i) = contains(all_files(i).name, '.nii');
        end
        nii = nii&(~BrainMask);
        for i  = 1:length(all_files)
            if nii(i)
                nii_length(i) = strlength(all_files(i).name);
            else
                nii_length(i) = 999;
            end
        end
        [~,nii_idx] = min(nii_length,[],2,'linear');
        copyfile([folder_Brainvisa,filesep,all_files(nii_idx).name], folder_Brainvisa_BR);
        display(sprintf('Copied txt file %s to %s', all_files(nii_idx).name, folder_Brainvisa_BR));
        
    end
else
end

%% get into the Brainvisa_BR and make the Bipola adjustment

cd(folder_Brainvisa_BR);
fname = dir;
if ~any(strcmp({fname.name}, 'ElecResults.csv'))
    waitfor(msgbox(['Choose electrode name file to make it bipola pattern ']));
    [FILENAME, PATHNAME] = uigetfile([folder_Brainvisa_BR ,'.txt'],'All Files (*.*)','MultiSelect','on');
    elecInfo=csv2Cell([PATHNAME,FILENAME]);
    
    % exception
    elecInfo_raw = elecInfo;
    disp('This is how the electrode looks like ------ ');
    disp(elecInfo_raw);
    prompt = (['You need to make it something like ''A8_L'''...
        '\nClick 1, if it is ''A_R8'''...
        '\nClick 2, if it is ''A8''\n']);
    
    ID = input(prompt,'s');
    
    if strcmp(ID, '1')
        tmp_chan_names = cell(length(elecInfo),1);
        tmp_chan_hem = cell(length(elecInfo),1);
        
        for i = 1:length(elecInfo)
            idcs = strfind(elecInfo{i},'_');
            tmp_chan_names{i,1} = [elecInfo{i,1}(1:idcs-1),elecInfo{i,1}(idcs+2:end)];
            tmp_chan_hem{i,1} = elecInfo{i,1}(idcs+1);
        end
        
        prompt = 'What are the is left channel? ';
        isleft = input(prompt);
        for i = 1:length(elecInfo)
            if ismember(i,isleft)
                tmp_chan_hem{i,1} = 'L';
            else
                tmp_chan_hem{i,1} = 'R';
            end
        end
    else
        tmp_chan_names = cell(length(elecInfo),1);
        tmp_chan_hem = cell(length(elecInfo),1);
        
        for i = 1:length(elecInfo)
            tmp_chan_names{i,1} = elecInfo{i,:};
        end
        
        prompt = 'What are the is left channel? ';
        isleft = input(prompt);
        for i = 1:length(elecInfo)
            if ismember(i,isleft)
                tmp_chan_hem{i,1} = 'L';
            else
                tmp_chan_hem{i,1} = 'R';
            end
        end
    end
    
    
    %sort name based on alphabet
    
    channame_zero = cell(length(elecInfo),1);
    capExpr = '[a-zA-Z]';
    numExpr = '\d*';
    for ei = 1:length(elecInfo)
        channame_cap = regexp(tmp_chan_names{ei,1},capExpr,'match');
        channame_num = regexp(tmp_chan_names{ei,1},numExpr,'match');
        channame_zero{ei} = strcat(channame_cap{:},sprintf('%02d', str2num(channame_num{:})));
    end
    
    
    [channame_zero_sort,izero,~] = unique(channame_zero);
    if length(channame_zero) > length(channame_zero_sort)
        error(['The name is duplicated, most likely that the previous electrode slot has not been closed!']);
    else
    end
    % deal with the elecInfo accordingly!!!
    elecInfo = tmp_chan_names(izero,:);
    
    %build the bipolar sequence
    fs_BR_idx = [];
    for bii=1:length(elecInfo)-1
        name1 = join(regexp(string(elecInfo{bii,1}),'[a-z]','Match','ignorecase'),'');
        name2 = join(regexp(string(elecInfo{bii+1,1}),'[a-z]','Match','ignorecase'),'');
        if strcmp(name1,name2)
            fs_BR_idx = [fs_BR_idx,bii];
        else
            continue
        end
    end
    
    elecInfo_BR = cell(length(fs_BR_idx),1);
    for i = 1:length(fs_BR_idx)
        ii = fs_BR_idx(i);
        elecInfo_BR{i,1} = [elecInfo{ii,1},'-',elecInfo{ii+1,1}];
    end
    
    chan_hem = tmp_chan_hem(izero);
    chan_hem = chan_hem(fs_BR_idx);
    
    fileID = fopen([sbj_name '_Name_BR_adjust.txt'],'w');
    fprintf(fileID,'%2s\n',elecInfo_BR{:});
    fclose(fileID);
    
    %% save Pos.txt
    
    waitfor(msgbox(['Choose Position name file to make it bipola pattern ']));
    [FILENAME, PATHNAME] = uigetfile([folder_Brainvisa_BR ,'.txt'],'All Files (*.*)','MultiSelect','on');
    PosInfo=csv2Cell([PATHNAME,FILENAME],' ');
    
    PosInfo = PosInfo(izero,:);
    
    PosInfo_BR = cell(length(fs_BR_idx),3);
    for i = 1:length(fs_BR_idx)
        ii = fs_BR_idx(i);
        PosInfo_BR{i,1} = num2str((str2double(PosInfo{ii,1})+str2double(PosInfo{ii+1,1}))/2);
        PosInfo_BR{i,2} = num2str((str2double(PosInfo{ii,2})+str2double(PosInfo{ii+1,2}))/2);
        PosInfo_BR{i,3} = num2str((str2double(PosInfo{ii,3})+str2double(PosInfo{ii+1,3}))/2);
    end
    
    PosInfo_BR_row = cell(3,length(fs_BR_idx));
    for ci = 1: length(fs_BR_idx)
        PosInfo_BR_row(:,ci) = PosInfo_BR(ci,:);
    end
    fileID = fopen([sbj_name  '_Pos_BR_adjust.txt'],'w');
    fprintf(fileID,'%2s %2s %2s\n',PosInfo_BR_row{:});
    fclose(fileID);
    
    %% save MNI_Pos.txt
    
    waitfor(msgbox(['Choose MNI Position name file to make it bipola pattern ']));
    [FILENAME, PATHNAME] = uigetfile([folder_Brainvisa_BR ,'.txt'],'All Files (*.*)','MultiSelect','on');
    MNIInfo=csv2Cell([PATHNAME,FILENAME],' ');
    
    MNIInfo = MNIInfo(izero,:);
    
    MNIInfo_BR = cell(length(fs_BR_idx),3);
    for i = 1:length(fs_BR_idx)
        ii = fs_BR_idx(i);
        MNIInfo_BR{i,1} = num2str((str2double(MNIInfo{ii,1})+str2double(MNIInfo{ii+1,1}))/2);
        MNIInfo_BR{i,2} = num2str((str2double(MNIInfo{ii,2})+str2double(MNIInfo{ii+1,2}))/2);
        MNIInfo_BR{i,3} = num2str((str2double(MNIInfo{ii,3})+str2double(MNIInfo{ii+1,3}))/2);
    end
    
    MNIInfo_BR_row = cell(3,length(fs_BR_idx));
    for ci = 1: length(fs_BR_idx)
        MNIInfo_BR_row(:,ci) = MNIInfo_BR(ci,:);
    end
    fileID = fopen([sbj_name  '_MNI_Pos_BR_adjust.txt'],'w');
    fprintf(fileID,'%2s %2s %2s\n',MNIInfo_BR_row{:});
    fclose(fileID);
    %% this is from Dello_label
    
    %get rid of the non_BR files
    waitfor(msgbox('Delete non Bipolar files'));
    [FILENAME, PATHNAME] = uigetfile([folder_Brainvisa_BR ,'.txt'],'All Files (*.*)','MultiSelect','on');
    for i = 1:length(FILENAME)
        delete ([PATHNAME,FILENAME{i}]);
    end
    save('BR_info.mat','elecInfo_BR','chan_hem');
    
else
    load('BR_info.mat');
end
fname = dir;
if ~any(strcmp({fname.name}, 'ElecResults.csv'))
    % run the SPM
    SegMRI
    
    % Atlas
    DepthEle = DepthElectrodes;
    
    % Parameters
    DepthEle.WhiteMatterPercentageThreshold = 0.9; % At least 90% of surrounding voxel should be grey matter
    AnatF     = dir('BNI*.nii');
    GreyF     = dir('c1*.nii');
    EleName   = dir('*_Name_BR_adjust.txt');
    ElePos    = dir('*_Pos_BR_adjust.txt');
    ElePos    = ElePos(~contains({ElePos.name},'MNI'));
    EleMNIPos = dir('*_MNI_Pos_BR_adjust.txt');
    
    DepthEle.AnatFile            = AnatF.name;
    DepthEle.ElectrodeNameFile   = EleName.name;
    DepthEle.ElectrodePosFile    = ElePos.name;
    DepthEle.ElectrodePosMNIFile = EleMNIPos.name;
    
    DepthEle.BrainMask = 'BrainMask.nii';
    
    DepthEle.GreyMask  = GreyF.name;
    
    DepthEle.AAL3Atlas  = '/Users/tony/Documents/Stanford/DELLO/DELLO/Atlas/AAL3v1_1mm.nii';
    DepthEle.Yeo7Atlas  = '/Users/tony/Documents/Stanford/DELLO/DELLO/Atlas/Yeo2011_7Networks_MNI152 (Yeo 2011).nii';
    DepthEle.AAL3LabelsFile = '/Users/tony/Documents/Stanford/DELLO/DELLO/Atlas/AAL3v1_1mm.nii.txt';
    DepthEle.Yeo7LabelsFile = '/Users/tony/Documents/Stanford/DELLO/DELLO/Atlas/Yeo2011_7Networks_MNI152 (Yeo 2011).txt';
    
    
    DepthEle.ReadData;
    DepthEle.LabelOutBrainElectrodes;
    DepthEle.LabelWhiteMatterElectrodes;
    DepthEle.LabelAAL3
    DepthEle.LabelYeo7
    DepthEle.ExportResultTable
    
    %save
    save('DepthEle.mat','DepthEle')
    load('ElecResults.mat')
else
    load('DepthEle.mat')
    load('ElecResults.mat')
end
%% adap to sheet

ID = cell(length(DepthEle.ElectrodeName),1);
bv_chan_names = elecInfo_BR;
bv_chan_hem = chan_hem;

for i = 1:length(DepthEle.ElectrodeName)
    ID{i,1} = sbj_name;
end

GreyPercentage = DepthEle.ElectrodeGreyPercentage;

AAL3idx = DepthEle.AAL3Index;

AAL3 = T.AAL3;

YEO7idx = DepthEle.Yeo7Index;

YEO7 = T.Yeo7;

RAS_coord = DepthEle.ElectrodePos;

MNI152_volume = DepthEle.ElectrodePosMNI;

MNI305_volume = zeros(length(DepthEle.ElectrodeName),3);%这个是想放在膨胀脑上的，目前有严重的问题

trans = [1.0022    0.0071   -0.0177    0.0528;...
    -0.0146    0.9990    0.0027   -1.5519;...
    0.0129    0.0094    1.0027   -1.2012];

for i = 1:length(DepthEle.ElectrodeName)
    coords = trans * [MNI152_volume(i,:),1]';
    coords = coords';
    MNI305_volume(i,:) = coords;
    coords = [];
end
% https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems
%% plot electrodes

visual_check = plotelec_brainvisa(bv_chan_names,bv_chan_hem,MNI152_volume,MNI305_volume, YEO7idx);




%% global and Dello reconcile
ppt_chan_names = globalVar.channame_BR';%this is where I changed Chao
ppt_chan_names = ppt_chan_names(~cellfun(@isempty, ppt_chan_names)); % remove empty cells
ppt_chan_names = cellfun(@(x) strrep(x, ' ', ''), ppt_chan_names, 'UniformOutput', false); % Remove eventual spaces
chan_comp = ppt_chan_names;

nchan_cmp = length(chan_comp);
nchan_bv = length(bv_chan_names);

% channels in EDF/TDT which are not in FS
% this could be those channels were not recorded, because someone forgot to
% close the channels

in_chan_cmp = false(1,nchan_bv);
for i = 1:nchan_bv
    in_chan_cmp(i) = ismember(bv_chan_names(i),chan_comp);%ismember(A,B),if element of A belong to element in B
end

in_bv = false(1,nchan_cmp);
for i = 1:nchan_cmp
    in_bv(i) = ismember(chan_comp(i),bv_chan_names);
end


if sum(in_chan_cmp) == length(in_chan_cmp) && sum(in_bv) == length(in_bv)
    disp('perfect match of brainvisa and channles from EDF')
elseif sum(in_chan_cmp) < length(in_chan_cmp) && sum(in_bv) == length(in_bv)
    disp('More channels in brainvisa, keep the channle names in EDF')
    ID = ID(in_chan_cmp);
    bv_chan_names = bv_chan_names(in_chan_cmp);
    bv_chan_hem = bv_chan_hem(in_chan_cmp);
    GreyPercentage  = GreyPercentage(in_chan_cmp);
    AAL3idx = AAL3idx(in_chan_cmp);
    AAL3 = AAL3(in_chan_cmp);
    YEO7idx = YEO7idx(in_chan_cmp);
    YEO7 = YEO7(in_chan_cmp);
    RAS_coord = RAS_coord(in_chan_cmp,:);
    MNI152_volume = MNI152_volume(in_chan_cmp,:);
    MNI305_volume = MNI305_volume(in_chan_cmp,:);
elseif sum(in_chan_cmp) == length(in_chan_cmp) && sum(in_bv) < length(in_bv)
    disp('More channels in EDF, keep the brainvisa channle names')
    
    ID_tmp = cell(nchan_cmp,1);
    ID_tmp(:) = {sbj_name};
    ID = ID_tmp;
    
    bv_chan_names_tmp = cell(nchan_cmp,1);
    bv_chan_names_tmp(in_bv) = bv_chan_names;
    bv_chan_names_tmp(in_bv==0) = chan_comp(in_bv==0);
    bv_chan_names = bv_chan_names_tmp;
    
    bv_chan_hem_tmp = cell(nchan_cmp,1);
    bv_chan_hem_tmp(in_bv) = bv_chan_hem;
    bv_chan_hem_tmp(in_bv==0) = {'NaN'};
    bv_chan_hem = bv_chan_hem_tmp;
    
    GreyPercentage_tmp = nan(nchan_cmp,1,1);
    GreyPercentage_tmp(in_bv,:) = GreyPercentage;
    GreyPercentage = GreyPercentage_tmp;
    
    AAL3idx_tmp = nan(nchan_cmp,1,1);
    AAL3idx_tmp(in_bv,:) = AAL3idx;
    AAL3idx = AAL3idx_tmp;
    
    AAL3_tmp = cell(nchan_cmp,1);
    AAL3_tmp(in_bv) = AAL3;
    AAL3_tmp(in_bv == 0) = {'NaN'};
    AAL3 = AAL3_tmp;
    
    
    YEO7idx_tmp = nan(nchan_cmp,1,1);
    YEO7idx_tmp(in_bv,:) = YEO7idx;
    YEO7idx = YEO7idx_tmp;
    
    YEO7_tmp = cell(nchan_cmp,1);
    YEO7_tmp(in_bv) = YEO7;
    YEO7_tmp(in_bv == 0) = {'NaN'};
    YEO7 = YEO7_tmp;
    
    RAS_coord_tmp = nan(nchan_cmp,3,1);
    RAS_coord_tmp(in_bv,:) = RAS_coord;
    RAS_coord = RAS_coord_tmp;
    
    MNI152_volume_tmp = nan(nchan_cmp,3,1);
    MNI152_volume_tmp(in_bv,:) = MNI152_volume;
    MNI152_volume = MNI152_volume_tmp;
    
    MNI305_volume_tmp = nan(nchan_cmp,3,1);
    MNI305_volume_tmp(in_bv,:) = MNI305_volume;
    MNI305_volume = MNI305_volume_tmp;
    
    % More in
elseif sum(in_chan_cmp) < length(in_chan_cmp) && sum(in_bv) < length(in_bv)
    
    disp(sbj_name)
    disp('channels in EDF which are not in bv')
    chan_comp(in_bv == 0)
    disp('channels in bv which are not in EDF/TDT')
    bv_chan_names(in_chan_cmp == 0)
    warning('this exception is not automatically fixable, please fix');
    mismatch_labels = 1;
else
end

%% save SubjVar_volume
if ~exist('mismatch_labels')
    %% Reorder and save in subjVar
    new_order = nan(1,nchan_cmp);
    for i = 1:nchan_cmp
        tmp = find(ismember(bv_chan_names, chan_comp(i)));
        if ~isempty(tmp)
            new_order(i) = tmp(1);
        end
    end
    
    RAS_coord = RAS_coord(new_order,:); %order sequency based on the EDF
    MNI152_volume = MNI152_volume(new_order,:);
    MNI305_volume = MNI305_volume(new_order,:);
    bv_chan_names = bv_chan_names(new_order,:);
    bv_chan_hem = bv_chan_hem(new_order,:);
    subjVar_volume.label_EDF = globalVar.channame;
    
    
    bv_names = bv_chan_names(new_order,:);
    bv_hem = bv_chan_hem(new_order,:);
    GreyPercentage  = GreyPercentage(new_order,:);
    AAL3idx = AAL3idx(new_order,:);
    AAL3 = AAL3(new_order,:);
    YEO7idx = YEO7idx(new_order,:);
    YEO7 = YEO7(new_order,:);
    RAS_coord = RAS_coord(new_order,:);
    MNI152_volume = MNI152_volume(new_order,:);
    MNI305_volume = MNI305_volume(new_order,:);
    
    subjVar_volume_BR.eleinfo = table(bv_names,bv_hem,GreyPercentage,AAL3idx,AAL3,YEO7idx,YEO7,RAS_coord,MNI152_volume,MNI305_volume);
    subjVar_volume_BR.label_EDF = globalVar.channame_BR;
    subjVar_volume_BR.visual_check = visual_check;
    
    
    cd(folder_Brainvisa_BR);
    
    %% Save subjVar
    if exist([dirs.original_data filesep sbj_name filesep 'subjVar_volume_BR_' sbj_name '.mat'], 'file')
        prompt = ['subjVar_volume_BR already exist for ' sbj_name ' . Replace it? (y or n):'] ;
        ID = input(prompt,'s');
        ID = 'y';
        if strcmp(ID, 'y')
            save([dirs.original_data filesep sbj_name filesep 'subjVar_volume_BR_' sbj_name '.mat'], 'subjVar_volume_BR')
            disp(['subjVar_volume_BR saved for ' sbj_name])
            subjVar_volume_BR_created = 1;
        else
            warning(['subjVar NOT saved for ' sbj_name])
        end
    else
        save([dirs.original_data filesep sbj_name filesep 'subjVar_volume_BR_' sbj_name '.mat'], 'subjVar_volume_BR')
        disp(['subjVar_volume_BR saved for ' sbj_name])
        subjVar_volume_BR_created = 1;
    end
    
else
    subjVar_volume_BR_created = 0;
end
end














