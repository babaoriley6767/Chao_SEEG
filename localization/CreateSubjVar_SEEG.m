function [subjVar, subjVar_created]  = CreateSubjVar_SEEG(sbj_name, comp_root, server_root, code_root,center,remark)
% function [subjVar, subjVar_created]  = CreateSubjVar(sbj_name, dirs, data_format)

%% Coordinate systems
% LEPTO_coord and native_coord the same
% so we use LEPTO_coord for when ploting in native space and
% MNI_coord for when ploting in fsaverage
% Sometimes the empty means cut channel (literaly physically cut from the grid)
% Maybe remove this after
%data_format = GetFSdataFormat(sbj_name, 'Stanford');

[fs_iEEG, data_format] = GetFSdataFormat(center);
dirs = InitializeDirs('SEEG_test', sbj_name, comp_root, server_root, code_root); % 'Pedro_NeuroSpin2T'
%sprintf('creating subjVar for subject %s', sbj_name)


if strcmp(data_format, 'edf')
    % Load a given globalVar
    gv_dir = dir(fullfile([dirs.data_root filesep 'originalData/' sbj_name]));
    gv_inds = arrayfun(@(x) contains(x.name, 'global') && ~contains(x.name, '._'), gv_dir);
    fn_tmp = gv_dir(gv_inds);
    load([fn_tmp(1).folder filesep fn_tmp(1).name])
else
end

%% remark the electrodes
%freesurfer location and filename
dirs.freesurfer

subDir = dirs.freesurfer;
subj = dir(subDir);
subj=subj(~ismember({subj.name},{'.','..', '.DS_Store'}) & horzcat(subj.isdir) == 1);
subj = subj.name;
subDir = [subDir subj];
avgDir=dirs.fsDir_local;

% Import electrode names
cd([subDir filesep 'elec_recon'])
elecFname=fullfile(subDir,'elec_recon',[subj '.electrodeNames']);
elecInfo=csv2Cell(elecFname,' ',2);


for ei = 1:length(elecInfo)
    elecInfo{ei,2} = 'D';
end

if remark
    
    %%%%%
    % statement is in here
    % remark %%%% working on this
    
    capExpr = '[a-zA-Z]';
    numExpr = '\d*';
    for ei = 1:length(elecInfo)
        ch_cap = regexp(elecInfo{ei,1},capExpr,'match');
        ch_num = regexp(elecInfo{ei,1},numExpr,'match');
        ch_type = regexp(elecInfo{ei,2},capExpr,'match');
        ch_hem = regexp(elecInfo{ei,3},capExpr,'match');
        ch_imread = strcat(ch_cap{:},'_',sprintf('%d', str2num(ch_num{:})));%sprintf('%02d', str2num(ch_num{:}))
        image_chan = imread([subDir filesep 'elec_recon' filesep 'PICS' filesep subj '_' ch_hem{:} ch_type{:} '_' ch_imread 'Slices.jpg']);
        
        %plot the raw electrodes position
        figureDim = [1 1 .8 .5];
        figure('units', 'normalized', 'outerposition', figureDim)
        image(image_chan);
        title([ch_cap{:}  ch_num{:} ' ' 'localization manully check'])
        set(gca,'fontsize',20)
        grid off;
        set(gca,'xticklabel',{[]})
        set(gca,'yticklabel',{[]})
        
        prompt = ['You need to identify if this electrodes belongs to the surface of depth?'...
            '\nClick y, this is a site belong to the surface'...
            '\nClick any other buttons this is a site belong to the depth structure\n'];
        
        ID = input(prompt,'s');
        
        if strcmp(ID, 'y')
            elecInfo{ei,2} = 'S';
        else
        end
        close all
    end
    
    
else
    
end
elecInfo_row = cell(3,length(elecInfo));
for ci = 1: length(elecInfo)
    elecInfo_row(:,ci) = elecInfo(ci,:);
end

fileID = fopen([subj '.electrodeNames'],'w');
fprintf(fileID,'%2s\n','This is for the SEEG EI/functional network Project');
fprintf(fileID,'%2s\n','Chao Zhang');
fprintf(fileID,'%2s%2s%2s\n',elecInfo_row{:});
fclose(fileID);

%% Cortex and channel label correction
subjVar = [];
cortex = getcort(dirs);
% native_coord = importCoordsFreesurfer(dirs);
% fs_chan_names = importElectNames(dirs);
[MNI_coord, chanInfo, RAS_coord] = sub2AvgBrainCustom_CSS_forBeijing([],dirs, sbj_name, dirs.fsDir_local);
[MGRID_coord, elect_names] = getmgrid(dirs);%elect_names pre with L/R Chao


%CLARA CHNAGE HERE
fs_chan_names = chanInfo.Name;
%fs_chan_names = chanInfo.elecNames

%introducing a chance for one patient here C18_29 %STILL NEEDS TO BE FIXED
% add ''' to the left side of incorrect PT029 fs_chan_namese
if strcmp(sbj_name, 'C18_29')||strcmp(sbj_name, 'C17_21')
    fs_chan_names_temp = chanInfo.elecNames;
    for i=1:length(chanInfo.elecNames)
        if strcmp(fs_chan_names_temp{i}(1),'L')
            digit_indx = isstrprop(chanInfo.Name{i,1},'digit');
            digit_indx = find(digit_indx,1);
            fs_chan_names_new{i,1} = [chanInfo.Name{i}(1:digit_indx-1),char(39),chanInfo.Name{i}(digit_indx:end)];%%chao modified here, plz thinking of the problem of a'1 and aO'1
        elseif strcmp(fs_chan_names_temp{i}(1),'R')
            fs_chan_names_new{i,1} = chanInfo.Name{i};
        end
        fs_chan_names = fs_chan_names_new;
    end
else
    fs_chan_names = chanInfo.Name;
end



close all
V = importVolumes(dirs);

subjVar.sbj_name = sbj_name;
subjVar.cortex = cortex;
subjVar.V = V;


%% Correct channel name
% Load naming from google sheet

%%CSS change here:
% 
% if isempty(sheet)
%     sheet = 'generic';
% else
% end


% DOCID = '1ukMue_vc69c5XwrrOKQ9ll59P_qhs3QHYT3fh0kmPVM';
% GID = '0';
% 
% if strcmp(data_format, 'edf') && ~strcmp(sbj_name, 'S17_118_TW')
%     [DOCID,GID] = getGoogleSheetInfo('chan_names_ppt', 'chan_names_fs_figures');
% else
%     [DOCID,GID] = getGoogleSheetInfo('chan_names_ppt', 'chan_names_ppt_log');
% end

% googleSheet = GetGoogleSpreadsheet(DOCID, GID);
% ppt_chan_names = googleSheet.(sbj_name);
ppt_chan_names = globalVar.channame';%this is where I changed Chao
ppt_chan_names = ppt_chan_names(~cellfun(@isempty, ppt_chan_names)); % remove empty cells
ppt_chan_names = cellfun(@(x) strrep(x, ' ', ''), ppt_chan_names, 'UniformOutput', false); % Remove eventual spaces

nchan_fs = length(fs_chan_names);
if strcmp(sbj_name, 'S17_117_MC')
    chan_comp = globalVar.channame;
    nchan_cmp = length(globalVar.channame);
else
    chan_comp = ppt_chan_names;
    nchan_cmp = length(ppt_chan_names);
end

% channels in EDF/TDT which are not in FS
% this could be because those channels were not recorded, but only implanted (EDF)
%               because somebody forgot to paste the correspondent .mat TDT file
%(this would indicate an error, and could potentially vary across block)
in_chan_cmp = false(1,nchan_fs);
for i = 1:nchan_fs
    in_chan_cmp(i) = ismember(fs_chan_names(i),chan_comp);%ismember(A,B),if element of A belong to element in B
end

in_fs = false(1,nchan_cmp);
for i = 1:nchan_cmp
    in_fs(i) = ismember(chan_comp(i),fs_chan_names);
end

    % perfect match, do nothing
if sum(in_chan_cmp) == length(in_chan_cmp) && sum(in_fs) == length(in_fs)
    disp('perfect match of Freesurfer and channles from googlesheet/EDF')
    % 1: More channels in freesurfer
elseif sum(in_chan_cmp) < length(in_chan_cmp) && sum(in_fs) == length(in_fs)
    disp('More channels in freesurfer, keep the channle names in googlesheet/EDF')
    fs_chan_names = fs_chan_names(in_chan_cmp);
    RAS_coord = RAS_coord(in_chan_cmp,:);
    MNI_coord = MNI_coord(in_chan_cmp,:);
    MGRID_coord = MGRID_coord(in_chan_cmp,:);%remove fs chan which was not shown in googlesheet
    % 2: More channels in EDF/TDT
elseif sum(in_chan_cmp) == length(in_chan_cmp) && sum(in_fs) < length(in_fs)
    disp('More channels in googlesheet/EDF, keep the freesurfer channle names')
    fs_chan_names_tmp = cell(nchan_cmp,1);
    fs_chan_names_tmp(in_fs) = fs_chan_names;
    fs_chan_names_tmp(in_fs==0) = chan_comp(in_fs==0);
    fs_chan_names = fs_chan_names_tmp;
    
    RAS_coord_tmp = nan(nchan_cmp,3,1);
    RAS_coord_tmp(in_fs,:) = RAS_coord;
    RAS_coord = RAS_coord_tmp;
    
    MNI_coord_tmp = nan(nchan_cmp,3,1);
    MNI_coord_tmp(in_fs,:) = MNI_coord;
    MNI_coord = MNI_coord_tmp;
    
    MGRID_coord_tmp = nan(nchan_cmp,3,1);
    MGRID_coord_tmp(in_fs,:) = MGRID_coord; %remove google chan which was not shown in fs
    MGRID_coord = MGRID_coord_tmp;

    
    % More in
elseif sum(in_chan_cmp) < length(in_chan_cmp) && sum(in_fs) < length(in_fs)
    
    disp(sbj_name)
    disp('channels in EDF/TDT which are not in FS')
    chan_comp(in_fs == 0)
    disp('channels in FS which are not in EDF/TDT')
    fs_chan_names(in_chan_cmp == 0)
    warning('this exception is not automatically fixable, please decide:')
    
    %     prompt = 'Do you want to remove the FS-only and add the EDF/TDT-only?';
    %     ID = input(prompt,'s');
    ID = 'y';
    if strcmp(ID, 'y')
        % First remove the FS which are not in EDF/TDT
        fs_chan_names = fs_chan_names(in_chan_cmp);
        RAS_coord = RAS_coord(in_chan_cmp,:);
        MNI_coord = MNI_coord(in_chan_cmp,:);
        MGRID_coord =  MGRID_coord(in_chan_cmp,:);
        
        % Second add the EDF/TDT which are not in FS
        fs_chan_names_tmp = cell(nchan_cmp,1);
        fs_chan_names_tmp(in_fs) = fs_chan_names;
        fs_chan_names_tmp(in_fs==0) = chan_comp(in_fs==0);
        fs_chan_names = fs_chan_names_tmp;
        
        %         native_coord_tmp = nan(size(native_coord,1),size(native_coord,2),1);
        RAS_coord_tmp = nan(size(RAS_coord,1),size(RAS_coord,2),1);
        MNI_coord_tmp = nan(size(MNI_coord,1),size(MNI_coord,2),1);
        MGRID_coord_tmp = nan(size(MGRID_coord,1),size(MGRID_coord,2),1);

        if in_fs(end) == 0
            RAS_coord_tmp(end+1,:) = nan;
            MNI_coord_tmp(end+1,:) = nan;
            MGRID_coord_tmp(end+1,:) = nan;
        else
        end
        
        RAS_coord_tmp(in_fs,:) = RAS_coord;
        RAS_coord = RAS_coord_tmp;
        
        MNI_coord_tmp(in_fs,:) = MNI_coord;
        MNI_coord = MNI_coord_tmp;
        
        MGRID_coord_tmp(in_fs,:) = MGRID_coord;
        MGRID_coord = MGRID_coord_tmp;        
    else
        warning('channel labels not fixed, please double check PPT/FS')
        mismatch_labels = 1;
    end
end

% CHANGE HERE CSS
if ~exist('mismatch_labels')
    %% Reorder and save in subjVar
    new_order = nan(1,nchan_cmp);
    for i = 1:nchan_cmp
        tmp = find(ismember(fs_chan_names, chan_comp(i)));
        if ~isempty(tmp)
            new_order(i) = tmp(1);
        end
    end
    
    subjVar.LEPTO_coord = RAS_coord(new_order,:); %order sequency based on the googlesheet/EDF
    subjVar.MNI_coord = MNI_coord(new_order,:);
    subjVar.MGRID_coord = MGRID_coord(new_order,:);
    
    % labels mean the corrected names
    if strcmp(data_format, 'TDT')
        subjVar.labels = chan_comp;
        subjVar.labels_EDF = [];
    else
        %     subjVar.elect_names = chan_comp;
        if isfield(globalVar, 'channame')
            subjVar.labels_EDF = globalVar.channame;
            subjVar.labels = chan_comp;
        else
            subjVar.labels = chan_comp;
            subjVar.labels_EDF = [];
        end
    end
    
    
    %% Demographics % chao changes here
    subjVar.demographics = [];
    
%     subjVar.demographics = GetDemographics(sbj_name);
%     if isempty(subjVar.demographics)
%         warning(['There is no demographic info for ' sbj_name '. Please add it to the google sheet.'])
%     else
%     end
%     
%     if ~exist([dirs.original_data filesep sbj_name], 'dir')
%         mkdir([dirs.original_data filesep sbj_name])
%     else
%     end
    
    %% Electrode labelling
    subjVar = localizeElect_CSS_forBeijing(subjVar,dirs,chanInfo);
    
    %% Save subjVar
    if exist([dirs.original_data filesep sbj_name filesep 'subjVar_' sbj_name '.mat'], 'file')
        %         prompt = ['subjVar already exist for ' sbj_name ' . Replace it? (y or n):'] ;
        %         ID = input(prompt,'s');
        ID = 'y';
        if strcmp(ID, 'y')
            save([dirs.original_data filesep sbj_name filesep 'subjVar_' sbj_name '.mat'], 'subjVar')
            disp(['subjVar saved for ' sbj_name])
            subjVar_created = 1;
        else
            warning(['subjVar NOT saved for ' sbj_name])
        end
    else
        save([dirs.original_data filesep sbj_name filesep 'subjVar_' sbj_name '.mat'], 'subjVar')
        disp(['subjVar saved for ' sbj_name])
        subjVar_created = 1;
    end
    
else
    subjVar_created = 0;
end

end

