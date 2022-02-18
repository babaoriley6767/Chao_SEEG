function [subjVar_BR, subjVar_created]  = CreateSubjVar_SEEG_BR(sbj_name, comp_root, server_root, code_root,center,remark)
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

%% copy the files and change it into BR format


dirs.freesurfer

subDir = dirs.freesurfer;
subj = dir(subDir);
subj=subj(~ismember({subj.name},{'.','..', '.DS_Store'}) & horzcat(subj.isdir) == 1);
subj = subj.name;
subDir = [subDir subj];
avgDir=dirs.fsDir_local;

%this code can copy files and 1 layer of folder into another folder
if ~exist([subDir filesep 'elec_recon_BR'])
    mkdir([subDir filesep 'elec_recon_BR']);
    cd([subDir filesep 'elec_recon'])
    filenames=dir;
    for i=3:length(filenames)
        if ~isfolder([filenames(i).folder,filesep, filenames(i).name])
            copyfile(filenames(i).name,[subDir filesep 'elec_recon_BR'])
        else
            mkdir([subDir filesep 'elec_recon_BR' filesep filenames(i).name]);
            cd([subDir filesep 'elec_recon' filesep filenames(i).name]);
            filenames_fold = dir;
            for j = 3:length(filenames_fold)
                copyfile(filenames_fold(j).name,[subDir filesep 'elec_recon_BR' filesep filenames(i).name]);
            end
            cd([subDir filesep 'elec_recon'])
        end
    end
else
end

%% take care of the .electrode .lepto .leptovox .pial .pialvox imploc.txt and .mgrid files
% Import electrode names
cd([subDir filesep 'elec_recon_BR'])
elecFname=fullfile(subDir,'elec_recon_BR',[subj '.electrodeNames']);
elecInfo=csv2Cell(elecFname,' ',2);


%%%%%
% statement is in here
% remark %%%% working on this

capExpr = '[a-zA-Z]';
numExpr = '\d*';
if remark
    for ei = 1:length(elecInfo)
    ch_cap = regexp(elecInfo{ei,1},capExpr,'match');
    ch_num = regexp(elecInfo{ei,1},numExpr,'match');
    ch_type = regexp(elecInfo{ei,2},capExpr,'match');
    ch_hem = regexp(elecInfo{ei,3},capExpr,'match');
    ch_imread = strcat(ch_cap{:},'_',sprintf('%d', str2num(ch_num{:})));%sprintf('%02d', str2num(ch_num{:}))
    image_chan = imread([subDir filesep 'elec_recon_BR' filesep 'PICS' filesep subj '_' ch_hem{:} ch_type{:} '_' ch_imread 'Slices.jpg']);
    
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


%sort name based on alphabet

channame_zero = cell(length(elecInfo),1);
capExpr = '[a-zA-Z]';
numExpr = '\d*';
for ei = 1:length(elecInfo)
    channame_cap = regexp(elecInfo{ei,1},capExpr,'match');
    channame_num = regexp(elecInfo{ei,1},numExpr,'match');
    channame_zero{ei} = strcat(channame_cap{:},sprintf('%02d', str2num(channame_num{:})));
end


[channame_zero_sort,izero,~] = unique(channame_zero);
if length(channame_zero) > length(channame_zero_sort)
    error(['The name is duplicated, most likely that the previous electrode slot has not been closed!',...
        ])
end
% deal with the elecInfo accordingly!!!
elecInfo = elecInfo(izero,:);

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

elecInfo_BR = cell(length(fs_BR_idx),3);
for i = 1:length(fs_BR_idx)
    ii = fs_BR_idx(i);
    elecInfo_BR{i,1} = [elecInfo{ii,1},'-',elecInfo{ii+1,1}];
    elecInfo_BR{i,2} = elecInfo{ii+1,2};%ii is strict for surfer and ii+1 is loose for surfer
    elecInfo_BR{i,3} = elecInfo{ii,3};
end



elecInfo_BR_row = cell(3,length(elecInfo_BR));
for ci = 1: length(elecInfo_BR)
    elecInfo_BR_row(:,ci) = elecInfo_BR(ci,:);
end
fileID = fopen([subj '.electrodeNames_BR'],'w');
fprintf(fileID,'%2s\n','This is for the SEEG EI/functional network Project');
fprintf(fileID,'%2s\n','Chao Zhang');
fprintf(fileID,'%2s%2s%2s\n',elecInfo_BR_row{:});
fclose(fileID);

% XFname=fullfile(subDir,'elec_recon_BR',[subj '.electrodeNames_BR']);
% XInfo=csv2Cell(XFname,' ');

%% save .CT
CTFname=fullfile(subDir,'elec_recon_BR',[subj '.CT']);
CTInfo=csv2Cell(CTFname,' ',2);

CTInfo = CTInfo(izero,:);

CTInfo_BR = cell(length(fs_BR_idx),3);
for i = 1:length(fs_BR_idx)
    ii = fs_BR_idx(i);
    CTInfo_BR{i,1} = num2str((str2double(CTInfo{ii,1})+str2double(CTInfo{ii+1,1}))/2);
    CTInfo_BR{i,2} = num2str((str2double(CTInfo{ii,2})+str2double(CTInfo{ii+1,2}))/2);
    CTInfo_BR{i,3} = num2str((str2double(CTInfo{ii,3})+str2double(CTInfo{ii+1,3}))/2);
end

CTInfo_BR_row = cell(3,length(fs_BR_idx));
for ci = 1: length(fs_BR_idx)
    CTInfo_BR_row(:,ci) = CTInfo_BR(ci,:);
end
fileID = fopen([subj '.CT_BR'],'w');
fprintf(fileID,'%2s\n','This is for the SEEG EI/functional network Project');
fprintf(fileID,'%2s\n','Chao Zhang');
fprintf(fileID,'%2s %2s %2s\n',CTInfo_BR_row{:});
fclose(fileID);
% 
% XFname=fullfile(subDir,'elec_recon_BR',[subj '.CT_BR']);
% XInfo=csv2Cell(XFname,' ');

%% save .LEPTO
LEPTOFname=fullfile(subDir,'elec_recon_BR',[subj '.LEPTO']);
LEPTOInfo=csv2Cell(LEPTOFname,' ',2);

LEPTOInfo = LEPTOInfo(izero,:);

LEPTOInfo_BR = cell(length(fs_BR_idx),3);
for i = 1:length(fs_BR_idx)
    ii = fs_BR_idx(i);
    LEPTOInfo_BR{i,1} = num2str((str2double(LEPTOInfo{ii,1})+str2double(LEPTOInfo{ii+1,1}))/2);
    LEPTOInfo_BR{i,2} = num2str((str2double(LEPTOInfo{ii,2})+str2double(LEPTOInfo{ii+1,2}))/2);
    LEPTOInfo_BR{i,3} = num2str((str2double(LEPTOInfo{ii,3})+str2double(LEPTOInfo{ii+1,3}))/2);
end

LEPTOInfo_BR_row = cell(3,length(fs_BR_idx));
for ci = 1: length(fs_BR_idx)
    LEPTOInfo_BR_row(:,ci) = LEPTOInfo_BR(ci,:);
end
fileID = fopen([subj '.LEPTO_BR'],'w');
fprintf(fileID,'%2s\n','This is for the SEEG EI/functional network Project');
fprintf(fileID,'%2s\n','Chao Zhang');
fprintf(fileID,'%2s %2s %2s\n',LEPTOInfo_BR_row{:});
fclose(fileID);

% XFname=fullfile(subDir,'elec_recon_BR',[subj '.LEPTO_BR']);
% XInfo=csv2Cell(XFname,' ');

%% save .LEPTOVOX
LEPTOVOXFname=fullfile(subDir,'elec_recon_BR',[subj '.LEPTOVOX']);
LEPTOVOXInfo=csv2Cell(LEPTOVOXFname,' ',2);

LEPTOVOXInfo = LEPTOVOXInfo(izero,:);

LEPTOVOXInfo_BR = cell(length(fs_BR_idx),3);
for i = 1:length(fs_BR_idx)
    ii = fs_BR_idx(i);
    LEPTOVOXInfo_BR{i,1} = num2str((str2double(LEPTOVOXInfo{ii,1})+str2double(LEPTOVOXInfo{ii+1,1}))/2);
    LEPTOVOXInfo_BR{i,2} = num2str((str2double(LEPTOVOXInfo{ii,2})+str2double(LEPTOVOXInfo{ii+1,2}))/2);
    LEPTOVOXInfo_BR{i,3} = num2str((str2double(LEPTOVOXInfo{ii,3})+str2double(LEPTOVOXInfo{ii+1,3}))/2);
end

LEPTOVOXInfo_BR_row = cell(3,length(fs_BR_idx));
for ci = 1: length(fs_BR_idx)
    LEPTOVOXInfo_BR_row(:,ci) = LEPTOVOXInfo_BR(ci,:);
end
fileID = fopen([subj '.LEPTOVOX_BR'],'w');
fprintf(fileID,'%2s\n','This is for the SEEG EI/functional network Project');
fprintf(fileID,'%2s\n','Chao Zhang');
fprintf(fileID,'%2s %2s %2s\n',LEPTOVOXInfo_BR_row{:});
fclose(fileID);

% XFname=fullfile(subDir,'elec_recon_BR',[subj '.LEPTOVOX_BR']);
% XInfo=csv2Cell(XFname,' ');

%% save .PIAL
PIALFname=fullfile(subDir,'elec_recon_BR',[subj '.PIAL']);
PIALInfo=csv2Cell(PIALFname,' ',2);

PIALInfo = PIALInfo(izero,:);

PIALInfo_BR = cell(length(fs_BR_idx),3);
for i = 1:length(fs_BR_idx)
    ii = fs_BR_idx(i);
    PIALInfo_BR{i,1} = num2str((str2double(PIALInfo{ii,1})+str2double(PIALInfo{ii+1,1}))/2);
    PIALInfo_BR{i,2} = num2str((str2double(PIALInfo{ii,2})+str2double(PIALInfo{ii+1,2}))/2);
    PIALInfo_BR{i,3} = num2str((str2double(PIALInfo{ii,3})+str2double(PIALInfo{ii+1,3}))/2);
end

PIALInfo_BR_row = cell(3,length(fs_BR_idx));
for ci = 1: length(fs_BR_idx)
    PIALInfo_BR_row(:,ci) = PIALInfo_BR(ci,:);
end
fileID = fopen([subj '.PIAL_BR'],'w');
fprintf(fileID,'%2s\n','This is for the SEEG EI/functional network Project');
fprintf(fileID,'%2s\n','Chao Zhang');
fprintf(fileID,'%2s %2s %2s\n',PIALInfo_BR_row{:});
fclose(fileID);
% 
% XFname=fullfile(subDir,'elec_recon_BR',[subj '.PIAL_BR']);
% XInfo=csv2Cell(XFname,' ');
%% save .PIALVOX
PIALVOXFname=fullfile(subDir,'elec_recon_BR',[subj '.PIALVOX']);
PIALVOXInfo=csv2Cell(PIALVOXFname,' ',2);

PIALVOXInfo = PIALVOXInfo(izero,:);

PIALVOXInfo_BR = cell(length(fs_BR_idx),3);
for i = 1:length(fs_BR_idx)
    ii = fs_BR_idx(i);
    PIALVOXInfo_BR{i,1} = num2str((str2double(PIALVOXInfo{ii,1})+str2double(PIALVOXInfo{ii+1,1}))/2);
    PIALVOXInfo_BR{i,2} = num2str((str2double(PIALVOXInfo{ii,2})+str2double(PIALVOXInfo{ii+1,2}))/2);
    PIALVOXInfo_BR{i,3} = num2str((str2double(PIALVOXInfo{ii,3})+str2double(PIALVOXInfo{ii+1,3}))/2);
end

PIALVOXInfo_BR_row = cell(3,length(fs_BR_idx));
for ci = 1: length(fs_BR_idx)
    PIALVOXInfo_BR_row(:,ci) = PIALVOXInfo_BR(ci,:);
end
fileID = fopen([subj '.PIALVOX_BR'],'w');
fprintf(fileID,'%2s\n','This is for the SEEG EI/functional network Project');
fprintf(fileID,'%2s\n','Chao Zhang');
fprintf(fileID,'%2s %2s %2s\n',PIALVOXInfo_BR_row{:});
fclose(fileID);

% XFname=fullfile(subDir,'elec_recon_BR',[subj '.PIALVOX_BR']);
% XInfo=csv2Cell(XFname,' ');

%% save PostimpLoc.txt  undone
% PostimpLocFname=fullfile(subDir,'elec_recon_BR',[subj 'PostimpLoc.txt']);
% PostimpLocInfo=csv2Cell(PostimpLocFname,' ');
% 
% PostimpLocInfo = PostimpLocInfo(izero,:);
% 
% PostimpLocInfo_BR = cell(length(fs_BR_idx),3);
% for i = 1:length(fs_BR_idx)
%     ii = fs_BR_idx(i);
%     elecInfo_BR{i,1} = [elecInfo{ii,1},'-',elecInfo{ii+1,1}];
%     elecInfo_BR{i,2} = elecInfo{ii,2};
%     elecInfo_BR{i,3} = elecInfo{ii,3};
%     PIALVOXInfo_BR{i,1} = num2str((str2double(PIALVOXInfo{ii,1})+str2double(PIALVOXInfo{ii+1,1}))/2);
%     PIALVOXInfo_BR{i,2} = num2str((str2double(PIALVOXInfo{ii,2})+str2double(PIALVOXInfo{ii+1,2}))/2);
%     PIALVOXInfo_BR{i,3} = num2str((str2double(PIALVOXInfo{ii,3})+str2double(PIALVOXInfo{ii+1,3}))/2);
% end
% 
% %% save .mgrid  undone
% ddddd
% 
% %%  not sure if this part is useable
% 
% 
% elecInfo_table = table(elecInfo(:,1), elecInfo(:,2), elecInfo(:,3));
% elecInfo_table.Properties.VariableNames = {'Name', 'Depth_Strip_Grid', 'Hem'};
% elecNames=elecInfo(:,1);
% nElec=size(elecInfo,1);
% isLeft=zeros(nElec,1);
% isSubdural=zeros(nElec,1);
% for a=1:nElec
%     if ~strcmpi(elecInfo{a,2},'D')
%         isSubdural(a)=1;
%     end
%     if strcmpi(elecInfo{a,3},'L')
%         isLeft(a)=1;
%     end
% end


%% Cortex and channel label correction
subjVar_BR = [];
cortex = getcort(dirs);  %
% native_coord = importCoordsFreesurfer(dirs);
% fs_chan_names = importElectNames(dirs);
[MNI_coord, chanInfo, RAS_coord] = sub2AvgBrain_Chao_BR([],dirs, sbj_name, dirs.fsDir_local);

% [MGRID_coord, elect_names] = getmgrid(dirs);%elect_names pre with L/R Chao

 %the following is BR-wise MGRID output
MGRID_coord = zeros(length(fs_BR_idx),3);
for i = 1:length(fs_BR_idx)
    MGRID_coord(i,1) = str2double(PIALVOXInfo_BR{i,1});
    MGRID_coord(i,2) = str2double(PIALVOXInfo_BR{i,2});
    MGRID_coord(i,3) = str2double(PIALVOXInfo_BR{i,3});
end


elect_names = cell(1,length(fs_BR_idx));
for i = 1:length(fs_BR_idx)
   elect_names{1,i} = append (elecInfo_BR{i,3},elecInfo_BR{i,2},elecInfo_BR{i,1});
end

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

subjVar_BR.sbj_name = sbj_name;
subjVar_BR.cortex = cortex;
subjVar_BR.V = V;


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
ppt_chan_names = globalVar.channame_BR;%this is where I changed Chao
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
    disp('More channels in freesurfer, keep the channle names in googlesheet/EDF,PLZ manully change the code')
    fs_chan_names = fs_chan_names(in_chan_cmp);
    RAS_coord = RAS_coord(in_chan_cmp,:);
    MNI_coord = MNI_coord(in_chan_cmp,:);
    MGRID_coord = MGRID_coord(in_chan_cmp,:);%remove fs chan which was not shown in googlesheet
    % 2: More channels in EDF/TDT
elseif sum(in_chan_cmp) == length(in_chan_cmp) && sum(in_fs) < length(in_fs)
    disp('More channels in googlesheet/EDF, keep the freesurfer channle names,PLZ manully change the code')
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
    disp('channels in EDF/TDT which are not in FS,PLZ manully change the code')
    chan_comp(in_fs == 0)
    disp('channels in FS which are not in EDF/TDT,PLZ manully change the code')
    fs_chan_names(in_chan_cmp == 0)
    warning('this exception is not automatically fixable, please decide:PLZ manully change the code')
    
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
    
    subjVar_BR.LEPTO_coord = RAS_coord(new_order,:); %order sequency based on the googlesheet/EDF
    subjVar_BR.MNI_coord = MNI_coord(new_order,:);
    subjVar_BR.MGRID_coord = MGRID_coord(new_order,:);
    
    % labels mean the corrected names
    if strcmp(data_format, 'TDT')
        subjVar_BR.labels = chan_comp;
        subjVar_BR.labels_EDF = [];
    else
        %     subjVar_BR.elect_names = chan_comp;
        if isfield(globalVar, 'channame')
            subjVar_BR.labels_EDF = globalVar.channame_BR;
            subjVar_BR.labels = chan_comp;
        else
            subjVar_BR.labels = chan_comp;
            subjVar_BR.labels_EDF = [];
        end
    end
    
    
    %% Demographics % chao changes here
    subjVar_BR.demographics = [];
    
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
    subjVar_BR = localizeElect_Chao_BR(subjVar_BR,dirs,chanInfo);
    
    %% Save subjVar
    if exist([dirs.original_data filesep sbj_name filesep 'subjVar_BR_' sbj_name '.mat'], 'file')
        %         prompt = ['subjVar already exist for ' sbj_name ' . Replace it? (y or n):'] ;
        %         ID = input(prompt,'s');
        ID = 'y';
        if strcmp(ID, 'y')
            save([dirs.original_data filesep sbj_name filesep 'subjVar_BR_' sbj_name '.mat'], 'subjVar_BR')
            disp(['subjVar_BR saved for ' sbj_name])
            subjVar_created = 1;
        else
            warning(['subjVar_BR NOT saved for ' sbj_name])
        end
    else
        save([dirs.original_data filesep sbj_name filesep 'subjVar_BR_' sbj_name '.mat'], 'subjVar_BR')
        disp(['subjVar_BR saved for ' sbj_name])
        subjVar_created = 1;
    end
    
else
    subjVar_created = 0;
end

end

