%% Branch 1 basic config
%set roots
addpath(genpath('/Users/tony/Documents/Stanford/Chao_SEEG/'))
addpath(genpath('/Users/tony/Documents/Stanford/DELLO'))
addpath(genpath('/Users/tony/Desktop/function_tools/iELVis/'))
addpath(genpath('/Users/tony/Desktop/function_tools/for_plot/CCEP_fMRI/'))
addpath('/Users/tony/Desktop/function_tools/spm12_7219')
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%home

%Initialize Directories
project_name = 'SEEG_test';

% Retrieve subject information
sbj_name = '2017_FT_SP_018_xiangxiang';
sbj_name = '2019_TT_SP_001_houwei';
sbj_name = 'group1_SJY';
sbj_name = 'group1_WXL';
sbj_name = 'group1_WYH';
sbj_name = 'group1_XX';
sbj_name = 'group2_CBB';
sbj_name = 'group2_YZZ';
sbj_name = 'group2_ZS'
sbj_name = 'test_bowen';

% center
if strcmp(sbj_name(6:7),'FT')
    center = 'fengtai';
else
    center = 'tiantan';
end

% Get block names
block_names = BlockBySubj(sbj_name,project_name)

% take care of the directory 
dirs = InitializeDirs(project_name, sbj_name, comp_root, server_root, code_root); 

% Get iEEG sampling rate and data format
[fs_iEEG, data_format] = GetFSdataFormat(center);

% Create subject folders and global_var
CreateFolders(sbj_name, project_name, block_names, center, dirs, 1);

% Copy the iEEG and behavioral files from server to local folders
parfor i = 1:length(block_names)
    CopyFilesServer(sbj_name,project_name,block_names{i},dirs) % CZ_important here
end

%% Branch 2 - data conversion(rawdata extraction, naming, decimate etc.)
ref_chan = [];
epi_chan = [];
empty_chan = [];

if strcmp(data_format, 'edf')
    SaveDataNihonKohden_sortData(sbj_name, project_name, block_names, dirs, ref_chan, epi_chan, empty_chan) 
else
    error('Data format has to be edf ')
end
%% Branch 3 - notch filtering & bipolar rereference

notch_and_BR(sbj_name, project_name, block_names, dirs)


%% localization method1

[subjVar_volume, subjVar_volume_created]  = CreateSubjVar_volume(sbj_name, comp_root, server_root, code_root, center);


remark = true; % Manually determine whether the electrode belongs to surface or depth

[subjVar, subjVar_created]  = CreateSubjVar_SEEG(sbj_name, comp_root, server_root, code_root,center,remark);

[subjVar_BR, subjVar_created]  = CreateSubjVar_SEEG_BR(sbj_name, comp_root, server_root, code_root,center,remark);

%% localization method2

[subjVar_volume_BR, subjVar_volume_created]  = CreateSubjVar_volume_BR(sbj_name, comp_root, server_root, code_root, center);



%% CCEP
%% H2
%% EI
%Step 1 
ER_pipeline_Chao(sbj_name, project_name, block_names, dirs , 'Band', 'EI') 
%Step 2
Nu = 1.2;
Lamda = 85;
Eta = 2;
EI_window = [65 110];%sjy;[50 150];%seconds
EI_pipeline_Chao(sbj_name, project_name, block_names(1), dirs, 'Band', 'EI',Nu,Lamda,Eta,EI_window)%test

Nu = 1.0;
Lamda = 70;
Eta = 4;
EI_window = [65 110];%sjy;[50 150];%seconds
EI_pipeline_Chao(sbj_name, project_name, block_names{2}, dirs, 'Band', 'EI',Nu,Lamda,Eta,EI_window)

%% Group analysis of the subjVar_volume_BR.eleinfo

y = zeros(length(0:0.2:20));
for i = 1:length(0:0.2:20)
        y(i) = 1/i;
end

figure
plot([0:0.2:20],y);
        