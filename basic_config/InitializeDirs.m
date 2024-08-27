function dirs = InitializeDirs(project_name,sbj_name, comp_root, server_root, code_root)
%% Initialize directories

dirs.server_root = server_root; 
dirs.comp_root = comp_root;
dirs.code_root = code_root;
dirs.data_root = sprintf('%s/neuralData',dirs.comp_root);
dirs.result_root = sprintf('%s/Results',dirs.comp_root);
dirs.project = sprintf('%s/Results/%s',dirs.comp_root,project_name);
%dirs.psych_root = sprintf('%s/psychData',dirs.comp_root);%
% dirs.elec = sprintf('%s/ECoG Patient Info/Electrodes/Native_elecs',dirs.comp_root);
% dirs.mni_elec = sprintf('%s/ECoG Patient Info/Electrodes/MNI_elecs',dirs.comp_root);
% dirs.mni_cortex = sprintf('%s/ECoG Patient Info/Cortex/ColinCortex',dirs.comp_root);
% dirs.cortex = sprintf('%s/ECoG Patient Info/Cortex/Native_cortex',dirs.comp_root);
% dirs.ROI = sprintf('%s/ECoG Patient Info/ROIs',dirs.comp_root);
dirs.original_data = [dirs.data_root '/originalData'];
% dirs.MVPAData = [dirs.comp_root,'/MVPAData'];

% Define freesurfer local folder
comp = computer();
if contains(comp, 'MACI64') % General for MAC
    dirs.fsDir_local = '/Applications/freesurfer/subjects/fsaverage';
else % If not MAC, user should choose the folder
    dirs.fsDir_local = uigetdir();
end


%Set freesurfer folder
all_folders = dir(fullfile('/Volumes/workstation/EPNETWORK/server/'));% /Volumes/CHAO_IRON_M/data_SEEG/server/
if isempty(all_folders)
    warning('You are not connected to the server, therefore no Fressurfer folder will be specified.')
else
    for i = 1:length(all_folders)
        tpm(i) = contains(all_folders(i).name, sbj_name);
    end
    sbj_folder_name = all_folders(find(tpm == 1)).name;

    all_folders_sbj = dir(fullfile(['/Volumes/workstation/EPNETWORK/server/' sbj_folder_name]));%/Volumes/CHAO_IRON_M/data_SEEG/server/
    for i = 1:length(all_folders_sbj)
        tpm_2(i) = contains(all_folders_sbj(i).name, 'surfer');
    end
    if sum(tpm_2) == 0
        warning('There is no Freesurfer folder')
        dirs.freesurfer = [];
    else
        dirs.freesurfer = ['/Volumes/workstation/EPNETWORK/server/' sbj_folder_name '/' all_folders_sbj(tpm_2).name '/'];%/Volumes/CHAO_IRON_M/data_SEEG/server/
    end
end

%Set brainvisa folder
all_folders = dir(fullfile('/Volumes/workstation/EPNETWORK/server/'));% /Volumes/CHAO_IRON_M/data_SEEG/server/
if isempty(all_folders)
    warning('You are not connected to the server, therefore no Brainvisa folder will be specified.')
else
    for i = 1:length(all_folders)
        tpm(i) = contains(all_folders(i).name, sbj_name);
    end
    sbj_folder_name = all_folders(find(tpm == 1)).name;

    all_folders_sbj = dir(fullfile(['/Volumes/workstation/EPNETWORK/server/' sbj_folder_name]));%/Volumes/CHAO_IRON_M/data_SEEG/server/
    for i = 1:length(all_folders_sbj)
        tpm_2(i) = contains(all_folders_sbj(i).name, 'visa');
    end
    if sum(tpm_2) == 0
        warning('There is no Brianvisa folder')
        dirs.brainvisa = [];
    else
        dirs.brainvisa = ['/Volumes/workstation/EPNETWORK/server/' sbj_folder_name '/' all_folders_sbj(tpm_2).name '/'];%/Volumes/CHAO_IRON_M/data_SEEG/server/
    end
end 
    
end