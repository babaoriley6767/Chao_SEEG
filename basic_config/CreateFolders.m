function CreateFolders(sbj_name, project_name, block_name, center, dirs, import_server)
%% Create folders

% Get generic name  to match the server
sbj_name_generic = sbj_name;

if import_server
    all_folders = dir(fullfile('/Volumes/workstation/EPNETWORK/server/'));%/Volumes/CHAO_IRON_M/data_SEEG/server/
    for i = 1:length(all_folders)
        tpm(i) = contains(all_folders(i).name, sbj_name_generic);
    end
    sbj_folder_name = all_folders(find(tpm == 1)).name;
end

%% more laayer under neuraldata

folder_sublayers={'SpecData', 'BandData'};%need to find out if I need these

for i = 1:length(folder_sublayers)
    foldersublayers.(folder_sublayers{i}) = sprintf('%s/%s/%s',dirs.data_root,folder_sublayers{i});
    if ~exist(foldersublayers.(folder_sublayers{i}))
        mkdir(foldersublayers.(folder_sublayers{i}));
    end
end


%% data per subject (original, CAR and results)

folder_names = {'originalData', 'CARData','BRData'};%need to find out how to adapt these 'CAR is common average reference'
% 'BR is bipolar reference'

for i = 1:length(folder_names)
    folders.(folder_names{i}) = sprintf('%s/%s/%s',dirs.data_root,folder_names{i},sbj_name);
end

% folders.psych_dir = sprintf('%s/%s',dirs.psych_root,sbj_name);
folders.result_dir = sprintf('%s/%s/%s',dirs.result_root,project_name,sbj_name);

fieldname_folders = fieldnames(folders);

for i = 1:length(fieldname_folders)
    if ~exist(folders.(fieldname_folders{i}))
        mkdir(folders.(fieldname_folders{i}));
    end
end

%% Check if the globalval.mat exist
globalfile= dir([folders.originalData,filesep,'global*.mat']);

for bn = 1:length(block_name)
   
    clear globalVar
    if ~isempty(globalfile)
        load([folders.originalData,filesep,globalfile(1).name])
    end
    
    %% Per block - create folders and globalVar
    globalVar.block_name = block_name{bn};
    globalVar.sbj_name = sbj_name;
    globalVar.project_name = project_name;
    globalVar.center = center;
    %%
    for i = 1:length(folder_sublayers)
        globalVar.(folder_sublayers{i})=foldersublayers.(folder_sublayers{i});
    end
    
    %%
    for i = 1:length(fieldname_folders)
        globalVar.(fieldname_folders{i}) = [folders.(fieldname_folders{i}) '/' block_name{bn}];
        if ~exist(globalVar.(fieldname_folders{i}))
            mkdir(globalVar.(fieldname_folders{i}));
        end
        if strcmp(fieldname_folders{i}, 'result_dir') || strcmp(fieldname_folders{i}, 'originalData')
        else
            if ~exist([globalVar.(fieldname_folders{i}) ])
                mkdir([globalVar.(fieldname_folders{i}) ]);
            end
        end
    end
    %% Original folders from the server
    % iEEG data
    if import_server
            waitfor(msgbox(['Choose server file for iEEG data of block ' block_name{bn}]));
            [FILENAME, PATHNAME] = uigetfile(['/Volumes/workstation/EPNETWORK/server/' sbj_folder_name,'.edf'],'All Files (*.*)','MultiSelect','on');%/Volumes/CHAO_IRON_M/data_SEEG/server/
            globalVar.iEEG_data_server_path = [PATHNAME, FILENAME];
    % Save globalVariable
    fn = [folders.originalData '/' sprintf('global_%s_%s_%s.mat',project_name,sbj_name,block_name{bn})];
    save(fn,'globalVar');
end

end
