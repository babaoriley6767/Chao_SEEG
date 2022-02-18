function CopyFilesServer(sbj_name,project_name, block_names,dirs)

% Load globalVar
fn = sprintf('%s/originalData/%s/global_%s_%s_%s.mat',dirs.data_root,sbj_name,project_name,sbj_name,block_names);
load(fn,'globalVar');


copyfile(globalVar.iEEG_data_server_path, globalVar.originalData)
display(sprintf('Copied EDF file %s to %s', globalVar.iEEG_data_server_path, globalVar.originalData));

end
