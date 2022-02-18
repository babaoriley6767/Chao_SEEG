function V = importVolumes(dirs)
%% Get folder name
subDir=dirs.freesurfer;
subj = dir(subDir);
subj = subj(end).name;
subDir = [subDir subj];

% Unzip nifti in case they are zipped
dirs.freesurfer

gunzip(fullfile(subDir,'elec_recon', 'brainmask.nii.gz'))
V = niftiread(fullfile(subDir,'elec_recon', 'brainmask.nii'));

end