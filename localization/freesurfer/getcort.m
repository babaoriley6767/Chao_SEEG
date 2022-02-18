function cortex = getcort(dirs)

cortex = [];
cortex.left = [];
cortex.right = [];

subDir=dirs.freesurfer;
subj = dir(dirs.freesurfer);
subj = subj(end).name;
subDir = [subDir subj];

left_path = fullfile(subDir,'surf', 'lh.PIAL'); % paths to lh/rh.pial
right_path = fullfile(subDir,'surf','rh.PIAL'); % paths to lh/rh.pial
[cortex.left.vert, cortex.left.tri] = read_surf(left_path);
[cortex.right.vert, cortex.right.tri] = read_surf(right_path);

cortex.left = cortex.left';
cortex.right = cortex.right';
end


