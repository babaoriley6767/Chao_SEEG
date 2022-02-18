function [mgrid_coord, elect_names] = getmgrid(dirs)

subDir=dirs.freesurfer;
subj = dir(dirs.freesurfer);
subj = subj(end).name;
subDir = [subDir subj];

[mgrid_coord, elect_names] = mgrid2matlab_custom(subj,0, dirs.freesurfer);

end

