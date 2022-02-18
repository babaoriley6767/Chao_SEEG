function [fs_iEEG, data_format] = GetFSdataFormat(center)


if strcmp(center, 'fengtai')
    
     fs_iEEG = 2000;% there might be problem
    data_format = 'edf';

elseif strcmp(center, 'tiantan')
    
    fs_iEEG = 2000;% there might be problem
    data_format = 'edf';
    
end

end