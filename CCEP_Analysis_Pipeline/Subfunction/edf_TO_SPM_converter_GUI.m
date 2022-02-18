function [D,IndexInUse] = edf_TO_SPM_converter_GUI(filename,index_to_keep,prefix)
% This function gives you a interactive way to select channels you want to
% keep. Convert edf file data to SPM object format
% Baotian @ Beijing 20180811
% 20190228 update Baotian
% 20190429 update by Baotian

% initialize spm EEG module
spm('defaults', 'EEG');

if nargin < 1 || isempty(filename)
    [filename, path] = uigetfile('*.edf','Please select edf file to convert');
    filename_full = [path filename];
end
if nargin < 2 || isempty(index_to_keep)
    % read edf header using function in field trip
    edf_header = ft_read_header(filename);
    [index_to_keep,~] = listdlg('ListString',edf_header.label,'PromptString','Select Channels to Keep' );
    IndexInUse = index_to_keep;
end
if nargin < 3 || isempty(prefix)
    % add prefix to the converted file, by default 'spmeeg_'
    prefix = 'spmeeg_';
end

IndexInUse = index_to_keep;
filename_full = filename;
edf_header = ft_read_header(filename_full);
% convert the data
S.dataset = filename_full;
S.mode  = 'continuous';
if strcmp(index_to_keep,'all')
    S.channels = 'all';
else
    S.channels = edf_header.label(index_to_keep);
end
S.outfile  = [prefix filename];
D = spm_eeg_convert(S);
end