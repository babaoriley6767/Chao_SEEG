function [server_root, comp_root, code_root] = AddPaths(user)

if strcmp(user,'Chao_iMAC')
    code_root = sprintf('/Users/tony/Documents/Stanford/Chao_SEEG');% location of the main code
    comp_root = sprintf('/Users/tony/Desktop/data_SEEG/data/');% location of analysis_SEEG folder /Volumes/CHAO_IRON_M/data_SEEG/data/
else
end
server_root = '/Users/tony/Desktop/data_SEEG/server';% location of the raw data /Volumes/CHAO_IRON_M/data_SEEG/server
end
