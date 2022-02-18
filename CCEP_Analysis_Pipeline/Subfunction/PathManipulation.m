function  [ScriptFolder,SubjectFolder] = PathManipulation(User)
%PATHMANIPULATION Summary of this function goes here
%   Detailed explanation goes here
switch User
    case 'BaotianMacmini'
        % Add SPM
        % addpath('C:\Users\THIENC\Desktop\spm12_7219');
        addpath('/Users/Baotian/Desktop/ToolBoxes/spm12');
        spm('defaults', 'EEG');
        
        % Add path
        % addpath(genpath('C:\Users\THIENC\Desktop\EmotionalFaces'))
        addpath(genpath('/Users/Baotian/Desktop/ToolBoxes/EmotionalFaces'))
        ScriptFolder = '/Users/Baotian/Desktop/ToolBoxes/EmotionalFaces';
        SubjectFolder = '/Users/Baotian/Desktop/Data/EmotionalFaces/PT050';
    case 'BaotianZ620'
        % Add SPM
        % addpath('C:\Users\THIENC\Desktop\spm12_7219');
        addpath('D:\spm12_7219');
        spm('defaults', 'EEG');
        
        % Add path
        % addpath(genpath('C:\Users\THIENC\Desktop\EmotionalFaces'))
        addpath(genpath('D:\CCEP_Pipeline'))
        ScriptFolder = 'D:\CCEP_Pipeline';
        SubjectFolder = 'E:\';
end

end

