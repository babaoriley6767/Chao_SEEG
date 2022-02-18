% Organize the generated preprocessed results
% Baotian Zhao 20190421
clear
cd('L:\')
SubjectsFolders = dir('PT*');
for i = 1:length(SubjectsFolders)
    cd(SubjectsFolders(i).name)
    Pairs = dir(pwd);
    Pairs = Pairs(3:end);
    for j = 1:length(Pairs)
        cd(Pairs(j).name)
        % Remove files not png
        delete *.edf
        delete *.mat
        delete *.dat
        cd ..
    end
    cd ..
end
        

