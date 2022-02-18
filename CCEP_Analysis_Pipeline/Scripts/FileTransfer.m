% Transfer the data 
PatientList = [2,3,5,8,9:16,38,39,41,47,48,49,51]; % for CCEP cohort
% Move the edf files to the analysis folder
for i = 1:5
    cd('T:\')
    MainDir = dir('*0*');
    for j = 1:length(MainDir)
        if contains(MainDir(j).name,['00' num2str(PatientList(i))])
            mkdir(['E:\CCEP\' 'PT00' num2str(PatientList(i))])
            cd(MainDir(j).name)
            CCEPDir = dir('CCEP_Raw*');
            cd(CCEPDir.name)
            copyfile('*.edf',['E:\CCEP\' 'PT00' num2str(PatientList(i))])
            cd ..
            cd ..
        end
    end
end

for i = 6:12
    cd('T:\')
    MainDir = dir('*0*');
    for j = 1:length(MainDir)
        if contains(MainDir(j).name,['0' num2str(PatientList(i))])
            mkdir(['E:\CCEP\' 'PT0' num2str(PatientList(i))])
            cd(MainDir(j).name)
            CCEPDir = dir('CCEP_Raw*');
            cd(CCEPDir.name)
            copyfile('*.edf',['E:\CCEP\' 'PT0' num2str(PatientList(i))])
            cd ..
            cd ..
        end
    end
end

for i = 13:19
    cd('T:\')
    MainDir = dir('*0*');
    for j = 1:length(MainDir)
        if contains(MainDir(j).name,['0' num2str(PatientList(i))])
            mkdir(['E:\CCEP\' 'PT0' num2str(PatientList(i))])
            cd(MainDir(j).name)
            ECoGDir = dir('*ECoG');
            cd(ECoGDir.name)
            CCEPDir = dir('CCEP_Raw*');
            cd(CCEPDir.name)
            copyfile('*.edf',['E:\CCEP\' 'PT0' num2str(PatientList(i))])
            cd ..
            cd ..
            cd ..
        end
    end
end
% % Move the mat files to the analysis folder
% for i = 1:length(PatientList)
%     cd('T:\')
%     MainDir = dir('*0*');
%     for j = 1:length(MainDir)
%         if contains(MainDir(j).name,['0' num2str(PatientList(i))])
%             cd(MainDir(j).name)
%             ECoGDir = dir('*Behavioral');
%             cd(ECoGDir.name)
%             EmotionalDir = dir('Emotional*');
%             cd(EmotionalDir.name)
%             copyfile('*.mat',['E:\EmotionalFaces\' 'PT0' num2str(PatientList(i))])
%             cd ..
%             cd ..
%             cd ..
%         end
%     end
% end