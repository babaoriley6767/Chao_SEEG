[n t s] = xlsread('test.xlsx');
Pattern = '-';
s(:,1) = Remove_Name_Pattern(s(:,1),Pattern);
for i =1:length(s)-1
    ss(i).EI = s{i+1,3};
    ss(i).chan = s{i+1,1};
    if ss(i).EI>0.3
        ss(i).tag = 1;
    end
end
    
    % denfine the ss.tag by EI
    
    for i =1:length(ss)
        temp = find(ismember(chanlabels_new,ss(i).chan));
        Channels(temp).EI = ss(i).EI;
        Channels(temp).tag = ss(i).tag;
    end
    
    num_EZ = 0; num_PZ = 0; num_NIZ = 0 ;EI =[];     %tag 1:EZ 2:PZ 3:NIZ
    for i =1:length(Channels)
        if Channels(i).tag == 1
            num_EZ = num_EZ+1;
            EI(1,num_EZ) = i;
        end
    end
    for i =1:length(Channels)
        if Channels(i).tag == 2
            num_PZ = num_PZ+1;
            EI(2,num_PZ) = i;
        end
    end
    for i =1:length(Channels)
        if Channels(i).tag == 3
            num_NIZ = num_NIZ+1;
            EI(3,num_NIZ) = i;
        end
    end
    len(1) = num_EZ;
    len(2) = num_PZ;
    len(3) = num_NIZ;
    for i =1:3
        for j =1:3
            temp = 0;
            temp_n = 0;
            for k =1:len(i)
                for l =1:len(j)
                    temp_s1 = chanlabels_new{EI(i,k)};
                    temp_s2 = chanlabels_new{EI(j,l)};
                    if strncmp(temp_s1,temp_s2,1) ~=1
                        temp_n =temp_n+1;                        
                        temp = temp+M(EI(i,k),EI(j,l));
                    end                    
                end
            end
            tempss(i,j) =temp_n;
            N(i,j) = temp/temp_n;
        end
    end

    
    save('PTXXX_EI.mat','EI');
    save('PTXXX_N.mat','N');
    save('PTXXX_ss.mat','ss');
    save('PTXXX_Channels.mat','Channels');
    
    