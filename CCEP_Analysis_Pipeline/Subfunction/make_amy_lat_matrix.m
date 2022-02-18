%to calc the amplitude and latency of CCEP N1
CCEPdir = dir();
CCEPdir = CCEPdir(3:end);
for i = 1:length(CCEPdir)
    cd(CCEPdir(i).name);
    edfFiles = dir('*.mat');
    temp = 0;max_num= 0 ;
    for j = 1:length(edfFiles)
        if length(edfFiles(j).name)>temp
            temp = length(edfFiles(j).name);
            max_num = j;
        end
    end
    edfFile = edfFiles(max_num);
    D = spm_eeg_load(edfFile.name);
    Pattern = '-';
    chanlabels_new = D.chanlabels;
    chanlabels_new = Remove_Name_Pattern(chanlabels_new,Pattern);
    Evoked = find(ismember(chanlabels_new,CCEPdir(i).name));
    for j= 1:length(D.chanlabels)
        datas = squeeze(D(j,:,:));
        data = mean(datas,2);
        [Amp(Evoked,j),Latency(Evoked,j)] = findmax(data,10,50);
    end
    cd ..
end

load('PTXXX_Channels.mat');
num = 0;
for i =1:length(Channels)
    if Channels(i).keep == 1 
        num = num+1;
        Chan_need(num) = i;
    end
end
for i = 1:length(Chan_need)
    for j = 1:length(Chan_need)
      Amp_grey(i,j) = Amp(Chan_need(i),Chan_need(j));
      Latency_grey(i,j) = Latency(Chan_need(i),Chan_need(j));
    end
end


load('PTXXX_ss.mat');
    for i =1:length(ss)
        temp = find(ismember(chanlabels_new,ss(i).chan));
        Channels(temp).tag = ss(i).tag;
    end
    
    num_EZ = 0; num_PZ = 0; num_NIZ = 0 ;     %tag 1:EZ 2:PZ 3:NIZ
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
            tempa = 0;
            templ = 0;
            temp_n = 0;
            for k =1:len(i)
                for l =1:len(j)
                    temp_s1 = chanlabels_new{EI(i,k)};
                    temp_s2 = chanlabels_new{EI(j,l)};
                    if strncmp(temp_s1,temp_s2,1) ~=1
                        temp_n =temp_n+1;                        
                        tempa = tempa+Amp(EI(i,k),EI(j,l));
                        templ = templ+Latency(EI(i,k),EI(j,l));
                    end                    
                end
            end
            tempss(i,j) =temp_n;
            NA(i,j) = tempa/temp_n;
            NL(i,j) = templ/temp_n;
        end
    end
    