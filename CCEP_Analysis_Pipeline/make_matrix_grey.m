load('PTXXX_Matrix.mat');
load('PTXXX_chanlabels.mat');
for i = 1:length(chanlabels_new)
    Channels(i).id = chanlabels_new{i};
    Channels(i).keep = 0;
end

%define keep or not

num = 0;
for i =1:length(Channels)
    if Channels(i).keep == 1 
        num = num+1;
        Chan_need(num) = i;
    end
end
for i = 1:length(Chan_need)
    for j = 1:length(Chan_need)
      M_grey(i,j) = M(Chan_need(i),Chan_need(j));
    end
end
save('PTXXX_Matrix_grey.mat','M_grey');
save('PTXXX_Channels.mat','Channels');

