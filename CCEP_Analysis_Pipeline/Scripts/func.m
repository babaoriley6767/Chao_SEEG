spm;
D = spm_eeg_load;                                    %读取
% figure
% plot(D(38,10000:180000))
Hippocampus = D(26:27,350000:550000);                %需要观察的两个channel，刺激的时间区间
RefHippo = Hippocampus(1,:) - Hippocampus(2,:);
% figure;
% plot(RefHippo)

% find time trigger
TriggerChannel = D(38,350000:550000);                % DC 10 找trigger

% Binarize
maxTrigger = max(TriggerChannel);
threshold = 0.3*maxTrigger;

BinTrigger = (TriggerChannel >= threshold);
timeStamps = diff(BinTrigger);

timeStameNew = find(timeStamps == 1);

EpochedData = zeros(40,401);
for i = 1:40
    EpochedData(i,:) = RefHippo(timeStameNew(i)-100:timeStameNew(i)+300);    %选取的观察时间区间
end


for i = 1:40
    a = EpochedData(i,:);
    EpochedData_New(i,:) = a;
    EpochedData_New(i,100:109) = remove_art(a);       %去除电刺激伪迹
end      

for i = 1:40
    a = EpochedData_New(i,:);                      %计算RMS
    RMS(i) = Calc_RMS(a);
end
yy=mean(RMS);
yy=median(RMS);
aa=mean(EpochedData_New(:,:));
yyy=Calc_RMS(aa);


figure
plot(EpochedData(:,:)')
axis tight

figure
plot(EpochedData_New(:,:)')
axis tight

figure
plot(mean(EpochedData_New(:,:)',2))
axis tight

