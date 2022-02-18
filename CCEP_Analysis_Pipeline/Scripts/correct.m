spm;
D = spm_eeg_load;
D.chanlabels
figure
plot(D(11,10000:200000))
Hippocampus = D(58:59,10000:200000);
RefHippo = Hippocampus(1,:) - Hippocampus(2,:);

figure;
plot(RefHippo)

% find time trigger
TriggerChannel = D(39,10000:200000);

% Binarize
maxTrigger = max(TriggerChannel);
threshold = 0.3*maxTrigger;

BinTrigger = (TriggerChannel >= threshold);
timeStamps = diff(BinTrigger);

timeStameNew = find(timeStamps == 1);


EpochedData = zeros(40,201);
for i = 1:40
    EpochedData(i,:) = RefHippo(timeStameNew(i)-100:timeStameNew(i)+100);
end

% for i = 1:40
%     figure
%     plot(EpochedData(i,:))
%     hold on
%     % remove artifact
%     a = EpochedData(i,:);
%     EpochedData_New(i,:) = a;
%     t1 = a(91:101);t2 = a(111:121);
%     r1 = flip(t1);r2 = flip(t2);
%     x1 = 1;x2 = 0;tt=1/10;
%       for j=1:11 
%           t1(j) = r1(j)*x1;
%           t2(j) = r2(j)*x2;
%           x1 = x1-tt;
%           
%           x2 = x2+tt;
%       end
%     y = t1+t2;
%     EpochedData_New(i,101:111) = y;
%     plot(EpochedData_New(i,:));
%     axis tight
%     grid on
% end

for i = 1:40
    a = EpochedData(i,:);
    EpochedData_New(i,:) = a;
    t1 = a(91:100);t2 = a(111:120);
    r1 = flip(t1);r2 = flip(t2);
    x1 = 1;x2 = 0;tt=1/10;
      for j=1:10 
          t1(j) = r1(j)*x1;
          t2(j) = r2(j)*x2;
          x1 = x1-tt;          
          x2 = x2+tt;
      end
    y = t1+t2;
    EpochedData_New(i,100:109) = y;
end          %È¥´Ì¼¤Î±¼£



figure
plot(EpochedData(:,:)')
axis tight

figure
plot(EpochedData_New(:,:)')
axis tight

        
figure
plot(mean(EpochedData(:,:)',2))
axis tight


% cc = zeros(1,201)
% for i =1:40
%     if EpochedData(i,101) > 200 
%        EpochedData(i,:)=cc;
%     end
% end