function CCEPDCTriggers = FindCCEPTriggers(DC)
% Function to extract CCEP triggers from DC channels
% DC triggers timestamps are in seconds, the onset is zero
% Baotian Zhao 20190404 @ Beijing


if nargin < 1
    DC = spm_eeg_load();
end

maxTrigger = max(DC(:));
threshold = 0.3*maxTrigger;
BinTrigger = (DC(:) >= threshold);
timeStamps = diff(BinTrigger);
timeStampNew = find(timeStamps == 1);

CCEPDCTriggers = timeStampNew/DC.fsample;

% plot EEG triggers
figure
hold on
for i = 1:length(CCEPDCTriggers)
    line([CCEPDCTriggers(i) CCEPDCTriggers(i)],[-500 500],'Color','r','LineWidth',1)
end
axis tight
grid on
title(['Trigger Num = ' num2str(length(CCEPDCTriggers))])
print(['TriggerNum=' num2str(length(CCEPDCTriggers))],'-dpng')
end

