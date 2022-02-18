%% load data
load('Data.mat')

%% Task one 
figure
plot(Data.Time,Data.EEG(3,:));
axis tight
grid on 
xlabel('Time(s)')
ylabel('Amp(uA)')
title('DC Channel Triggers')

%% Task two
% Bipolar rereference
BipolarRef = Data.EEG(1,:) - Data.EEG(2,:);
figure
subplot(3,1,1)
plot(Data.Time(2001:5000),Data.EEG(1,2001:5000) - mean(Data.EEG(1,2001:5000)))
title(Data.ChannelLabel(1))
FigLegend
subplot(3,1,2)
plot(Data.Time(2001:5000),Data.EEG(2,2001:5000) - mean(Data.EEG(2,2001:5000)))
title(Data.ChannelLabel(2))
FigLegend
subplot(3,1,3)
plot(Data.Time(2001:5000),BipolarRef(2001:5000) - mean(BipolarRef(2001:5000)))
title('Bipolar A1-A2')
FigLegend
%% Task three
Delta     = [4];
Theta     = [4 7];
Alpha     = [8 15];
Beta      = [16 31];
Gamma     = [33 70];
HighGamma = [70 170];
% filter data using Butterworth 5 order section filter.
DeltaTrace     = ButterworthFilter(BipolarRef,'lowpassiir',10,Delta,Data.SamplingRate);
ThetaTrace     = ButterworthFilter(BipolarRef,'bandpassiir',10,Theta,Data.SamplingRate);
AlphaTrace     = ButterworthFilter(BipolarRef,'bandpassiir',10,Alpha,Data.SamplingRate);
BetaTrace      = ButterworthFilter(BipolarRef,'bandpassiir',10,Beta,Data.SamplingRate);
GammaTrace     = ButterworthFilter(BipolarRef,'bandpassiir',10,Gamma,Data.SamplingRate);
HighGammaTrace = ButterworthFilter(BipolarRef,'bandpassiir',10,HighGamma,Data.SamplingRate);

figure
subplot(6,1,1)
plot(Data.Time(2001:5000),DeltaTrace(2001:5000) - mean(DeltaTrace(2001:5000)))
hold on
plot(Data.Time(2001:5000),BipolarRef(2001:5000) - mean(BipolarRef(2001:5000)))
FigLegend
title('DeltaTrace')
subplot(6,1,2)
plot(Data.Time(2001:5000),ThetaTrace(2001:5000) - mean(ThetaTrace(2001:5000)))
hold on
plot(Data.Time(2001:5000),BipolarRef(2001:5000) - mean(BipolarRef(2001:5000)))
FigLegend
title('ThetaTrace')
subplot(6,1,3)
plot(Data.Time(2001:5000),AlphaTrace(2001:5000) - mean(ThetaTrace(2001:5000)))
hold on
plot(Data.Time(2001:5000),BipolarRef(2001:5000) - mean(BipolarRef(2001:5000)))
FigLegend
title('AlphaTrace')
subplot(6,1,4)
plot(Data.Time(2001:5000),BetaTrace(2001:5000) - mean(BetaTrace(2001:5000)))
hold on
plot(Data.Time(2001:5000),BipolarRef(2001:5000) - mean(BipolarRef(2001:5000)))
FigLegend
title('BetaTrace')
subplot(6,1,5)
plot(Data.Time(2001:5000),GammaTrace(2001:5000) - mean(GammaTrace(2001:5000)))
hold on
plot(Data.Time(2001:5000),BipolarRef(2001:5000) - mean(BipolarRef(2001:5000)))
FigLegend
title('GammaTrace')
subplot(6,1,6)
plot(Data.Time(2001:5000),HighGammaTrace(2001:5000) - mean(HighGammaTrace(2001:5000)))
hold on
plot(Data.Time(2001:5000),BipolarRef(2001:5000) - mean(BipolarRef(2001:5000)))
FigLegend
title('HighGammaTrace')
%% epoch data
timeStameResults = GetTimeStamp(Data.EEG(3,:));
% visually theck the onsets from DC channel
figure
plot(Data.EEG(3,:))
hold on
for i = 1:length(timeStameResults)
    plot(timeStameResults(i),0,'ro')
end
timeStameResults = timeStameResults(2:end);
EpochedData = zeros(length(timeStameResults),2201);
for i = 1:length(timeStameResults)
    EpochedData(i,:) = BipolarRef(timeStameResults(i)-600:timeStameResults(i)+1600) - ...
        mean(BipolarRef(timeStameResults(i)-600:timeStameResults(i)+1600));
end

figure
plot(EpochedData(:,:)')
xticks(0:200:2201)
xticklabels(-0.6:0.2:1.6)
FigLegend
title('Epoched Pile Plot')
%% Time frequency transform
% TF transform take some time
FourierValues = zeros(length(timeStameResults),501,2074);
for i = 1:length(timeStameResults)
[FourierValues(i,:,:),f,t] = spectrogram(EpochedData(i,:),128,127,1000,Data.SamplingRate,'yaxis');
end
FourierValuesCut = FourierValues(:,1:201,37:2037);

figure
for i = 1:6
subplot(2,3,i)
imagesc(abs(squeeze(FourierValuesCut(i,:,:))))
set(gca,'YDir','normal')
ylabel('Hz')
xlabel('Time(s)')
title(['Trial','-',num2str(i)])
xticks(0:500:2001)
xticklabels(-0.5:0.5:1.5)
end

figure
imagesc(squeeze(mean(abs(FourierValuesCut(:,:,:)))))
set(gca,'YDir','normal')
ylabel('Hz')
xlabel('Time(s)')
title('AllTrialMean')
xticks(0:500:2001)
xticklabels(-0.5:0.5:1.5)
%% Baseline normaliziton
% Z-score normalizaiton for each trial
Z_scoreTrials = zeros(size(FourierValuesCut));
for i = 1:size(FourierValuesCut,1)
    Z_scoreTrials(i,:,:) = abs(squeeze(FourierValuesCut(i,:,:)))./std(abs(squeeze(FourierValuesCut(i,:,1:500)))')';
end
figure
imagesc(squeeze(mean(Z_scoreTrials)))
set(gca,'YDir','normal')
ylabel('Hz')
xlabel('Time(s)')
title('AllTrialMean')
xticks(0:500:2001)
xticklabels(-0.5:0.5:1.5)
colormap('jet')
colorbar

%% High Gamma traces
% Mean the high gamma power
HighGamma = squeeze(mean(Z_scoreTrials));
HighGamma = mean(HighGamma(70:170,:)) - mean(mean(HighGamma(70:170,1:500)));
figure
plot(HighGamma)
ylabel('Power(Z-score)')
xlabel('Time(s)')
title('HighGamma')
xticks(0:500:2001)
xticklabels(-0.5:0.5:1.5)
axis tight
grid on

%% Subfunctions
function FigLegend
axis tight
grid on
xlabel('Time(s)')
ylabel('Amp(uA)')
end

function FilteredSignal = ButterworthFilter(Sig,Type,Order,CutoffFrequency,SamplingRate)
if length(CutoffFrequency) == 1
    Filter = designfilt(Type,'FilterOrder',Order, ...
        'HalfPowerFrequency',CutoffFrequency,'DesignMethod','butter','SampleRate',SamplingRate);
    FilteredSignal = filtfilt(Filter,Sig);
elseif length(CutoffFrequency) == 2
    Filter = designfilt(Type,'FilterOrder',Order, ...
        'HalfPowerFrequency1',CutoffFrequency(1),'HalfPowerFrequency2',CutoffFrequency(2),...
        'DesignMethod','butter','SampleRate',SamplingRate);
    FilteredSignal = filtfilt(Filter,Sig);
end
end

function timeStameResults = GetTimeStamp(TriggerChannel)
% Binarize
maxTrigger = max(TriggerChannel);
threshold = 0.3*maxTrigger;

BinTrigger = (TriggerChannel >= threshold);
timeStamps = diff(BinTrigger);
timeStameResults = find(timeStamps == 1);
end

