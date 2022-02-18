function EI_pipeline(sbj_name, project_name, block_names, dirs) %unfinished 21_3_23 Chao
%% INPUTS
%   sbj_name:               subject name
%   project_name:           name of task
%   block_names:                     names of blocks to be analyed (cell of strings)
%   dirs:                   directories pointing to files of interest (generated by InitializeDirs)
for bi = 1:length(block_names)

%%
% Load globalVar
bn = block_names{bi};
    
%% Load globalVar
fn = sprintf('%s/originalData/%s/global_%s_%s_%s.mat',dirs.data_root,sbj_name,project_name,sbj_name,bn);
load(fn,'globalVar');

for el = 1:length(globalVar.channame_BR)

%need to figure out a output direction??? el wise? our upper logcial wise?
% dir_out = [data_root,freq_band,'Data',filesep,sbj_name,filesep,bn];c
dir_out = [data_root,freq_band,filesep,sbj_name,filesep,bn];
if ~exist(dir_out, 'dir')
    mkdir(dir_out)
end


%% Per electrode
%loda each BR data
load(sprintf('%s/BRiEEG%s_%.2d.mat',globalVar.BRData,bn,el));
data.label = globalVar.channame_BR{el};

%% step 1 calculate the ER
data_iEEG = data.wave;
fs = data.fsample; % frequency sampling is defined here
data_length = length(data.wave); % how many dots are in this data
data_time_length = fix(data_length / fs); % data time in seconds(s)
time_series = 0:1/fs:(data_length-1)/fs;
% defination of deferent frequency band.
% attention, this is slightly different from the article, part of the
% reason is we need to consider the frequency resoluton after FFT, by using 1000 fs,
% it is better to use integers only, however you can use zero padding to
% increase the frequency resoluton, but I think that is not necessary.
theta_band = [4:7];
alpha_band = [8:12];
beta_band = [13:24];    
gamma_band = [24:97];

% define sliding window duration D
D = 1; % time window to perform FFT in seconds (s), time window in article is probably 0.5 according to their map % this could be changed 
window_length = D * fs;%(points)
% define step length for the sliding window
step_length = 0.01; % step_length in seconds(s), 0.01 acturally means 10 dot with 1000 fs, step length affact time consuming seriously % this could be changed 

%% fft transform & frequency sub-band energy extraction & ER calculation

indx_matrix = buffer(1:data_length,round(window_length),round((1-step_length)*fs),'nodelay');
for i = 1:size(indx_matrix,2)-1
    data_in_window = data_iEEG(indx_matrix(:,i));
% for i = 1:(step_length * fs):(data_length - D * fs) % slide the window
%     data_in_window = data_iEEG(i:(window_length+i-1)); % extract time series data in one window %chao do we need minus 1????
    freq_resolution = fs / window_length; % define freq_resolution %chao ???啥
    X_data = periodogram(data_in_window); % I need to read more books
    freq_range = 0:(length(X_data)-1)*freq_resolution; % define the frequency axis>??????
    % extract energy of different band  %chao need to make sure the
    % frequency sampling and energy 
    E_theta = sum(X_data(theta_band + 1));% chao, there is engergy result of 0 Hz in periodogram, this is why +1
    E_alpha = sum(X_data(alpha_band + 1));
    E_beta = sum(X_data(beta_band + 1));
    E_gamma = sum(X_data(gamma_band + 1));
    % calculate time-varying ER
%     ER_number = fix(i/(step_length*fs))+1;
    ER(i) = (E_beta + E_gamma)/(E_theta + E_alpha); 
end

%% step 2: optimal detection of rapid discharges
% calculate the vector of ER_n
for i = 1:length(ER)
    ER_n(i) = mean(ER(1:i));
end
v = 0.5;%%%??
for i = 1:length(ER)
    UN_temp(i) = ER(i) - ER_n(i) - v;
    UN(i) = sum(UN_temp(1:i));
end

% find the local minimum of UN and determine when lamda is reached
lamda = 300;%%%??
% calculate local minimum

u_N = zeros(1,length(ER));
for i = 1:length(UN)
    u_N(i) = min(UN(1:i));
end
% calculate threshold
N_a = [];
for i = 1:length(UN)
    if UN(i) - u_N(i) >= lamda
        N_a = (i - 1) * step_length; % alarm time in seconds
        break
    end
end

% calculate N_d, Another particular case
% is when structure Si does not generate fast activity. In this case,
% Ndi was set to be equal to the time corresponding to seizure
% termination. This operation leads to very low value of EIi, as
% expected.
if isempty(N_a)
    N_d = length(UN);% I need to adjust this value
else
    for i = 1:length(UN)
        if u_N(i) == u_N(round(N_a/step_length))
            N_d = (i - 1) * step_length; % detection time in seconds
        break
        end
    end
end

Nd = N_d;
% ER = ER;









for i = 1:EEG.Channels
    [Nd(i),ER(i,:)] = EI_calc(EEG,EEG.data(i,:),i);
end

% calculate N_0 and EI for each channel
N_0 = min(Nd);
H = 5;
for i = 1:EEG.Channels
    ER_temp = ER(i,:);
       if Nd(i) <= length(ER_temp)/100 - H % acturally it is H * fs / step_length * fs
           EI(i) = (1 / (Nd(i) - N_0 + 1)) * sum(ER_temp(round(Nd(i)*100):(round(Nd(i)*100) + H*100)));
       else
           EI(i) = (1 / (Nd(i) - N_0 + 1));
       end
end
% normalization EI index
for i = 1:EEG.Channels
    EI(i) = EI(i)/max(EI);
end



%% calculate time-frequency parameters
EEG.times = time_series;
EEG.srate = fs;
EEG.pnts = data_length;

% definitions, selections...

min_freq =  2;
max_freq = 90;
num_frex = 30;

% define wavelet parameters
time = -1:1/EEG.srate:1;
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s    =  3./(2*pi*frex); % this line is for figure 13.14
% s    = 10./(2*pi*frex); % this line is for figure 13.14

% definte convolution parameters
n_wavelet            = length(time);
n_data               = EEG.pnts;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% get FFT of data
eegfft = fft(EEG_channel_data,n_conv_pow2);

% initialize
eegpower = zeros(num_frex,EEG.pnts); % frequencies X time X trials


% loop through frequencies and compute synchronization
for fi=1:num_frex
    
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    eegpower(fi,:) = 10*log10(abs(eegconv).^2);
end

%% plot figure of each channel
% figure('color','w');
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)],'color','w')
subplot(4,1,1);
plot(t,EEG_channel_data);
ylabel('EEG');
axis tight
grid on
hold on 
plot([N_d N_d],[min(EEG_channel_data) max(EEG_channel_data)],'r');
channel_title = char(EEG.labels(channel_index));
title(channel_title);

subplot(4,1,2);
contourf(EEG.times,frex,eegpower,40,'linecolor','none');
ylabel('Time-frequency map');
axis tight
colormap jet
hold on 
plot([N_d N_d],[1 89],'b');

subplot(4,1,3);
plot(t_2,ER);
ylabel('Energy ratio ER[n]');
axis tight
grid on
hold on 
plot([N_d N_d],[min(ER) max(ER)],'r');

subplot(4,1,4);
plot(t_2,UN);
ylabel('Cumulative sumU_N');
axis tight
grid on
hold on 
plot([N_d N_d],[min(UN) max(UN)],'r'); 

set(gcf,'PaperPositionMode','auto')
print(gcf,'-dpng','-r0',channel_title)
% print(gcf,'-dpng',channel_title);
close all;
 




































fn_out = sprintf('%s%s%s%s%s%s%s%siEEG%s_%.2d.mat',globalVar.([datatype,'Data']),freq_band,filesep,sbj_name,filesep,bn,filesep,freq_band,bn,el);
% fn_out = [globalVar.([datatype,'Data']),filesep,freq_band,filesep,

save(fn_out,'data')
disp(['Wavelet filtering: Block ', bn,', Elec ',num2str(el)])

% EI caluculations for SEEG data ref Bartolomei et al. Brain 2008
% Sheng Jingwei 2015.11.14
% input data
% Updated by Baotian 20170609, Beijing


end
end