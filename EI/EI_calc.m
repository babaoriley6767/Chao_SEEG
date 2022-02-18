%% EI calculator Ver. 1.1.20170611 by ZC,WY,BT
% This is the first editon of EI calculation from 
% Brain, 2008, Fabrice Bartolomei
%
% Written by BT 20170609 in Beijing
%
% This is the matlab code version, matlab GUI version is under develepment
%
% All the variable symbols are from the article, you may refer to the
% article to have a better understanding about the whole procedure and make
% modification yourself
function [Nd,ER] = EI_calc(EEG,EEG_channel_data,channel_index)
%% step 1: definition of a statistic that abruptly increases as fast oscillations appear in the signal
%
% data preparation and define variables in this section
% tic;
% clear;
% clc;
% load data_getaiming; % data should be a 1*n vector
data = EEG_channel_data;
fs = EEG.srate; % frequency sampling is defined here
data_length = EEG.TimePoints; % how many dots are in this data
% data_time_length = fix(data_length / fs); % data time in seconds
t = 0:1/fs:(data_length-1)/fs;
% defination of deferent frequency band.
% attention, this is slightly different from the article, part of the
% reason is we need to consider the frequency resoluton after FFT, by using 1000 fs,
% it is better to use integers only, however you can use zero padding to
% increase the frequency resoluton, but I think that is not necessary.
theta_band = 4:7;
alpha_band = 8:12;
beta_band = 13:24;
% gamma_band = 24:97; % this is the raw parameter in the article
gamma_band = 24:97; % revised according to clinical phenomenon

% define sliding window duration D
D = 1; % time window to perform FFT in seconds, time window in article is probably 0.5 according to their map
window_length = D * fs;
% define step length for the sliding window
step_length = 0.01; % step_length in seconds, 0.01 acturally means 10 dot with 1000 fs, step length affect time consuming seriously
t_2 = 0:step_length:(data_length-1)/fs;
% zero paddin the data
zero_padding = zeros(1,window_length - 1);
data = [data zero_padding];
%% fft transform & frequency sub-band energy extraction & ER calculation
ER = zeros(1,fix(data_length/(step_length*fs))); % initialize
for i = 1:(step_length * fs):data_length % slide the window
    data_in_window = data(i:(window_length+i-1)); % extract time series data in one window
%     freq_resolution = fs / window_length; % define freq_resolution
    X_data = periodogram(data_in_window,[],1000); % I need to read more books
%     freq_range = 0:(length(X_data)-1)*freq_resolution; % define the frequency axis
    % extract energy of different band
    E_theta = sum(X_data(theta_band + 1));
    E_alpha = sum(X_data(alpha_band + 1));
    E_beta = sum(X_data(beta_band + 1));
    E_gamma = sum(X_data(gamma_band + 1));
    % calculate time-varying ER
    ER_number = fix(i/(step_length*fs))+1;
    ER(ER_number) = (E_beta + E_gamma)/(E_theta + E_alpha); 
end

%% step 2: optimal detection of rapid discharges
% calculate the vector of ER_n
ER_n = zeros(1,length(ER));
for i = 1:length(ER)
    ER_n(i) = mean(ER(1:i));
end

UN_temp = zeros(1,length(ER));
UN = zeros(1,length(ER));
v = 0.5;
for i = 1:length(ER)
    UN_temp(i) = ER(i) - ER_n(i) - v;
    UN(i) = sum(UN_temp(1:i));
end
% find the local minimum of UN and determine when lamda is reached
lamda = 300;
% calculate local minimum

u_N = zeros(1,length(ER));
for i = 1:length(UN);
    u_N(i) = min(UN(1:i));
end
% calculate threshold
N_a = [];
for i = 1:length(UN);
    if UN(i) - u_N(i) >= lamda;
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
    N_d = length(UN)/100;
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


% toc;
%% calculate time-frequency parameters
EEG.times = t;
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
n_conv_pow2          = pow2(nextpow2(n_convolution));%for the speed consern
half_of_wavelet_size = (n_wavelet-1)/2;

% get FFT of data
eegfft = fft(EEG_channel_data,n_conv_pow2);%%chao fft?

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
end
    
    











