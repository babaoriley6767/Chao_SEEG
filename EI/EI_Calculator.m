%% EI calculator Ver. 1.0.20170609 by ZC,WY,BT
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

%% step 1: definition of a statistic that abruptly increases as fast oscillations appear in the signal
%
% data preparation and define variables in this section
clear;
clc;
load data_getaiming; % data should be a 1*n vector
fs = 1000; % frequency sampling is defined here
data_length = length(data); % how many dots are in this data
data_time_length = fix(data_length / fs); % data time in seconds
t = 0:1/fs:(data_length-1)/fs;
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
D = 1; % time window to perform FFT in seconds, time window in article is probably 0.5 according to their map % this could be changed 
window_length = D * fs;
% define step length for the sliding window
step_length = 0.01; % step_length in seconds, 0.01 acturally means 10 dot with 1000 fs, step length affact time consuming seriously % this could be changed 

%% fft transform & frequency sub-band energy extraction & ER calculation
for i = 1:(step_length * fs):(data_length - D * fs) % slide the window
    data_in_window = data(i:(window_length+i-1)); % extract time series data in one window
    freq_resolution = fs / window_length; % define freq_resolution
    X_data = periodogram(data_in_window); % I need to read more books
    freq_range = 0:(length(X_data)-1)*freq_resolution; % define the frequency axis
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
for i = 1:length(ER)
    ER_n(i) = mean(ER(1:i));
end
v = 0.5;
for i = 1:length(ER)
    UN_temp(i) = ER(i) - ER_n(i) - v;
    UN(i) = sum(UN_temp(1:i));
end

    
    
    
    











