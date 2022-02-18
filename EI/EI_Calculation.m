% EI caluculations for SEEG data ref Bartolomei et al. Brain 2008
% Sheng Jingwei 2015.11.14
% input data
% Updated by Baotian 20170609, Beijing

%% step 1 calculate the ER
% load 'data.mat' % timeseries of single channel
data_length = length(data);
fs = 1000; % Sampling rate in Hz 
time_length = data_length/fs;
D = 1; % sliding window in seconds
n = time_length/D; % number of intervals
n = fix(n);
N = D*fs; % number of points in each interval
% calculating ER for n intervals
for i = 1:n
x = data(1+N*(i-1):N*i); %extracting time series of each interval
df = fs/(N-1); % resolution of frequency
X = fft(x);% fourior transform to x
X = X(1:N/2);
f = (1:N/2)*df;
% frequency boundary in Hz
theta_boundary = 3.5; 
alpha_boundary = 7.5;
beta_boundary = 12.4;
gama_boundary = 24;
band_boundary = 97;
% location in the corresponding matrix
theta_down = find(abs(f-theta_boundary) == min(abs(f-theta_boundary)));
theta_up = find(abs(f-alpha_boundary) == min(abs(f-alpha_boundary))) - 1;
alpha_down = find(abs(f-alpha_boundary) == min(abs(f-alpha_boundary)));
alpha_up = find(abs(f-beta_boundary) == min(abs(f-beta_boundary))) - 1;
beta_down = find(abs(f-beta_boundary) == min(abs(f-beta_boundary)));
beta_up = find(abs(f-gama_boundary) == min(abs(f-gama_boundary))) - 1;
gama_down = find(abs(f-gama_boundary) == min(abs(f-gama_boundary)));
gama_up = find(abs(f-band_boundary) == min(abs(f-band_boundary)))-1;

X_theta = X(theta_down:theta_up);
X_alpha = X(alpha_down:alpha_up);
X_beta = X(beta_down:beta_up);
X_gama = X(gama_down:gama_up);

E_theta = sum(X_theta.*conj(X_theta))/(2*pi);
E_alpha = sum(X_alpha.*conj(X_alpha))/(2*pi);
E_beta = sum(X_beta.*conj(X_beta))/(2*pi);
E_gama = sum(X_gama.*conj(X_gama))/(2*pi);

ER(i) = (E_beta + E_gama)/(E_alpha + E_theta);
end

%% step 2
% by BT,20170606
v = 0.2;
for i = 1:n
    ERn(i) = sum(ER(1:i))/i;
end
for i = 1:n
    UN(i) = sum(ER(i) - ERn(i) - v);
end











