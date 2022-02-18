% tic;
load data_getaiming;
EEG.times = 0:0.001:51.999;
EEG.srate = 1000;
EEG.pnts = length(data);
EEG.trials = 1;
%% time-frequency calculation

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
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% get FFT of data
eegfft = fft(data,n_conv_pow2);

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
% toc;

% figure
% subplot(121)
% contourf(EEG.times,frex,eegpower,40,'linecolor','none')
% set(gca,'clim',[-3 3],'xlim',[0 52000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
% title('Logarithmic frequency scaling')
% 
% subplot(122)
% contourf(EEG.times,frex,eegpower,40,'linecolor','none')
% set(gca,'clim',[-3 3],'xlim',[0 52000])
% title('Linear frequency scaling')