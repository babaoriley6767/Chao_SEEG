function wave = sine_wave(freq,duration,amp)
% this function is used to produce a vector of sin wave, the sampling rate
% is fixed to 1000, you need to put the frequency and duration parameter in
% the input of this function
srate = 1000;
t = 0:1/srate:duration;
wave = amp.*sin(freq*2*pi*t);
% plot(t,wave);
xlabel('t');
ylabel('amp');
end