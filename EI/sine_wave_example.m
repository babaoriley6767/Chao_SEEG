wave_5  = sine_wave(5,0.5,1);
wave_5_big  = sine_wave(5,0.5,2);
wave_10 = sine_wave(10,0.5,1);
wave_18 = sine_wave(18,0.5,1);
wave_50 = sine_wave(50,0.5,1);
t = 0:0.001:0.5;
figure('color','w');
subplot(5,1,1);
plot(t,wave_5);
title('5 Hz');
axis tight
grid on

subplot(5,1,2);
plot(t,wave_5_big);
title('5 Hz');
axis tight
grid on

subplot(5,1,3);
plot(t,wave_10);
title('10 Hz');
axis tight
grid on

subplot(5,1,4);
plot(t,wave_18);
title('18 Hz');
axis tight
grid on

subplot(5,1,5);
plot(t,wave_50);
title('50 Hz');
axis tight
grid on

wave_example = wave_5 + wave_10 + wave_18 + wave_50
plot(wave_example)
t = 0:2:500;
fft = fft(wave_example);
fft = fft(1:251);
plot3(t,real(fft),imag(fft))