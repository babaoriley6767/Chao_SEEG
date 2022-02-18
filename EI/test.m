
t = 0:0.001:60.004;
plot(t,data)
fft_r = fft(data,length(data));
y = abs(fft_r);
plot(t,y)

