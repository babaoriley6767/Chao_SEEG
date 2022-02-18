clear all;
a = sine_wave(10,30,1);
b = sine_wave(50,5,1);
b2 = sine_wave(30,5,1);
b3 = sine_wave(20,5,1);
b4 = sine_wave(10,15,1);
c = [a b b2 b3 b4];
t = 0:0.001:((length(c)-1)/1000);
plot(t,c);
axis tight