a = EpochedData(10,:);
figure
plot(a);
t1 = a(91:101);t2 = a(111:121);
r1 = flip(t1);r2 = flip(t2);
x1 = 1;x2 = 0;tt=1/10;
for i=1:11 
    t1(i) = r1(i)*x1;
    t2(i) = r2(i)*x2;
    x1 = x1-tt;
    x2 = x2-tt;
end
y = t1+t2;
a(101:111) = y;
figure
plot(a);