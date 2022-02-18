figure('color','w');
subplot(4,1,1);
plot(t,data_raw);
ylabel('EEG');
axis tight
grid on
hold on 
plot([N_d N_d],[-850 400],'r');

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
plot([N_d N_d],[0 6.8],'r');

subplot(4,1,4);
plot(t_2,UN);
ylabel('Cumulative sumU_N');
axis tight
grid on
hold on 
plot([N_d N_d],[-1250 -0.5],'r');

