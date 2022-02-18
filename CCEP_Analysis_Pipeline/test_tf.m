files =dir();
for ii =1:length(files)-2
    cd(files(ii+2).name);
    clear t_files;
    t_files =dir();
    for jj =1:length(t_files)-2
        cd(t_files(jj+2).name);
        clear D;clear lab;clear trial;
        if (ii==1|ii==5|ii==6)
          load(strcat('BipM_fffffbemeeg_',t_files(jj+2).name,'.mat'))
        else load(strcat('BipM_fffffbedmeeg_',t_files(jj+2).name,'.mat')) 
        end
        load('seq.mat');
        lab = D.channels; trial = D.trials; D = meeg(D);
        for i =1:length(seq)
            mkdir(lab(seq(i)).label);
            cd(lab(seq(i)).label);
            for k =1:length(trial)
                    clear DD;clear DDD;clear tDD;clear P;clear PP;clear l_PP;clear l_PP_baseline;clear tPP;
                    DD = D(seq(i),:,k);
                    tDD =DD;
                    t = 0:1/1000:1.5;
                    freq = 1:300;
                    P = morlet_transform(tDD, t, freq);
                    PP = squeeze(P);
                    l_PP = PP;
                    l_PP_baseline = l_PP(:,100:450);
                    tPP = (l_PP- mean(l_PP_baseline,2))./std(l_PP_baseline,0,2);
                    baselined_PP(:,:,k) = tPP;
                    imagesc(tPP); axis xy; set(gca,'CLim',[0 10]); print(num2str(k),'-dpng'); close;
                    gamma_n1(k) = sum(sum(tPP(30:100,510:550)))/71/41;
                    gamma_n2(k) = sum(sum(tPP(30:100,600:700)))/71/101;
                    ripple_n1(k) = sum(sum(tPP(100:200,510:550)))/101/41;
                    ripple_n2(k) = sum(sum(tPP(100:200,600:700)))/101/101;
             end
                save('gamma_n1.mat','gamma_n1');
                save('gamma_n2.mat','gamma_n2');
                save('ripple_n1.mat','ripple_n1');
                save('ripple_n2.mat','ripple_n2');
                save('baselined_PP.mat','baselined_PP');
                clear gamma_n1;clear gamma_n2;clear ripple_n1;clear ripple_n2;clear baselined_PP;        
                                       
            cd ..;        
                       
        end                       
        cd ..;    
    end    
    cd ..;
end

load('BipM_fffffbemeeg_A1A2.mat');
load('seq.mat');
lab = D.channels;trial = D.trials;
seq = [23 24 25 34 35]; save('seq.mat','seq');

D = meeg(D);
for i =1:length(seq)
    mkdir(lab(seq(i)).label);
    cd(lab(seq(i)).label);
    for j =1:length(seq)
        mkdir(lab(seq(j)).label);
        cd(lab(seq(j)).label);
        for k =1:length(trial)
            clear DD;clear DDD;clear tDD;clear P;clear PP;clear l_PP;clear l_PP_baseline;clear tPP;
            DD = D(seq(j),:,k);

        
    end
            
        
        
        
        
        cd ..;
    end
    cd ..;
end

        
        
    




lab = D.chanlabels;
DD = D(13,:,:);
squeeze(DD);
figure;
DD = D(13,:,1);
plot(DD);
DDD = squeeze(DD);
for j = 12:13      % see this!!     1:length(D.chanlabels)41
    mkdir(D.chanlabels{j});
    cd(D.chanlabels{j});
    for i =1:length(D.events)
        clear DD;clear DDD;clear tDD;clear P;clear PP;clear l_PP;clear l_PP_baseline;clear tPP;
        DD = D(j,:,i);
%       plot(DD);grid on; figure;plot(tDD); grid on;
%       tDD = remove_art(DD,500);
        tDD =DD;
        t = 0:1/1000:1.5;
        freq = 1:300;
        P = morlet_transform(tDD, t, freq);
        PP = squeeze(P);
        o_PP{i} = PP;
        %imagesc(PP);
        %imagesc(log10(PP));   axis xy;    set(gca,'CLim',[0 10]);    print(num2str(i),'-dpng');        close;
        %l_PP =log10(PP);
        l_PP = PP;
        l_PP_baseline = l_PP(:,100:450);
        tPP = (l_PP- mean(l_PP_baseline,2))./std(l_PP_baseline,0,2);
        %tPP = (l_PP - mean(l_PP_baseline,2)) ./mean(l_PP_baseline,2);
        imagesc(tPP);        axis xy;        set(gca,'CLim',[0 10]); print(num2str(i),'-dpng');  close;
        baselined_PP(:,:,i) = tPP;
        %tPP =l_PP;
        gamma_n1(i) = sum(sum(tPP(30:100,510:550)))/71/41;
        gamma_n2(i) = sum(sum(tPP(30:100,600:700)))/71/101;
        ripple_n1(i) = sum(sum(tPP(100:200,510:550)))/101/41;
        ripple_n2(i) = sum(sum(tPP(100:200,600:700)))/101/101;
        
    end
    save('gamma_n1.mat','gamma_n1');
    save('gamma_n2.mat','gamma_n2');
    save('ripple_n1.mat','ripple_n1');
    save('ripple_n2.mat','ripple_n2');
    save('baselined_PP.mat','baselined_PP');
    cd ..;
    clear gamma_n1;clear gamma_n2;clear ripple_n1;clear ripple_n2;clear baselined_PP;
end


imagesc(mean(baselined_PP,3));
axis xy;
set(gca,'CLim',[0 10]);
print('averaged','-dpng');


l_PP = log10(PP);
imagesc(l_PP);
axis xy;
set(gca,'CLim',[-1 1]);
l_PP_baseline = l_PP(:,200:300);
imagesc(l_PP_baseline);
tPP = (l_PP - mean(l_PP_baseline,2)) ./mean(l_PP_baseline,2);
imagesc(tPP);        
axis xy;
set(gca,'CLim',[-1 1]);

DD(12,:,:);
squeeze(ans);
aans = mean(ans,2);
plot(DDD);
plot(aans);
tt =aans;

taans = remove_art(tt,500);
plot(taans);
t = -0.5:1/1000:1;
freq = 1:300;
P = morlet_transform(taans', t, freq);

PP = squeeze(P);
imagesc(log10(PP));
axis xy;
set(gca,'CLim',[-1.5 2]);

imagesc(lPP);

PP_baseline = PP(:,200:300);
imagesc(log10(PP_baseline));

plot(gamma_n1,'.');
mean(gamma_n1); median(gamma_n1);

Q1=prctile(gamma_n1,25);
Q3=prctile(gamma_n1,75);
IQR = Q3-Q1;
B = gamma_n1;
B(B>1.5*IQR+Q3|B<Q1-1.5*IQR)=[];
plot(B,'*');mean(B);median(B);
for  i =1:length(gamma_n1)
    if gamma_n1(i)

plot(gamma_n2,'.');
mean(gamma_n2); median(gamma_n2);


plot(ripple_n1,'.');
mean(ripple_n1); median(ripple_n1);

plot(ripple_n2,'.');
mean(ripple_n2); median(ripple_n2);



