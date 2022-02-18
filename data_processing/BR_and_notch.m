function BadChanRejectCAR(sbj_name, project_name, bns, dirs)

for i = 1:length(bns)
    bn = bns{i};
    
    % Load globalVar
    fn = sprintf('%s/originalData/%s/global_%s_%s_%s.mat',dirs.data_root,sbj_name,project_name,sbj_name,bn);
    load(fn,'globalVar');

    
    %% Step 1: Bad channel detection based on raw power
    %This part is for initialization
    
    % Format the data
    sub_info.data.name = sprintf('%s/iEEG%s_',globalVar.originalData,bn);
    
    load([sub_info.data.name '01.mat'])
    ttot=size(wave,2); % ttot = total time in units of samples
    clear wave
    
    data=zeros(ttot,globalVar.nchan); % initialize data variable
    
    %% loop across channels
    for n=1:globalVar.nchan
        if n<10
            nl=['0' num2str(n)];
        else
            nl=num2str(n);
        end % lame hack-around for channel names
        % add notch here.
        load([sub_info.data.name nl '.mat'])
        disp(['notch filtering electrode ' nl ' out of ' num2str(globalVar.nchan)])
        
        %Filtering 60 Hz line noise & harmonics
        wave = noiseFiltData(globalVar, wave);
        wave = wave.'; %because has incorrect ordering for efficiency, should be time x 1
        data(:,n)=-wave; % invert: data are recorded inverted . this might be changed later
        clear wave nl%??????
    end
    data_all = data;
    data=single(data);
    
    %% Run algorithm
       
    a=var(data);
    b=find(a>(10*median(a))); % 10 * greated than median. default 5
    c=find(a<(median(a)/50)); %1/50 smaller than the var. default 1/5
    
    globalVar.too_N_or_S = [b,c];

    channame_spa = globalVar.channame;
    for i = 1: length(globalVar.channame)
        channame_spa{i} = [channame_spa{i},' '];
    end
    if ~isempty([b c])
        disp(['channels that are too noisy or too silent: ' channame_spa{globalVar.too_N_or_S}]);
    end
    
    
    % %Plot bad channels
    figureDim = [0 0 .5 .5];
    figure('units', 'normalized', 'outerposition', figureDim)
    
    for ii = 1:length(globalVar.channame)
        hold on
        if ismember(ii, globalVar.too_N_or_S)
            plot(zscore(data(:,ii))+5*ii,'r');
        else
            plot(zscore(data(:,ii))+5*ii,'b')
        end
    end
    
    ylim([0 length(globalVar.channame)*5+5])
    xlim([0 length(data)])
    
    x = zeros(1,length(globalVar.channame));
    y=(1:length(globalVar.channame))*5;
    
    plot(x,y,'--');
    set(gca,'ytick',y,'yticklabel',globalVar.channame(:))
    
    
    title('Bad electrodes based on raw power (variance)')
    xlabel('Time (s)')
    ylabel('Electrode name')
       
    % Update the globalVar.badChan
    globalVar.badChan = [globalVar.too_N_or_S];
    
    % remove bad channels
    chan_lbls=1:size(data,2);
    data(:,globalVar.badChan)=[];
    chan_lbls(globalVar.badChan)=[];
    
    
    %% Bad channel detection step 3: Bad channel detection based on powerspectra
    set_ov=0; % overlap
    f = 0:250; %
    data_pxx=zeros(length(f),size(data,2));
    
    for k = 1:size(data,2)
        %[Pxx,f] = pwelch(data(1:100*globalVar.iEEG_rate,k),floor(globalVar.iEEG_rate),set_ov,f,floor(globalVar.iEEG_rate));
        [Pxx,f] = pwelch(data(:,k),floor(globalVar.iEEG_rate),set_ov,f,floor(globalVar.iEEG_rate));
        data_pxx(:,k)=Pxx;
    end
    
    % Plot chnas in different colors:
    %     figureDim = [0 0 .5 1];
    %     figure('units', 'normalized', 'outerposition', figureDim)
    plotthis=log(data_pxx);
    figureDim = [0 0 .5 1];
    figure('units', 'normalized', 'outerposition', figureDim)
    for fi = 1:size(plotthis,2)
        hold on
        plot(f,plotthis(:,fi))
        text(0:20:size(plotthis,1), plotthis(1:20:size(plotthis,1),fi), num2str(chan_lbls(fi)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
    end
    xlim([1 size(plotthis,1)])
    ylim([min(plotthis(:)) max(plotthis(:))])
    xlabel('Frequency')
    ylabel('Power')
    
    % Prompt for bad channels
    globalVar.bad_chan_spec = promptBadChanSpec; % of the remaining channels %?command window???????
    
    % Update globalVar.badChan
    globalVar.badChan = [globalVar.badChan globalVar.bad_chan_spec];
    

    %% Re-referencing data to the common average reference CAR - and save
    clear data
    data_all = data_all';
%     data_good = data_all;
%     data_good(globalVar.badChan,:) = [];
%     % Demean good electrodes
%     for ci = 1:size(data_good,1)
%         data_good(ci,:) = data_good(ci,:) - mean(data_good(ci,:));
%     end
%     CAR = mean(data_good,1); % common average reference
%     % Subtract CAR and save
    BR_idx = [];
    for bii=1:size(data_all,1)-1
        name1 = join(regexp(string(globalVar.channame{bii}),'[a-z]','Match','ignorecase'),'');
        name2 = join(regexp(string(globalVar.channame{bii+1}),'[a-z]','Match','ignorecase'),'');
        if strcmp(name1,name2)
            data.wave=data_all(bii,:)-data_all(bii+1,:);
            data.fsample = globalVar.iEEG_rate;
            disp(['Writing: ' sprintf('%s/BRiEEG%s_%.2d.mat',globalVar.BRData,bn, bii)])
        save(sprintf('%s/BRiEEG%s_%.2d.mat',globalVar.BRData,bn, bii),'data');
        clear data
        BR_idx = [BR_idx,bii];
        else
            continue
        end
    end
   disp(['there''re ',num2str(length(BR_idx)),' electrodes in terms of bipolar reference'])
    
    
        % %Plot bad channels
    figureDim = [0 0 .5 .5];
    figure('units', 'normalized', 'outerposition', figureDim)
    hold on
    
    for bii=BR_idx
        idx = find(BR_idx==bii);
        data = data_all(bii,:)-data_all(bii+1,:);
        if ismember(bii, globalVar.too_N_or_S)
            plot(zscore(data)+5*idx,'r');
        else
            plot(zscore(data)+5*idx,'b');
        end
    end
    

    ylim([0 length(BR_idx)*5+5])
    xlim([0 length(data)])
    
    channame_BR = cell(1,length(BR_idx));
    for i = 1:length(BR_idx)
        ii = BR_idx(i);
        channame_BR{i} = [globalVar.channame{ii},'-',globalVar.channame{ii+1}];
    end
    
    
    x = zeros(1,length(channame_BR));
    y=(1:length(channame_BR))*5;
    
    plot(x,y,'--');
    set(gca,'ytick',y,'yticklabel',channame_BR(:))
    
    
    title('Bad electrodes based on raw power (BR)')
    xlabel('Time (s)')
    ylabel('Electrode name (BR)')

    %% Save globalVar (For the last time; don't re-write after this point!)
    fn = sprintf('%s/originalData/%s/global_%s_%s_%s.mat',dirs.data_root,sbj_name,project_name,sbj_name,bn);
    save(fn,'globalVar');
end

end


