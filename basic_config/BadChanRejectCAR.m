function BadChanRejectCAR(sbj_name, project_name, bns, dirs)

for i = 1:length(bns)
    bn = bns{i};
    
    % Load globalVar
    fn = sprintf('%s/originalData/%s/global_%s_%s_%s.mat',dirs.data_root,sbj_name,project_name,sbj_name,bn);
    load(fn,'globalVar');
    
    %Chan info PROMPT
    globalVar.badChan = [globalVar.refChan globalVar.epiChan];
    
    %% Step 1: Bad channel detection based on raw power
    %This part saves data.mat
    
    % Format the data
    sub_info.data.name = sprintf('%s/iEEG%s_',globalVar.originalData,bn);
    sub_info.Pdiode.name = sprintf('%s/Pdio%s_02',globalVar.originalData,bn);% ?????
    
    data_nr=1;
    load([sub_info.data(data_nr).name '01.mat'])
    ttot=size(wave,2); % ttot = total time in units of samples  %????globalVar?????????
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
        load([sub_info.data(data_nr).name nl '.mat'])
        disp(['notch filtering electrode ' nl ' out of ' num2str(globalVar.nchan)])
        
        %Filtering 60 Hz line noise & harmonics
        wave = noiseFiltData(globalVar, wave);
        wave = wave.'; %because has incorrect ordering for efficiency, should be time x 1
        data(:,n)=-wave; % invert: data are recorded inverted
        clear wave nl%??????
    end
    data_all = data;
    data=single(data);
    
    %% Run algorithm
    bad_chans=[globalVar.refChan globalVar.badChan globalVar.epiChan globalVar.emptyChan];
    a=var(data);
    b=find(a>(5*median(a))); % 5 * greated than median. ?????????  %???????????
    c=find(a<(median(a)/5)); %1/5 smaller than the var. ?????????
    if ~isempty([b c])
        disp(['additional bad channels: ' int2str(setdiff([b c],bad_chans))]);
    end
    disp('add additional bad channel to subj_info-bad_channels2 if you want to cut them')
    
    % %Plot bad channels
    bad_cha_tmp = setdiff([b c],[globalVar.badChan globalVar.epiChan]);
    figureDim = [0 0 .5 .5];
    figure('units', 'normalized', 'outerposition', figureDim)
    subplot(2,1,1)
    for ii = bad_cha_tmp
        hold on
        plot(zscore(data(:,ii))+ii);%?2?ii???????overlap????
    end
    title('Bad electrodes based on raw power')
    xlabel('Time (s)')
    ylabel('Electrode number')
    % Update the globalVar.badChan
    globalVar.badChan = [bad_cha_tmp globalVar.badChan];
    
    % remove bad channels
    chan_lbls=1:size(data,2);
    data(:,globalVar.badChan)=[];
    chan_lbls(globalVar.badChan)=[];
    
    %% Step 2: Bad channel detection based on spikes in the raw data% ?????????80uv??????????????
    nr_jumps=zeros(1,size(data,2));
    for k=1:size(data,2)
        nr_jumps(k)=length(find(diff(data(:,k))>80)); % find jumps>80uV, calculate the sum of the point that larger than 80
    end
    
    subplot(2,1,2)
    plot(nr_jumps);
    ylabel('Number of spikes')
    xlabel('Electrode number')
    ej= floor(globalVar.chanLength/globalVar.iEEG_rate);% 1 jump per second in average
    jm=find(nr_jumps>ej);% Find channels with more than ... jumps
    clcl=chan_lbls(jm);% Real channel numbers of outliers
    disp(['spiky channels: ' int2str(clcl)]);
    
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
    bad_chan_spec = promptBadChanSpec; % of the remaining channels %?command window???????
    
    % Update globalVar.badChan
    globalVar.badChan = [bad_chan_spec globalVar.badChan];
    
    %% Bad channel detection step 4: Bad channel detection based on HFOs
    [pathological_chan_id,pathological_event] = find_paChan(data_all,globalVar.channame,globalVar.iEEG_rate, 1.5);
    % pathological_event are in bipolar montage
    
    %% Inspect all bad channels
    % blue channels - detected in steps 1,2,3
    % red channels - detected in step 4
    % green channels - detected in both steps (1,2,3) and 4
    figureDim = [0 0 .5 .5];
    figure('units', 'normalized', 'outerposition', figureDim)
    for ii = unique([globalVar.badChan pathological_chan_id])
        hold on
        if ~isempty(find(pathological_chan_id == ii)) && ~isempty(find(globalVar.badChan == ii))
            color_plot = [0 1 0];
        elseif ~isempty(find(pathological_chan_id == ii)) && isempty(find(globalVar.badChan == ii))
            color_plot = [1 0 0];
        elseif isempty(find(pathological_chan_id == ii)) && ~isempty(find(globalVar.badChan == ii))
            color_plot = [0 0 1];
        end
        plot(zscore(data_all(:,ii))+ii, 'Color', color_plot);
    end
    title('All bad electrodes')
    xlabel('Time (s)')
    ylabel('Electrode number')
    
    
    % Update globalVar.badChan
    globalVar.badChan = unique([pathological_chan_id globalVar.badChan]);
    globalVar.pathological_event_bipolar_montage = pathological_event;
    
    %% Eyeball the rereferenced data after removing the bad channels
    % This should be interactive - ask Su to help creating a gui.
%     data_all(:,globalVar.badChan)=[];
%     data_down_car = car(data);
    
    %% Plot CAR data for eyeballing
%     for iii = 1:size(data_down_car,2)
%         hold on
%         plot(zscore(data_down_car(1:round(globalVar.iEEG_rate*20),iii))+iii);
%     end
    
    %% Re-referencing data to the common average reference CAR - and save
    clear data
    data_all = data_all';
    data_good = data_all;
    data_good(globalVar.badChan,:) = [];
    % Demean good electrodes
    for ci = 1:size(data_good,1)
        data_good(ci,:) = data_good(ci,:) - mean(data_good(ci,:));
    end
    CAR = mean(data_good,1); % common average reference
    % Subtract CAR and save   
    for cii = 1:size(data_all,1)
        %%% WHY DO I HAVE TO SWITCH SIGN TO MATCH PREVIOUS (with separate FiltData)?!!! %%%
        data.wave = -single(data_all(cii,:) - CAR); % subtract CAR 
        %%% WHY DO I HAVE TO SWITCH SIGN TO MATCH PREVIOUS? (with separate FiltData)!!! %%%
        data.fsample = globalVar.iEEG_rate;
        disp(['Writing: ' sprintf('%s/CARiEEG%s_%.2d.mat',globalVar.CARData,bn, cii)])
        save(sprintf('%s/CARiEEG%s_%.2d.mat',globalVar.CARData,bn, cii),'data');
        clear data
    end
    
    %% Save globalVar (For the last time; don't re-write after this point!)
    fn = sprintf('%s/originalData/%s/global_%s_%s_%s.mat',dirs.data_root,sbj_name,project_name,sbj_name,bn);
    save(fn,'globalVar');
end

end



