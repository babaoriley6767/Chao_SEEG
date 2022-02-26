function EI_pipeline_Chao(sbj_name,project_name,block_names,dirs,datatype,freq_band,Nu,Lamda,Eta,EI_window)
%% INPUTS
%   sbj_name:               subject name
%   project_name:           name of task
%   block_names:                     names of blocks to be analyed (cell of strings)
%   dirs:                   directories pointing to files of interest (generated by InitializeDirs)

load([dirs.original_data filesep sbj_name filesep 'subjVar_volume_BR_' sbj_name '.mat'], 'subjVar_volume_BR');
anat_idx = ~strcmp(subjVar_volume_BR.eleinfo.AAL3,'NotAvailable');
YEO7_idx = ~strcmp(subjVar_volume_BR.eleinfo.YEO7,'NotAvailable');
gray_idx = YEO7_idx & anat_idx;


for bi = 1%:length(block_names)
    bn = block_names{bi};
    
    %% Load globalVar
    fn = sprintf('%s/originalData/%s/global_%s_%s_%s.mat',dirs.data_root,sbj_name,project_name,sbj_name,bn);
    load(fn,'globalVar');
    
    if ~isempty([dirs.original_data filesep sbj_name filesep 'sz_loc_EI_' [sbj_name,'_',bn] '.mat'])
        prompt = ['sz_loc_EI already exist for ' [sbj_name,'_',bn] ' . Do you want to use the previous parameters?'] ;
        ID = input(prompt,'s');
        if strcmp(ID, 'y')
            load([dirs.original_data filesep sbj_name filesep 'sz_loc_EI_' [sbj_name,'_',bn] '.mat'])
            Nu = sz_loc_EI.Nu;
            Lamda = sz_loc_EI.Lamda;
            Eta = sz_loc_EI.Eta;
            EI_window = sz_loc_EI.EI_window;
            disp('we''re now using the previous parameters~')
        else
            warning('even though the sz_loc_EI already exist, we''re now using the new parameter!!!')
        end
    else
    end
    
    if strcmp(datatype,'Band')
        data_root=globalVar.BandData;
    else
        data_root=globalVar.SpecData;
    end
    
    time_series = 0:1/globalVar.iEEG_rate:(1/globalVar.iEEG_rate)*(globalVar.chanLength-1);
    
    if ~isempty(EI_window)
        EEG_win = dsearchn(time_series',EI_window');
        ER_win  = dsearchn(globalVar.ER_series',EI_window');
        globalVar.chanLength = EEG_win(end)-EEG_win(1)+1;
        globalVar.ER_length = ER_win(end)-ER_win(1)+1;
        time_series = time_series(EEG_win(1):EEG_win(end));
        globalVar.ER_series = globalVar.ER_series(ER_win(1):ER_win(end));
    else
        EI_window = zeros(1,2);
        EI_window(1) = round(globalVar.ER_series(1));
        EI_window(2) = round(globalVar.ER_series(end));
        EEG_win = dsearchn(time_series',EI_window');
        ER_win  = dsearchn(globalVar.ER_series',EI_window');
        globalVar.chanLength = EEG_win(end)-EEG_win(1)+1;
        globalVar.ER_length = ER_win(end)-ER_win(1)+1;
        time_series = time_series(EEG_win(1):EEG_win(end));
        globalVar.ER_series = globalVar.ER_series(ER_win(1):ER_win(end));
    end
    
    %data initialization
    BR_matrix = zeros(length(globalVar.channame_BR(gray_idx)),globalVar.chanLength);
    ER_matrix = zeros(length(globalVar.channame_BR(gray_idx)),globalVar.ER_length);
    
    
    for i = 1:length(globalVar.channame_BR(gray_idx))
        channame_idx = find(gray_idx,i);
        el = channame_idx(end);
        %load ER
        fn_out = sprintf('%s%s%s%s%s%s%s%siEEG%s_%.2d.mat',globalVar.([datatype,'Data']),freq_band,filesep,sbj_name,filesep,bn,filesep,freq_band,bn,el);
        load(fn_out);
        BR_matrix(i,:) = data.wave(EEG_win(1):EEG_win(end));
        ER_matrix(i,:) = data.ER(ER_win(1):ER_win(end));
    end
    
    ER_step = globalVar.ER_series(2) - globalVar.ER_series(1);
    
    %Nu = 1;%%%%
    
    %Lamda = 400;%%%%%
    %channel sellection
    %decimate shoud be 250Hz
    
    ER_n = zeros(size(ER_matrix));
    UN_temp = zeros(size(ER_matrix));
    UN = zeros(size(ER_matrix));
    u_N = zeros(size(ER_matrix));
    
    for i = 1:length(ER_matrix)
        ER_n(:,i) = mean(ER_matrix(:,1:i),2);
    end
    
    %     for i = 1:length(ER_matrix)
    %         UN_temp(:,i) = ER_matrix(:,i) - ER_n(:,i)-Nu;
    %         UN(:,i) = sum(UN_temp(:,1:i),2);
    %     end
    for i = 1:length(ER_matrix)
        UN_temp(:,i) = ER_matrix(:,i) - ER_n(:,i)-Nu;
    end
    
    
    for i = 1:length(ER_matrix)
        UN(:,i) = sum(UN_temp(:,1:i),2);
    end
    
    for i = 1:length(ER_matrix)
        u_N(:,i ) = min(UN(:,1:i),[],2);
    end
    
    N_a = cell(length(globalVar.channame_BR(gray_idx)),1);%alarm time
    N_d = cell(length(globalVar.channame_BR(gray_idx)),1);%detection time
    
    for i = 1:length(N_a)
        N_a{i,1} = find((UN(i,:) - u_N(i,:)) >= Lamda,1);
        if isempty(N_a{i,1})
            N_a{i,1} = globalVar.ER_length - 1;
            N_d{i,1} = globalVar.ER_length - 1;
        else
            N_d{i,1} = find(u_N(i,:) == u_N(i,N_a{i,1}),1);
        end
    end
    
    
    
    %      = reshape(normalize(reshape(ER_matrix,1,[]),'range'),length(globalVar.channame_BR),[]);
    ER_matrix_norm = normalize(ER_matrix,2,'range');
    
    
    amp = 1;
    % pattern selection?
    %% Plot bad channels
    figureDim = [1 1 1 1];%[0 0 .5 .5]
    figure('units', 'normalized', 'outerposition', figureDim)
    hold on
    
    BR_matrix_plot = bsxfun(@plus, zscore(BR_matrix,0,2)*amp, [0:-5:-5*(length(globalVar.channame_BR(gray_idx))-1)]');
    
    x_channame = zeros(1,length(globalVar.channame_BR(gray_idx)));
    
    y_axis = (length(globalVar.channame_BR(gray_idx))-1)*-5 :5 :0;
    
    y_axis_all = length(globalVar.channame_BR(gray_idx))*-5 :5 :5;
    
    
    ER_plot = contourf(globalVar.ER_series,y_axis_all,flip([zeros(1,length(globalVar.ER_series));ER_matrix;zeros(1,length(globalVar.ER_series))]),100,'linecolor','none');
    set(gca,'clim',[0 25]) %ER plot
    BR_plot = plot(time_series,BR_matrix_plot,'k');%iEEG_plot
    
    for i = 1:length(globalVar.channame_BR(gray_idx))
        N_a_plot = plot(globalVar.ER_series,repmat(0-5*(i-1),1,globalVar.ER_length),'y','LineStyle', 'none' ,'Marker','x','MarkerSize', 12,'LineWidth',2,'MarkerIndices',N_a{i});
    end
    
    for i = 1:length(globalVar.channame_BR(gray_idx))
        N_d_plot = plot(globalVar.ER_series,repmat(0-5*(i-1),1,globalVar.ER_length),'y','LineStyle', 'none' ,'Marker','s','MarkerSize', 12,'LineWidth',2,'MarkerIndices',N_d{i});
    end
    
    
    % plot(x_channame,y_axis,'--');%channame
    
    channame_BR_id = arrayfun(@num2str,(1:length(globalVar.channame_BR)),'UniformOutput',0);
    
    space_connect   = cell(1, length(globalVar.channame_BR));
    space_connect(:) = {'   '};
    
    EI_y_label = strcat(globalVar.channame_BR,space_connect,channame_BR_id);
    EI_y_label = EI_y_label(gray_idx);
    EI_y_label = flip(EI_y_label);
    
    set(gca,'ytick',y_axis,'yticklabel',EI_y_label);
    
    %     ylim([-5 length(globalVar.channame_BR)*5+5])
    ylim([length(globalVar.channame_BR(gray_idx))*-5 5])
    %     xlim([0 length(BR_matrix)])
    xlim([0,round(globalVar.ER_series(end))])
    
    
    if ~isempty(EI_window)
        xlim(EI_window)
    else
        xlim([0,round(globalVar.ER_series(end))])
    end
    
    
    
    title(['iEEG raw ',bn])
    xlabel('Time (s)')
    ylabel('Electrode name (BR)')
    
    %% identify the bad channel mannuly
    
    if ~isempty([dirs.original_data filesep sbj_name filesep 'sz_loc_EI_' [sbj_name,'_',bn] '.mat'])
        prompt = ['sz_loc_EI.bad_chan already exist for ' [sbj_name,'_',bn] ' . Do you want to use the previous parameters?'] ;
        ID = input(prompt,'s');
        if strcmp(ID, 'y')
            load([dirs.original_data filesep sbj_name filesep 'sz_loc_EI_' [sbj_name,'_',bn] '.mat'])
            bad_chan = sz_loc_EI.bad_chan;
            disp('we''re now using the previous parameters~')
        else
            [bad_chan] = promptBadChanSpec;
            sz_loc_EI.bad_chan = bad_chan;
            warning('even though the sz_loc_EI.bad_chan already exist, we''re now using the new parameter!!!')
        end
    else
        [bad_chan] = promptBadChanSpec;
        sz_loc_EI.bad_chan = bad_chan;
    end
    
    gray_idx(bad_chan) = 0;
    %% plot again
    close all
    
    %data initialization
    BR_matrix = zeros(length(globalVar.channame_BR(gray_idx)),globalVar.chanLength);
    ER_matrix = zeros(length(globalVar.channame_BR(gray_idx)),globalVar.ER_length);
    
    
    for i = 1:length(globalVar.channame_BR(gray_idx))
        channame_idx = find(gray_idx,i);
        el = channame_idx(end);
        %load ER
        fn_out = sprintf('%s%s%s%s%s%s%s%siEEG%s_%.2d.mat',globalVar.([datatype,'Data']),freq_band,filesep,sbj_name,filesep,bn,filesep,freq_band,bn,el);
        load(fn_out);
        if ~isempty(EI_window)
            BR_matrix(i,:) = data.wave(EEG_win(1):EEG_win(end));
            ER_matrix(i,:) = data.ER(ER_win(1):ER_win(end));
        else
            BR_matrix(i,:) = data.wave;
            ER_matrix(i,:) = data.ER;
        end
    end
    
    ER_step = globalVar.ER_series(2) - globalVar.ER_series(1);
    
    %Nu = 1;%%%%
    
    %Lamda = 400;%%%%%
    %channel sellection
    %decimate shoud be 250Hz
    
    ER_n = zeros(size(ER_matrix));
    UN_temp = zeros(size(ER_matrix));
    UN = zeros(size(ER_matrix));
    u_N = zeros(size(ER_matrix));
    
    for i = 1:length(ER_matrix)
        ER_n(:,i) = mean(ER_matrix(:,1:i),2);
    end
    
    %     for i = 1:length(ER_matrix)
    %         UN_temp(:,i) = ER_matrix(:,i) - ER_n(:,i)-Nu;
    %         UN(:,i) = sum(UN_temp(:,1:i),2);
    %     end
    for i = 1:length(ER_matrix)
        UN_temp(:,i) = ER_matrix(:,i) - ER_n(:,i)-Nu;
    end
    
    
    for i = 1:length(ER_matrix)
        UN(:,i) = sum(UN_temp(:,1:i),2);
    end
    
    for i = 1:length(ER_matrix)
        u_N(:,i ) = min(UN(:,1:i),[],2);
    end
    
    N_a = cell(length(globalVar.channame_BR(gray_idx)),1);%alarm time
    N_d = cell(length(globalVar.channame_BR(gray_idx)),1);%detection time
    
    for i = 1:length(N_a)
        N_a{i,1} = find((UN(i,:) - u_N(i,:)) >= Lamda,1);
        if isempty(N_a{i,1})
            N_a{i,1} = globalVar.ER_length - 1;
            N_d{i,1} = globalVar.ER_length - 1;
        else
            N_d{i,1} = find(u_N(i,:) == u_N(i,N_a{i,1}),1);
        end
    end
    
    
    figureDim = [1 1 1 1];%[0 0 .5 .5]
    figure('units', 'normalized', 'outerposition', figureDim)
    hold on
    
    BR_matrix_plot = bsxfun(@plus, zscore(BR_matrix,0,2)*amp, [0:-5:-5*(length(globalVar.channame_BR(gray_idx))-1)]');
    
    x_channame = zeros(1,length(globalVar.channame_BR(gray_idx)));
    
    y_axis = (length(globalVar.channame_BR(gray_idx))-1)*-5 :5 :0;
    
    y_axis_all = length(globalVar.channame_BR(gray_idx))*-5 :5 :5;
    
    
    ER_plot = contourf(globalVar.ER_series,y_axis_all,flip([zeros(1,length(globalVar.ER_series));ER_matrix;zeros(1,length(globalVar.ER_series))]),100,'linecolor','none');
    set(gca,'clim',[0 25]) %ER plot
    BR_plot = plot(time_series,BR_matrix_plot,'k');%iEEG_plot
    
    for i = 1:length(globalVar.channame_BR(gray_idx))
        N_a_plot = plot(globalVar.ER_series,repmat(0-5*(i-1),1,globalVar.ER_length),'y','LineStyle', 'none' ,'Marker','x','MarkerSize', 12,'LineWidth',2,'MarkerIndices',N_a{i});
    end
    
    for i = 1:length(globalVar.channame_BR(gray_idx))
        N_d_plot = plot(globalVar.ER_series,repmat(0-5*(i-1),1,globalVar.ER_length),'y','LineStyle', 'none' ,'Marker','s','MarkerSize', 12,'LineWidth',2,'MarkerIndices',N_d{i});
    end
    
    
    % plot(x_channame,y_axis,'--');%channame
    
    channame_BR_id = arrayfun(@num2str,(1:length(globalVar.channame_BR)),'UniformOutput',0);
    
    space_connect   = cell(1, length(globalVar.channame_BR));
    space_connect(:) = {'   '};
    
    EI_y_label = strcat(globalVar.channame_BR,space_connect,channame_BR_id);
    EI_y_label = EI_y_label(gray_idx);
    EI_y_label = flip(EI_y_label);
    
    set(gca,'ytick',y_axis,'yticklabel',EI_y_label);
    
    %     ylim([-5 length(globalVar.channame_BR)*5+5])
    ylim([length(globalVar.channame_BR(gray_idx))*-5 5])
    %     xlim([0 length(BR_matrix)])
    xlim([0,round(globalVar.ER_series(end))])
    
    
    if ~isempty(EI_window)
        xlim(EI_window)
    else
        xlim([0,round(globalVar.ER_series(end))])
    end
    
    
    
    title(['iEEG raw ',bn])
    xlabel('Time (s)')
    ylabel('Electrode name (BR)')
    
    %% calculate the EI index
    N_0 = min([N_d{:}]);
    
    EI = cell(length(globalVar.channame_BR(gray_idx)),1);
    EI_index = cell(length(globalVar.channame_BR(gray_idx)),1);
    
    N_d_time = cell(length(globalVar.channame_BR(gray_idx)),1);
    N_d_time_index = cell(length(globalVar.channame_BR(gray_idx)),1);
    
    N_d_ER = cell(length(globalVar.channame_BR(gray_idx)),1);
    N_d_ER_index = cell(length(globalVar.channame_BR(gray_idx)),1);
    
    for i = 1:length(globalVar.channame_BR(gray_idx))
        if N_d{i,1}==globalVar.ER_length - 1
            EI{i,1} = nan;
            N_d_time{i,1} = nan;
            N_d_ER{i,1} = nan;
        else
            if globalVar.ER_series(N_d{i,1}) <= round(globalVar.ER_series(end)) - Eta
                EI{i,1} = (1/(globalVar.ER_series(N_d{i,1})-globalVar.ER_series(N_0)+1))*sum(ER_matrix(i,N_d{i,1}:(N_d{i,1}+round(Eta/ER_step))));
                N_d_time{i,1} = globalVar.ER_series(N_d{i,1})-globalVar.ER_series(N_0);
                N_d_ER{i,1} = sum(ER_matrix(i,N_d{i,1}:(N_d{i,1}+round(Eta/ER_step))));
            else
                %                 EI{i,1} = 1/(globalVar.ER_series(N_d{i,1})-globalVar.ER_series(N_0)+1);
                EI{i,1} = nan;
                N_d_time{i,1} = nan;
                N_d_ER{i,1} = nan;
            end
        end
    end
    
    EI_max = max([EI{:}]);
    N_d_ER_max = max([N_d_ER{:}]);
    
    for i = 1:length(globalVar.channame_BR(gray_idx))
        EI_index{i,1} = EI{i,1}./EI_max;
    end
    
    for i = 1:length(globalVar.channame_BR(gray_idx))
        N_d_ER_index{i,1} = N_d_ER{i,1}./N_d_ER_max;
    end
    
    sz_loc_EI.eleinfo = subjVar_volume_BR.eleinfo;
    
    % creat the  EI data in each table element
    name = 'EI';
    EI_table = cell(length(gray_idx),1);
    EI_table(gray_idx) = EI_index(:);
    EI_table(~gray_idx) = {NaN};
    sz_loc_EI.eleinfo.(name) = EI_table;%%%
    
    % creat the time interval between the EI =1 and other electrodes in
    % each table element
    name = 'time_interval';
    N_d_time_table = cell(length(gray_idx),1);
    N_d_time_table(gray_idx) = N_d_time(:);
    N_d_time_table(~gray_idx) = {NaN};
    sz_loc_EI.eleinfo.(name) = N_d_time_table;
    
    % creat the time interval index in
    % each table element
    TII_seq = [0:0.2:20]';%1/f parameter
    name = 'time_interval_index';
    N_d_time_index_table = cell(length(gray_idx),1);
    for j = 1:length(N_d_time )
        if isnan(N_d_time{j})
            N_d_time_index{j} = NaN;
        else
            TIIid = dsearchn(TII_seq,N_d_time{j});
            N_d_time_index{j} = 1/TIIid;
        end
    end
    N_d_time_index_table(gray_idx) = N_d_time_index(:);
    N_d_time_index_table(~gray_idx) = {NaN};
    sz_loc_EI.eleinfo.(name) = N_d_time_index_table;
    
    % creat the ER in each table element
    name = 'ER';
    N_d_ER_table = cell(length(gray_idx),1);
    N_d_ER_table(gray_idx) = N_d_ER(:);
    N_d_ER_table(~gray_idx) = {NaN} ;
    sz_loc_EI.eleinfo.(name) = N_d_ER_table;
    
    % creat the ER index in each table element
    name = 'ER_index';
    N_d_ER_index_table = cell(length(gray_idx),1);
    N_d_ER_index_table(gray_idx) = N_d_ER_index(:);
    N_d_ER_index_table(~gray_idx) = {NaN} ;
    sz_loc_EI.eleinfo.(name) = N_d_ER_index_table;
    
    %creat the Euclidean distance in each table element
    
    
    name = 'Eu_dis';
    coords = sz_loc_EI.eleinfo.MNI305_volume(gray_idx,:);
    Eu_dis_table = cell(length(gray_idx),1);
    Eu_dis = cell(length(EI_index),1);
    for j = 1:length(EI_index)
        if isnan(EI_index{j})
            Eu_dis{j} = NaN;
        else
            Eu_dis{j}= norm(coords(cat(1,EI_index{:})==1,:)-coords(j,:));
        end
    end
    Eu_dis_table(gray_idx) = Eu_dis;
    Eu_dis_table(~gray_idx) = {NaN};
    sz_loc_EI.eleinfo.(name) = Eu_dis_table;
    
    
    %% plot the EI on inflated brain
    
    %     addpath(genpath('/Users/tony/Desktop/function_tools/for_plot/iELVis-master/'))
    %
    %     global globalFsDir;
    %     globalFsDir ='/Users/tony/Desktop/iELVis/Plot_Elelctrodes/';
    fsDir='/Users/tony/Desktop/iELVis/Plot_Elelctrodes/';%%% seems useful
    cd([fsDir]);
    
    % plotting parameters
    
    bv_chan_hem = subjVar_volume_BR.eleinfo.bv_hem(gray_idx);
    isLeft = [];
    for j = 1:length(bv_chan_hem)
        if strcmp(bv_chan_hem{j},'R')
            isLeft(j,1) = 0;
        else
            isLeft(j,1) = 1;
        end
    end
    
    bv_chan_names = sz_loc_EI.eleinfo.bv_names(gray_idx,:);
    elecCoord152 = sz_loc_EI.eleinfo.MNI152_volume(gray_idx,:);
    elecCoord305 = sz_loc_EI.eleinfo.MNI305_volume(gray_idx,:);
    %     EI = cell2mat(sz_loc_EI.eleinfo.EI(:));
    EI_color = cell2mat(EI_index);
    EI_color(isnan(EI_color)) = 0;
    
    
    view_side = lower(bv_chan_hem{1});
    cfg=[];
    cfg.view=view_side;
    cfg.clickElec = 'y';
    cfg.elecSize=15;
    cfg.surfType='inflated';
    cfg.opaqueness=0.1;
    cfg.ignoreDepthElec='n';
    cfg.elecNames = bv_chan_names;
    cfg.showLabels='y';
    cfg.elecCoord=[elecCoord305 isLeft];
    cfg.elecColors = EI_color;
    cfg.elecColorScale=[0 1];
    cfg.title = [];
    cfg.backgroundColor = [1,1,1];
    %     cfg.overlayParcellation='Y7';
    cfgOut = plotPialSurf_v2('fsaverage',cfg);
    
    
    disp(['there are ', num2str(sum(EI_color>.3)), ' electrodes that have EI > .3']);
    disp(['the ' bv_chan_names{EI_color==1,1}, '''EI is 1']);
    
    if sz_loc_EI.eleinfo.time_interval{cell2mat(sz_loc_EI.eleinfo.EI)==1} ~= 0
        warning('The detection time of max EI electrode is not the first one!')
    else
    end
    
    prompt = ('Does the EI value meet expectations?y/n') ;
    ID = input(prompt,'s');
    if strcmp(ID, 'y')
    else
        warning('reset the parameters % you could put a break point here and change the EI mannuly in sz_loc_EI.eleinfo')
        return
    end
    
    %% Save subjVar
    sz_loc_EI.Nu = Nu;
    sz_loc_EI.Lamda = Lamda;
    sz_loc_EI.Eta = Eta;
    sz_loc_EI.EI_window = EI_window;
    sz_loc_EI.name = [sbj_name,'_',bn];
    
    sz_loc_EI_idx = find(cat(1,sz_loc_EI.eleinfo.EI{:})==1);
    
    %     sz_loc_EI.AAL3idx_Group = sz_loc_EI.eleinfo.AAL3idx(sz_loc_EI_idx);
    sz_loc_EI.YEO7idx_Group = sz_loc_EI.eleinfo.YEO7idx(sz_loc_EI_idx);
    
    %% stats
    %statement
    eleinfo = sz_loc_EI.eleinfo;
    data_EI = cell2mat(eleinfo.EI);
    data_TII = cell2mat(eleinfo.time_interval_index);
    data_ER_idx = cell2mat(eleinfo.ER_index);
    data_Eu_dis = cell2mat(eleinfo.Eu_dis);
    %     idx_AAL3 = eleinfo.AAL3idx;
    idx_YEO7 = eleinfo.YEO7idx;
    YEO7_v = [5,7,1,2,3,4,6];% the sequence if from PNAS
    %     AAL3_v = [41,42,45,46];%left hippo, right hippo, left amy, right amy
    %     EI_inclu = ismember(idx_AAL3,AAL3_v) | idx_YEO7~=0 ;% all the meaningful EI
    idx_EZ = data_EI>=0.3 & data_EI < 1;
    idx_PZ = data_EI>0 & data_EI<0.3;
    
    %     if ~ismember(sz_loc_EI.YEO7idx_Group,YEO7_v)&&~ismember(sz_loc_EI.AAL3idx_Group,AAL3_v)
    %         warning('we need to fix the EI electrode parcellation')
    %     else
    
    %EI_EZ
    EI_EZ = cell(1,length(YEO7_v));
    EI_PZ = cell(1,length(YEO7_v));
    EI_EZPZ = cell(1,length(YEO7_v));
    
    TII_EZ = cell(1,length(YEO7_v));
    TII_PZ = cell(1,length(YEO7_v));
    TII_EZPZ = cell(1,length(YEO7_v));
    
    ER_idx_EZ = cell(1,length(YEO7_v));
    ER_idx_PZ = cell(1,length(YEO7_v));
    ER_idx_EZPZ = cell(1,length(YEO7_v));
    
    Eu_dis_EZ = cell(1,length(YEO7_v));
    Eu_dis_PZ = cell(1,length(YEO7_v));
    Eu_dis_EZPZ = cell(1,length(YEO7_v));
    
    
    EI_EZ_mat = nan(1,length(YEO7_v));
    for i = 1:length(YEO7_v)
        EI_EZ{i} = data_EI(idx_EZ & ismember(idx_YEO7,YEO7_v(i)));
        EI_EZ_mat(i) = nanmean(data_EI(idx_EZ & ismember(idx_YEO7,YEO7_v(i))));
    end
    
    mat_EI_EZ = nan(7,7);
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    mat_EI_EZ(row_idx,:) = EI_EZ_mat;
    
    
    EI_PZ_mat = nan(1,length(YEO7_v));
    for i = 1:length(YEO7_v)
        EI_PZ{i} = data_EI(idx_PZ & ismember(idx_YEO7,YEO7_v(i)));
        EI_PZ_mat(i) = nanmean(data_EI(idx_PZ & ismember(idx_YEO7,YEO7_v(i))));
    end
    
    mat_EI_PZ = nan(7,7);
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    mat_EI_PZ(row_idx,:) = EI_PZ_mat;
    
    EI_EZPZ_mat = nan(1,length(YEO7_v));
    for i = 1:length(YEO7_v)
        EI_EZPZ{i} = data_EI((idx_EZ|idx_PZ) & ismember(idx_YEO7,YEO7_v(i)));
        EI_EZPZ_mat(i) = nanmean(data_EI((idx_EZ|idx_PZ) & ismember(idx_YEO7,YEO7_v(i))));
    end
    
    mat_EI_EZPZ = nan(7,7);
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    mat_EI_EZPZ(row_idx,:) = EI_EZPZ_mat;
    
    TII_EZ_mat = nan(1,length(YEO7_v));
    for i = 1:length(YEO7_v)
        TII_EZ{i} = data_TII(idx_EZ & ismember(idx_YEO7,YEO7_v(i)));
        TII_EZ_mat(i) = nanmean(data_TII(idx_EZ & ismember(idx_YEO7,YEO7_v(i))));
    end
    
    mat_TII_EZ = nan(7,7);
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    mat_TII_EZ(row_idx,:) = TII_EZ_mat;
    
    TII_PZ_mat = nan(1,length(YEO7_v));
    for i = 1:length(YEO7_v)
        TII_PZ{i} = data_TII(idx_PZ & ismember(idx_YEO7,YEO7_v(i)));
        TII_PZ_mat(i) = nanmean(data_TII(idx_PZ & ismember(idx_YEO7,YEO7_v(i))));
    end
    
    mat_TII_PZ = nan(7,7);
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    mat_TII_PZ(row_idx,:) = TII_PZ_mat;
    
    TII_EZPZ_mat = nan(1,length(YEO7_v));
    for i = 1:length(YEO7_v)
        TII_EZPZ{i} = data_TII((idx_EZ|idx_PZ) & ismember(idx_YEO7,YEO7_v(i)));
        TII_EZPZ_mat(i) = nanmean(data_TII((idx_EZ|idx_PZ) & ismember(idx_YEO7,YEO7_v(i))));
    end
    
    mat_TII_EZPZ = nan(7,7);
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    mat_TII_EZPZ(row_idx,:) = TII_EZPZ_mat;
    
    ER_idx_EZ_mat = nan(1,length(YEO7_v));
    for i = 1:length(YEO7_v)
        ER_idx_EZ{i} = data_ER_idx(idx_EZ & ismember(idx_YEO7,YEO7_v(i)));
        ER_idx_EZ_mat(i) = nanmean(data_ER_idx(idx_EZ & ismember(idx_YEO7,YEO7_v(i))));
    end
    
    mat_ER_idx_EZ = nan(7,7);
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    mat_ER_idx_EZ(row_idx,:) = ER_idx_EZ_mat;
    
    ER_idx_PZ_mat = nan(1,length(YEO7_v));
    for i = 1:length(YEO7_v)
        ER_idx_PZ{i} = data_ER_idx(idx_PZ & ismember(idx_YEO7,YEO7_v(i)));
        ER_idx_PZ_mat(i) = nanmean(data_ER_idx(idx_PZ & ismember(idx_YEO7,YEO7_v(i))));
    end
    
    mat_ER_idx_PZ = nan(7,7);
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    mat_ER_idx_PZ(row_idx,:) = ER_idx_PZ_mat;
    
    ER_idx_EZPZ_mat = nan(1,length(YEO7_v));
    for i = 1:length(YEO7_v)
        ER_idx_EZPZ{i} = data_ER_idx((idx_EZ|idx_PZ) & ismember(idx_YEO7,YEO7_v(i)));
        ER_idx_EZPZ_mat(i) = nanmean(data_ER_idx((idx_EZ|idx_PZ) & ismember(idx_YEO7,YEO7_v(i))));
    end
    
    mat_ER_idx_EZPZ = nan(7,7);
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    mat_ER_idx_EZPZ(row_idx,:) = ER_idx_EZPZ_mat;
    
    for i = 1:length(YEO7_v)
        Eu_dis_EZ{i} = data_Eu_dis(idx_EZ & ismember(idx_YEO7,YEO7_v(i)));
    end
    
    for i = 1:length(YEO7_v)
        Eu_dis_PZ{i} = data_Eu_dis(idx_PZ & ismember(idx_YEO7,YEO7_v(i)));
    end
    
    for i = 1:length(YEO7_v)
        Eu_dis_EZPZ{i} = data_Eu_dis((idx_EZ|idx_PZ) & ismember(idx_YEO7,YEO7_v(i)));
    end
    
    data_stats = cat(1,EI_EZ,EI_PZ,EI_EZPZ,TII_EZ,TII_PZ,TII_EZPZ,ER_idx_EZ,ER_idx_PZ,ER_idx_EZPZ,Eu_dis_EZ,Eu_dis_PZ,Eu_dis_EZPZ);
    
    T = cell2table(data_stats,...
        'VariableNames',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'},...
        'RowNames',{'EI_EZ','EI_PZ','EI_EZPZ','TII_EZ','TII_PZ','TII_EZPZ','ER_idx_EZ','ER_idx_PZ','ER_idx_EZPZ','Eu_dis_EZ','Eu_dis_PZ','Eu_dis_EZPZ'}');
    
    %% figure1 EI
    figure(1)
    subplot(3,3,1)
    Spider_plot(EI_EZ_mat,'AxesLabels', { 'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'},...
        'AxesLimits', [zeros(1,7); ones(1,7)],...
        'AxesInterval', 2,...
        'FillOption', 'on',...
        'LineStyle', {'none'},...
        'FillTransparency', 0.1);
    legend('EI EZ')
    
    subplot(3,3,2)
    Spider_plot(EI_PZ_mat,'AxesLabels', { 'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'},...
        'AxesLimits', [zeros(1,7); ones(1,7)],...
        'AxesInterval', 2,...
        'FillOption', 'on',...
        'LineStyle', {'none'},...
        'FillTransparency', 0.1);
    legend('EI PZ')
    
    subplot(3,3,3)
    Spider_plot(EI_EZPZ_mat,'AxesLabels', { 'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'},...
        'AxesLimits', [zeros(1,7); ones(1,7)],...
        'AxesInterval', 2,...
        'FillOption', 'on',...
        'LineStyle', {'none'},...
        'FillTransparency', 0.1);
    legend('EI EZPZ')
    
    
    subplot(3,3,4)
    imagesc(mat_EI_EZ)
    colormap(gca,'parula');
    colorbar() ; % Add color bar and make sure the color ranges from 0:1
    caxis([-1,1]);
    set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
    set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
    set(gca,'XTickLabelRotation',45)
    
    subplot(3,3,5)
    imagesc(mat_EI_PZ)
    colormap(gca,'parula');
    colorbar() ; % Add color bar and make sure the color ranges from 0:1
    caxis([-1,1]);
    set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
    set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
    set(gca,'XTickLabelRotation',45)
    
    subplot(3,3,6)
    imagesc(mat_EI_EZPZ)
    colormap(gca,'parula');
    colorbar() ; % Add color bar and make sure the color ranges from 0:1
    caxis([-1,1]);
    set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
    set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
    set(gca,'XTickLabelRotation',45)
    
    data2 = cell2mat(T{'EI_EZ',:}');
    data1 = cell2mat(T{'Eu_dis_EZ',:}');
    subplot(3,3,7)
    scatter(data1, data2, '+', 'MarkerFaceColor', 'k');
    ylabel('EI EZ');
    xlabel ('Eu dis EZ');
    box 'on'
    axis square;
    r = corrcoef(data1, data2);
    %     disp(r(1,2));
    tmp=corrcoef(data1,data2);
    str=sprintf('r= %1.2f',tmp(1,2));
    Tx = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(Tx, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    data2 = cell2mat(T{'EI_PZ',:}');
    data1 = cell2mat(T{'Eu_dis_PZ',:}');
    subplot(3,3,8)
    scatter(data1, data2, '+', 'MarkerFaceColor', 'k');
    ylabel('EI PZ');
    xlabel ('Eu dis PZ');
    box 'on'
    axis square;
    r = corrcoef(data1, data2);
    %     disp(r(1,2));
    tmp=corrcoef(data1,data2);
    str=sprintf('r= %1.2f',tmp(1,2));
    Tx = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(Tx, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    data2 = cell2mat(T{'EI_EZPZ',:}');
    data1 = cell2mat(T{'Eu_dis_EZPZ',:}');
    subplot(3,3,9)
    scatter(data1, data2, '+', 'MarkerFaceColor', 'k');
    ylabel('EI EZPZ');
    xlabel ('Eu dis EZPZ');
    box 'on'
    axis square;
    r = corrcoef(data1, data2);
    %     disp(r(1,2));
    tmp=corrcoef(data1,data2);
    str=sprintf('r= %1.2f',tmp(1,2));
    Tx = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(Tx, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    %% figure TII
    figure(2)
    subplot(3,3,1)
    Spider_plot(TII_EZ_mat,'AxesLabels', { 'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'},...
        'AxesLimits', [zeros(1,7); ones(1,7)],...
        'AxesInterval', 2,...
        'FillOption', 'on',...
        'LineStyle', {'none'},...
        'FillTransparency', 0.1);
    legend('TII EZ')
    
    subplot(3,3,2)
    Spider_plot(TII_PZ_mat,'AxesLabels', { 'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'},...
        'AxesLimits', [zeros(1,7); ones(1,7)],...
        'AxesInterval', 2,...
        'FillOption', 'on',...
        'LineStyle', {'none'},...
        'FillTransparency', 0.1);
    legend('TII PZ')
    
    subplot(3,3,3)
    Spider_plot(TII_EZPZ_mat,'AxesLabels', { 'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'},...
        'AxesLimits', [zeros(1,7); ones(1,7)],...
        'AxesInterval', 2,...
        'FillOption', 'on',...
        'LineStyle', {'none'},...
        'FillTransparency', 0.1);
    legend('TII EZPZ')
    
    
    subplot(3,3,4)
    imagesc(mat_TII_EZ)
    colormap(gca,'parula');
    colorbar() ; % Add color bar and make sure the color ranges from 0:1
    caxis([-1,1]);
    set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
    set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
    set(gca,'XTickLabelRotation',45)
    
    subplot(3,3,5)
    imagesc(mat_TII_PZ)
    colormap(gca,'parula');
    colorbar() ; % Add color bar and make sure the color ranges from 0:1
    caxis([-1,1]);
    set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
    set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
    set(gca,'XTickLabelRotation',45)
    
    subplot(3,3,6)
    imagesc(mat_TII_EZPZ)
    colormap(gca,'parula');
    colorbar() ; % Add color bar and make sure the color ranges from 0:1
    caxis([-1,1]);
    set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
    set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
    set(gca,'XTickLabelRotation',45)
    
    data2 = cell2mat(T{'TII_EZ',:}');
    data1 = cell2mat(T{'Eu_dis_EZ',:}');
    subplot(3,3,7)
    scatter(data1, data2, '+', 'MarkerFaceColor', 'k');
    ylabel('TII EZ');
    xlabel ('Eu dis EZ');
    box 'on'
    axis square;
    r = corrcoef(data1, data2);
    %     disp(r(1,2));
    tmp=corrcoef(data1,data2);
    str=sprintf('r= %1.2f',tmp(1,2));
    Tx = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(Tx, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    data2 = cell2mat(T{'TII_PZ',:}');
    data1 = cell2mat(T{'Eu_dis_PZ',:}');
    subplot(3,3,8)
    scatter(data1, data2, '+', 'MarkerFaceColor', 'k');
    ylabel('TII PZ');
    xlabel ('Eu dis PZ');
    box 'on'
    axis square;
    r = corrcoef(data1, data2);
    %     disp(r(1,2));
    tmp=corrcoef(data1,data2);
    str=sprintf('r= %1.2f',tmp(1,2));
    Tx = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(Tx, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    data2 = cell2mat(T{'TII_EZPZ',:}');
    data1 = cell2mat(T{'Eu_dis_EZPZ',:}');
    subplot(3,3,9)
    scatter(data1, data2, '+', 'MarkerFaceColor', 'k');
    ylabel('TII EZPZ');
    xlabel ('Eu dis EZPZ');
    box 'on'
    axis square;
    r = corrcoef(data1, data2);
    %     disp(r(1,2));
    tmp=corrcoef(data1,data2);
    str=sprintf('r= %1.2f',tmp(1,2));
    Tx = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(Tx, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    %% figure ER idx
    figure(3)
    subplot(3,3,1)
    Spider_plot(ER_idx_EZ_mat,'AxesLabels', { 'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'},...
        'AxesLimits', [zeros(1,7); ones(1,7)],...
        'AxesInterval', 2,...
        'FillOption', 'on',...
        'LineStyle', {'none'},...
        'FillTransparency', 0.1);
    legend('ER idx EZ')
    
    subplot(3,3,2)
    Spider_plot(ER_idx_PZ_mat,'AxesLabels', { 'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'},...
        'AxesLimits', [zeros(1,7); ones(1,7)],...
        'AxesInterval', 2,...
        'FillOption', 'on',...
        'LineStyle', {'none'},...
        'FillTransparency', 0.1);
    legend('ER idx PZ')
    
    subplot(3,3,3)
    Spider_plot(ER_idx_EZPZ_mat,'AxesLabels', { 'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'},...
        'AxesLimits', [zeros(1,7); ones(1,7)],...
        'AxesInterval', 2,...
        'FillOption', 'on',...
        'LineStyle', {'none'},...
        'FillTransparency', 0.1);
    legend('ER idx EZPZ')
    
    
    subplot(3,3,4)
    imagesc(mat_ER_idx_EZ)
    colormap(gca,'parula');
    colorbar() ; % Add color bar and make sure the color ranges from 0:1
    caxis([-1,1]);
    set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
    set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
    set(gca,'XTickLabelRotation',45)
    
    subplot(3,3,5)
    imagesc(mat_ER_idx_PZ)
    colormap(gca,'parula');
    colorbar() ; % Add color bar and make sure the color ranges from 0:1
    caxis([-1,1]);
    set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
    set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
    set(gca,'XTickLabelRotation',45)
    
    subplot(3,3,6)
    imagesc(mat_ER_idx_EZPZ)
    colormap(gca,'parula');
    colorbar() ; % Add color bar and make sure the color ranges from 0:1
    caxis([-1,1]);
    set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
    set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
    set(gca,'XTickLabelRotation',45)
    
    data2 = cell2mat(T{'ER_idx_EZ',:}');
    data1 = cell2mat(T{'Eu_dis_EZ',:}');
    subplot(3,3,7)
    scatter(data1, data2, '+', 'MarkerFaceColor', 'k');
    ylabel('ER idx EZ');
    xlabel ('Eu dis EZ');
    box 'on'
    axis square;
    r = corrcoef(data1, data2);
    %     disp(r(1,2));
    tmp=corrcoef(data1,data2);
    str=sprintf('r= %1.2f',tmp(1,2));
    Tx = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(Tx, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    data2 = cell2mat(T{'ER_idx_PZ',:}');
    data1 = cell2mat(T{'Eu_dis_PZ',:}');
    subplot(3,3,8)
    scatter(data1, data2, '+', 'MarkerFaceColor', 'k');
    xlabel('ER idx PZ');
    ylabel ('Eu dis PZ');
    box 'on'
    axis square;
    r = corrcoef(data1, data2);
    %     disp(r(1,2));
    tmp=corrcoef(data1,data2);
    str=sprintf('r= %1.2f',tmp(1,2));
    Tx = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(Tx, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    data2 = cell2mat(T{'ER_idx_EZPZ',:}');
    data1 = cell2mat(T{'Eu_dis_EZPZ',:}');
    subplot(3,3,9)
    scatter(data1, data2, '+', 'MarkerFaceColor', 'k');
    ylabel('ER idx EZPZ');
    xlabel ('Eu dis EZPZ');
    box 'on'
    axis square;
    r = corrcoef(data1, data2);
    %     disp(r(1,2));
    tmp=corrcoef(data1,data2);
    str=sprintf('r= %1.2f',tmp(1,2));
    Tx = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(Tx, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    
    %% save data
    sz_loc_EI.mat_EI_EZ = mat_EI_EZ;
    sz_loc_EI.mat_EI_PZ = mat_EI_PZ;
    sz_loc_EI.mat_EI_EZPZ = mat_EI_EZPZ;
    sz_loc_EI.mat_TII_EZ = mat_EI_EZ;
    sz_loc_EI.mat_TII_PZ = mat_EI_PZ;
    sz_loc_EI.mat_TII_EZPZ = mat_EI_EZPZ;
    sz_loc_EI.T = T;
    sz_loc_EI.TII_seq = TII_seq;
    
    
    
    
    if exist([dirs.original_data filesep sbj_name filesep 'sz_loc_EI_' [sbj_name,'_',bn] '.mat'], 'file')
        prompt = ['sz_loc_EI already exist for ' [sbj_name,'_',bn] ' . Replace it? (y or n):'] ;
        ID = input(prompt,'s');
        if strcmp(ID, 'y')
            save([dirs.original_data filesep sbj_name filesep 'sz_loc_EI_' [sbj_name,'_',bn] '.mat'], 'sz_loc_EI')
            disp(['sz_loc_EI saved for ' [sbj_name,'_',bn]])
        else
            warning(['sz_loc_EI NOT saved for ' [sbj_name,'_',bn]])
        end
    else
        save([dirs.original_data filesep sbj_name filesep 'sz_loc_EI_' [sbj_name,'_',bn] '.mat'], 'sz_loc_EI')
        disp(['sz_loc_EI saved for ' [sbj_name,'_',bn]])
    end
    
    
    
end
end

