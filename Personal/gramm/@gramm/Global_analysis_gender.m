%get ready with the addpath,plz adjust accordingly
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_personal-master/'))
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_preproc-master/'))
%addpath(genpath('/Users/chao/Desktop/function_tools/gramm-master'))% this is a matlab based graph toolbox which is similar to the ggplot2
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%home

% define the cohort
sbj_names_all = {'C17_20';'C17_21';'C18_22';'C18_23';'C18_24';'C18_25';'C18_26';'C18_27';'C18_28';'C18_29';'C18_30'...
    ;'C18_31';'C18_32';'C18_33';'C18_34';'C18_35';'C18_37';'C18_38';'C18_39';'C18_40';'C18_41';'C18_42';'C18_43';'C18_44'...
    ;'C18_45';'C18_46';'C18_47';'C18_49';'C19_50';'C19_51';'C19_52';'C19_53';'C19_55';'C19_58';'C19_60';'C19_62';'S17_114_EB'...
    ;'S17_116_AA';'S17_118_TW';'S20_148_SM';'S20_149_DR';'S20_150_CM';'S20_152_HT'};
sbj_genders = {'M';'F';'F';'F';'M';'F';'F';'F';'M';'F';'F';'F';'F';'F';'F';'F';'F';'F';'M';'F';'F';'F';'M';'M';'M';'F';'F';...
    'M';'F';'F';'M';'F';'F';'M';'F';'M';'F';'M';'F';'F';'M';'F';'F'};

%make a specific selection of cohort
sbj_names = sbj_names_all(1:36);%China
% sbj_names = sbj_names_all(37:end);%Stanford
% sbj_names = sbj_names_all;%all

indx_female = strcmp(sbj_genders(1:36),'F');%this part is the gender of the patient
indx_male  = strcmp(sbj_genders(1:36),'M');%this part is the gender of the patient

sbj_names=sbj_names(indx_female)%this part is the gender of the patient
sbj_names=sbj_names(indx_male)%this part is the gender of the patient

sbj_names = {'C18_49'};

%sbj_names = sbj_names_all([1:36,40,43]);% asian
%sbj_names = sbj_names_all([37:39,41,42]);% white

% define the the abbreviations of kinds of brian structures
anat_all = {'SFG','SFS','MFG','IFS','IFG','OFC','MPG','SMA','VMPFC','ACC','MCC','PCC','STG','STS','MTG','ITS','ITG','AMY','HIPPO A','HIPPO M','HIPPO P'...
    ,'TP','OTS','FG','CS','PHG','PRECENTRAL G','POSTCENTRAL G','SPL','IPL','IPS','PCG','CG','POF','CF','LG','SOG','MOG','IOG','WM','OUT','EC',...
    'FOP','POP','TOP','EMPTY','PARL','LESION','INSULA','BASAL'};

%make a specific selection of anatomical structures
anat = {'INSULA'};anat_name = 'INSULA';
anat = {'ACC','MCC'};anat_name = 'ACC MCC';
anat = {'ACC'};anat_name = 'ACC';
anat = {'MCC'};anat_name = 'MCC';
anat = {'FG'};anat_name = 'FG';
anat = {'FG','OTS','CS'};anat_name = 'FG OTS CS';
anat = {'IFS','IFG'};anat_name = 'IFS IFG';%
anat = {'SFS','SFG'};anat_name = 'SFS SFG';
anat = {'OFC'};anat_name = 'OFC';
anat = {'AMY'};anat_name = 'amygdala';
anat = {'SMA','PARL'};anat_name = 'SMA and PARL';
anat = {'SOG','MOG','IOG'}; anat_name= 'occipital';
anat = {'FOP','POP','TOP'}; anat_name= 'operculum';
anat = {'PCG','CG','POF','CF','PCC'}; anat_name='posterior medial';
anat = {'SPL','IPL'};anat_name='parieltal area';
anat = {'POSTCENTRAL G'};anat_name='central area';
anat = {'PRECENTRAL G'};anat_name='central area';
anat = {'STG','STS','MTG','ITS','ITG'};anat_name='temporal area lateral';
anat = {'PCG'};anat_name='central area';

side = 'none';%'L','R','none'


anat = {'STS'};anat_name='STS';
% Check on the meaning of the abbreviations
anat_displ =  importdata('/Users/chao/Documents/Stanford/code/lbcn_personal-master/Chao/anat_abbreviation.txt');%pls select a directory to store the 
disp(anat_displ);


%% Visit each excel table, add a name column, and concatenate them into a cell
T = cell(size(sbj_names,1), 1);
for i = 1:length(sbj_names)
    %cd(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/' sbj_names{i}])%plz adjust accordingly to your ecosystem
    T{i} = readtable(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/' sbj_names{i} '/' sbj_names{i} '_stats.xlsx']);
    sbj_name_channame = cell(size(T{i},1),1);
    for j = 1:size(T{i},1)
        sbj_name_channame{j} = [sbj_names{i},'-',T{i}.glv_channame{j}];
    end
    T{i}.sbj_name_channame = sbj_name_channame;
end
%Creat another table with rows of specific cohorts and column of specific anatomical
%structures
sz = [size(sbj_names,1) size(anat,2)];
varTypes = cell(1,size(anat,2));
varTypes(:) = {'cell'};
T2 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',anat,'RowNames',sbj_names);
%put the glv_index into each space of the table as a vector
if isempty(side)||strcmp(side,'none')
    for i = 1:length(sbj_names)
        for j = 1:length(anat)
            idx1 = strcmp(T{i}.label,anat{j});
            idx2 = T{i}.any_activation;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
            idx = idx1 & idx2;
            T2{sbj_names{i},anat{j}} = {T{i}.glv_index(idx)'};
        end
    end
else
    for i = 1:length(sbj_names)
        for j = 1:length(anat)
            idx1 = strcmp(T{i}.label,anat{j});
            idx2 = T{i}.any_activation;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
            idx3 = strcmp(T{i}.LvsR,side);
            idx = idx1 & idx2 & idx3;
            T2{sbj_names{i},anat{j}} = {T{i}.glv_index(idx)'};
        end
    end
end
%Since there may be empty electrodes in the stats sheet, Glv_index may be str. Here, all Glv_index in T2 will be transformed into vetor
for i = 1:length(sbj_names)
    for j = 1:length(anat)
        if iscell(T2{sbj_names{i},anat{j}}{:})
            T2{sbj_names{i},anat{j}}{:} = str2double(T2{sbj_names{i},anat{j}}{:});
        end
    end
end
%Creat a third table that horizontally concatenate all the specific
%anatomical structures in together and get rid of the empty rows in T3
sz = [size(sbj_names,1) 1];
varTypes = cell(1,1);
varTypes(:) = {'cell'};
T3 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',{'anat'},'RowNames',sbj_names);
for i = 1:size(T3,1)
    T3{i,:}{:} = horzcat(T2{i,:}{:});
end
loc=cellfun('isempty', T3{:,'anat'} );% 
T3(loc,:)=[];
%%
%define the plot and stats parameters first
project_name ='race_encoding_simple'% 'race_encoding_simple'or'Grad_CPT'
plot_params = genPlotParams(project_name,'timecourse');
plot_params.single_trial_replot = true;
plot_params.single_trial_thr = 20;%the threshold of HFB it could be like 10 15 20 ...
stats_params = genStatsParams(project_name);
%%
%make a specific selection of conditions and coloring

conditions = {'female','male'}; column = 'condNames3';% this part is the gender of the test
plot_params.col = [cdcol.raspberry_red;%this part is the gender of the test
cdcol.blue];%this part is the gender of the test


 conditions = {'asian','black','white'}; column = 'condNames';

%%
%concatenate  data for conditions
plot_data = cell(1,length(conditions));
plot_data_all = cell(1,length(conditions));
stats_data = cell(1,length(conditions));
stats_data_all = cell(1,length(conditions));
for i = 1:length(T3.Properties.RowNames)
    if ~isempty(T3.anat{i})
        indx = i;
        sbj_name = T3.Properties.RowNames{indx};%set basic pipeline parameters
        if contains(sbj_name,'C17')||contains(sbj_name,'C18')||contains(sbj_name,'C19')
            center = 'China';
        else
            center = 'Stanford';
        end
        block_names = BlockBySubj(sbj_name,project_name);
        dirs = InitializeDirs(project_name, sbj_name, comp_root, server_root, code_root);
        load([dirs.data_root,filesep,'OriginalData',filesep,sbj_name,filesep,'global_',project_name,'_',sbj_name,'_',block_names{1},'.mat'])
        for j = 1:length(T3.anat{i})
            data_all = concatBlocks(sbj_name, block_names,dirs,T3.anat{i}(j),'HFB','Band',{'wave'},['stimlock_bl_corr']);%'stimlock_bl_corr'
            
            %smooth is in the trial level
            winSize = floor(data_all.fsample*plot_params.sm);%this part is smooth1
            gusWin = gausswin(winSize)/sum(gausswin(winSize));%this part is smooth2
            data_all.wave_sm = convn(data_all.wave,gusWin','same');%this part is smooth3
            
            [grouped_trials_all,~] = groupConds(conditions,data_all.trialinfo,column,'none',[],false);
            [grouped_trials,cond_names] = groupConds(conditions,data_all.trialinfo,column,stats_params.noise_method,stats_params.noise_fields_trials,false);
            % this part is to exclude HFB over a fixed threshold
            if plot_params.single_trial_replot
                thr_raw =[];
                for di = 1:size(data_all.wave,1)
                    if ~isempty(find(data_all.wave(di,:)>=plot_params.single_trial_thr));
                        fprintf('You have deleted the data over threshold %d from the data \n',plot_params.single_trial_thr);
                    else
                    end
                end
                [thr_raw,thr_column] = find(data_all.wave >= plot_params.single_trial_thr);
                thr_raw = unique(thr_raw);
            end

                for ci = 1:length(conditions)
                    grouped_trials{ci} = setdiff(grouped_trials{ci},thr_raw);% 
                    plot_data{ci} = [plot_data{ci};nanmean(data_all.wave_sm(grouped_trials{ci},:),1)];
                    plot_data_all{ci} = [plot_data_all{ci};nanmean(data_all.wave_sm(grouped_trials_all{ci},:),1)];
                    stats_data{ci} = [plot_data{ci};nanmean(data_all.wave(grouped_trials{ci},:),1)]; % this part of data were prepared for further comparison (original non-smoothed data)
                    stats_data_all{ci} = [plot_data_all{ci};nanmean(data_all.wave(grouped_trials_all{ci},:),1)];
                end
        end
    else
    end
end
% randomly pick half of the other_races with fixed rng
% if strcmp(column,'condNames9')
%     rng('default');
%     half_index = randsample(size(plot_data{2},1),round(size(plot_data{2},1)/2));
%     plot_data{2} = plot_data{2}(half_index,:);
% else
% end
%% plot figure based on aboving data
clear h
load('cdcol.mat')
figureDim = [100 100 .23 .35 ];
figure('units', 'normalized', 'outerposition', figureDim)
for ci = 1:length(conditions)
    lineprops.col{1} = plot_params.col(ci,:);
%     if plot_params.single_trial_replot
%         for di = 1:size(plot_data{ci},1)
%             if ~isempty(find(plot_data{ci}(di,:)>=plot_params.single_trial_thr))
%                 fprintf('You have deleted the data over threshold %d from the condition %d \n;',plot_params.single_trial_thr,ci);
%             else
%             end
%         end
%         [thr_raw,thr_column] = find(plot_data{ci} >= plot_params.single_trial_thr);
%         thr_raw = unique(thr_raw);
%         plot_data{ci}(thr_raw,:)=[];
        %plot_data_stac{ci} = plot_data{ci}(:,(find(abs(data.time-plot_params.clust_per_win(1))<0.001):find(abs(data.time-plot_params.clust_per_win(2))<0.001)));%%%%%Chao the time window to do the permutation
%       cond_names{ci} = [cond_names{ci},' (',num2str(size(plot_data{ci},1)),' of ',num2str(size(plot_data_all{ci},1)), ' trials)'];
%     end
    if ~strcmp(plot_params.eb,'none')
        lineprops.style= '-';
        lineprops.width = plot_params.lw;
        lineprops.edgestyle = '-';
        if strcmp(plot_params.eb,'ste')
            mseb(data_all.time,nanmean(plot_data{ci}),nanstd(plot_data{ci})/sqrt(size(plot_data{ci},1)),lineprops,1);
            hold on
        else %'std'
            mseb(data_all.time,nanmean(plot_data{ci}),nanstd(plot_data{ci}),lineprops,1);
            hold on
        end
    else
    end
    h(ci)=plot(data_all.time,nanmean(plot_data{ci}),'LineWidth',plot_params.lw,'Color',plot_params.col(ci,:));
    hold on
end
xlim(plot_params.xlim)


if contains(sbj_names{1},'C17_20') && contains(sbj_names{end},'S20_152_HT')
    xlabel(plot_params.xlabel);
else
end


if strcmp(anat_name,'INSULA')
    ylabel(plot_params.ylabel);
else
end

%xlabel(plot_params.xlabel);
%ylabel(plot_params.ylabel);
set(gca,'fontsize',plot_params.textsize)
box off

% if strcmp(anat_name, 'INSULA')
%     y_lim = [-0.2 0.7];
% elseif strcmp(anat_name, 'ACC MCC')
%     y_lim = [-.1 .7];
% elseif strcmp(anat_name, 'FG')
%     y_lim = [-.5 4.2];
% elseif strcmp(anat_name, 'IFS IFG')
%     y_lim = [-.2 2];
% elseif strcmp(anat_name, 'SFS SFG')
%     y_lim = [-.2 .9];
% elseif strcmp(anat_name, 'OFC')
%     y_lim = [-.4 1.5];
% elseif strcmp(anat_name, 'amygdala')
%     y_lim = [-.2 0.5];
% else
    y_lim = ylim;
% end

if size(data_all.trialinfo.allonsets,2) > 1
    time_events = cumsum(nanmean(diff(data_all.trialinfo.allonsets,1,2)));
    for i = 1:length(time_events)
        plot([time_events(i) time_events(i)],y_lim,'Color', [.5 .5 .5], 'LineWidth',1)
    end
else
end
plot([0 0],y_lim, 'Color', [0 0 0], 'LineWidth',2)
plot(xlim,[0 0], 'Color', [.5 .5 .5], 'LineWidth',1)
ylim(y_lim)
box on 
leg = legend(h,cond_names,'Location','Northeast', 'AutoUpdate','off');%cond_names has the trial infomation(default), and cond_names2 is about the category
legend boxoff
set(leg,'fontsize',plot_params.legendfontsize, 'Interpreter', 'none')

%legend off
% set(gca,'XLabel','Time(S)');%chao
xlabel('Time(S)') 

sites_num = sum(cellfun(@numel, T3{:,'anat'} ));
sbj_names_num = size(T3,1);
title([num2str(sites_num),' sites in ' anat_name ' from ',num2str(sbj_names_num),' Subjects'])

%cd('/Users/chao/Desktop/Project_in_Stanford/RACE/4_working_data/globe_analysis_figures');%plz adjust accordingly

%% stats
[H,P,CI,STATS] = ttest2(stats_data{1},stats_data{2}); STATS.H = H; STATS.P = P; STATS.CI = CI;

%% plot the distribution of sites among cases
figureDim = [100 100 .23 .35 ];
figure('units', 'normalized', 'outerposition', figureDim)
%x = 1:sbj_names_num;
x = cell(1,size(T3,1));
for i = 1:size(T3,1)
    x{i} = T3.Row{i};
end
x=categorical(x);
y = cellfun(@numel, T3{:,'anat'} );
barh(x,y)
set(gca,'fontsize',plot_params.textsize)
sbj_names_num = size(T3,1);
title(['distribution of sites in ' anat_name ' from ',num2str(sbj_names_num),' Subjects'])

%cd('/Users/chao/Desktop/Project_in_Stanford/RACE/4_working_data/globe_analysis_figures');%plz adjust accordingly
%% plot the sites in MNI space 
%This part can work but the figure is not pretty, I'm still working on this part, if you know some way that we could plot a nicer brain
%plz tell me know, chao
%PlotCoverageElect
%PlotCoverageElectComb

%define the figure parameters first
cfg = [];
% cfg.highlight_col = [1 0 0];
% cfg.chan_highlight = [];
% cfg.correction_factor = 0;
% cfg.MarkerSize = 20;
% cfg.alpha = 0.3;
% cfg.lobe = 'medial';
% cfg.correction_factor = 0;
cfg.flip_right = true;% this is to flipped all the sites

% get the MNI coords and L/R information from variable T
coords = struct;
coords.MNI_coord = [];
coords.LvsR = [];
coords.isleft = [];
coords.channame = [];

if isempty(side)||strcmp(side,'none')
    for i = 1:length(sbj_names)
        for j = 1:length(anat)
            idx1 = strcmp(T{i}.label,anat{j});
            idx2 = T{i}.any_activation;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
            idx = idx1 & idx2;
            coords_in_T = [T{i}.MNI_coord_1 T{i}.MNI_coord_2 T{i}.MNI_coord_3];
            coords_in_T = coords_in_T(idx,:);
            LvsR_in_T = T{i}.LvsR(idx,:);
            channame_in_T = T{i}.sbj_name_channame(idx,:);
            coords.MNI_coord = [coords.MNI_coord;coords_in_T];
            coords.LvsR = [coords.LvsR;LvsR_in_T];
            coords.channame = [coords.channame; channame_in_T];
        end
    end
else
    for i = 1:length(sbj_names)
        for j = 1:length(anat)
            idx1 = strcmp(T{i}.label,anat{j});
            idx2 = T{i}.any_activation;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
            idx3 = strcmp(T{i}.LvsR,side);
            idx = idx1 & idx2 & idx3;
            coords_in_T = [T{i}.MNI_coord_1 T{i}.MNI_coord_2 T{i}.MNI_coord_3];
            coords_in_T = coords_in_T(idx,:);
            LvsR_in_T = T{i}.LvsR(idx,:);
            channame_in_T = T{i}.sbj_name_channame(idx,:);
            coords.MNI_coord = [coords.MNI_coord;coords_in_T];
            coords.LvsR = [coords.LvsR;LvsR_in_T];
            coords.channame = [coords.channame; channame_in_T];
        end
    end
end


coords.isleft = [];
for i = 1:length(coords.LvsR)
    if strcmp(coords.LvsR{i},'L')
        coords.isleft(i,1) = 1;
    else
        coords.isleft(i,1) = 0;
    end
end


% this is for flip
if cfg.flip_right
for i = 1:length(coords.MNI_coord)
    if coords.MNI_coord(i,1) <= 0
        coords.MNI_coord(i,1) = coords.MNI_coord(i,1)+2*abs(coords.MNI_coord(i,1));
    end
end
else
end
coords.isleft = zeros(length(coords.isleft),1);

% plot the brain and the sites

% PlotCoverageGroup(coords,cfg)
%% plot the sites in MNI space by using iELVis
addpath(genpath('/Users/chao/Desktop/function_tools/for_plot/CCEP_fMRI/'))
addpath(genpath('/Users/chao/Desktop/function_tools/for_plot/iELVis-master/'))
global globalFsDir;
globalFsDir ='/Volumes/CHAO_IRON_M/iELVis_working/Plot_Elelctrodes';
fsDir='/Volumes/CHAO_IRON_M/iELVis_working/Plot_Elelctrodes';
cd([fsDir]);
%plotting=0;
load('DDS_parc_colors.mat');
load('cdcol_stanford.mat');
cdcol.loc_yellow = hex2rgb('#fff301');


% asian race orange #ff7c01     [1 0.4863 0.0039]
% Black race purple #612696     [0.3804 0.1490 0.5882]
% white race green  #02a755     [0.0078 0.6549 0.3333]
% localization yellow #fff301   [1 0.9529 0.0039]




eleColors =[];

for i = 1:size(coords.channame)
%      if contains(coods.channame{i},'C17')|| contains(coods.channame{i},'C18')||contains(coods.channame{i},'C19')
%         eleColors(i,:) = [0.9216,0.5686,0.3059];
% %     elseif flag(i) == 2
% %         eleColors(i,:) = [0.9196,0.1412,0.3451]
%     else
        eleColors(i,:) =[cdcol.loc_yellow];
end

cfg=[];
cfg.view='r';
cfg.elecSize=12;
cfg.surfType='inflated';    
cfg.opaqueness=1;
cfg.ignoreDepthElec='n';
cfg.elecCoord=[coords.MNI_coord coords.isleft];
cfg.elecNames=coords.channame;
cfg.elecColors = eleColors;
cfg.elecColorScale=[0 1];
cfg.title = [];
cfg.backgroundColor = [1,1,1];
cfgOut=plotPialSurf_v2('fsaverage',cfg);




%% Check whether the left and right coordinates are correct in the the subjvar 
for i = 1:length(sbj_names_all)
    load(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/',sbj_names_all{i},'/','subjVar_',sbj_names_all{i},'.mat'])
    elinfo_link = subjVar.elinfo;
    
    for j = 1:size(elinfo_link,1)
        if strcmp(elinfo_link.LvsR{j},'L')&&elinfo_link.MNI_coord(j,1)<0
            disp('good to go')
        elseif strcmp(elinfo_link.LvsR{j},'R')&&elinfo_link.MNI_coord(j,1)>0
            disp('good to go')
        elseif strcmp(elinfo_link.LvsR{j},'L')&&elinfo_link.MNI_coord(j,1)>0
            warning(['please check the subjvar' sbj_names{i}])
            return
        elseif strcmp(elinfo_link.LvsR{j},'R')&& elinfo_link.MNI_coord(j,1)<0
            warning(['please check the subjvar' sbj_names{i}])
            return
        else
            disp('catched empty site')
        end
    end
end

%% pedro method
task = 'MMR'
cond_names = {'math', 'autobio'};
column = 'condNames';
data = concatenate_multiple_elect(elect_list, task, dirs, 'Band', 'HFB', 'stim');
plot_group_elect(data, task, cond_names, column);
%% about the sheet
% take care of the str and num of FS_ind and glv_index
for i = 1:size(T,1)
    if ~iscell(T{i}.FS_ind)
       T{i}.FS_ind = num2cell(T{i}.FS_ind);
       for j = 1:size(T{i},1)
           T{i}.FS_ind{j} = num2str(T{i}.FS_ind{j});
       end
    else
    end 
end
for i = 1:size(T,1)
    if ~iscell(T{i}.glv_index)
       T{i}.glv_index = num2cell(T{i}.glv_index);
       for j = 1:size(T{i},1)
           T{i}.glv_index{j} = num2str(T{i}.glv_index{j});
       end
    else
    end
end

for i = 1:size(T,1)
    if ~iscell(T{i}.label)
       T{i}.label = num2cell(T{i}.label);
       for j = 1:size(T{i},1)
           T{i}.label{j} = num2str(T{i}.label{j});
       end
    else
    end
end
% concatenate table vertically
bigtable = [T{1};T{2}];
for tidx = 3:numel(T)
   bigtable = [bigtable; T{tidx}];
end

%% Pedro's code 

% % Group analyses
% 
% 1. Select which electrodes to include
% -Statistical
%     Compare all trails agains baseline
%     permutation test between avg baseline period whithin trial vs. avg 1s period within trial
% -Anatomical
%     ROIS
% After these steps you should have the electrode_list
% 
% 2. Concatenate all electrodes of interest
chan_num = [];
sbj_name_p = [];
for i = 1:size(T3,1)
   for j = 1:size(T3.anat{i},2)
    a = T3.anat{i}(j);
    chan_num = [chan_num;a];
    b = T3.Row(i);
    sbj_name_p = [sbj_name_p;b];
   end
end

project_name = 'race_encoding_simple';
electrode_list = table(chan_num,sbj_name_p);
dirs = InitializeDirs(project_name, sbj_name_p{1}, comp_root, server_root, code_root);
data_all = concatenate_multiple_elect(electrode_list, project_name, dirs, 'Band', 'HFB', 'stim');
cond_names = {'asian','black','white'};column = 'condNames';
plot_group_elect(data_all, project_name, cond_names, column);
% 3. Compare two conditions across electrodes
% -Define conditions and time window
%     cond_names = {'autobio', 'math'};
%     column = 'condNames';
stats_params = genStatsParams(project_name);
% -Average each electrode per condition within the 1s period
%     (now you have a single value per electrode per condition)
%     Ready to compare electrodes with independent sample t-test
STATS = stats_group_elect(data_all,data_all.time, project_name,cond_names, column, stats_params);

