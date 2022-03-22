% this is the EP network group analysis

%% set roots
addpath(genpath('/Users/tony/Documents/Stanford/Chao_SEEG/'))
addpath(genpath('/Users/tony/Documents/Stanford/DELLO'))
addpath(genpath('/Users/tony/Desktop/function_tools/iELVis/'))
addpath(genpath('/Users/tony/Desktop/function_tools/for_plot/CCEP_fMRI/'))
addpath('/Users/tony/Desktop/function_tools/spm12_7219')
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%home

project_name = 'EPnetwork';

sbj_names = {'2015_FT_hankunpeng';
    '2015_FT_huangzutu';
    '2015_FT_songziwei';
    '2015_FT_wangmingming';
    '2015_FT_wangzheng';
    '2015_FT_xumengjiao';
    '2015_FT_yangbihong';
    '2015_FT_yinfengting';
    '2015_FT_yuanjinrui';
    '2015_FT_zhaoxichun';
    '2017_FT_chengang';
    '2017_FT_genghaowen';
    '2017_FT_guyuxuan';
    '2017_FT_jielei';
    '2017_FT_lijingjing';
    '2017_FT_lilechun';
    '2017_FT_linzixuan';
    '2017_FT_liuboqiang';
    '2017_FT_liuruilin';
    '2017_FT_masinan';
    '2017_FT_panzhiyong';
    '2017_FT_wangmengru';
    '2017_FT_wangyan';
    '2017_FT_xionghuihui';
    '2017_FT_yangchenglei';
    '2017_FT_yangrui';
    '2017_FT_yaodongyuan';
    '2017_FT_yena';
    '2017_FT_yuanye';
    '2017_FT_zhanggenhong'};

dirs = InitializeDirs(project_name, sbj_names{1}, comp_root, server_root, code_root);
work_path = '/Users/tony/Desktop/4_working_data/Group_analysis/data';%%% this folder is for the EI ER TII plot
work_sbjs_path = '/Users/tony/Desktop/4_working_data/Group_analysis/data_sbjs';%%% this folder is for the electrodes pot
%% copy the EI mat file to work folder

if ~exist(work_path)
    mkdir(work_path);
end

for i = 1:length(sbj_names)
    block_names = BlockBySubj(sbj_names{i},project_name);
    for j = 1:length(block_names)
        fn = sprintf('%s/%s/sz_loc_EI_%s_%s.mat',dirs.original_data,sbj_names{i},sbj_names{i},block_names{j});
        copyfile(fn, work_path);
        fprintf('EI files Copied EI mat file %s to %s \n', fn, work_path);
    end
end

if ~exist(work_sbjs_path)
    mkdir(work_sbjs_path);
end

for i = 1:length(sbj_names)
    block_names = BlockBySubj(sbj_names{i},project_name);
    for j = 1
        fn = sprintf('%s/%s/sz_loc_EI_%s_%s.mat',dirs.original_data,sbj_names{i},sbj_names{i},block_names{j});
        copyfile(fn, work_sbjs_path);
        fprintf('subjs electrodes filses Copied EI mat file %s to %s \n', fn, work_path);
    end
end
%% figure plot the electrodes
% We have here the logic of several layers of different plot electrodes
%The first layer is all the electrode contacts of SEEG                                                                  layer1_SEEG_idx;
%The second layer is the electrode contact that belongs to the YEO7 network and is manually removed as the bad site     layer2_YEO7_idx;
%The third layer is the contact of all epilepsy spread (EI>0)                                                           layer3_EP_idx;
%The fourth layer is where all EI = 1 belong                                                                            layer4_EI_idx;

%% figure SEEG and YEO7 electrodes
folders = dir(work_sbjs_path);
folders = string({folders.name});
folders = folders(~startsWith(folders,"."))' ;

elecCoord152_all_layer1_SEEG = [];
elecCoord152_all_layer2_YEO7 = [];
elecCoord152_all_layer3_EP   = [];
elecCoord152_all_layer4_EI   = [];

elecCoord305_all_layer1_SEEG = [];
elecCoord305_all_layer2_YEO7 = [];
elecCoord305_all_layer3_EP   = [];
elecCoord305_all_layer4_EI   = [];

%In the first stage, it is very likely that errors will occur
%when judging left and right in the preprocessing, so it is highly
%recommended to determine the left and right according to the coordinates
%of MNI, left <0 right>0

isleft_all_layer1_SEEG = [];
isleft_all_layer2_YEO7 = [];
isleft_all_layer3_EP   = [];
isleft_all_layer4_EI   = [];

chan_names_all_layer1_SEEG = [];
chan_names_all_layer2_YEO7 = [];
chan_names_all_layer3_EP = [];
chan_names_all_layer4_EI = [];


for i = 1:length(folders) %%% access each mat file
    fn = sprintf('%s/%s',work_sbjs_path,folders(i,1));
    load(fn);
    bad_chan = sz_loc_EI.bad_chan;
    eleinfo = sz_loc_EI.eleinfo;
    %     Nu = sz_loc_EI.Nu;
    %     Lamda = sz_loc_EI.Lamda;
    %     Eta = sz_loc_EI.Eta;
    %     EI_window = sz_loc_EI.EI_window;
    %     mat_EI_EZ = sz_loc_EI.mat_EI_EZ;
    %     mat_EI_PZ = sz_loc_EI.mat_EI_PZ;
    %     mat_EI_EZPZ = sz_loc_EI.mat_EI_EZPZ;
    %     mat_TII_EZ = sz_loc_EI.mat_TII_EZ;
    %     mat_TII_PZ = sz_loc_EI.mat_TII_PZ;
    %     mat_TII_EZPZ = sz_loc_EI.mat_TII_EZPZ;
    %     T = sz_loc_EI.T;
    %     TII_seq = sz_loc_EI.T;
    name = sz_loc_EI.name;
    BR_chan_length = size(eleinfo,1);
    bad_chan_idx = ones(size(eleinfo,1),1);
    bad_chan_idx(bad_chan',1) = 0;
    anat_idx = ~strcmp(eleinfo.AAL3,'NotAvailable');
    YEO7_idx = ~strcmp(eleinfo.YEO7,'NotAvailable');
    
    
    
    layer1_SEEG_idx = ~strcmp(eleinfo.YEO7,'NaN');
    layer2_YEO7_idx = layer1_SEEG_idx & anat_idx & YEO7_idx & bad_chan_idx;
    layer3_EP_idx = layer2_YEO7_idx & ~isnan(cell2mat(eleinfo.EI));
    layer4_EI_idx = layer3_EP_idx & cell2mat(eleinfo.EI)==1;
    
    if sum(layer3_EP_idx) ~= sum(~isnan(cell2mat(eleinfo.EI)))
        warning(['idx cracked',name]);
    elseif sum(layer4_EI_idx) ~= sum(cell2mat(eleinfo.EI)==1)
        warning(['idx cracked',name]);
    else
    end
    
    
    elecCoord152_layer1_SEEG = eleinfo.MNI152_volume(layer1_SEEG_idx,:);
    elecCoord152_all_layer1_SEEG = [elecCoord152_all_layer1_SEEG;elecCoord152_layer1_SEEG];
    elecCoord152_layer2_YEO7 = eleinfo.MNI152_volume(layer2_YEO7_idx,:);
    elecCoord152_all_layer2_YEO7 = [elecCoord152_all_layer2_YEO7;elecCoord152_layer2_YEO7];
    elecCoord152_layer3_EP   = eleinfo.MNI152_volume(layer3_EP_idx,:);
    elecCoord152_all_layer3_EP = [elecCoord152_all_layer3_EP;elecCoord152_layer3_EP];
    elecCoord152_layer4_EI   = eleinfo.MNI152_volume(layer4_EI_idx,:);
    elecCoord152_all_layer4_EI = [elecCoord152_all_layer4_EI;elecCoord152_layer4_EI];
    
    elecCoord305_layer1_SEEG = eleinfo.MNI305_volume(layer1_SEEG_idx,:);
    elecCoord305_all_layer1_SEEG = [elecCoord305_all_layer1_SEEG;elecCoord305_layer1_SEEG];
    elecCoord305_layer2_YEO7 = eleinfo.MNI305_volume(layer2_YEO7_idx,:);
    elecCoord305_all_layer2_YEO7 = [elecCoord305_all_layer2_YEO7;elecCoord305_layer2_YEO7];
    elecCoord305_layer3_EP   = eleinfo.MNI305_volume(layer3_EP_idx,:);
    elecCoord305_all_layer3_EP = [elecCoord305_all_layer3_EP;elecCoord305_layer3_EP];
    elecCoord305_layer4_EI   = eleinfo.MNI305_volume(layer4_EI_idx,:);
    elecCoord305_all_layer4_EI = [elecCoord305_all_layer4_EI;elecCoord305_layer4_EI];
    
    
    
    isleft_vector = ones(size(eleinfo,1),1);
    for j = 1:BR_chan_length
        if eleinfo.MNI152_volume(j,1)<=0
            isleft_vector(j,1) = 1;
        elseif eleinfo.MNI152_volume(j,1)>0
            isleft_vector(j,1) = 0;
        elseif strcmp(eleinfo.bv_hem{j,:},'NaN')
            isleft_vector(j,1) = NaN;
        else
        end
    end
    
    isleft_layer1_SEEG = isleft_vector(layer1_SEEG_idx,:);
    isleft_all_layer1_SEEG = [isleft_all_layer1_SEEG;isleft_layer1_SEEG];
    isleft_layer2_YEO7 = isleft_vector(layer2_YEO7_idx,:);
    isleft_all_layer2_YEO7 = [isleft_all_layer2_YEO7;isleft_layer2_YEO7];
    isleft_layer3_EP   = isleft_vector(layer3_EP_idx,:);
    isleft_all_layer3_EP = [isleft_all_layer3_EP;isleft_layer3_EP];
    isleft_layer4_EI   = isleft_vector(layer4_EI_idx,:);
    isleft_all_layer4_EI = [isleft_all_layer4_EI;isleft_layer4_EI];
    
    
    chan_names = cell(length(eleinfo.bv_names),1);
    for k = 1:length(eleinfo.bv_names)
        chan_names{k} = [name,'_',eleinfo.bv_names{k}] ;
    end
    
    
    chan_names_layer1_SEEG = chan_names(layer1_SEEG_idx);
    chan_names_all_layer1_SEEG = [chan_names_all_layer1_SEEG;chan_names_layer1_SEEG];
    chan_names_layer2_YEO7 = chan_names(layer2_YEO7_idx);
    chan_names_all_layer2_YEO7 = [chan_names_all_layer2_YEO7;chan_names_layer2_YEO7];
    chan_names_layer3_EP = chan_names(layer3_EP_idx);
    chan_names_all_layer3_EP =[chan_names_all_layer3_EP;chan_names_layer3_EP];
    chan_names_layer4_EI = chan_names(layer4_EI_idx);
    chan_names_all_layer4_EI = [chan_names_all_layer4_EI;chan_names_layer4_EI];
    
    
end

%figure part

load /Users/tony/Documents/Stanford/Chao_SEEG/visualization/colormaps/cdcol_2018.mat

global globalFsDir;
globalFsDir ='/Users/tony/Desktop/iELVis/Plot_Elelctrodes/';
fsDir='/Users/tony/Desktop/iELVis/Plot_Elelctrodes/';%%% seems useful
cd([fsDir]);

SEEG_YEO7_idx = ismember(elecCoord305_all_layer1_SEEG,elecCoord305_all_layer2_YEO7,'rows');
if sum(SEEG_YEO7_idx)~= length(elecCoord305_all_layer2_YEO7)
    warning('plz checking the position of the electrodes')
else
end

eleColors_SEEG_YEO7 = ones(length(SEEG_YEO7_idx),3);

for i = 1:length(SEEG_YEO7_idx)
    if SEEG_YEO7_idx(i)
        eleColors_SEEG_YEO7(i,:) = cdcol.gold_cadmium_yellow;
    else
        eleColors_SEEG_YEO7(i,:) = cdcol.grey;
    end
end

cfg=[];
cfg.view='lm';
cfg.elecSize=15;
cfg.surfType='inflated';
cfg.opaqueness=1;
cfg.ignoreDepthElec='n';
cfg.elecNames = chan_names_all_layer1_SEEG;
cfg.elecCoord=[elecCoord305_all_layer1_SEEG isleft_all_layer1_SEEG];
%cfg. ignoreChans = {'PT049-X7'};
cfg.elecColors = eleColors_SEEG_YEO7;
cfg.elecColorScale=[0 1];
cfg.title = [];
cfg.backgroundColor = [1,1,1];
% cfg.overlayParcellation='Y7';
cfgOut = plotPialSurf_v2('fsaverage',cfg);

%% figure EP and EI electrodes

folders = dir(work_path);
folders = string({folders.name});
folders = folders(~startsWith(folders,"."))' ;


data_cells = cell(length(folders),1);
for i = 1:length(folders)
    fn = sprintf('%s/%s',work_path,folders(i,1));
    load(fn);
    bad_chan = sz_loc_EI.bad_chan;
    eleinfo = sz_loc_EI.eleinfo;
    name = sz_loc_EI.name;
    
    BR_chan_length = size(eleinfo,1);
    
    bad_chan_idx = ones(size(eleinfo,1),1);
    bad_chan_idx(bad_chan',1) = 0;
    eleinfo.bad_chan = bad_chan_idx;
    
    sbjs_bv_names = cell(BR_chan_length,1);
    for j = 1:BR_chan_length
        sbjs_bv_names{j} = [name,'_',eleinfo.bv_names{j}];
    end
    
    eleinfo.sbjs_bv_names = sbjs_bv_names;
    eleinfo.YEO7idx_Group = ones(length(eleinfo.bv_names),1).*sz_loc_EI.YEO7idx_Group;
    data_cells{i,1} = eleinfo;
end

G_sheet  = vertcat(data_cells{:});

anat_idx = ~strcmp(G_sheet.AAL3,'NotAvailable');
YEO7_idx = ~strcmp(G_sheet.YEO7,'NotAvailable');

layer1_SEEG_idx = ~strcmp(G_sheet.YEO7,'NaN');
layer2_YEO7_idx = layer1_SEEG_idx & anat_idx & YEO7_idx & G_sheet.bad_chan;
layer3_EP_idx = layer2_YEO7_idx & ~isnan(cell2mat(G_sheet.EI));
layer4_EI_idx = layer3_EP_idx & cell2mat(G_sheet.EI)==1;

EP_coords_305_all = G_sheet.MNI305_volume(layer3_EP_idx,:);
[EP_coords_305_UNQ,InEP_all,InEP_UNQ] = unique(EP_coords_305_all,'rows');
EI_coords_305_all = G_sheet.MNI305_volume(layer4_EI_idx,:);
[EI_coords_305_UNQ,InEI_all,InEI_UNQ] = unique(EI_coords_305_all,'rows');

EP_EI_idx = ismember(EP_coords_305_UNQ,EI_coords_305_UNQ,'rows');
if sum(EP_EI_idx)~= length(EI_coords_305_UNQ)
    warning('plz checking the position of the electrodes')
else
end

chan_names_EP_all = G_sheet.sbjs_bv_names(layer3_EP_idx );
chan_names_EP_UNQ = chan_names_EP_all(InEP_all);
chan_names_EI_all = G_sheet.sbjs_bv_names(layer4_EI_idx );
chan_names_EI_UNQ = chan_names_EI_all(InEI_all);

isleft_EP_UNQ = ones(size(chan_names_EP_UNQ,1),1);
for j = 1:length(isleft_EP_UNQ)
    if EP_coords_305_UNQ(j,1)<=0
        isleft_EP_UNQ(j,1) = 1;
    elseif EP_coords_305_UNQ(j,1)>0
        isleft_EP_UNQ(j,1) = 0;
    else
    end
end




%figure part

load /Users/tony/Documents/Stanford/Chao_SEEG/visualization/colormaps/cdcol_2018.mat

global globalFsDir;
globalFsDir ='/Users/tony/Desktop/iELVis/Plot_Elelctrodes/';
fsDir='/Users/tony/Desktop/iELVis/Plot_Elelctrodes/';%%% seems useful
cd([fsDir]);



eleColors_EP_EI = ones(length(EP_EI_idx),3);

for i = 1:length(EP_EI_idx)
    if EP_EI_idx(i)
        eleColors_EP_EI(i,:) = cdcol.lemon_yellow;
    else
        eleColors_EP_EI(i,:) = cdcol.permanent_blue;
    end
end

cfg=[];
cfg.view='lm';
cfg.elecSize=15;
cfg.surfType='inflated';
cfg.opaqueness=1;
cfg.ignoreDepthElec='n';
cfg.elecNames = chan_names_EP_UNQ;
cfg.elecCoord=[EP_coords_305_UNQ isleft_EP_UNQ];
%cfg. ignoreChans = {'PT049-X7'};
cfg.elecColors = eleColors_EP_EI;
cfg.elecColorScale=[0 1];
cfg.title = [];
cfg.backgroundColor = [1,1,1];
% cfg.overlayParcellation='Y7';
cfgOut = plotPialSurf_v2('fsaverage',cfg);

%% figure matrix

folders = dir(work_path);
folders = string({folders.name});
folders = folders(~startsWith(folders,"."))' ;


data_cells = cell(length(folders),1);
data_EI_cells = cell(length(folders),3);
data_TII_cells = cell(length(folders),3);
data_ER_cells = cell(length(folders),3);
for i = 1:length(folders)
    fn = sprintf('%s/%s',work_path,folders(i,1));
    load(fn);
    bad_chan = sz_loc_EI.bad_chan;
    eleinfo = sz_loc_EI.eleinfo;
    name = sz_loc_EI.name;
    
    BR_chan_length = size(eleinfo,1);
    
    bad_chan_idx = ones(size(eleinfo,1),1);
    bad_chan_idx(bad_chan',1) = 0;
    eleinfo.bad_chan = bad_chan_idx;
    
    sbjs_bv_names = cell(BR_chan_length,1);
    for j = 1:BR_chan_length
        sbjs_bv_names{j} = [name,'_',eleinfo.bv_names{j}];
    end
    
    eleinfo.sbjs_bv_names = sbjs_bv_names;
    eleinfo.YEO7idx_Group = ones(length(eleinfo.bv_names),1).*sz_loc_EI.YEO7idx_Group;
    
    data_cells{i,1} = eleinfo;
    
    
    % we could extract EI mat though the easy way like this
    data_EI_cells{i,1} = sz_loc_EI.mat_EI_EZ;
    data_EI_cells{i,2} = sz_loc_EI.mat_EI_PZ;
    data_EI_cells{i,3} = sz_loc_EI.mat_EI_EZPZ;

    
    T = sz_loc_EI.T;
    
    % ER_idx extraction
    
    YEO7_v = [5,7,1,2,3,4,6];
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    
    mat_ER_idx_EZ = nan(7,7);
    for j=1:length(YEO7_v)
        mat_ER_idx_EZ(row_idx,j) = nanmean(T{'ER_idx_EZ',j}{:});
    end
    data_ER_cells{i,1} = mat_ER_idx_EZ;
    
    mat_ER_idx_PZ = nan(7,7);
    for j=1:length(YEO7_v)
        mat_ER_idx_PZ(row_idx,j) = nanmean(T{'ER_idx_PZ',j}{:});
    end
    data_ER_cells{i,2} = mat_ER_idx_PZ;
    
    mat_ER_idx_EZPZ = nan(7,7);
    for j=1:length(YEO7_v)
        mat_ER_idx_EZPZ(row_idx,j) = nanmean(T{'ER_idx_EZPZ',j}{:});
    end
    data_ER_cells{i,3} = mat_ER_idx_EZPZ;
    
    % TII extraction
    
    data_TII_cells{i,1} = sz_loc_EI.mat_TII_EZ;
    data_TII_cells{i,2} = sz_loc_EI.mat_TII_PZ;
    data_TII_cells{i,3} = sz_loc_EI.mat_TII_EZPZ;
    
    mat_TII_EZ = nan(7,7);
    for j=1:length(YEO7_v)
        mat_TII_EZ(row_idx,j) = nanmean(T{'TII_EZ',j}{:});
    end
    data_TII_cells{i,1} = mat_TII_EZ;
    
    mat_TII_PZ = nan(7,7);
    for j=1:length(YEO7_v)
        mat_TII_PZ(row_idx,j) = nanmean(T{'TII_PZ',j}{:});
    end
    data_TII_cells{i,2} = mat_TII_PZ;
    
    mat_TII_EZPZ = nan(7,7);
    for j=1:length(YEO7_v)
        mat_TII_EZPZ(row_idx,j) = nanmean(T{'TII_EZPZ',j}{:});
    end
    data_TII_cells{i,3} = mat_TII_EZPZ;
    
end

G_sheet  = vertcat(data_cells{:});

EI_EZ_matrix = nanmean(cat(3,data_EI_cells{:,1}),3);
EI_PZ_matrix = nanmean(cat(3,data_EI_cells{:,2}),3);
EI_EZPZ_matrix = nanmean(cat(3,data_EI_cells{:,3}),3);
ER_EZ_matrix = nanmean(cat(3,data_ER_cells{:,1}),3);
ER_PZ_matrix = nanmean(cat(3,data_ER_cells{:,2}),3);
ER_EZPZ_matrix = nanmean(cat(3,data_ER_cells{:,3}),3);
TII_EZ_matrix = nanmean(cat(3,data_TII_cells{:,1}),3);
TII_PZ_matrix = nanmean(cat(3,data_TII_cells{:,2}),3);
TII_EZPZ_matrix = nanmean(cat(3,data_TII_cells{:,3}),3);

figure
imagesc(TII_EZPZ_matrix)
colormap(gca,'parula');
colorbar() ; % Add color bar and make sure the color ranges from 0:1
caxis([-1,1]);
set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
set(gca,'XTickLabelRotation',45)

%% inside and outside the YEO7
folders = dir(work_path);
folders = string({folders.name});
folders = folders(~startsWith(folders,"."))' ;


data_cells = cell(length(folders),1);
data_EI_cells = cell(length(folders),3);
data_TII_cells = cell(length(folders),3);
data_ER_cells = cell(length(folders),3);
for i = 1:length(folders)
    fn = sprintf('%s/%s',work_path,folders(i,1));
    load(fn);
    bad_chan = sz_loc_EI.bad_chan;
    eleinfo = sz_loc_EI.eleinfo;
    name = sz_loc_EI.name;
    
    BR_chan_length = size(eleinfo,1);
    
    bad_chan_idx = ones(size(eleinfo,1),1);
    bad_chan_idx(bad_chan',1) = 0;
    eleinfo.bad_chan = bad_chan_idx;
    
    sbjs_bv_names = cell(BR_chan_length,1);
    for j = 1:BR_chan_length
        sbjs_bv_names{j} = [name,'_',eleinfo.bv_names{j}];
    end
    
    eleinfo.sbjs_bv_names = sbjs_bv_names;
    eleinfo.YEO7idx_Group = ones(length(eleinfo.bv_names),1).*sz_loc_EI.YEO7idx_Group;
    
    data_cells{i,1} = eleinfo;
    
    data_EI_cells{i,1} = sz_loc_EI.mat_EI_EZ;
    data_EI_cells{i,2} = sz_loc_EI.mat_EI_PZ;
    data_EI_cells{i,3} = sz_loc_EI.mat_EI_EZPZ;
    
    
    data_TII_cells{i,1} = sz_loc_EI.mat_TII_EZ;
    data_TII_cells{i,2} = sz_loc_EI.mat_TII_PZ;
    data_TII_cells{i,3} = sz_loc_EI.mat_TII_EZPZ;
    
    
    T = sz_loc_EI.T;
    
    YEO7_v = [5,7,1,2,3,4,6];
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    
    mat_ER_idx_EZ = nan(7,7);
    for j=1:length(YEO7_v)
        mat_ER_idx_EZ(row_idx,j) = nanmean(T{'ER_idx_EZ',j}{:});
    end
    data_ER_cells{i,1} = mat_ER_idx_EZ;
    
    mat_ER_idx_PZ = nan(7,7);
    for j=1:length(YEO7_v)
        mat_ER_idx_PZ(row_idx,j) = nanmean(T{'ER_idx_PZ',j}{:});
    end
    data_ER_cells{i,2} = mat_ER_idx_PZ;
    
    mat_ER_idx_EZPZ = nan(7,7);
    for j=1:length(YEO7_v)
        mat_ER_idx_EZPZ(row_idx,j) = nanmean(T{'ER_idx_EZPZ',j}{:});
    end
    data_ER_cells{i,3} = mat_ER_idx_EZPZ;
    
    
end

G_sheet  = vertcat(data_cells{:});

anat_idx = ~strcmp(G_sheet.AAL3,'NotAvailable');
YEO7_idx = ~strcmp(G_sheet.YEO7,'NotAvailable');

layer1_SEEG_idx = ~strcmp(G_sheet.YEO7,'NaN');
layer2_YEO7_idx = layer1_SEEG_idx & anat_idx & YEO7_idx & G_sheet.bad_chan;
layer3_EP_idx = layer2_YEO7_idx & ~isnan(cell2mat(G_sheet.EI)) & cell2mat(G_sheet.EI)>0;
layer4_EI_idx = layer3_EP_idx & cell2mat(G_sheet.EI)==1;

idx_sameYEO = ones(size(G_sheet,1),1);
for i = 1:size(G_sheet,1)
    if G_sheet.YEO7idx(i)==G_sheet.YEO7idx_Group(i)
        idx_sameYEO(i) = 1;
    else
        idx_sameYEO(i) = 0;
    end
end

EI_vec = cell2mat(G_sheet.EI);
Eu_dis_vec = cell2mat(G_sheet.Eu_dis);

idx1 = layer3_EP_idx&idx_sameYEO;
dataEI1 = EI_vec(idx1);
sum(dataEI1 > 0.3)/length(dataEI1)
datadis1 = Eu_dis_vec(idx1);
idx2 = layer3_EP_idx&~idx_sameYEO;
dataEI2 = EI_vec(idx2);
sum(dataEI2 > 0.3)/length(dataEI2)
datadis2 = Eu_dis_vec(idx2);



data2 = dataEI2;
data1 = datadis2;
figure
scatter(data1, data2, '+', 'MarkerFaceColor', 'k');
ylabel('EI PZ');
xlabel ('Eu dis PZ');
box 'on'
axis square;

%     disp(r(1,2));
tmp=corrcoef(data1,data2);
if size(tmp,2) == 2
    str=sprintf('r= %1.2f',tmp(1,2));
else
    str = 'notavailable';
end
Tx = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(Tx, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%% BOX PLOT OF each network
folders = dir(work_path);
folders = string({folders.name});
folders = folders(~startsWith(folders,"."))' ;


data_cells = cell(length(folders),1);
data_EI_cells = cell(length(folders),3);
data_TII_cells = cell(length(folders),3);
data_ER_cells = cell(length(folders),3);
data_weight_EI_dis = cell(7,2);
for i = 1:length(folders)
    fn = sprintf('%s/%s',work_path,folders(i,1));
    load(fn);
    bad_chan = sz_loc_EI.bad_chan;
    eleinfo = sz_loc_EI.eleinfo;
    name = sz_loc_EI.name;
    
    BR_chan_length = size(eleinfo,1);
    
    bad_chan_idx = ones(size(eleinfo,1),1);
    bad_chan_idx(bad_chan',1) = 0;
    eleinfo.bad_chan = bad_chan_idx;
    
    sbjs_bv_names = cell(BR_chan_length,1);
    for j = 1:BR_chan_length
        sbjs_bv_names{j} = [name,'_',eleinfo.bv_names{j}];
    end
    
    eleinfo.sbjs_bv_names = sbjs_bv_names;
    eleinfo.YEO7idx_Group = ones(length(eleinfo.bv_names),1).*sz_loc_EI.YEO7idx_Group;
    
    data_cells{i,1} = eleinfo;
    
    data_EI_cells{i,1} = sz_loc_EI.mat_EI_EZ;
    data_EI_cells{i,2} = sz_loc_EI.mat_EI_PZ;
    data_EI_cells{i,3} = sz_loc_EI.mat_EI_EZPZ;
    
    
    data_TII_cells{i,1} = sz_loc_EI.mat_TII_EZ;
    data_TII_cells{i,2} = sz_loc_EI.mat_TII_PZ;
    data_TII_cells{i,3} = sz_loc_EI.mat_TII_EZPZ;
    
    
    T = sz_loc_EI.T;
    
    YEO7_v = [5,7,1,2,3,4,6];
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    
    mat_ER_idx_EZ = nan(7,7);
    for j=1:length(YEO7_v)
        mat_ER_idx_EZ(row_idx,j) = nanmean(T{'ER_idx_EZ',j}{:});
    end
    data_ER_cells{i,1} = mat_ER_idx_EZ;
    
    mat_ER_idx_PZ = nan(7,7);
    for j=1:length(YEO7_v)
        mat_ER_idx_PZ(row_idx,j) = nanmean(T{'ER_idx_PZ',j}{:});
    end
    data_ER_cells{i,2} = mat_ER_idx_PZ;
    
    mat_ER_idx_EZPZ = nan(7,7);
    for j=1:length(YEO7_v)
        mat_ER_idx_EZPZ(row_idx,j) = nanmean(T{'ER_idx_EZPZ',j}{:});
    end
    data_ER_cells{i,3} = mat_ER_idx_EZPZ;
    
    
end

G_sheet  = vertcat(data_cells{:});

anat_idx = ~strcmp(G_sheet.AAL3,'NotAvailable');
YEO7_idx = ~strcmp(G_sheet.YEO7,'NotAvailable');

layer1_SEEG_idx = ~strcmp(G_sheet.YEO7,'NaN');
layer2_YEO7_idx = layer1_SEEG_idx & anat_idx & YEO7_idx & G_sheet.bad_chan;
layer3_EP_idx = layer2_YEO7_idx & ~isnan(cell2mat(G_sheet.EI)) & cell2mat(G_sheet.EI)>0;
layer4_EI_idx = layer3_EP_idx & cell2mat(G_sheet.EI)==1;

idx_sameYEO = ones(size(G_sheet,1),1);
for i = 1:size(G_sheet,1)
    if G_sheet.YEO7idx(i)==G_sheet.YEO7idx_Group(i)
        idx_sameYEO(i) = 1;
    else
        idx_sameYEO(i) = 0;
    end
end

EI_vec = cell2mat(G_sheet.EI);
Eu_dis_vec = cell2mat(G_sheet.Eu_dis);

idx_YEO7_group = ones(length(layer3_EP_idx),7);
idx_YEO7_group(:,1) = layer3_EP_idx & G_sheet.YEO7idx_Group==5;
idx_YEO7_group(:,2) = layer3_EP_idx & G_sheet.YEO7idx_Group==7;
idx_YEO7_group(:,3) = layer3_EP_idx & G_sheet.YEO7idx_Group==1;
idx_YEO7_group(:,4) = layer3_EP_idx & G_sheet.YEO7idx_Group==2;
idx_YEO7_group(:,5) = layer3_EP_idx & G_sheet.YEO7idx_Group==3;
idx_YEO7_group(:,6) = layer3_EP_idx & G_sheet.YEO7idx_Group==4;
idx_YEO7_group(:,7) = layer3_EP_idx & G_sheet.YEO7idx_Group==6;


for i = YEO7_v
    data_weight_EI_dis{i,1} = EI_vec(idx_YEO7_group(:,i)&idx_sameYEO);
    data_weight_EI_dis{i,2} = EI_vec(idx_YEO7_group(:,i)&~idx_sameYEO);
end


figure
boxplot([data_weight_EI_dis{7,1};data_weight_EI_dis{7,2}],[ones(length(data_weight_EI_dis{7,1}),1);ones(length(data_weight_EI_dis{7,2}),1).*2])

data2 = dataEI2;
data1 = datadis2;
figure
scatter(data1, data2, '+', 'MarkerFaceColor', 'k');
ylabel('EI PZ');
xlabel ('Eu dis PZ');
box 'on'
axis square;

%     disp(r(1,2));
tmp=corrcoef(data1,data2);
if size(tmp,2) == 2
    str=sprintf('r= %1.2f',tmp(1,2));
else
    str = 'notavailable';
end
Tx = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(Tx, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%%
















for i = 1:length(folders)
    fn = sprintf('%s/%s',work_sbjs_path,folders(i,1));
    load(fn);
    bad_chan = sz_loc_EI.bad_chan;
    eleinfo = sz_loc_EI.eleinfo;
    eleinfo.sbjs_bv_names = sbjs_bv_names;
    data_cells{i,1} = eleinfo;
    
    
    
    Nu = sz_loc_EI.Nu;
    Lamda = sz_loc_EI.Lamda;
    Eta = sz_loc_EI.Eta;
    EI_window = sz_loc_EI.EI_window;
    mat_EI_EZ = sz_loc_EI.mat_EI_EZ;
    mat_EI_PZ = sz_loc_EI.mat_EI_PZ;
    mat_EI_EZPZ = sz_loc_EI.mat_EI_EZPZ;
    mat_TII_EZ = sz_loc_EI.mat_TII_EZ;
    mat_TII_PZ = sz_loc_EI.mat_TII_PZ;
    mat_TII_EZPZ = sz_loc_EI.mat_TII_EZPZ;
    T = sz_loc_EI.T;
    TII_seq = sz_loc_EI.T;
    name = sz_loc_EI.name;
    BR_chan_length = size(eleinfo,1);
    bad_chan_idx = ones(size(eleinfo,1),1);
    bad_chan_idx(bad_chan',1) = 0;
    anat_idx = ~strcmp(eleinfo.AAL3,'NotAvailable');
    YEO7_idx = ~strcmp(eleinfo.YEO7,'NotAvailable');
    
    
    
    
    data_cells = cell(length(folders),1);
    for i = 1:length(folders)
        fn = sprintf('%s/%s',work_path,folders(i,1));
        load(fn);
        bad_chan = sz_loc_EI.bad_chan;
        eleinfo = sz_loc_EI.eleinfo;
        name = sz_loc_EI.name;
        
        BR_chan_length = size(eleinfo,1);
        
        bad_chan_idx = ones(size(eleinfo,1),1);
        bad_chan_idx(bad_chan',1) = 0;
        eleinfo.bad_chan = bad_chan_idx;
        
        sbjs_bv_names = cell(BR_chan_length,1);
        for j = 1:BR_chan_length
            sbjs_bv_names{j} = [name,'_',eleinfo.bv_names{j}];
        end
        
        eleinfo.sbjs_bv_names = sbjs_bv_names;
        data_cells{i,1} = eleinfo;
    end
    
    %%
    
    
    eleColors_EP_EI = ones(length(EP_EI_idx),3); %% we have 2 SZ for each patient so this may be different;
    
    for i = 1:length(EP_EI_idx)
        if EP_EI_idx(i)
            eleColors_SEEG_YEO7(i,:) = cdcol.gold_cadmium_yellow;
        else
            eleColors_SEEG_YEO7(i,:) = cdcol.grey;
        end
    end
    
    
    
    cfg=[];
    cfg.view=view_side;
    cfg.elecSize=15;
    cfg.surfType='inflated';
    cfg.opaqueness=1;
    cfg.ignoreDepthElec='n';
    cfg.elecNames = chan_names_all_layer1_SEEG;
    cfg.elecCoord=[elecCoord305 isLeft];
    %cfg. ignoreChans = {'PT049-X7'};
    cfg.elecColors = eleColorsYEO7;
    cfg.elecColorScale=[0 1];
    cfg.title = [];
    cfg.backgroundColor = [1,1,1];
    cfg.overlayParcellation='Y7';
    cfgOut = plotPialSurf_v2('fsaverage',cfg);
    
    
end
