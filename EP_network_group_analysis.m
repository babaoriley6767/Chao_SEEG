% this is the EP network group analysis

%% set roots
addpath(genpath('/Users/tony/Documents/Stanford/Chao_SEEG/'))
addpath(genpath('/Users/tony/Documents/Stanford/DELLO'))
addpath(genpath('/Users/tony/Desktop/function_tools/iELVis/'))
addpath(genpath('/Users/tony/Desktop/function_tools/for_plot/CCEP_fMRI/'))
addpath('/Users/tony/Desktop/function_tools/spm12_7219')
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%home

project_name = 'EPnetwork';

% sbj_names = {'2015_FT_hankunpeng';
%     '2015_FT_huangzutu';
%     '2015_FT_songziwei';
%     '2015_FT_wangmingming';
%     '2015_FT_wangzheng';
%     '2015_FT_xumengjiao';
%     '2015_FT_yangbihong';
%     '2015_FT_yinfengting';
%     '2015_FT_yuanjinrui';
%     '2015_FT_zhaoxichun';
%     '2017_FT_chengang';
%     '2017_FT_genghaowen';
%     '2017_FT_guyuxuan';
%     '2017_FT_jielei';
%     '2017_FT_lijingjing';
%     '2017_FT_lilechun';
%     '2017_FT_linzixuan';
%     '2017_FT_liuboqiang';
%     '2017_FT_liuruilin';
%     '2017_FT_masinan';
%     '2017_FT_panzhiyong';
%     '2017_FT_wangmengru';
%     '2017_FT_wangyan';
%     '2017_FT_xionghuihui';
%     '2017_FT_yangchenglei';
%     '2017_FT_yangrui';
%     '2017_FT_yaodongyuan';
%     '2017_FT_yena';
%     '2017_FT_yuanye';
%     '2017_FT_zhanggenhong'};

sbj_names = {'2015_FT_Yangzhizhong';
    '2016_FT_Chenbinbin';
    '2016_FT_fanglei';
    '2016_FT_liujiafeng';
    '2016_FT_Wuyuhua';
    '2017_FT_Wuxiaolong';
    '2017_FT_Xiangxiang';
    '2019_TT_wudi';
    '2020_TT_jiawenqiao';
    '2020_TT_liubingtong';
    '2020_TT_yurongqi';
    '2020_TT_duguoyin';
    '2019_TT_wangjingfan';
    '2015_FT_hankunpeng';
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
work_path = '/Users/tony/Desktop/4_working_data/Group_analysis/data_2';%%% this folder is for the EI ER TII plot
work_sbjs_path = '/Users/tony/Desktop/4_working_data/Group_analysis/data_sbjs';%%% this folder is for the electrodes pot
R_path = '/Users/tony/Desktop/4_working_data/Group_analysis/R_PART';
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
        eleColors_EP_EI(i,:) = cdcol.bismuth_yellow;
    else
        eleColors_EP_EI(i,:) = cdcol.cyan;
    end
end

cfg=[];
cfg.view='lm';
cfg.elecSize=20;
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

%% YEO7,EP and EI electrodes
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

YEO7_coords_305_all = G_sheet.MNI305_volume(layer2_YEO7_idx,:);
[YEO7_coords_305_UNQ,InSEEG_all,InSEEG_UNQ] = unique(YEO7_coords_305_all,'rows');
EP_coords_305_all = G_sheet.MNI305_volume(layer3_EP_idx,:);
[EP_coords_305_UNQ,InEP_all,InEP_UNQ] = unique(EP_coords_305_all,'rows');
EI_coords_305_all = G_sheet.MNI305_volume(layer4_EI_idx,:);
[EI_coords_305_UNQ,InEI_all,InEI_UNQ] = unique(EI_coords_305_all,'rows');

YEO7_EP_idx = ismember(YEO7_coords_305_UNQ,EP_coords_305_UNQ,'rows');
if sum(YEO7_EP_idx)~= length(EP_coords_305_UNQ)
    warning('plz checking the position of the electrodes  SEEG EP')
end

YEO7_EI_idx = ismember(YEO7_coords_305_UNQ,EI_coords_305_UNQ,'row');
if sum(YEO7_EI_idx) ~= length(EI_coords_305_UNQ)
    warning('plz checking the position of the electrodes  SEEG EI')
end

YEO7_EP_EI_idx = YEO7_EP_idx+YEO7_EI_idx;

chan_names_YEO7_all = G_sheet.sbjs_bv_names(layer2_YEO7_idx);
chan_names_YEO7_UNQ = chan_names_YEO7_all(InSEEG_all);
chan_names_EP_all = G_sheet.sbjs_bv_names(layer3_EP_idx );
chan_names_EP_UNQ = chan_names_EP_all(InEP_all);
chan_names_EI_all = G_sheet.sbjs_bv_names(layer4_EI_idx );
chan_names_EI_UNQ = chan_names_EI_all(InEI_all);

isleft_YEO7_UNQ = ones(size(chan_names_YEO7_UNQ,1),1);
for j = 1:length(isleft_YEO7_UNQ)
    if YEO7_coords_305_UNQ(j,1)<=0
        isleft_YEO7_UNQ(j,1) = 1;
    elseif YEO7_coords_305_UNQ(j,1)>0
        isleft_YEO7_UNQ(j,1) = 0;
    else
    end
end




%figure part

load /Users/tony/Documents/Stanford/Chao_SEEG/visualization/colormaps/cdcol_2018.mat

global globalFsDir;
globalFsDir ='/Users/tony/Desktop/iELVis/Plot_Elelctrodes/';
fsDir='/Users/tony/Desktop/iELVis/Plot_Elelctrodes/';%%% seems useful
cd([fsDir]);



eleColors_SEEG_EP_EI = ones(length(YEO7_EP_EI_idx),3);

for i = 1:length(YEO7_EP_EI_idx)
    if YEO7_EP_EI_idx(i) == 2
        eleColors_SEEG_EP_EI(i,:) = cdcol.bismuth_yellow;
    elseif YEO7_EP_EI_idx(i) == 1
        eleColors_SEEG_EP_EI(i,:) = cdcol.cyan;
    else
        eleColors_SEEG_EP_EI(i,:) = cdcol.grey;
    end
end

cfg=[];
cfg.view='rm';
cfg.elecSize=24;
cfg.surfType='inflated';
cfg.opaqueness=1;
cfg.ignoreDepthElec='n';
cfg.elecNames = chan_names_SEEG_UNQ;
cfg.elecCoord=[YEO7_coords_305_UNQ isleft_YEO7_UNQ];
%cfg. ignoreChans = {'PT049-X7'};
cfg.elecColors = eleColors_SEEG_EP_EI;
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
    
end

G_sheet  = vertcat(data_cells{:});

EI_EZ_matrix = nanmean(cat(3,data_EI_cells{:,1}),3);
EI_PZ_matrix = nanmean(cat(3,data_EI_cells{:,2}),3);
EI_EZPZ_matrix = nanmean(cat(3,data_EI_cells{:,3}),3);
ER_EZ_matrix = nanmean(cat(3,data_ER_cells{:,1}),3);
ER_PZ_matrix = nanmean(cat(3,data_ER_cells{:,2}),3);
ER_EZPZ_matrix = nanmean(cat(3,data_ER_cells{:,3}),3);

% we need to adjust the TII because sometimes the EI =1 electrodes is not
% the earliest one
folders = dir(work_path);
folders = string({folders.name});
folders = folders(~startsWith(folders,"."))' ;

data_cells = cell(length(folders),1);
data_TII_cells = cell(length(folders),3);

exTII_no_all = [];
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
    
    YEO7_v = [5,7,1,2,3,4,6];
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    
    if eleinfo.time_interval{cell2mat(eleinfo.EI)==1}~=0
        disp([name ', this SZ EI=1 is not the earliest one'])
    else
    end
    
    eleinfo.time_interval_adjust = cell2mat(eleinfo.time_interval) - eleinfo.time_interval{cell2mat(eleinfo.EI)==1};
    exTII_idx =  eleinfo.time_interval_adjust <0;
    
    
    
    TII_seq = sz_loc_EI.TII_seq;
    
    time_interval_index_adjust = nan(size(eleinfo,1),1);
    for j = 1:size(eleinfo,1)
        if isnan(eleinfo.time_interval_adjust(j))
            time_interval_index_adjust(j) = NaN;
        elseif exTII_idx(j)
            time_interval_index_adjust(j) = NaN;
        elseif ~exTII_idx(j)
            TIIid = dsearchn(TII_seq,eleinfo.time_interval_adjust(j));
            time_interval_index_adjust(j) = 1/TIIid;
        end
    end
    eleinfo.time_interval_index_adjust = time_interval_index_adjust;
    exTII_no = sum(exTII_idx);
    exTII_no_all = [exTII_no_all;exTII_no];
    
    
    anat_idx = ~strcmp(eleinfo.AAL3,'NotAvailable');
    YEO7_idx = ~strcmp(eleinfo.YEO7,'NotAvailable');
    layer1_SEEG_idx = ~strcmp(eleinfo.YEO7,'NaN');
    layer2_YEO7_idx = layer1_SEEG_idx & anat_idx & YEO7_idx & bad_chan_idx;
    layer3_EP_idx = layer2_YEO7_idx & ~isnan(cell2mat(eleinfo.EI));
    layer4_EI_idx = layer3_EP_idx & cell2mat(eleinfo.EI)==1;
    layer4_EZ_idx = layer3_EP_idx & cell2mat(eleinfo.EI)>=0.3 & cell2mat(eleinfo.EI) < 1;
    layer4_PZ_idx = layer3_EP_idx & cell2mat(eleinfo.EI)>0 & cell2mat(eleinfo.EI)<0.3;
    layer4_EZPZ_idx = layer3_EP_idx & cell2mat(eleinfo.EI)>0& cell2mat(eleinfo.EI) < 1;
    
    
    
    TII_adjust_EZ_mat = nan(1,length(YEO7_v));
    TII_adjust_PZ_mat = nan(1,length(YEO7_v));
    TII_adjust_EPPZ_mat = nan(1,length(YEO7_v));
    for j = 1:length(YEO7_v)
        TII_adjust_EZ_mat(j) = nanmean(eleinfo.time_interval_index_adjust(layer4_EZ_idx & ismember(eleinfo.YEO7idx,YEO7_v(j))));
        TII_adjust_PZ_mat(j) = nanmean(eleinfo.time_interval_index_adjust(layer4_PZ_idx & ismember(eleinfo.YEO7idx,YEO7_v(j))));
        TII_adjust_EZPZ_mat(j) = nanmean(eleinfo.time_interval_index_adjust(layer4_EZPZ_idx & ismember(eleinfo.YEO7idx,YEO7_v(j))));
    end
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    mat_TII_adjust_EZ = nan(7,7);
    mat_TII_adjust_PZ = nan(7,7);
    mat_TII_adjust_EZPZ = nan(7,7);
    
    mat_TII_adjust_EZ(row_idx,:) = TII_adjust_EZ_mat;
    mat_TII_adjust_PZ(row_idx,:) =TII_adjust_PZ_mat;
    mat_TII_adjust_EZPZ(row_idx,:) = TII_adjust_EZPZ_mat;
    
    data_TII_cells{i,1} = mat_TII_adjust_EZ;
    data_TII_cells{i,2} = mat_TII_adjust_PZ;
    data_TII_cells{i,3} = mat_TII_adjust_EZPZ;
    
    data_cells{i,1} = eleinfo;
    
end
disp(['we have exclude ' mat2str(sum(exTII_no_all)) ' electrodes']);
G_adjust_sheet  = vertcat(data_cells{:});

TII_adjust_EZ_matrix = nanmean(cat(3,data_TII_cells{:,1}),3);
TII_adjust_PZ_matrix = nanmean(cat(3,data_TII_cells{:,2}),3);
TII_adjust_EZPZ_matrix = nanmean(cat(3,data_TII_cells{:,3}),3);


% all of the 9 matrix
figure
subplot(3,3,1)
imagesc(EI_EZ_matrix)
colormap(gca,'parula');
colorbar() ; % Add color bar and make sure the color ranges from 0:1
caxis([-1,1]);
set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
set(gca,'XTickLabelRotation',45)
title('EI EZ')

subplot(3,3,2)
imagesc(EI_PZ_matrix)
colormap(gca,'parula');
colorbar() ; % Add color bar and make sure the color ranges from 0:1
caxis([-1,1]);
set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
set(gca,'XTickLabelRotation',45)
title('EI PZ')

subplot(3,3,3)
imagesc(EI_EZPZ_matrix)
colormap(gca,'parula');
colorbar() ; % Add color bar and make sure the color ranges from 0:1
caxis([-1,1]);
set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
set(gca,'XTickLabelRotation',45)
title('EI EZPZ')

subplot(3,3,4)
imagesc(TII_adjust_EZ_matrix)
colormap(gca,'parula');
colorbar() ; % Add color bar and make sure the color ranges from 0:1
caxis([-1,1]);
set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
set(gca,'XTickLabelRotation',45)
title('TII EZ')

subplot(3,3,5)
imagesc(TII_adjust_PZ_matrix)
colormap(gca,'parula');
colorbar() ; % Add color bar and make sure the color ranges from 0:1
caxis([-1,1]);
set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
set(gca,'XTickLabelRotation',45)
title('TII PZ')

subplot(3,3,6)
imagesc(TII_adjust_EZPZ_matrix)
colormap(gca,'parula');
colorbar() ; % Add color bar and make sure the color ranges from 0:1
caxis([-1,1]);
set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
set(gca,'XTickLabelRotation',45)
title('TII EZPZ')

subplot(3,3,7)
imagesc(ER_EZ_matrix)
colormap(gca,'parula');
colorbar() ; % Add color bar and make sure the color ranges from 0:1
caxis([-1,1]);
set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
set(gca,'XTickLabelRotation',45)
title('ER EZ')

subplot(3,3,8)
imagesc(ER_PZ_matrix)
colormap(gca,'parula');
colorbar() ; % Add color bar and make sure the color ranges from 0:1
caxis([-1,1]);
set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
set(gca,'XTickLabelRotation',45)
title('ER PZ')

subplot(3,3,9)
imagesc(ER_EZPZ_matrix)
colormap(gca,'parula');
colorbar() ; % Add color bar and make sure the color ranges from 0:1
caxis([-1,1]);
set(gca,'XTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
set(gca,'YTicklabel',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'})
set(gca,'XTickLabelRotation',45)
title('ER EZPZ')
% figure in the paper
%% adapt to R folder, upper results
EI_EZPZ_matrix;
TII_adjust_EZPZ_matrix;
ER_EZPZ_matrix;

EI_EZPZ_matrix_T = array2table(EI_EZPZ_matrix,...
    'VariableNames',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'},...
    'RowNames',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
EI_EZPZ_matrix_T_flip = flip(EI_EZPZ_matrix_T);

TII_adjust_EZPZ_matrix_T = array2table(TII_adjust_EZPZ_matrix,...
    'VariableNames',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'},...
    'RowNames',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
TII_adjust_EZPZ_matrix_T_flip = flip(TII_adjust_EZPZ_matrix_T);

ER_EZPZ_matrix_T = array2table(ER_EZPZ_matrix,...
    'VariableNames',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'},...
    'RowNames',{'Limbic', 'Default', 'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Frontalparietal'});
ER_EZPZ_matrix_T_flip = flip(ER_EZPZ_matrix_T);

YEO7_name_v = {'Limbic'; 'Default'; 'Visual';'Somatomotor';'Dorsal Attention';'Ventral Attention';'Frontoparietal'};
YEO7_name_v_flip = flip(YEO7_name_v);

cd(R_path);
writetable(EI_EZPZ_matrix_T_flip,'EI_EZPZ_matrix_T_flip.csv')
writetable(TII_adjust_EZPZ_matrix_T_flip,'TII_adjust_EZPZ_matrix_T_flip.csv')
writetable(ER_EZPZ_matrix_T_flip,'ER_EZPZ_matrix_T_flip.csv')



%% inside and outside the YEO7, making the G_adjust_sheet_T file for R (quick)
folders = dir(work_path);
folders = string({folders.name});
folders = folders(~startsWith(folders,"."))' ;

data_cells = cell(length(folders),1);
data_TII_cells = cell(length(folders),3);

exTII_no_all = [];
sz_id_all = [];
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
    
    YEO7_v = [5,7,1,2,3,4,6];
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    
    if eleinfo.time_interval{cell2mat(eleinfo.EI)==1}~=0
        disp([name ', this SZ EI=1 is not the earliest one'])
    else
    end
    
    eleinfo.time_interval_adjust = cell2mat(eleinfo.time_interval) - eleinfo.time_interval{cell2mat(eleinfo.EI)==1};
    exTII_idx =  eleinfo.time_interval_adjust <0;
    
    for j = 1:size(eleinfo,1)
        if  exTII_idx(j)
            eleinfo.time_interval_adjust(j) = NaN;
        else
        end
    end
    
    
    TII_seq = sz_loc_EI.TII_seq;
    
    time_interval_index_adjust = nan(size(eleinfo,1),1);
    for j = 1:size(eleinfo,1)
        if isnan(eleinfo.time_interval_adjust(j))
            time_interval_index_adjust(j) = NaN;
        elseif exTII_idx(j)
            time_interval_index_adjust(j) = NaN;
        elseif ~exTII_idx(j)
            TIIid = dsearchn(TII_seq,eleinfo.time_interval_adjust(j));
            time_interval_index_adjust(j) = 1/TIIid;
        end
    end
    eleinfo.time_interval_index_adjust = time_interval_index_adjust;
    exTII_no = sum(exTII_idx);
    exTII_no_all = [exTII_no_all;exTII_no];
    
    
    anat_idx = ~strcmp(eleinfo.AAL3,'NotAvailable');
    YEO7_idx = ~strcmp(eleinfo.YEO7,'NotAvailable');
    layer1_SEEG_idx = ~strcmp(eleinfo.YEO7,'NaN');
    layer2_YEO7_idx = layer1_SEEG_idx & anat_idx & YEO7_idx & bad_chan_idx;
    layer3_EP_idx = layer2_YEO7_idx & ~isnan(cell2mat(eleinfo.EI));
    layer4_EI_idx = layer3_EP_idx & cell2mat(eleinfo.EI)==1;
    layer4_EZ_idx = layer3_EP_idx & cell2mat(eleinfo.EI)>=0.3 & cell2mat(eleinfo.EI) < 1;
    layer4_PZ_idx = layer3_EP_idx & cell2mat(eleinfo.EI)>0 & cell2mat(eleinfo.EI)<0.3;
    layer4_EZPZ_idx = layer3_EP_idx & cell2mat(eleinfo.EI)>0& cell2mat(eleinfo.EI) < 1;
    
    
    
    TII_adjust_EZ_mat = nan(1,length(YEO7_v));
    TII_adjust_PZ_mat = nan(1,length(YEO7_v));
    TII_adjust_EPPZ_mat = nan(1,length(YEO7_v));
    for j = 1:length(YEO7_v)
        TII_adjust_EZ_mat(j) = nanmean(eleinfo.time_interval_index_adjust(layer4_EZ_idx & ismember(eleinfo.YEO7idx,YEO7_v(j))));
        TII_adjust_PZ_mat(j) = nanmean(eleinfo.time_interval_index_adjust(layer4_PZ_idx & ismember(eleinfo.YEO7idx,YEO7_v(j))));
        TII_adjust_EZPZ_mat(j) = nanmean(eleinfo.time_interval_index_adjust(layer4_EZPZ_idx & ismember(eleinfo.YEO7idx,YEO7_v(j))));
    end
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    mat_TII_adjust_EZ = nan(7,7);
    mat_TII_adjust_PZ = nan(7,7);
    mat_TII_adjust_EZPZ = nan(7,7);
    
    mat_TII_adjust_EZ(row_idx,:) = TII_adjust_EZ_mat;
    mat_TII_adjust_PZ(row_idx,:) =TII_adjust_PZ_mat;
    mat_TII_adjust_EZPZ(row_idx,:) = TII_adjust_EZPZ_mat;
    
    data_TII_cells{i,1} = mat_TII_adjust_EZ;
    data_TII_cells{i,2} = mat_TII_adjust_PZ;
    data_TII_cells{i,3} = mat_TII_adjust_EZPZ;
    
    data_cells{i,1} = eleinfo;
    
    sz_id = ones(size(data_cells{i,1},1),1).*i;
    sz_id_all = [sz_id_all;sz_id];
    
    
end
disp(['we have exclude ' num2str(sum(exTII_no_all)) ' electrodes']);
G_adjust_sheet  = vertcat(data_cells{:});

sbj_id_idx = ones(size(folders));
for i = 2:size(folders,1)
    curr_id = convertStringsToChars(folders(i));
    prev_id = convertStringsToChars(folders(i-1));
    if strcmp(curr_id(1:end-8),prev_id(1:end-8))
        sbj_id_idx(i,1) = sbj_id_idx(i-1,1);
    else
        sbj_id_idx(i,1) = sbj_id_idx(i-1,1)+1;
    end
end

sbj_numid_all = [];
for i = 1:size(data_cells,1)
    sbj_numid = ones(size(data_cells{i,1},1),1).*sbj_id_idx(i);
    sbj_numid_all = [sbj_numid_all;sbj_numid];
end
    





anat_idx = ~strcmp(G_adjust_sheet.AAL3,'NotAvailable');
YEO7_idx = ~strcmp(G_adjust_sheet.YEO7,'NotAvailable');
layer1_SEEG_idx = ~strcmp(G_adjust_sheet.YEO7,'NaN');
layer2_YEO7_idx = layer1_SEEG_idx & anat_idx & YEO7_idx & G_adjust_sheet.bad_chan;
layer3_EP_idx = layer2_YEO7_idx & ~isnan(cell2mat(G_adjust_sheet.EI)) & cell2mat(G_adjust_sheet.EI)>0;
layer4_EI_idx = layer3_EP_idx & cell2mat(G_adjust_sheet.EI)==1;
layer4_EZPZ_idx = layer3_EP_idx & cell2mat(G_adjust_sheet.EI)>0& cell2mat(G_adjust_sheet.EI) < 1;


idx_sameYEO = ones(size(G_adjust_sheet,1),1);
for i = 1:size(G_adjust_sheet,1)
    if G_adjust_sheet.YEO7idx(i)==G_adjust_sheet.YEO7idx_Group(i)
        idx_sameYEO(i) = 1;
    else
        idx_sameYEO(i) = 0;
    end
end

EI_vec = cell2mat(G_adjust_sheet.EI);
TII_adjust_vec = G_adjust_sheet.time_interval_index_adjust;
TI_adjust_vec = G_adjust_sheet.time_interval_adjust;
ER_vec = cell2mat(G_adjust_sheet.ER_index);
Eu_dis_vec = cell2mat(G_adjust_sheet.Eu_dis);

intra_idx = layer4_EZPZ_idx &  idx_sameYEO;
inter_idx = layer4_EZPZ_idx & ~idx_sameYEO;

intra_EI = EI_vec(intra_idx);
intra_TI_adjust = TI_adjust_vec(intra_idx);
intra_TII_adjust = TII_adjust_vec(intra_idx);
intra_ER = ER_vec(intra_idx);
intra_Eu_dis = Eu_dis_vec(intra_idx);

inter_EI = EI_vec(inter_idx);
inter_TI_adjust = TI_adjust_vec(inter_idx);
inter_TII_adjust = TII_adjust_vec(inter_idx);
inter_ER = ER_vec(inter_idx);
inter_Eu_dis = Eu_dis_vec(inter_idx);


G_adjust_sheet_T = table(G_adjust_sheet.YEO7idx_Group,...
    G_adjust_sheet.YEO7idx,...
    layer4_EZPZ_idx,...
    idx_sameYEO,...
    EI_vec,...
    TI_adjust_vec,...
    TII_adjust_vec,...
    ER_vec,...
    Eu_dis_vec,...
    sz_id_all,...
    sbj_numid_all,...
    'VariableNames',{'YEO7idx_Group',...
    'YEO7idx',...
    'layer4_EZPZ_idx',...
    'idx_sameYEO',...
    'EI_vec',...
    'TI_adjust_vec',...
    'TII_adjust_vec',...
    'ER_vec',...
    'Eu_dis_vec',...
    'sz_id_all',...
    'sbj_numid_all'});

data2 = inter_TI_adjust;
data1 = inter_Eu_dis;
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

%% adapt to R folder, upper results
cd(R_path);
writetable(G_adjust_sheet_T,'G_adjust_sheet_T.csv')


%% Bar PLOT OF each network unfinish
folders = dir(work_path);
folders = string({folders.name});
folders = folders(~startsWith(folders,"."))' ;

data_cells = cell(length(folders),1);
data_TII_cells = cell(length(folders),3);

exTII_no_all = [];
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
    
    YEO7_v = [5,7,1,2,3,4,6];
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    
    if eleinfo.time_interval{cell2mat(eleinfo.EI)==1}~=0
        disp([name ', this SZ EI=1 is not the earliest one'])
    else
    end
    
    eleinfo.time_interval_adjust = cell2mat(eleinfo.time_interval) - eleinfo.time_interval{cell2mat(eleinfo.EI)==1};
    exTII_idx =  eleinfo.time_interval_adjust <0;
    
    for j = 1:size(eleinfo,1)
        if  exTII_idx(j)
            eleinfo.time_interval_adjust(j) = NaN;
        else
        end
    end
    
    
    TII_seq = sz_loc_EI.TII_seq;
    
    time_interval_index_adjust = nan(size(eleinfo,1),1);
    for j = 1:size(eleinfo,1)
        if isnan(eleinfo.time_interval_adjust(j))
            time_interval_index_adjust(j) = NaN;
        elseif exTII_idx(j)
            time_interval_index_adjust(j) = NaN;
        elseif ~exTII_idx(j)
            TIIid = dsearchn(TII_seq,eleinfo.time_interval_adjust(j));
            time_interval_index_adjust(j) = 1/TIIid;
        end
    end
    eleinfo.time_interval_index_adjust = time_interval_index_adjust;
    exTII_no = sum(exTII_idx);
    exTII_no_all = [exTII_no_all;exTII_no];
    
    
    anat_idx = ~strcmp(eleinfo.AAL3,'NotAvailable');
    YEO7_idx = ~strcmp(eleinfo.YEO7,'NotAvailable');
    layer1_SEEG_idx = ~strcmp(eleinfo.YEO7,'NaN');
    layer2_YEO7_idx = layer1_SEEG_idx & anat_idx & YEO7_idx & bad_chan_idx;
    layer3_EP_idx = layer2_YEO7_idx & ~isnan(cell2mat(eleinfo.EI));
    layer4_EI_idx = layer3_EP_idx & cell2mat(eleinfo.EI)==1;
    layer4_EZ_idx = layer3_EP_idx & cell2mat(eleinfo.EI)>=0.3 & cell2mat(eleinfo.EI) < 1;
    layer4_PZ_idx = layer3_EP_idx & cell2mat(eleinfo.EI)>0 & cell2mat(eleinfo.EI)<0.3;
    layer4_EZPZ_idx = layer3_EP_idx & cell2mat(eleinfo.EI)>0& cell2mat(eleinfo.EI) < 1;
    
    
    
    TII_adjust_EZ_mat = nan(1,length(YEO7_v));
    TII_adjust_PZ_mat = nan(1,length(YEO7_v));
    TII_adjust_EPPZ_mat = nan(1,length(YEO7_v));
    for j = 1:length(YEO7_v)
        TII_adjust_EZ_mat(j) = nanmean(eleinfo.time_interval_index_adjust(layer4_EZ_idx & ismember(eleinfo.YEO7idx,YEO7_v(j))));
        TII_adjust_PZ_mat(j) = nanmean(eleinfo.time_interval_index_adjust(layer4_PZ_idx & ismember(eleinfo.YEO7idx,YEO7_v(j))));
        TII_adjust_EZPZ_mat(j) = nanmean(eleinfo.time_interval_index_adjust(layer4_EZPZ_idx & ismember(eleinfo.YEO7idx,YEO7_v(j))));
    end
    row_idx = dsearchn(YEO7_v',sz_loc_EI.YEO7idx_Group);
    mat_TII_adjust_EZ = nan(7,7);
    mat_TII_adjust_PZ = nan(7,7);
    mat_TII_adjust_EZPZ = nan(7,7);
    
    mat_TII_adjust_EZ(row_idx,:) = TII_adjust_EZ_mat;
    mat_TII_adjust_PZ(row_idx,:) =TII_adjust_PZ_mat;
    mat_TII_adjust_EZPZ(row_idx,:) = TII_adjust_EZPZ_mat;
    
    data_TII_cells{i,1} = mat_TII_adjust_EZ;
    data_TII_cells{i,2} = mat_TII_adjust_PZ;
    data_TII_cells{i,3} = mat_TII_adjust_EZPZ;
    
    data_cells{i,1} = eleinfo;
    
end
disp(['we have exclude ' num2str(sum(exTII_no_all)) ' electrodes']);
G_adjust_sheet  = vertcat(data_cells{:});



anat_idx = ~strcmp(G_adjust_sheet.AAL3,'NotAvailable');
YEO7_idx = ~strcmp(G_adjust_sheet.YEO7,'NotAvailable');
layer1_SEEG_idx = ~strcmp(G_adjust_sheet.YEO7,'NaN');
layer2_YEO7_idx = layer1_SEEG_idx & anat_idx & YEO7_idx & G_adjust_sheet.bad_chan;
layer3_EP_idx = layer2_YEO7_idx & ~isnan(cell2mat(G_adjust_sheet.EI)) & cell2mat(G_adjust_sheet.EI)>0;
layer4_EI_idx = layer3_EP_idx & cell2mat(G_adjust_sheet.EI)==1;
layer4_EZPZ_idx = layer3_EP_idx & cell2mat(G_adjust_sheet.EI)>0& cell2mat(G_adjust_sheet.EI) < 1;

YEO7_sheet = G_adjust_sheet(layer4_EI_idx,:);
YEO7_data = YEO7_sheet.YEO7idx_Group;
YEO7_group = cell(size(YEO7_data,1),1);
for i = 1:size(YEO7_data,1)
    if contains(YEO7_sheet(i,:).YEO7{:},'Ventral_Attention')
        YEO7_group{i,1} = 'Ventral Attention';
    elseif contains(YEO7_sheet(i,:).YEO7{:},'Default')
        YEO7_group{i,1} = 'Default';
    elseif contains(YEO7_sheet(i,:).YEO7{:},'Frontoparietal')
        YEO7_group{i,1} = 'Frontoparietal';
    elseif contains(YEO7_sheet(i,:).YEO7{:},'Limbic')
        YEO7_group{i,1} = 'Limbic';
    elseif contains(YEO7_sheet(i,:).YEO7{:},'Dorsal_Attention')
        YEO7_group{i,1} = 'Dorsal Attention';
    elseif contains(YEO7_sheet(i,:).YEO7{:},'Visual')
        YEO7_group{i,1} = 'Visual';
    elseif contains(YEO7_sheet(i,:).YEO7{:},'Somatomotor')
        YEO7_group{i,1} = 'Somatomotor';
    else
    end
end 
no_group = ones(7,1);
YEO7_name_v = {'Limbic'; 'Default'; 'Visual';'Somatomotor';'Dorsal Attention';'Ventral Attention';'Frontoparietal'};
YEO7_name_v_flip = flip(YEO7_name_v);
    
YEO7_name_v_flip = {'Frontoparietal';'Ventral Attention';'Dorsal Attention';'Somatomotor';'Visual';'Default';'Limbic'};
for i = 1:7
    no_group(i) = sum(strcmp(YEO7_group(:),YEO7_name_v{i}));
end
YEO7_group_table = table(no_group,YEO7_name_v,'VariableNames',{'no_group','YEO7_name'});
YEO7_group_table_flip = flip(YEO7_group_table);
cd(R_path);
writetable(YEO7_group_table,'YEO7_group.csv')
writetable(YEO7_group_table_flip,'YEO7_group_flip.csv')

%% modify the G_adjust_sheet_to make the line plot (quick)
win_length = 31;
stepping_length = 1;
vec_sheer = ~(isnan(G_adjust_sheet_T.EI_vec) | isnan(G_adjust_sheet_T.TII_adjust_vec) | isnan(G_adjust_sheet_T.TI_adjust_vec)...
    | isnan(G_adjust_sheet_T.ER_vec) | (G_adjust_sheet_T.EI_vec == 1)|(G_adjust_sheet_T.layer4_EZPZ_idx)==0|(G_adjust_sheet_T.TI_adjust_vec)==0) ;
T_sheer = G_adjust_sheet_T(vec_sheer,:);
time_wins = buffer(0:round(max(T_sheer.Eu_dis_vec)),win_length,win_length-stepping_length,'nodelay');
T_sheer_cells = cell(size(time_wins,2),1);
% for i  =1:length(T_sheer_cells)
%     idx_buffer = time_wins(:,i);
%     idx_T_sheer = T_sheer.Eu_dis_vec>=idx_buffer(1)&T_sheer.Eu_dis_vec<idx_buffer(end);
%     data_sheer = T_sheer(idx_T_sheer,:);
%     T_sheer_cells{i,1} = data_sheer;
%     T_sheer_cells{i,2} = size(data_sheer,1);% elec number at this window
%     T_sheer_cells{i,3} = sum(data_sheer.idx_sameYEO);% elec number of intra at this window
%     T_sheer_cells{i,4} = size(data_sheer,1)-sum(data_sheer.idx_sameYEO);% elec number of inter at this window
%     
%     mean_intra = nanmean(data_sheer.EI_vec(data_sheer.idx_sameYEO==1,:));
%     T_sheer_cells{i,5} = mean_intra;%mean EI_vec intra
%     mean_inter = nanmean(data_sheer.EI_vec(data_sheer.idx_sameYEO==0,:));
%     T_sheer_cells{i,6} = mean_inter;%mean EI_vec inter
%     lme = fitlme(data_sheer,'EI_vec~idx_sameYEO +(1|sbj_numid_all)');
%     T_sheer_cells{i,7} = lme.Coefficients.pValue(1);%lme p value for EI vec
%     [h,p,ci,stats] = ttest2(data_sheer.EI_vec(data_sheer.idx_sameYEO==1,:),data_sheer.EI_vec(data_sheer.idx_sameYEO==0,:));
%     T_sheer_cells{i,8} = p;
%     
%     
%     mean_intra = nanmean(data_sheer.TI_adjust_vec(data_sheer.idx_sameYEO==1,:));
%     T_sheer_cells{i,9} = mean_intra;%mean TI_vec intra
%     mean_inter = nanmean(data_sheer.TI_adjust_vec(data_sheer.idx_sameYEO==0,:));
%     T_sheer_cells{i,10} = mean_inter;%mean TI_vec inter
%     lme = fitlme(data_sheer,'TI_adjust_vec~idx_sameYEO +(1|sbj_numid_all)');
%     T_sheer_cells{i,11} = lme.Coefficients.pValue(1);%lme p value for EI vec
%     [h,p,ci,stats] = ttest2(data_sheer.TI_adjust_vec(data_sheer.idx_sameYEO==1,:),data_sheer.TI_adjust_vec(data_sheer.idx_sameYEO==0,:));
%     T_sheer_cells{i,12} = p;
%     
%     
%     mean_intra = nanmean(data_sheer.ER_vec(data_sheer.idx_sameYEO==1,:));
%     T_sheer_cells{i,13} = mean_intra;%mean ER_vec intra
%     mean_inter = nanmean(data_sheer.ER_vec(data_sheer.idx_sameYEO==0,:));
%     T_sheer_cells{i,14} = mean_inter;%mean ER_vec inter
%     lme = fitlme(data_sheer,'ER_vec~idx_sameYEO +(1|sbj_numid_all)');
%     T_sheer_cells{i,15} = lme.Coefficients.pValue(1);%lme p value for EI vec
%     [h,p,ci,stats] = ttest2(data_sheer.ER_vec(data_sheer.idx_sameYEO==1,:),data_sheer.ER_vec(data_sheer.idx_sameYEO==0,:));
%     T_sheer_cells{i,16} = p;
% 
% end

cd(R_path);
%writetable(T_sheer,'G_adjust_sheet_T2.csv')
writetable(T_sheer,'T_sheer.csv')
% cd(R_path);
% writetable(G_adjust_sheet_T,'G_adjust_sheet_T.csv')

%% cingulate
            smidx = G_adjust_sheet.YEO7idx_Group==2;
            smidx2 = cell2mat(G_adjust_sheet.EI)==1;
            
            anat_idx = ~strcmp(G_adjust_sheet.AAL3,'NotAvailable');
            YEO7_idx = ~strcmp(G_adjust_sheet.YEO7,'NotAvailable');
            layer1_SEEG_idx = ~strcmp(G_adjust_sheet.YEO7,'NaN');
            layer2_YEO7_idx = layer1_SEEG_idx & anat_idx & YEO7_idx & G_adjust_sheet.bad_chan;
            layer3_EP_idx = layer2_YEO7_idx & ~isnan(cell2mat(G_adjust_sheet.EI)) & cell2mat(G_adjust_sheet.EI)>0;
            layer4_EI_idx = layer3_EP_idx & cell2mat(G_adjust_sheet.EI)==1;
            layer4_EZPZ_idx = layer3_EP_idx & cell2mat(G_adjust_sheet.EI)>0& cell2mat(G_adjust_sheet.EI) < 1;
            
            sm_sheet = G_adjust_sheet(smidx&layer3_EP_idx ,:);
            
            
            addpath(genpath('/Users/tony/Desktop/function_tools/for_plot/iELVis-master/'))
            %
            global globalFsDir;
            globalFsDir ='/Users/tony/Desktop/iELVis/Plot_Elelctrodes/';
            fsDir='/Users/tony/Desktop/iELVis/Plot_Elelctrodes/';%%% seems useful
            cd([fsDir]);
            
            % plotting parameters
            
            sm_coords = sm_sheet.MNI305_volume;
            
            isleft_vector = zeros(size(sm_sheet,1),1);
            for j = 1:length(isleft_vector)
                if sm_sheet.MNI305_volume(j,1)<0
                    sm_coords(j,1) = sm_coords(j,1)+2*abs(sm_coords(j,1));
                else
                end
            end
            
            eleColors_sm = cell2mat(sm_sheet.EI);
            
            
            
            cfg=[];
            cfg.view='rm';
            cfg.clickElec = 'y';
            cfg.elecSize=15;
            cfg.surfType='inflated';
            cfg.opaqueness=0.1;
            cfg.ignoreDepthElec='n';
            cfg.elecNames = sm_sheet.sbjs_bv_names;
            cfg.showLabels='n';
            cfg.elecCoord=[sm_coords isleft_vector];
            cfg.elecColors = eleColors_sm;
            cfg.elecColorScale=[0 1];
            cfg.title = [];
            cfg.backgroundColor = [1,1,1];
            %     cfg.overlayParcellation='Y7';
            cfgOut = plotPialSurf_v2('fsaverage',cfg);
            
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
            
