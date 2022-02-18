function[visual_check] = plotelec_brainvisa(bv_chan_names,bv_chan_hem,MNI152_volume,MNI305_volume, YEO7idx)
%% addpath
% I am using this
addpath(genpath('/Users/tony/Desktop/function_tools/for_plot/iELVis-master/'))
% global globalFsDir;
% globalFsDir ='/Volumes/CHAO_IRON_M/iELVis_working/Plot_Elelctrodes';
% fsDir='/Volumes/CHAO_IRON_M/iELVis_working/Plot_Elelctrodes';%%% seems useful
% cd([fsDir]);



global globalFsDir;
globalFsDir ='/Users/tony/Desktop/iELVis/Plot_Elelctrodes/';
fsDir='/Users/tony/Desktop/iELVis/Plot_Elelctrodes/';%%% seems useful
cd([fsDir]);

%% plotting parameters


isLeft = [];
for j = 1:length(bv_chan_hem)
    if strcmp(bv_chan_hem{j},'R')
        isLeft(j,1) = 0;
    else
        isLeft(j,1) = 1;
    end
end

%%


elecCoord152 = MNI152_volume;
elecCoord305 = MNI305_volume;

eleColors =[];
for j = 1:length(bv_chan_hem)
    eleColors(j,:) = [0.3137,0.4510,0.4980];
end

eleColorsYEO7 = [];
for j = 1:length(bv_chan_hem)
    if YEO7idx(j) == 0
        eleColorsYEO7(j,:) = [0	0	0];
    elseif YEO7idx(j) == 1
        eleColorsYEO7(j,:) = [0.470588235294118	0.0705882352941177	0.525490196078431];
    elseif YEO7idx(j) == 2
        eleColorsYEO7(j,:) = [0.274509803921569	0.509803921568627	0.705882352941177];
    elseif YEO7idx(j) == 3
        eleColorsYEO7(j,:) = [0	0.462745098039216	0.0549019607843137];
    elseif YEO7idx(j) == 4
        eleColorsYEO7(j,:) = [0.768627450980392	0.227450980392157	0.980392156862745];
    elseif YEO7idx(j) == 5
        eleColorsYEO7(j,:) = [0.862745098039216	0.972549019607843	0.643137254901961];
    elseif YEO7idx(j) == 6
        eleColorsYEO7(j,:) = [0.901960784313726	0.580392156862745	0.133333333333333];
    elseif YEO7idx(j) == 7
        eleColorsYEO7(j,:) = [0.803921568627451	0.243137254901961	0.305882352941177];
    else
    end
end


% plot the figure of contacts of insula and others in Native space

view_side = lower(bv_chan_hem{1});

cfg=[];
cfg.view=view_side;
cfg.elecSize=15;
cfg.surfType='inflated';
cfg.opaqueness=1;
cfg.ignoreDepthElec='n';
cfg.elecNames = bv_chan_names;
cfg.elecCoord=[elecCoord305 isLeft];
%cfg. ignoreChans = {'PT049-X7'};
cfg.elecColors = eleColorsYEO7;
cfg.elecColorScale=[0 1];
cfg.title = [];
cfg.backgroundColor = [1,1,1];
cfg.overlayParcellation='Y7';
cfgOut = plotPialSurf_v2('fsaverage',cfg);

visual_check = false;

prompt = ['\n You need to mannuly identify if these electrodes localization'...
    '\n y, this is FINE'...
    '\n Click any other buttons that there is errors of the electrodes localization\n'];

ID = input(prompt,'s');

if strcmp(ID, 'y')
    visual_check = true;
else
end
close all
end