function mannuly_Epi_onset_time

load(sprintf('%s/CARiEEG%s_%.2d.mat',globalVar.CARData,bn,el));

data = WaveletFilter(data.wave,data.fsample,fs_targ,freqs,span,norm,avgfreq);
data.label = globalVar.channame{el};

load(sprintf('%s/originalData/%s/global_%s_%s_%s.mat',dirs.data_root,sbj_name,project_name,sbj_name,block_names{1}),'globalVar');

figureDim = [0 0 .5 .7];%
figure('units', 'normalized', 'outerposition', figureDim)%
hold on%
axes('XLim', [0 globalVar.chanLength], 'units','pixels', ...
     'position',[0 0 .5 .5], 'NextPlot', 'add');%
 
 
 
    for bii=BR_idx
        idx = find(BR_idx==bii);
        data = data_all(bii,:)-data_all(bii+1,:);
        if ismember(bii, globalVar.too_N_or_S)
            plot(zscore(data)+5*idx,'r');
        else
            plot(zscore(data)+5*idx,'b');
        end
    end
    

    ylim([-5 length(BR_idx)*5+5])
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
    
    
    title(['Bad electrodes based on raw power (BR) ',bn])
    xlabel('Time (s)')
    ylabel('Electrode name (BR)')



    x = zeros(1,length(channame_BR));
    y=(1:length(channame_BR))*5;
    
    plot(x,y,'--');
    set(gca,'ytick',y,'yticklabel',channame_BR(:))
    
    
    title(['Bad electrodes based on raw power (BR) ',bn])
    xlabel('Time (s)')
    ylabel('Electrode name (BR)')
    
    ylim([-5 length(BR_idx)*5+5])
    xlim([0 length(data)])






FigH = figure('position',[360 500 400 400]);
axes('XLim', [0 4*pi], 'units','pixels', ...
     'position',[100 50 200 200], 'NextPlot', 'add');
x     = linspace(0, 4*pi, 400);
y     = sin(x);
LineH = plot(x,y);
TextH = uicontrol('style','text',...
    'position',[170 340 40 15]);
SliderH = uicontrol('style','slider','position',[100 280 200 20],...
    'min', 0, 'max', 4*pi);
addlistener(SliderH, 'Value', 'PostSet', @callbackfn);
movegui(FigH, 'center')
    function callbackfn(source, eventdata)
    num          = get(eventdata.AffectedObject, 'Value');
    LineH.YData  = sin(num * x);
    TextH.String = num2str(num);
    end
  end