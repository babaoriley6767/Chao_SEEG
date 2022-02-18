classdef Orienter < handle
%Zoom and pan plotted data (referencePlot) using keypresses (arrow keys, shift, and command).
%Keypresses adjust the limits of the axes (Ax2Adjust) where the data are plotted. 
%A separate axes, TrackerAx, shows the user which portion of the referencePlot is visible in Ax2Adjust. 
%TrackerAx contains copy of referencePlot (miniPlot) such that all its X and Y data are visible. 
%A semi-opaque yellow square (TrackerPatch) highlights the portion of data currently visible in Ax2Adjust.
%See MATLAB documentation on Patch to can change the color etc. of TrackerPatch

%Orienter automatically responds to events.
%If referencePlot data or Ax2Adjust limits changes, listeners automatically 
%adjust miniPlot and TrackerAx accordingly. These allows you link axes with Ax2Adjust and to 
%edit referencePlot data without having to alter Orienter.

%User can adjust the speed/precision (aka moveSpeed) of zooming and panning.
%More specifically, user can adjust how much the axes limits changes with
%each key press by clicking the patch object.

%Orienter should function in standalone applications/GUIs.

%TrackerPatch appearance can be altered by accessing Patch properties 
%(e.g. ?FaceColor? for color, ?FaceAlpha? for transparency).

%% limitations
% Only changes in referencePlot XData and YData properties will cause miniPlot to update.
% If properties of referencePlot change BESIDES X-Data and Y-Data (e.g. color)
% these changes may not transfer automically transfer to miniPlot. 
% you can always ensure miniplot matches referencePlot by using
% update_miniPlot method. When referencePlot X/Y Data do change, ALL
% properties of referencePlot are transferred to miniPlot.

%Ellen Zakreski (2016), efzakreski@mac.com
    %% KEY PRESS 
    % || Key       || Description
    % ||(+modifier)||
    %----------------------------------------------------------------
 %PANING-------------------------------------------------------------
    % || >         || pan right  (moves both left and right sides of patch)
    % || <         || pan left
    % || ^         || pan up     (moves both top and bottom sides of patch)
    % || v         || pan down
 %ZOOMING (horizontally)---------------------------------------------
    % || > shift   || zoom in horizontally (moves only left side of patch right)
    % || < shift   || zoom out horizontally (moves only left side left)
    % || > command || zoom out horizontally  (moves right side right)
    % || < command || zoom in horizontally (moves right side left) 
 %ZOOMING (vertically)-----------------------------------------------
    % || v shift   || zoom in vertically (moves only top side of patch down)
    % || ^ shift   || zoom out vertically (moves top side of patch up)
    % || v command || zoom out vertically  (move bottom side of patch down)
    % || ^ command || zoom in vertically (moves bottom side up) 

    properties
        referencePlot %plot handle or any graphics object with XData and YData
        Ax2Adjust %will become the axes whose limits will be adjusted (by default, is the first axes that is ancestor to the reference plot).
        
        keypressCatcher %figure that when clicked will detect keypresses that will move the patch (by default, is the parent of the TrackerAx).
        
        listeners = struct('plotXData',[],'plotYData',[],'axXLim',[],'axYLim',[]); %see addAllListeners method below
        
        XmoveSpeed = 2; %percentage by which X-axis limits change with each keypress.
        YmoveSpeed = 10; %percentage by which Y-axis limits change with each keypress.
        
        %% TRACKER
        miniPlot %copy of referencePlot (plotted in TrackerAx);
        %properties below are initially structs but become objects 
        %TrackerFig = struct('MenuBar','none','ToolBar','none','Position',[100 100 600 200]); 
        TrackerFig = struct('Position',[100 100 600 200]); 
        TrackerAx = struct('NextPlot','add',...%TrackerAx (by default plotted on TrackerFig) shows the entire reference plot
            'XLimMode','manual','XLim',[1 2],...
            'YLimMode','manual','YLim',[1 2]);
        TrackerPatch 
    end %end properties
    
    properties (Dependent)
        totalXLim % X limits of reference plot 
        totalYLim % Y limits of reference plot
    end
    
    methods
        function obj = Orienter(referencePlot)
            %referencePlot can be plot or any graphics object with XData and YData
            if nargin ~=0
                obj.referencePlot = referencePlot;
                %get axes to adjust
                p = referencePlot.Parent;
                while ~strcmpi(p.Type,'axes')
                    p = p.Parent;
                end
                obj.Ax2Adjust = p; clear p
                % set to default Ax2Adjust axes limits 
                defaultview(obj);
                
                % Add key press function to tracker fig(so that key presses trigger patch movement when trackerfig is in focus.
                obj.TrackerFig = figure(obj.TrackerFig);
               
                % TrackerAx shows the entire reference plot
                obj.TrackerAx = axes(obj.TrackerAx,'parent',obj.TrackerFig);
                
                obj.TrackerPatch = patch('parent',obj.TrackerAx,...
                    'XData',[1 1 2 2 1],...%TrackerPatch is square showing the portion of the total data in  view.
                    'YData',[1 2 2 1 1],...
                    'FaceColor',[1 1 0],...%yellow
                    'FaceAlpha',0.5,... %opacity = 50%.
                    'EdgeColor','none',...
                    'buttondownfcn',@(src,ev)adjust_moveSpeed(obj));
                
                
                % make miniPlot if reference plot changes
                obj.miniPlot = copyobj(obj.referencePlot,obj.TrackerAx);
                
                addAllListeners(obj)
                
                
                 %make miniPlot reflect referencePlot
                 update_plotData(obj,'X');  update_plotData(obj,'Y');
                 %make patch reflect portion of data currently in view
                 update_axLim(obj,'X');  update_axLim(obj,'Y');
                
            end
        end
        %default view
        function defaultview(obj)
            %by default show 20% of the total x data,
            obj.Ax2Adjust.XLim = [obj.totalXLim(1), obj.totalXLim(1) + diff(obj.totalXLim)*0.20];
            % and 80% of the total y data
            obj.Ax2Adjust.YLim = [obj.totalYLim(1), obj.totalYLim(1) + diff(obj.totalYLim)*0.80];
        end
        
        
        %% listeners
        function addAllListeners(obj)
            structfun(@delete,obj.listeners);
            % listen to changes in reference plot X/YData
            
            obj.listeners.plotXData = addlistener(obj.referencePlot,'XData','PostSet',@(src,ev)update_plotData(obj,'X'));
            obj.listeners.plotYData  = addlistener(obj.referencePlot,'YData','PostSet',@(src,ev)update_plotData(obj,'Y'));
           
            %listen to changes in Ax2Adjust
            obj.listeners.axXLim = addlistener(obj.Ax2Adjust,'XLim','PostSet',@(src,ev)update_axLim(obj,'X'));
            obj.listeners.axYLim = addlistener(obj.Ax2Adjust,'YLim','PostSet',@(src,ev)update_axLim(obj,'Y'));
        end
        
        % listern callbacks
        function update_plotData(obj,xory)
            %if referencePlot X/YData changes
            xory=upper(xory);
            %% update minplot only X/YData
            %set(obj.miniPlot,[xory,'Data'],obj.referencePlot.([xory,'Data']));
            %% just copy miniplot again.
            update_miniPlot(obj);
            
            %change the axes limits of TrackerAx
            newlim =obj.(['total',xory,'Lim']);
            if isempty(newlim) || any(isnan(newlim))
                warning('Could not reset tracker axes limits because plot %s data are all empty or NaN''s.', xory)
            else
                set(obj.TrackerAx,[xory,'Lim'],newlim);
            end
        end
        
        function update_axLim(obj,xory) 
            %if Ax2Adjust axes limits change update patch X/YData
            %accordingly so it overlaps with the portion of referencePlot that is currently
            %visible in Ax2Adjust.
            xory=upper(xory);
            % get indexes for converting axes limits to patch X/YData/
            switch xory
                case 'Y'
                    ind = [1 2 2 1 1]; %B T T B B (B=bottom, T=top)
                case 'X'
                    ind = [1 1 2 2 1]; %L L R R L  (L=left, R=right)
            end
            L = obj.Ax2Adjust.([xory,'Lim']);
            %update patch position.
            obj.TrackerPatch.([xory,'Data']) = L(ind);
        end
        
        %% Get
        function l=get.totalXLim(obj)
            l=feval(@(d)[nanmin(d),nanmax(d)],obj.referencePlot.XData);
        end
        function l=get.totalYLim(obj)
            l=feval(@(d)[nanmin(d),nanmax(d)],obj.referencePlot.YData);
        end
        
        
        
        %% Set
        function set.referencePlot(obj,newval)
            %validate
            if ~ishandle(newval) || ~isvalid(newval)
                error('referencePlot must be valid handle')
            end
            assert(all(ismember({'XData','YData'},fieldnames(newval))),'referencePlot must have properties called XData and YData');
            
            obj.referencePlot=newval;
        end
        function set.Ax2Adjust(obj,newval)
            if isstruct(newval)
                 newval = axes(newval);
             end
            if ~isa(newval,'matlab.graphics.axis.Axes') || ~isvalid(newval)
                error('Ax2Adjust must be valid axes handle');
            end
            set(newval,'NextPlot','add','XLimMode','manual','YLimMode','manual');
            obj.Ax2Adjust = newval;
            
        end
        function set.TrackerFig(obj,newval)
            if isstruct(newval)
                newval = figure(newval);
            end
            if ~isa(newval,'matlab.ui.Figure') || ~isvalid(newval)
                error('TrackerFig must be valid figure handle');
            end
            % Add key press function to tracker fig(so that key presses trigger patch movement when trackerfig is in focus.
            newval.KeyPressFcn=@(src,ev)keypress(obj,ev);
            obj.TrackerFig = newval;
        end
        function set.TrackerAx(obj,newval)
             if isstruct(newval)
                 newval = axes(newval);
             end
             if ~isa(newval,'matlab.graphics.axis.Axes') || ~isvalid(newval)
                error('TrackerAx must be valid axes handle');
            end
            set(newval,'NextPlot','add',...
                'XLimMode','manual','XLim',[1 2],'YLimMode','manual','YLim',[1 2]);
            obj.TrackerAx = newval;
        end
        
        function update_miniPlot(obj)
            delete(obj.miniPlot);
            obj.miniPlot = copyobj(obj.referencePlot,obj.TrackerAx);
            uistack(obj.TrackerPatch,'top'); %put TrackerPatch on top of all other plots in TrackerAx.
        end
        
        
        %% moveSpeed
        function adjust_moveSpeed(obj)
            xyls = {'X','Y'};
            prompt = cellfun(@(xy)sprintf('%% of %s-axis that changes per key press:',xy),xyls,'uniformoutput',0);
            defval = cellfun(@(xy)num2str(obj.([xy,'moveSpeed'])),xyls,'uniformoutput',0);
            answer = inputdlg(prompt,...
                'Adjust pan/zoom speed/precision.',...%title
                1,...%numlines
                defval);
            if isempty(answer); return; end
            for n=1:numel(xyls)
                val = str2double(answer{n});
                xy = xyls{n};
                if isnan(val)
                    uiwait(errordlg('%input for must be numeric. %s-axis was not changed.',xy));
                    continue
                end
                obj.([xy,'moveSpeed']) = val;
            end
        end
        
        function keypress(obj,ev)
            TL = ismember('shift',ev.Modifier); %shift moves top/left side only
            BR = ismember('command',ev.Modifier); %command moves bottom/right side only
            %This is because on the left side of the keyboard,
            %shift is above and to the left of the command key.
            
            ispan = ~xor(TL,BR); %ispan if either TL or BR, or not TL nor BR.
            %panning moves both sides move together, if zoom, only one side moves.
            if ispan; TL = true; BR = true; end;
            coef = [0,0]; % determines how much and in what direction Ax2Adjust axes limits change by.
            switch ev.Key
                case 'uparrow'
                    xory = 'Y'; 
                    coef([BR,TL]) = 1; %since for Y-axis limits [ymin,ymax], ymin is Bottom, ymax is Top [B,T]
                case 'downarrow'
                    xory = 'Y';
                    coef([BR,TL]) = -1; %since for Y-axis limits [ymin,ymax], ymin is Bottom, ymax is Top [B,T]
                case 'leftarrow'
                    xory  ='X';
                    coef([TL,BR]) = -1; %since for X-axis limits [xmin,xmax], xmin is Left, xmax is Right [L,R]
                case 'rightarrow'
                    xory  ='X'; 
                    coef([TL,BR]) = 1; %since for X-axis limits [xmin,xmax], xmin is Left, xmax is Right [L,R]
                otherwise
                    return
            end
            movepatch(obj,xory,coef);
        end
        
        function movepatch(obj,xory,coef)
            xory = upper(xory);
            mspeed = obj.([xory,'moveSpeed']);
            span = diff(obj.(sprintf('total%sLim',xory)));
            moveby = span*coef*mspeed/100;
            %coeff elements indicate which sides of axis limits (min, max or both) and whether to increase or decrease
            %is -1 (decrease along axis) or 1 (increase along axis)
            %set new axes position
            obj.Ax2Adjust.([xory,'Lim']) = obj.Ax2Adjust.([xory,'Lim']) + moveby; %will trigger listener to update patch position
        end

    end

    
end

