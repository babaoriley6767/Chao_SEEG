function OrienterObj = example
    %example of Orienter. 
    XX= unique(randi(300,[100,1]),'sorted');
    YY = rand(size(XX));
    f=figure; a=axes('parent',f); refPlot = plot(XX,YY,'-r*','parent',a);
    title(a,'Ax2Adjust: Pan using arrowkeys. Use shift/command to zoom in and out.');
    xlabel(a,'Visible X-Data');
    ylabel(a,'Visible Y-Data');
    legend({'referencePlot'});
    
    OrienterObj = Orienter(refPlot);
    
    title(OrienterObj.TrackerAx,'TrackerAx: Click patch(yellow square) to adjust pan/zoom speed');
    xlabel(OrienterObj.TrackerAx,'Total X-Data');
    ylabel(OrienterObj.TrackerAx,'Total Y-Data');
end