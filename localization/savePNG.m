function savePNG(f, res, title)
    figpos=getpixelposition(f); 
    resolution=get(0,'ScreenPixelsPerInch'); 
    set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
    print(f,title,'-dpng',['-r',num2str(res)],'-opengl') 
end

