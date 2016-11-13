%
%  作者：王元恺  日期：2016-10-22
%  输出图像
%

function outPlot(value,matrix,plotname,yname)
    deltaX = (value(2)-value(1))/200;
    plot(value(1):deltaX:value(2),matrix);
    hold on;
    str = 'unknown';
    if value(3) == 1
        str = 'vanLeer';
    elseif value(3) == 2
        str = 'vanAlbada';
    elseif value(3) == 3
        str = 'minmod';
    elseif value(3) == 4
        str = 'superbee';
    end
    title({plotname,['Limiter = ',str,'  κ = ',num2str(value(4))]},'FontSize',12,'FontWeight','bold','FontName','phong');
    xlabel('Position','FontWeight','bold');
    ylabel(yname,'FontWeight','bold');
end
