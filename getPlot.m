%
%  作者：王元恺  日期：2016-10-22
%  读入文件数据
%

function [value, RVEP]=getPlot(filename)
    %%  获取图像数据
    fid=fopen([filename,'.txt'],'r');
    value = [];
    tline = fgetl(fid);
    tline = str2num(tline);
    value = tline;
    RVEP = [];
    for i = 1:4
        tline = fgetl(fid);
        tline = str2num(tline);
        RVEP = [RVEP;tline];
    end
    fclose(fid);
end
