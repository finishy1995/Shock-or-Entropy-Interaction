%
%  ���ߣ���Ԫ��  ���ڣ�2016-10-22
%  �����ļ�����
%

function [value, RVEP]=getPlot(filename)
    %%  ��ȡͼ������
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
