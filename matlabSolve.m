%
%  ���ߣ���Ԫ��  ���ڣ�2016-11-13
%  matlab����
%

%% ���ó���
n = 201;

%% �Ա�����MUSCL Limiter
[valueN, nonMuscl]=getPlot('nonMuscl');
[valueH, hasMuscl]=getPlot('hasMuscl');
deltaX = (valueN(2)-valueN(1))/(n-1);
deltaY = (valueH(2)-valueH(1))/(n-1);
figure(1)
hold on;
plot(valueN(1):deltaX:valueN(2),nonMuscl(1,:),'ro','Markersize',4);
plot(valueH(1):deltaY:valueH(2),hasMuscl(1,:),'b','LineWidth',1.5);
title('�ܶȷֲ�ͼ','FontSize',12,'FontWeight','bold','FontName','phong');
xlabel('Position','FontWeight','bold');
ylabel('Density','FontWeight','bold');
legend('nonMuscl','hasMuscl');
print('-dpng','ifMusclDensity');
hold off;
close(1);

figure(2)
hold on;
plot(valueN(1):deltaX:valueN(2),nonMuscl(2,:),'ro','Markersize',4);
plot(valueH(1):deltaY:valueH(2),hasMuscl(2,:),'b','LineWidth',1.5);
title('�ٶȷֲ�ͼ','FontSize',12,'FontWeight','bold','FontName','phong');
xlabel('Position','FontWeight','bold');
ylabel('Velocity','FontWeight','bold');
legend('nonMuscl','hasMuscl');
print('-dpng','ifMusclVelocity');
hold off;
close(2);

figure(3)
hold on;
plot(valueN(1):deltaX:valueN(2),nonMuscl(3,:),'ro','Markersize',4);
plot(valueH(1):deltaY:valueH(2),hasMuscl(3,:),'b','LineWidth',1.5);
title('���ֲܷ�ͼ','FontSize',12,'FontWeight','bold','FontName','phong');
xlabel('Position','FontWeight','bold');
ylabel('Energy','FontWeight','bold');
legend('nonMuscl','hasMuscl');
print('-dpng','ifMusclEnergy');
hold off;
close(3);

figure(4)
hold on;
plot(valueN(1):deltaX:valueN(2),nonMuscl(4,:),'ro','Markersize',4);
plot(valueH(1):deltaY:valueH(2),hasMuscl(4,:),'b','LineWidth',1.5);
title('ѹ���ֲ�ͼ','FontSize',12,'FontWeight','bold','FontName','phong');
xlabel('Position','FontWeight','bold');
ylabel('Pressure','FontWeight','bold');
legend('nonMuscl','hasMuscl');
print('-dpng','ifMusclPressure');
hold off;
close(4);

%% �Ա�kappa���ʣ���Ӱ��
[value1, kappa1]=getPlot('kappa13');
[value2, kappa2]=getPlot('kappa1');
[value3, kappa3]=getPlot('kappa0');
[value4, kappa4]=getPlot('kappa-1');

figure(5)
subplot(2,2,1);
outPlot(value1,kappa1(1,:),'�ܶȷֲ�ͼ','Density');
subplot(2,2,2);
outPlot(value2,kappa2(1,:),'�ܶȷֲ�ͼ','Density');
subplot(2,2,3);
outPlot(value3,kappa3(1,:),'�ܶȷֲ�ͼ','Density');
subplot(2,2,4);
outPlot(value4,kappa4(1,:),'�ܶȷֲ�ͼ','Density');
print('-dpng','kappaDensity');
close(5);

figure(6)
subplot(2,2,1);
outPlot(value1,kappa1(2,:),'�ٶȷֲ�ͼ','Velocity');
subplot(2,2,2);
outPlot(value2,kappa2(2,:),'�ٶȷֲ�ͼ','Velocity');
subplot(2,2,3);
outPlot(value3,kappa3(2,:),'�ٶȷֲ�ͼ','Velocity');
subplot(2,2,4);
outPlot(value4,kappa4(2,:),'�ٶȷֲ�ͼ','Velocity');
print('-dpng','kappaVelocity');
close(6);

figure(7)
subplot(2,2,1);
outPlot(value1,kappa1(3,:),'���ֲܷ�ͼ','Energy');
subplot(2,2,2);
outPlot(value2,kappa2(3,:),'���ֲܷ�ͼ','Energy');
subplot(2,2,3);
outPlot(value3,kappa3(3,:),'���ֲܷ�ͼ','Energy');
subplot(2,2,4);
outPlot(value4,kappa4(3,:),'���ֲܷ�ͼ','Energy');
print('-dpng','kappaEnergy');
close(7);

figure(8)
subplot(2,2,1);
outPlot(value1,kappa1(4,:),'ѹ���ֲ�ͼ','Pressure');
subplot(2,2,2);
outPlot(value2,kappa2(4,:),'ѹ���ֲ�ͼ','Pressure');
subplot(2,2,3);
outPlot(value3,kappa3(4,:),'ѹ���ֲ�ͼ','Pressure');
subplot(2,2,4);
outPlot(value4,kappa4(4,:),'ѹ���ֲ�ͼ','Pressure');
print('-dpng','kappaPressure');
close(8);

%% �ԱȲ�ͬ��������Ӱ��
[value1, limiter1]=getPlot('vanLeer');
[value2, limiter2]=getPlot('vanAlbada');
[value3, limiter3]=getPlot('superbee');
[value4, limiter4]=getPlot('minmod');

figure(9)
subplot(2,2,1);
outPlot(value1,limiter1(1,:),'�ܶȷֲ�ͼ','Density');
subplot(2,2,2);
outPlot(value2,limiter2(1,:),'�ܶȷֲ�ͼ','Density');
subplot(2,2,3);
outPlot(value3,limiter3(1,:),'�ܶȷֲ�ͼ','Density');
subplot(2,2,4);
outPlot(value4,limiter4(1,:),'�ܶȷֲ�ͼ','Density');
print('-dpng','limiterDensity');
close(9);

figure(10)
subplot(2,2,1);
outPlot(value1,limiter1(2,:),'�ٶȷֲ�ͼ','Velocity');
subplot(2,2,2);
outPlot(value2,limiter2(2,:),'�ٶȷֲ�ͼ','Velocity');
subplot(2,2,3);
outPlot(value3,limiter3(2,:),'�ٶȷֲ�ͼ','Velocity');
subplot(2,2,4);
outPlot(value4,limiter4(2,:),'�ٶȷֲ�ͼ','Velocity');
print('-dpng','limiterVelocity');
close(10);

figure(11)
subplot(2,2,1);
outPlot(value1,limiter1(3,:),'���ֲܷ�ͼ','Energy');
subplot(2,2,2);
outPlot(value2,limiter2(3,:),'���ֲܷ�ͼ','Energy');
subplot(2,2,3);
outPlot(value3,limiter3(3,:),'���ֲܷ�ͼ','Energy');
subplot(2,2,4);
outPlot(value4,limiter4(3,:),'���ֲܷ�ͼ','Energy');
print('-dpng','limiterEnergy');
close(11);

figure(12)
subplot(2,2,1);
outPlot(value1,limiter1(4,:),'ѹ���ֲ�ͼ','Pressure');
subplot(2,2,2);
outPlot(value2,limiter2(4,:),'ѹ���ֲ�ͼ','Pressure');
subplot(2,2,3);
outPlot(value3,limiter3(4,:),'ѹ���ֲ�ͼ','Pressure');
subplot(2,2,4);
outPlot(value4,limiter4(4,:),'ѹ���ֲ�ͼ','Pressure');
print('-dpng','limiterPressure');
close(12);
