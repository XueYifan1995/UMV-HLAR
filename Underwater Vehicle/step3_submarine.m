
data = npsAUV_zigzag_2505;
set(0,'defaultfigurecolor','w')

u = data(1,:);
v = data(2,:);
r = data(3,:);
d = data(4,:);
%Data Mapping
figure(1)
subplot(311),plot(T,u,'linewidth',2),xlabel('time (s)'),grid on;hold on
title('result of submarine')
subplot(312),plot(T,v,'linewidth',2),xlabel('time (s)'),grid on;hold on
subplot(313),plot(T,r,'linewidth',2),xlabel('time (s)'),grid on;hold on
