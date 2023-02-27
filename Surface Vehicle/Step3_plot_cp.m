%% Step 1: Load data
load HSVACPMCKVLCC2Z3005 HSVACPMCKVLCC2Z3005
load HSVACPMCKVLCC2Z2505 HSVACPMCKVLCC2Z2505
load HSVACPMCKVLCC2Z1505 HSVACPMCKVLCC2Z1505
load HSVACPMCKVLCC2Z3505 HSVACPMCKVLCC2Z3505
load HSVACPMCKVLCC2Z1010P HSVACPMCKVLCC2Z1010P 
load HSVACPMCKVLCC2Z1001 HSVACPMCKVLCC2Z1001
data = HSVACPMCKVLCC2Z3505;
set(0,'defaultfigurecolor','w')
t = data(:,1);
x = data(:,2);
y = data(:,3);
psi = data(:,4);
u = data(:,5);
v = data(:,6);
r = data(:,7);
phi = data(:,8);
d =  data(:,9);
%Data Mapping
figure(3)
figure_FontSize=20;
set(findobj('FontSize',30),'FontSize',figure_FontSize);set(gcf,'Position',[100,100,400,480]);
subplot(414),plot(t,psi,'linewidth',2);hold on
subplot(411),plot(t,u,'linewidth',2),xlabel('time (s)'),grid on;hold on
subplot(412),plot(t,v,'linewidth',2),xlabel('time (s)'),grid on;hold on
subplot(413),plot(t,r,'linewidth',2),xlabel('time (s)'),grid on;hold on

%Calculation errors
rmse_u = sqrt(mean((U_pre-data(1:size(U_pre,1),5)).^2))
rmse_v = sqrt(mean((V_pre-data(1:size(V_pre,1),6)).^2))
rmse_r = sqrt(mean((R_pre-data(1:size(R_pre,1),7)).^2))