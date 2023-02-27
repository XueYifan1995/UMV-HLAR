clc,clear,close all

%% Importing Data
load HSVACPMCKVLCC2Z1005 HSVACPMCKVLCC2Z1005    %Data->Variable Name
load HSVACPMCKVLCC2Z1001 HSVACPMCKVLCC2Z1001
load HSVACPMCKVLCC2Z1505 HSVACPMCKVLCC2Z1505
load HSVACPMCKVLCC2Z2005 HSVACPMCKVLCC2Z2005
load HSVACPMCKVLCC2Z2505 HSVACPMCKVLCC2Z2505
load HSVACPMCKVLCC2Z3005 HSVACPMCKVLCC2Z3005
load HSVACPMCKVLCC2Z3505 HSVACPMCKVLCC2Z3505
load HSVACPMCKVLCC2Z1010P HSVACPMCKVLCC2Z1010P 
load HSVACPMCKVLCC2Z2010P HSVACPMCKVLCC2Z2010P 
load HSVACPMCKVLCC2Z1001 HSVACPMCKVLCC2Z1001
% %The meaning of each column in the data
% t = data(:,1);
% x = data(:,2);
% y = data(:,3);
% psi = data(:,4);   %Bow angle
% u = data(:,5);
% v = data(:,6);
% r = data(:,7);
% phi = data(:,8);   %Transverse roll angle
% d =  data(:,9);   %Rudder angle

x = [HSVACPMCKVLCC2Z1005(1:3100,:);HSVACPMCKVLCC2Z2005(1:3100,:);HSVACPMCKVLCC2Z3005(1:3100,:)];
u = x(:,9)*pi/180;
va = x(:,5)-1.179*ones(size(x(:,5)));
x = [va x(:,6) x(:,7)*pi/180]; 
xaug = [x(3:end-3,:) u(3:end-3,:)];

%% Plotting raw data
figure
plot(xaug(:,1),'linewidth',6,'color',[0,0.45,0.74])  %velocity of u
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp1u','-r600','-dpng');  %print figures

figure
plot(xaug(:,2),'linewidth',6,'color',[0.93,0.69,0.13])  %velocity of v
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp2v','-r600','-dpng');  

figure
plot(xaug(:,2),'linewidth',6,'color',[0.47,0.67,0.19])  %velocity of v
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp22v','-r600','-dpng');  

figure
plot(xaug(:,3),'linewidth',6,'color',[0.47,0.67,0.19])  %velocity of r
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp3r','-r600','-dpng');  

figure
plot(xaug(:,3),'linewidth',6,'color',[0.93,0.69,0.13])  %velocity of r
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp33r','-r600','-dpng');  

figure
plot(xaug(:,4),'linewidth',6,'color',[0,0,0])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp4d','-r600','-dpng');  

figure
plot(xaug(:,4),'linewidth',6,'color',[1.00,0.41,0.16])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp44d','-r600','-dpng');  

figure
plot(xaug(:,1).*xaug(:,1),'linewidth',6,'color',[0,0.45,0.74])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp5uu','-r600','-dpng'); 

figure
plot(xaug(:,1).*xaug(:,3),'linewidth',6,'color',[1.00,0.41,0.16])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp6ur','-r600','-dpng');

figure
plot(xaug(:,1).*xaug(:,4),'linewidth',6,'color',[0.47,0.67,0.19])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp7ud','-r600','-dpng');  

figure
plot(xaug(:,2).*xaug(:,2),'linewidth',6,'color',[0.84,0.84,0.84])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp8vv','-r600','-dpng');  

figure
plot(xaug(:,2).*xaug(:,3),'linewidth',6,'color',[0,0.45,0.74])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp9vr','-r600','-dpng');  

figure
plot(xaug(:,3).*xaug(:,3),'linewidth',6,'color',[0,0.45,0.74])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp10rr','-r600','-dpng');  

figure
plot(xaug(:,4).*xaug(:,4),'linewidth',6,'color',[0,0.45,0.74])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp11dd','-r600','-dpng');  

figure
plot(xaug(:,2).*xaug(:,2).*xaug(:,2),'linewidth',6,'color',[0.84,0.84,0.84])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp12vvv','-r600','-dpng');  

figure
plot(xaug(:,3).*xaug(:,3).*xaug(:,3),'linewidth',6,'color',[0.47,0.67,0.19])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp13rrr','-r600','-dpng');  

figure
plot(xaug(:,2).*xaug(:,2).*xaug(:,3),'linewidth',6,'color',[0.84,0.84,0.84])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp14vvr','-r600','-dpng');  

figure
plot(xaug(:,2).*xaug(:,3).*xaug(:,3),'linewidth',6,'color',[0.47,0.67,0.19])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Exp15vrr','-r600','-dpng');  