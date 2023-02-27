clc,clear,close all

%% Importing Data
load HSVACPMCKVLCC2Z1005 HSVACPMCKVLCC2Z1005    
load HSVACPMCKVLCC2Z1001 HSVACPMCKVLCC2Z1001
load HSVACPMCKVLCC2Z1505 HSVACPMCKVLCC2Z1505
load HSVACPMCKVLCC2Z2005 HSVACPMCKVLCC2Z2005
load HSVACPMCKVLCC2Z2505 HSVACPMCKVLCC2Z2505
load HSVACPMCKVLCC2Z3005 HSVACPMCKVLCC2Z3005
load HSVACPMCKVLCC2Z3505 HSVACPMCKVLCC2Z3505
load HSVACPMCKVLCC2Z1010P HSVACPMCKVLCC2Z1010P 
load HSVACPMCKVLCC2Z2010P HSVACPMCKVLCC2Z2010P 
load HSVACPMCKVLCC2Z1001 HSVACPMCKVLCC2Z1001
% t = data(:,1);
% x = data(:,2);
% y = data(:,3);
% psi = data(:,4);   
% u = data(:,5);
% v = data(:,6);
% r = data(:,7);
% phi = data(:,8);   
% d =  data(:,9);   

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
print('velocity of u','-r600','-dpng');  
figure
plot(xaug(:,2),'linewidth',6,'color',[0.93,0.69,0.13])  %velocity of v
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('velocity of v','-r600','-dpng');  
figure
plot(xaug(:,3),'linewidth',6,'color',[0.47,0.67,0.19])  %velocity of r
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('velocity of r','-r600','-dpng'); 
figure
plot(xaug(:,4),'linewidth',6,'color',[0,0,0])  %excitation signal
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('excitation signal','-r600','-dpng');  

%% Plot the derivative of the data
dt=0.05; %In order to keep the prediction process consistent with the time in the actual data
dx = zeros(length(x)-5,3);
for i=3:length(x)-3
        for k=1:size(x,2)
            dx(i-2,k) = (1/(12*dt))*(-x(i+2,k)+8*x(i+1,k)-8*x(i-1,k)+x(i-2,k));   
        end
    end
figure
plot(dx(:,1),'linewidth',1.5,'color',[0,0.45,0.74])  %velocity of u
box off
axis([0 9500 -0.025 0.025])
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Acceleration of u','-r600','-dpng'); 

figure
plot(dx(:,2),'linewidth',1.5,'color',[0.93,0.69,0.13])  %velocity of v
box off
axis([0 9500 -0.1 0.1])
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Acceleration of v','-r600','-dpng');  

figure
plot(dx(:,3),'linewidth',1.5,'color',[0.47,0.67,0.19])  %velocity of r
box off
axis([0 9500 -0.015 0.015])
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Acceleration of r','-r600','-dpng');  