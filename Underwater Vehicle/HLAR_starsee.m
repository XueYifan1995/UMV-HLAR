clc,clear,close all
%% 导入数据
load ZigZag10 ZigZag10
load ZigZag20 ZigZag20
load D2ZigZag25 D2ZigZag25
%% 最优计算
% lamda = [0.5124, 1.954];
lamda = [0.05124, 0.1954];

% lamda = [0.05, 0.1];
x_train = [ZigZag10;ZigZag20];
x_vali = D2ZigZag25;
x = x_train;
x = wdenoise(x);
u = x(:,9)*pi/180;
va = x(:,7);
x = [va x(:,6) x(:,1)*pi/180]; 

% dt=0.02; %In order to align the prediction process with the time in the actual data
dt=0.2;
dx = zeros(length(x)-5,3);
for i=3:length(x)-3
        for k=1:size(x,2)
            dx(i-2,k) = (1/(12*dt))*(-x(i+2,k)+8*x(i+1,k)-8*x(i-1,k)+x(i-2,k)); 
        end
    end
% concatenate

xaug = [x(3:end-3,:) u(3:end-3,:)];
dx(:,size(x,2)+1) = 0*dx(:,size(x,2));


lamda2=lamda(1);
lamda3=lamda(2);
n = 4;
% Sparse regression
LibraryType = 5; %1 is the initial dictionary base, 2 is the dictionary base constructed according to the Fossen model
Theta = selfpooldata(xaug,LibraryType);
Theta_norm = zeros(size(Theta,2),1); %zeros(size(Theta,2),1);ones
for i = 1:size(Theta,2)
   Theta_norm(i) = norm(Theta(:,i));
   Theta(:,i) = Theta(:,i)./Theta_norm(i);
end
m = size(Theta,2);

% lambda_vec = [8.2,1.4,2];   %Set sparse knob Adjusts according to Bayesian optimization results
lambda_vec = [10,lamda2,lamda3];

if exist('lambda_vec') == 1
    Xi = sparsifyDynamicsIndependent(Theta,dx,lambda_vec,n-1);
else
    Xi = sparsifyDynamics(Theta,dx,lambda,n-1);
end

for i = 1:size(Theta,2)  %Normalise the coefficient matrix
   Xi(i,:) = Xi(i,:)./Theta_norm(i);
end

yout = selfpooldatalist(Xi,LibraryType);
%% Make a forecast

n=4;
Nvar = 3;
datapre = x_vali;
data_pre = wdenoise(datapre);
dtt = 0.2;

tspan=[0];
for k=1:size(data_pre,1)-1
    t =dtt*k;
    tspan = [tspan,t];
end

xv = [data_pre(1,7),data_pre(1,6),data_pre(1,1)*pi/180];
psi_pre = data_pre(1,10);  %存储艏向角信息
x_p(1,:)=xv;
u_p = data_pre(:,9)*pi/180; 

for k=1:size(data_pre,1)-1   %Prediction using Eulerian dispersion methods
    y=[x_p(k,:) u_p(k)];
    xPool = selfpooldata(y,LibraryType);
    dxPool = xPool*Xi(:,1:Nvar);
    x_p(k+1,:) = x_p(k,:)+(dt*dxPool) ;      %Update next status
    psi_pre(k+1) = psi_pre(k)+x_p(k+1,3)*dtt*180/pi;
end

x_p(:,3) = x_p(:,3)*180/pi;

%% 绘制预测验证图-----------------------------------------------------------
load r_KT r_KT
load psi_KT psi_KT

figure
subplot(411)
h1 = plot(tspan,x_p(:,2),'--','linewidth',1.5,'color',[0.93,0.69,0.13]);
hold on
h2 = plot(tspan,data_pre(:,6),'linewidth',1.5,'color',[0.15,0.15,0.15]);
grid on,box off
xlabel('time(s)','Fontsize',14)
ylabel('v(m/s)','Fontsize',14)

subplot(412)
h1 = plot(tspan,x_p(:,3),'--','linewidth',1.5,'color',[0.93,0.69,0.13]);
hold on
h2 = plot(tspan,data_pre(:,1),'linewidth',1.5,'color',[0.15,0.15,0.15]);
hold on
h3 = plot(tspan,r_KT*180/pi,'-.','linewidth',1.5,'color',[0.47,0.67,0.19]);
grid on,box off
xlabel('time(s)','Fontsize',14)
ylabel('r(°/s)','Fontsize',14)

subplot(413)
h1 = plot(tspan,psi_pre,'--','linewidth',1.5,'color',[0.93,0.69,0.13]);
hold on
h2 = plot(tspan,data_pre(:,10),'linewidth',1.5,'color',[0.15,0.15,0.15]);
hold on
h3 = plot(tspan,psi_KT*180/pi,'-.','linewidth',1.5,'color',[0.47,0.67,0.19]);
grid on,box off
xlabel('time(s)','Fontsize',14)
ylabel('psi(°)','Fontsize',14)

subplot(414)
h4 = plot(tspan,data_pre(:,9),'linewidth',1.5,'color',[0.8500 0.3250 0.0980]);
axis([0 inf -30 30])
grid on,box off
xlabel('time(s)','Fontsize',14)
ylabel('\delta(°/s)','Fontsize',14)
legend([h1,h2,h3,h4],'Proposed HLAR method','Experiment','semi- K-T equation','\delta','Fontsize',12)
set(gcf, 'Position', [100, 100, 800, 600]); % 设置Figure窗口的位置和大小
% print('2KT and sindy graph','-r600','-dpng')

%% --------------------------------------------------------------------------
%计算均方根误差
rmse_v = sqrt(mean((x_p(:,2)-data_pre(:,6)).^2))
rmse_r = sqrt(mean((x_p(:,3)-data_pre(:,1)).^2))*pi/180
rmse_psi = sqrt(mean((psi_pre'-data_pre(:,10)).^2))*pi/180
