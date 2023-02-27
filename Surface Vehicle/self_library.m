%% Import data
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
% The meaning of each column in the data
% t = data(:,1);
% x = data(:,2);
% y = data(:,3);
% psi = data(:,4);   %Bow angle
% u = data(:,5);
% v = data(:,6);
% r = data(:,7);
% phi = data(:,8);   %Roll angle
% d =  data(:,9);   %Rudder angle

%% Calculating the SINDy model
x = [HSVACPMCKVLCC2Z1005(1:3100,:);HSVACPMCKVLCC2Z2005(1:3100,:);HSVACPMCKVLCC2Z3005(1:3100,:)];
u = x(:,9)*pi/180;
va = x(:,5)-1.179*ones(size(x(:,5)));
x = [va x(:,6) x(:,7)*pi/180]; 

dt=0.05; %In order to align the prediction process with the time in the actual data
dx = zeros(length(x)-5,3);
for i=3:length(x)-3
        for k=1:size(x,2)
            dx(i-2,k) = (1/(12*dt))*(-x(i+2,k)+8*x(i+1,k)-8*x(i-1,k)+x(i-2,k)); 
        end
    end
% concatenate

xaug = [x(3:end-3,:) u(3:end-3,:)];
dx(:,size(x,2)+1) = 0*dx(:,size(x,2));

n = size(dx,2);

% Sparse regression
clear Theta Xi
LibraryType = 1; %1 is the initial dictionary base, 2 is the dictionary base constructed according to the Fossen model
Theta = selfpooldata(xaug,LibraryType);
Theta_norm = zeros(size(Theta,2),1); %zeros(size(Theta,2),1);ones
for i = 1:size(Theta,2)
   Theta_norm(i) = norm(Theta(:,i));
   Theta(:,i) = Theta(:,i)./Theta_norm(i);
end
m = size(Theta,2);

lambda_vec = [0.15668,0.34913,0.096613];   %Set sparse knob Adjusts according to Bayesian optimization results


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
tspan=[0];
data_pre =  HSVACPMCKVLCC2Z3505;
xv = [data_pre(1,5)-1.179,data_pre(1,6),data_pre(1,7)*pi/180];
x_p(1,:)=xv;
u_p = data_pre(:,9)*pi/180; 

for k=1:size(data_pre,1)-1
    t =dt*k;
    tspan = [tspan,t];
end

tic
for k=1:size(data_pre,1)-1   %Prediction using Eulerian dispersion methods
    y=[x_p(k,:) u_p(k)];
    xPool = selfpooldata(y,LibraryType);
    dxPool = xPool*Xi(:,1:Nvar);
    x_p(k+1,:) = x_p(k,:)+(dt*dxPool) ;      %Update next status
end
t_sindy = toc;
x_p(:,1) = x_p(:,1)+1.179*ones(size(x_p(:,1)));
x_p(:,3) = x_p(:,3)*180/pi;

%% Calculation of non-parametric Gaussian models
GuassPreShip


%% Calculation of parametric models

% Step1_dataconstruct3_new
% Step2_main_new

load U_pre_35 U_pre_35
load V_pre_35 V_pre_35
load R_pre_35 R_pre_35
U_pre_semi = U_pre_35;
V_pre_semi = V_pre_35;
R_pre_semi = R_pre_35;



%% Drawing
figure
subplot(3,1,1)
patch([T', fliplr(T')], [u_lower', fliplr(u_upper')], 1, 'FaceColor', [0.85,0.85,1], 'EdgeColor', 'none');hold on;
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
h1 = plot(T,U_pre,'LineWidth',1.5,'color',[0,0.45,0.74]);xlabel('time (s)'),ylabel('u (m/s)');grid on;hold on;
h4 = plot(tspan(1:3600),data_pre(1:3600,5),'linewidth',1.5,'color',[0.15,0.15,0.15]);
hold on
h2 = plot(tspan(1:3600),x_p(1:3600,1),'--','linewidth',1.5,'color',[0.93,0.69,0.13]);
hold on
h3 = plot(tspan(1:3600),U_pre_semi,'-.','linewidth',1.5,'color',[0.47,0.67,0.19]);
legend([h1,h2,h3,h4],'Gaussian process','Proposed method','Semi-Abkowitz','Experiment')
axis([0 180 -inf inf])


subplot(3,1,2)
patch([T', fliplr(T')], [v_lower', fliplr(v_upper')], 1, 'FaceColor', [0.85,0.85,1], 'EdgeColor', 'none');hold on
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(T,V_pre,'linewidth',1.5,'color',[0,0.45,0.74]),xlabel('time (s)'),ylabel('v (m/s)');grid on;hold on
plot(tspan(1:3600),data_pre(1:3600,6),'linewidth',1.5,'color',[0.15,0.15,0.15])
hold on
plot(tspan(1:3600),x_p(1:3600,2),'--','linewidth',1.5,'color',[0.93,0.69,0.13])
hold on
plot(tspan(1:3600),V_pre_semi,'-.','linewidth',1.5,'color',[0.47,0.67,0.19])
axis([0 180 -inf inf])

subplot(3,1,3)
patch([T', fliplr(T')], [r_lower', fliplr(r_upper')], 1, 'FaceColor', [0.85,0.85,1], 'EdgeColor', 'none');hold on
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(T,R_pre*pi/180,'linewidth',1.5,'color',[0,0.45,0.74]),xlabel('time (s)'),ylabel('r (rad/s)');grid on;hold on
plot(tspan(1:3600),data_pre(1:3600,7)*pi/180,'linewidth',1.5,'color',[0.15,0.15,0.15])
hold on
plot(tspan(1:3600),x_p(1:3600,3)*pi/180,'--','linewidth',1.5,'color',[0.93,0.69,0.13])
hold on
plot(tspan(1:3600),R_pre_semi*pi/180,'-.','linewidth',1.5,'color',[0.47,0.67,0.19])
axis([0 180 -0.1 0.1])

%% Calculate relative errors
rmse_u = sqrt(mean((x_p(:,1)-data_pre(:,5)).^2))
rmse_v = sqrt(mean((x_p(:,2)-data_pre(:,6)).^2))
rmse_r = sqrt(mean((x_p(:,3)*pi/180-data_pre(:,7)*pi/180).^2))

rmse_u_semi = sqrt(mean((U_pre_semi-data_pre(1:3600,5)).^2))
rmse_v_semi = sqrt(mean((V_pre_semi-data_pre(1:3600,6)).^2))
rmse_r_semi = sqrt(mean((R_pre_semi*pi/180-data_pre(1:3600,7)*pi/180).^2))

rmse_u_guass = sqrt(mean((U_pre-data_pre(13:12:3601,5)).^2))
rmse_v_guass = sqrt(mean((V_pre-data_pre(13:12:3601,6)).^2))
rmse_r_guass = sqrt(mean((R_pre*pi/180-data_pre(13:12:3601,7)*pi/180).^2))