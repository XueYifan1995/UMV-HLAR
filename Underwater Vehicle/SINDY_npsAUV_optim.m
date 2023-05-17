clc,close all, clear
%For each of the three methods, it is necessary to confirm that each method 
% becomes the correct set of predictions when changing the prediction set
load npsAUV_zigzag_1010_005p npsAUV_zigzag_1010_005p
load npsAUV_zigzag_2020_005p npsAUV_zigzag_2020_005p
load npsAUV_zigzag_2505_005 npsAUV_zigzag_2505_005
load npsAUV_zigzag_2510_005 npsAUV_zigzag_2510_005
%% Get test data (ZigZag movement)
load npsAUV_zigzag_test npsAUV_zigzag_test
x_1 = npsAUV_zigzag_test;
%% Calculating the SINDY model
% GetModelZigzag  %Get the best lamda
xinput = [npsAUV_zigzag_1010_005p,npsAUV_zigzag_2020_005p];
x_1_self = xinput(1:3,:);
u = xinput(4,:)*pi/180;

%% Calculate sindy initial parameters
x_1_self = x_1_self';
dt = 0.05;
dx = zeros(length(x_1_self)-5,3);
for i=3:length(x_1_self)-3
        for k=1:size(x_1_self,2)
            dx(i-2,k) = (1/(12*dt))*(-x_1_self(i+2,k)+8*x_1_self(i+1,k)-8*x_1_self(i-1,k)+x_1_self(i-2,k));   
        end
    end
% concatenate
u=u';
xaug = [x_1_self(3:end-3,:) u(3:end-3,:)];
dx(:,size(x_1_self,2)+1) = 0*dx(:,size(x_1_self,2));

n = size(dx,2);

% Sparse regression
polyorder=3;   
usesine = 0;   
LibraryType = 4;
Theta = selfpooldata(xaug,LibraryType);   
Theta_norm = zeros(size(Theta,2),1); 
for i = 1:size(Theta,2)
   Theta_norm(i) = norm(Theta(:,i));
   Theta(:,i) = Theta(:,i)./Theta_norm(i);
end
m = size(Theta,2);
lambda_vec = [0.39217,0.10735,0.10254];

if exist('lambda_vec') == 1
    Xi = sparsifyDynamicsIndependent(Theta,dx,lambda_vec,n-1);
else
    Xi = sparsifyDynamics(Theta,dx,lambda,n-1);
end


for i = 1:size(Theta,2)
   Xi(i,:) = Xi(i,:)./Theta_norm(i);
end
yout = selfpooldatalist(Xi,LibraryType);
Xi_sphs = Xi;

%% Forecasting with Eulerian dispersion - SINDY
GetModel_validation

%% Using the semi-conjugate plus PhD thesis model
load U_pre_doc U_pre_doc
load V_pre_doc V_pre_doc
load R_pre_doc R_pre_doc
rmse_u_doc = sqrt(mean((x_1(1,:)-U_pre_doc').^2))
rmse_v_doc = sqrt(mean((x_1(2,:)-V_pre_doc').^2))
rmse_r_doc = sqrt(mean((x_1(3,:)-R_pre_doc').^2))


%% Using the semi-conjugate Gashisenda model
load U_pre_submarine U_pre_submarine
load V_pre_submarine V_pre_submarine
load R_pre_submarine R_pre_submarine

rmse_u_submarine = sqrt(mean((x_1(1,:)-U_pre_submarine').^2))
rmse_v_submarine = sqrt(mean((x_1(2,:)-V_pre_submarine').^2))
rmse_r_submarine = sqrt(mean((x_1(3,:)-R_pre_submarine').^2))


%% Final results show
figure
subplot(3,1,1)
h3 = plot(tspan,x_1(1,:),'linewidth',1.5,'color',[0.15,0.15,0.15]);
hold on
h1 = plot(tspan,xp(1,:),'--','linewidth',1.5,'color',[0.93,0.69,0.13]);
xlabel('time (s)'),ylabel('u (m/s)');
hold on
h2 = plot(tspan,U_pre_submarine,'-.','linewidth',1.5,'color',[0.47,0.67,0.19]);
h4 = plot(tspan,U_pre_doc,'linewidth',1.5,'color',[0,0.45,0.74]);
legend([h1,h2,h3,h4],{'Proposed HLAR Method','Semi-Silvestre','Experiment','Semi-Gertler'})
axis([0 inf 1 1.6])
grid on
box off

subplot(3,1,2)
plot(tspan,x_1(2,:),'linewidth',1.5,'color',[0.15,0.15,0.15])
hold on
plot(tspan,xp(2,:),'--','linewidth',1.5,'color',[0.93,0.69,0.13])
hold on
plot(tspan,V_pre_submarine,'-.','linewidth',1.5,'color',[0.47,0.67,0.19])
hold on
plot(tspan,V_pre_doc,'linewidth',1.5,'color',[0,0.45,0.74]);
xlabel('time (s)'),ylabel('v (m/s)');
grid on
box off


subplot(3,1,3)
plot(tspan,x_1(3,:),'linewidth',1.5,'color',[0.15,0.15,0.15])
hold on
plot(tspan,xp(3,:),'--','linewidth',1.5,'color',[0.93,0.69,0.13])
hold on
plot(tspan,R_pre_submarine,'-.','linewidth',1.5,'color',[0.47,0.67,0.19])
hold on
plot(tspan,R_pre_doc,'linewidth',1.5,'color',[0,0.45,0.74]);
xlabel('time (s)'),ylabel('r (rad/s)');
grid on
box off
axis([0 inf -0.04 0.04])