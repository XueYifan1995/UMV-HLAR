echo on
% Use the container simulation noisy data to trian Noise Input Gaussian process 
%
% Calls:       Infante_NMTGP_3dof , Adjust angel
%
% Author:      Yifan Xue
% Date:        2020-01-15
% Revisions: 
echo off
set(0,'defaultfigurecolor','w')
%% Step1:load data 

load HSVACPMCKVLCC2Z1005 HSVACPMCKVLCC2Z1005
load HSVACPMCKVLCC2Z2005 HSVACPMCKVLCC2Z2005
load HSVACPMCKVLCC2Z3005 HSVACPMCKVLCC2Z3005

h=  0.05; sample= 9;
data_raw = [HSVACPMCKVLCC2Z1005(1:3200,:);HSVACPMCKVLCC2Z2005(1:3200,:);HSVACPMCKVLCC2Z3005(1:3200,:)];
% data = wdenoise(data_raw);
data = data_raw;
num_tr = size(data,1);
t=linspace(0,num_tr,num_tr+1).*h;
data(:,1)=t(1:end-1);
scatter(data(:,1),[data_raw(:,6),data(:,6)]);

data_train = data;
t =data_train(:,1);
psi = data_train(:,4)*pi/180;
u = data_train(:,5);
v = data_train(:,6);
r = data_train(:,7)*pi/180;
phi = data_train(:,8)*pi/180;
d =  data_train(:,9)*pi/180;

%% Step2 Construct data
tic
u_x = u(1:end-1);  u_y = u(2:end);
v_x = v(1:end-1);  v_y = v(2:end);
r_x = r(1:end-1);  r_y = r(2:end);
d_x = d(1:end-1);
Xm = [u_x,v_x,r_x,d_x];
Ym = [u_y,v_y,r_y];
% noise test
Xm2=[u_x,v_x,r_x];
Am = (Ym-Xm2)/h;
Am_d = wdenoise(Am);
scatter(t(1:end-1),[Am(:,2),Am_d(:,2)]);axis([0 490 -0.05 0.05]);
%间隔
t_t =t(1:sample:end-1);
Xm_t= Xm(1:sample:num_tr,:);
Am_t = Am(1:sample:num_tr,:);
Am_d_t = Am_d(1:sample:num_tr,:);
scatter(t_t,[Am_t(:,2),Am_d_t(:,2)]);axis([0 490 -0.05 0.05]);

SVM_u = fitrsvm(Xm_t,Am_d_t(:,1),'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus','MaxObjectiveEvaluations',3));

SVM_v = fitrsvm(Xm_t,Am_d_t(:,2),'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus','MaxObjectiveEvaluations',5));
SVM_r = fitrsvm(Xm_t,Am_d_t(:,3),'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus','MaxObjectiveEvaluations',3));
t5 = toc;  %time for training
%% Step3 : Apply SVM to predict
tic
dt= h*sample;
total_time = 160;
m= ceil(total_time /dt);    %节拍 
x = zeros(3,1); %临时状态变量
U_a = zeros(3,1);
TEMP_a = zeros(4,1); 
T = zeros(m,1);  %时间
Y = zeros(m,8); %状态变量
Uci = zeros(m,1);Vci = zeros(m,1);Rci = zeros(m,1);
%Inatialize the ship 
u0 = 1.175; v0 = 0;  r0=0;
x0 = 0; y0 = 0; psi0 = 0;
d0 = 0;flag_lr0= -1;
Initial_input = [u0;v0;r0];
Initial_ob =  [x0; y0; psi0;d0];

U = Initial_input;
x = Initial_ob;
sy = flag_lr0;

for i=1:1:m
    t = dt*i;
    T(i,1)=t;
    time = t;
    [TEMP_a,TEMP_sm,U_a]=Infante_SVM_3dof_HVSA_ac(U,x,sy,SVM_u,SVM_v,SVM_r)  ;
    %Euler
    U = U + dt.*U_a;
    x= x + dt.*TEMP_a;
    sy =TEMP_sm;
    
    %保存数据
    Y(i,1) = U(1);%u
    Y(i,2) = U(2);%v
    Y(i,3) = U(3);%r
    Y(i,5) = x(1);%x
    Y(i,6) = x(2);%y
    Y(i,7) = x(3);%psi
    Y(i,9) = x(4);%舵角   
    Y(i,10) =sy;%舵角状态
end

U_pre  = Y(:,1);
V_pre  = Y(:,2);
R_pre  = Y(:,3)*180/pi;
Xp = Y(:,5);
Yp = Y(:,6);
psi = Y(:,7);
duo = Y(:,9);

duo = duo*180/pi;
psi = psi *180/pi;

t8 = toc; % time for prediction

HVSA_z30_5_BSVM = [T,U_pre,V_pre,R_pre];
save HVSA_z30_5_BSVM HVSA_z30_5_BSVM ;


figure(3)
subplot(411),plot(T,U_pre,'linewidth',1.5),xlabel('time (s)'),ylabel('u (m/s)'),title('speed U (m/s)'),grid on;hold on
subplot(412),plot(T,V_pre,'linewidth',1.5),xlabel('time (s)'),ylabel('v (m/s)'),grid on;hold on
subplot(413),plot(T,R_pre,'linewidth',1.5),xlabel('time (s)'),ylabel('r (rad/s)'),grid on;hold on
subplot(414),plot(T,psi,'linewidth',1.6);hold on
% plot(T,duo);hold off
xlabel('time (s)'),title('yaw angle \psi (deg)'),grid on;hold on
% legend('\psi','\delta_c');
