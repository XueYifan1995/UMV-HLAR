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

h=  0.05; sample= 15;
data_raw = [HSVACPMCKVLCC2Z1005(1:3100,:);HSVACPMCKVLCC2Z2005(1:3100,:);HSVACPMCKVLCC2Z3005(1:3100,:)];
% data = wdenoise(data_raw);
data = data_raw;
num_tr = size(data,1);
t=linspace(0,num_tr,num_tr+1).*h;
data(:,1)=t(1:end-1);
% scatter(data(:,1),[data_raw(:,6),data(:,6)]);

% data_train = data(1:sample:num_tr,:);
data_train = data;
t =data_train(:,1);
psi = data_train(:,4)*pi/180;
u = data_train(:,5);
v = data_train(:,6);
r = data_train(:,7)*pi/180;
phi = data_train(:,8)*pi/180;
d =  data_train(:,9)*pi/180;

%add noise
rng(0,'twister');
nm= size(u);
umax = max(u);vmax = max(v); rmax = max(r); dmax = max(d);
umin = min(u);vmin= min(v); rmin = min(v);  dmin = min(d);
sfm = 0.023; % noise level 0.06/10
u_n = u + umax*0.07*sfm*randn(nm);
v_n = v + 1.2*vmax*sfm*randn(nm);
r_n = r + 1.3*rmax*sfm*randn(nm);
d_n = d + 0.2*sfm*randn(nm);
u_xn = u(1:end-1);  u_yn = u(2:end);
v_xn = v(1:end-1);  v_yn = v(2:end);
r_xn = r(1:end-1);  r_yn = r(2:end);
d_xn = d(1:end-1);
%% Step2 Construct data
%construct the training points
tic
u_x = u(1:end-1);  u_y = u(2:end);
v_x = v(1:end-1);  v_y = v(2:end);
r_x = r(1:end-1);  r_y = r(2:end);
d_x = d(1:end-1);
Xm = [u_x,v_x,r_x,d_x];
Xm_n = [u_xn,v_xn,r_xn,d_xn];
Ym = [u_y,v_y,r_y];
% noise test
Xm2=[u_x,v_x,r_x];
Am = (Ym-Xm2)/h;
Am_d = wdenoise(Am,4);
% scatter(t(1:end-1),[Am(:,2),Am_d(:,2)]);axis([0 490 -0.05 0.05]);
%间隔
t_t =t(1:sample:end-1);
Xm_t= Xm(1:sample:num_tr,:);
Xm_tn= Xm_n(1:sample:num_tr,:);
Am_t = Am(1:sample:num_tr,:);
Am_d_t = Am_d(1:sample:num_tr,:);
% scatter(t_t,[Am_t(:,2),Am_d_t(:,2)]);axis([0 490 -0.05 0.05]);

[model,f] = trainNIGP(Xm_t,Am_t,-200,1); % train the multi output u-NARX GP
% [model,f] = trainNIGP(Xm_t,Am_d_t,-200,1); % train the multi output u-NARX GP
t1 = toc;  %time for training
%% Step3 : Apply NARX GP to predict
load HSVACPMCKVLCC2Z1505 HSVACPMCKVLCC2Z1505 
st=20;
d_pre = HSVACPMCKVLCC2Z1505(st:end,9)*pi/180;
delta_ex= d_pre(1:sample:end,:);
tic
dt= h*sample;
total_time = 160;
m= ceil(total_time /dt);    %节拍 
U_a = zeros(3,1);
x = zeros(3,1); %临时状态变量
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
    delta_exin=delta_ex(i);
    % one Multi-output NARX
    [TEMP_a,TEMP_sm,U_a,uci,vci,rci]=Infante_NIMTGP_3dof_HVSA_ac(time,U,x,sy,model,delta_exin)  ;
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
    Y(i,9) = x(4);%rudder   
    Y(i,10) =sy;%舵角状态
    Uci(i) = uci; 
    Vci(i) = vci; 
    Rci(i) = rci*180/pi; 
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

t2 = toc; % time for prediction

% HVSA_z30_5_NIGP = [T,U_pre,V_pre,R_pre];
% save HVSA_z30_5_NIGP HVSA_z30_5_NIGP ;

figure(3)
kk= 2;
subplot(311),
u_upper = U_pre + kk*Uci;  % upper boundry
u_lower = U_pre - kk*Uci;  % lower boundrr
patch([T', fliplr(T')], [u_lower', fliplr(u_upper')], 1, 'FaceColor', [0.85,0.85,1], 'EdgeColor', 'none');hold on;
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(T,U_pre,'LineWidth',2),xlabel('time (s)');grid on;hold on;

subplot(312),
v_upper = V_pre + kk*Vci;  % upper boundry
v_lower = V_pre - kk*Vci;  % lower boundrr
patch([T', fliplr(T')], [v_lower', fliplr(v_upper')], 1, 'FaceColor', [0.85,0.85,1], 'EdgeColor', 'none');hold on
% set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(T,V_pre,'linewidth',2),xlabel('time (s)'),grid on;hold on

subplot(313),
r_upper = R_pre + 2*Rci;  % upper boundry
r_lower = R_pre - 2*Rci;  % lower boundrr
patch([T', fliplr(T')], [r_lower', fliplr(r_upper')], 1, 'FaceColor', [0.85,0.85,1], 'EdgeColor', 'none');hold on
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(T,R_pre,'linewidth',2),xlabel('time (s)'),grid on;hold on

% subplot(414),plot(T,psi,'linewidth',2);hold on
% % plot(T,duo);hold off
% xlabel('time (s)'),grid on;hold on
% % legend('\psi','\delta_c');
