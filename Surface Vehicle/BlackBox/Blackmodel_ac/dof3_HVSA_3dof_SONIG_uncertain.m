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

load HSVACPMCKVLCC2Z2005 HSVACPMCKVLCC2Z2005
load HSVACPMCKVLCC2Z3505 HSVACPMCKVLCC2Z3505
data = [HSVACPMCKVLCC2Z2005;HSVACPMCKVLCC2Z3505];
[n2,~] = size (HSVACPMCKVLCC2Z3505);
num_tr = size(data,1);
h=  0.05; sample= 25;
data_train = data(1:sample:num_tr,:);

t1     = HSVACPMCKVLCC2Z2005(:,1);
t1_last = ones(n2,1)*t1(end);
t2     = t1_last + HSVACPMCKVLCC2Z3505(:,1);
t = [t1;t2]; t = t(1:sample:num_tr,:);

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
sfm = 0.001; % noise level 0.06/10
u_n = u + umax*0.07*sfm*randn(nm);
v_n = v + 1.2*vmax*sfm*randn(nm);
r_n = r + 1.3*rmax*sfm*randn(nm);
d_n = d + 0.2*sfm*randn(nm);

%plot the polluted identification data
figure(1)
subplot(414),plot(t,psi*180/pi,'linewidth',2);hold on
plot(t,d*180/pi);hold off
xlabel('time (s)'),title('yaw angle \psi (deg)'),grid;hold on
legend('\psi','\delta_c')
subplot(411),plot(t,[u,u_n],'linewidth',1.5),xlabel('time (s)'),title('speed U (m/s)'),grid on;hold on;
subplot(412),plot(t,[v,v_n],'linewidth',1.5),xlabel('time (s)'),title('speed V (m/s)'),grid on;hold on
subplot(413),plot(t,[r*180/pi,r_n*180/pi],'linewidth',1.5),xlabel('time (s)'),title('speed R (rad/s)'),grid on;hold on;
% scatter(t,r_n*180/pi,15);

figure(2)
t=t*2;
y1 = r*180/pi; y3 = r_n*180/pi; y2 = d_n*180/pi;
[AX,H1,H2]=plotyy(t,[y3,y1],t,y2,'plot');grid on;%双轴
set(AX(1),'XColor','k','YColor','k'); %X轴和第一个Y轴的颜色
set(AX(2),'XColor','k','YColor','k'); %X轴和第二个Y轴的颜色
set(AX(1),'ylim',[-4,4],'ytick',[-4,-2,0,2,4]); %定第一个Y轴的范围
set(AX(2),'ylim',[-40,40],'ytick',[-40,-20,0,20,40]); %定第二个Y轴的范围
set(H1(1),'LineStyle','-','linewidth',1.4,'color','r'); %y1的线型设置
set(H1(2),'LineStyle','-','linewidth',2.5);%y3的线型设置
set(H2,'LineStyle','-.','linewidth',2.3,'color','k'); %y2的线型设置
legend([H1(1),H1(2),H2],{'polluted data';'orginal data';'rudder angle'});
xlabel('time(s)');

%% Step2 Train and Tune the parameters for NARX GP model
%construct the training points
tic
xm = [u_n,v_n,r_n]';
u= d_n';
nm=310;
xnigp = [u(1:nm-1);xm(:,1:nm-1)];
fmh = xm(:,2:nm);
[model, nigp] = trainNIGP(xnigp',fmh',-500,1);
hyp = NIGPModelToHyperparameters(model);
sonig = createSONIG(hyp); 
sonig.addIIPDistance = 0.2; % This is the distance (normalized with respect to the length scales) above which new inducing input points are added.
t1 = toc;  %time for training
%% Add new points
su = hyp.sx(1);
sx= [hyp.sx(2);hyp.sx(3);hyp.sx(4)];
% We now start to implement measurements.
tic;
lastDisplay = toc;
disp('Starting to implement measurements.');
xpo = xm(:,1:nm); % This will contain the posterior mean of x.
upo = u(:,1:nm); % This will contain the posterior mean of u.
xstd = repmat(sx,[1,size(xpo,2)]); % This will contain the poster standard deviation of x.
ustd = repmat(su,[1,size(upo,2)]); % This will contain the poster standard deviation of u.
jointMean = [0*su;0*sx;xm(:,1)]; % This will contain the mean vector of the SONIG input which we're currently applying.
jointCov = [su.^2,zeros(1,6);zeros(6,1),blkdiag(diag(sx.^2),diag(sx.^2))]; % This will contain the covariance matrix of the SONIG input which we're currently applying.
for i = 2:nm
	% We display regularly timed updates.
	if toc > lastDisplay + 5
		disp(['Time passed is ',num2str(toc),' seconds. We are currently at measurement ',num2str(i),' of ',num2str(nm),', with ',num2str(sonig.nu),' IIPs.']);
		lastDisplay = toc;
	end
	% We set up the input and the output distributions, taking into account all the covariances of the parameters.
	jointMean = [u(i-1);jointMean(5:7);xm(:,i)]; % We shift the mean matrix one further.
	jointCov = [su.^2,zeros(1,6);zeros(6,1),[jointCov(5:7,5:7),zeros(3,3);zeros(3,3),diag(sx.^2)]]; % We shift the covariance matrix one further.
	jointDist = createDistribution(jointMean, jointCov);
	inputDist = getSubDistribution(jointDist, 1:4);
	outputDist = getSubDistribution(jointDist, 5:7);
	% We implement the measurement into the SONIG algorithm.
	[sonig, inputPost, outputPost, jointPost] = implementMeasurement(sonig, inputDist, outputDist, jointDist);
	% We update the distributions of all our points.
	jointMean = jointPost.mean;
	jointCov = jointPost.cov;
	upo(:,i-1) = jointPost.mean(1);
	xpo(:,i-1:i) = reshape(jointPost.mean(2:7),[3,2]);
	stds = sqrt(diag(jointCov)');
	ustd(:,i-1) = stds(1);
	xstd(:,i-1:i) = reshape(stds(2:7),[3,2]);
end
disp(['Finished implementing ',num2str(sonig.nm),' measurements in ',num2str(toc),' seconds, using ',num2str(sonig.nu),' IIPs.']);

%% Step3 : Apply SONIG to predict
tic
dt=h*sample;
total_time = 160;
m= ceil(total_time /dt);    %节拍 
x = zeros(3,1); %临时状态变量
TEMP_a = zeros(4,1); 
T = zeros(m,1);  %时间
Y = zeros(m,8); %状态变量
Uci = zeros(m,1);Vci = zeros(m,1);Rci = zeros(m,1);
%Inatialize the ship 
u0 = 1.175; v0 = 0;  r0=0;
x0 = 0; y0 = 0; psi0 = 0;
d0 = 0;flag_lr0= 1;
Initial_input = [u0;v0;r0];
Initial_ob =  [x0; y0; psi0;d0];
Spost = [0;0;0];
U = Initial_input;
x = Initial_ob;
sy = flag_lr0;

for i=1:1:m
    t = dt*i;
    T(i,1)=t;
    time = t;
    % NARX-SONIG model
    [TEMP_a,TEMP_sm,Temp_U,uci,vci,rci]=Infante_SONIG_3dof_HVSA(time,U,Spost,x,sy,sonig,i)  ;
    Spost = [uci;vci;rci];
    %Euler
    x= x + dt.*TEMP_a;
    sy =TEMP_sm;
    U = Temp_U;
    
    %保存数据
    Y(i,1) = Temp_U(1);%u
    Y(i,2) = Temp_U(2);%v
    Y(i,3) = Temp_U(3);%r
    Y(i,5) = x(1);%x
    Y(i,6) = x(2);%y
    Y(i,7) = x(3);%psi
    Y(i,9) = x(4);%舵角   
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

t3 = toc; % time for prediction

% container_z20_20_NIGP = [T,U_pre,V_pre,R_pre];
% save container_z20_20_NIGP container_z20_20_NIGP ;
% container_t21_NIGP = [T,U_pre,V_pre,R_pre];
% save container_t21_NIGP container_t21_NIGP ;

figure(3)
kk= 2;
subplot(411),
u_upper = U_pre + kk*Uci;  % upper boundry
u_lower = U_pre - kk*Uci;  % lower boundrr
patch([T', fliplr(T')], [u_lower', fliplr(u_upper')], 1, 'FaceColor', [0.85,0.85,1], 'EdgeColor', 'none');hold on;
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(T,U_pre,'LineWidth',1.5),xlabel('time (s)'),title('speed U (m/s)');grid on;hold on;

subplot(412),
v_upper = V_pre + kk*Vci;  % upper boundry
v_lower = V_pre - kk*Vci;  % lower boundrr
patch([T', fliplr(T')], [v_lower', fliplr(v_upper')], 1, 'FaceColor', [0.85,0.85,1], 'EdgeColor', 'none');hold on
% set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(T,V_pre,'linewidth',1.5),xlabel('time (s)'),title('speed V (m/s)'),grid on;hold on

subplot(413),
r_upper = R_pre + 2*Rci;  % upper boundry
r_lower = R_pre - 2*Rci;  % lower boundrr
patch([T', fliplr(T')], [r_lower', fliplr(r_upper')], 1, 'FaceColor', [0.85,0.85,1], 'EdgeColor', 'none');hold on
set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
plot(T,R_pre,'linewidth',1.5),xlabel('time (s)'),title('speed R (rad/s)'),grid on;hold on

subplot(414),plot(T,psi,'linewidth',1.6);hold on
% plot(T,duo);hold off
xlabel('time (s)'),title('yaw angle \psi (deg)'),grid on;hold on
% legend('\psi','\delta_c');

% %turning circle
% figure(4)
% subplot(311),
% u_upper = U_pre + 2*Uci;  % upper boundry
% u_lower = U_pre - 2*Uci;  % lower boundrr
% patch([T', fliplr(T')], [u_lower', fliplr(u_upper')], 1, 'FaceColor', [0.85,0.85,1], 'EdgeColor', 'none');hold on;
% set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
% plot(T,U_pre,'LineWidth',1.5),xlabel('time (s)'),title('speed U (m/s)');grid on;hold on;
% 
% subplot(312),
% v_upper = V_pre + 2*Vci;  % upper boundry
% v_lower = V_pre - 2*Vci;  % lower boundrr
% patch([T', fliplr(T')], [v_lower', fliplr(v_upper')], 1, 'FaceColor', [0.85,0.85,1], 'EdgeColor', 'none');hold on
% % set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
% plot(T,V_pre,'linewidth',1.5),xlabel('time (s)'),title('speed V (m/s)'),grid on;hold on
% 
% subplot(313),
% r_upper = R_pre + 2*Rci;  % upper boundry
% r_lower = R_pre - 2*Rci;  % lower boundrr
% patch([T', fliplr(T')], [r_lower', fliplr(r_upper')], 1, 'FaceColor', [0.9,0.9,1], 'EdgeColor', 'none');hold on
% set(gca, 'layer', 'top'); % We make sure that the grid lines and axes are above the grey area.
% plot(T,R_pre,'linewidth',1.5),xlabel('time (s)'),title('speed R (rad/s)'),grid on;hold on
% 
% figure(5)
% plot(Xp,Yp,'linewidth',1.6),xlabel('x position(m)'),ylabel('x position(m)'),title('Motion trajectory)'),grid on;hold on