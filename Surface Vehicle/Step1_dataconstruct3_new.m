set(0,'defaultfigurecolor','w')
%% Step1:load data 

load HSVACPMCKVLCC2Z1005 HSVACPMCKVLCC2Z1005    
load HSVACPMCKVLCC2Z1505 HSVACPMCKVLCC2Z1505
load HSVACPMCKVLCC2Z2005 HSVACPMCKVLCC2Z2005
load HSVACPMCKVLCC2Z3005 HSVACPMCKVLCC2Z3005
load HSVACPMCKVLCC2Z3505 HSVACPMCKVLCC2Z3505
load HSVACPMCKVLCC2Z1010P HSVACPMCKVLCC2Z1010P 
load HSVACPMCKVLCC2Z2010P HSVACPMCKVLCC2Z2010P 
load HSVACPMCKVLCC2Z1001 HSVACPMCKVLCC2Z1001

h=  0.05; sample= 2;  %Sampling step
%h is the set time interval

data_raw = [HSVACPMCKVLCC2Z1005(1:3100,:);HSVACPMCKVLCC2Z2005(1:3100,:);HSVACPMCKVLCC2Z3005(1:3100,:)];

data = data_raw;
num_tr = size(data,1);
t=linspace(0,num_tr,num_tr+1).*h;
data(:,1)=t(1:end-1);


data_train = data;
t =data_train(:,1);
psi = data_train(:,4)*pi/180;    %Bow angle
u = data_train(:,5);
v = data_train(:,6);
r = data_train(:,7)*pi/180;
phi = data_train(:,8)*pi/180;   %Transverse roll angle
d =  data_train(:,9)*pi/180;     %Rudder angle

%Data Mapping
figure(1)
subplot(414),plot(t,[psi*180/pi,-d*180/pi],'linewidth',2);hold on
xlabel('time (s)'),ylabel('\psi & \delta (deg)'),grid on;hold on
legend('\psi','\delta')
subplot(411),plot(t,u,'linewidth',2),xlabel('time (s)'),ylabel('u (m/s)');grid on;hold on
subplot(412),plot(t,v,'linewidth',2),xlabel('time (s)'),ylabel('v (m/s)');grid on;hold on
subplot(413),plot(t,r*180/pi,'linewidth',2),xlabel('time (s)'),ylabel('r (deg/s)');grid on;hold on
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
t_t =t(1:sample:end-1);
Xm_t= Xm(1:sample:num_tr,:);
Am_t = Am(1:sample:num_tr,:);
Am_d_t = Am_d(1:sample:num_tr,:);

Len1= Am_d_t(1:end-1,1);
Len2= Am_d_t(1:end-1,2);
Len3=Am_d_t(1:end-1,3);

u_x=u_x(1:2:end);
v_x=v_x(1:2:end);
r_x=r_x(1:2:end);
d_x=d_x(1:2:end);

figure(2)
sam=10;
scatter(t_t(1:sam:end-1),Am_t(1:sam:end,2));axis([0 320 -0.05 0.05]);hold on;
scatter(t_t(1:sam:end-1),Am_d_t(1:sam:end,2));axis([0 320 -0.05 0.05]);grid on;hold on;

for i=1:num_tr/2-1
 ua(i)=u_x(i)-1.179;
 [X1(i,:),X2(i,:),X3(i,:)]=Construct3(ua,v_x,r_x,d_x,i);  
end

theta1 = (X1'*X1)\X1'*Len1;   %Coefficients of the model
theta2 = (X2'*X2)\X2'*Len2;
theta3 = (X3'*X3)\X3'*Len3;   

save theta1 theta1   %Variables->Data files
save theta2 theta2
save theta3 theta3