set(0,'defaultfigurecolor','w')
%% Step1:load data 

%Use zigzag for the training set
load npsAUV_zigzag_1005 npsAUV_zigzag_1005  
load npsAUV_zigzag_2005 npsAUV_zigzag_2005
load npsAUV_zigzag_3005 npsAUV_zigzag_3005
load npsAUV_zigzag_1010_005p npsAUV_zigzag_1010_005p
load npsAUV_zigzag_2020_005p npsAUV_zigzag_2020_005p

h=  0.05; sample= 2;  %Sampling step
%h is the set time interval

data_raw = [npsAUV_zigzag_1010_005p,npsAUV_zigzag_2020_005p]';
data = data_raw;

num_tr = size(data,1);  
t=linspace(0,num_tr,num_tr+1).*h;

data_train = data;
u = data_train(:,1);
v = data_train(:,2);
r = data_train(:,3);
d =  data_train(:,4)*pi/180;     %Rudder angle

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

for i=1:num_tr/2-1
 [X1(i,:),X2(i,:),X3(i,:)]=ConstructAUV_submarine(u_x,v_x,r_x,d_x,i);  
end


theta1 = pinv(X1'*X1)*X1'*Len1;   %Coefficients of the model
theta2 = pinv(X2'*X2)*X2'*Len2;
theta3 = pinv(X3'*X3)*X3'*Len3;   

t1 = toc;  %Record training time
