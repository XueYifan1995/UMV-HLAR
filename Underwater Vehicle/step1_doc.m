set(0,'defaultfigurecolor','w')
%% Step1:load data 

%Use zigzag for the training set
load npsAUV_zigzag_1010_005p npsAUV_zigzag_1010_005p 
load npsAUV_zigzag_2020_005p npsAUV_zigzag_2020_005p
load npsAUV_zigzag_2505_005 npsAUV_zigzag_2505_005
load npsAUV_zigzag_2510_005 npsAUV_zigzag_2510_005

h=  0.05; sample= 2;  %sampling step
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
u_x_doc = u(1:end-1);  u_y = u(2:end);
v_x_doc = v(1:end-1);  v_y = v(2:end);
r_x_doc = r(1:end-1);  r_y = r(2:end);
d_x_doc = d(1:end-1);
Xm = [u_x_doc,v_x_doc,r_x_doc,d_x_doc];
Ym = [u_y,v_y,r_y];
% noise test
Xm2=[u_x_doc,v_x_doc,r_x_doc];
Am = (Ym-Xm2)/h;
Am_d = wdenoise(Am);  
t_t =t(1:sample:end-1);
Xm_t= Xm(1:sample:num_tr,:);
Am_t = Am(1:sample:num_tr,:);
Am_d_t = Am_d(1:sample:num_tr,:);



Len1= Am_d_t(1:end-1,1);
Len2= Am_d_t(1:end-1,2);
Len3=Am_d_t(1:end-1,3);

u_x_doc=u_x_doc(1:2:end);
v_x_doc=v_x_doc(1:2:end);
r_x_doc=r_x_doc(1:2:end);
d_x_doc=d_x_doc(1:2:end);

for i=1:num_tr/2-1
 [X1_doc(i,:),X2_doc(i,:),X3_doc(i,:)]=ConstructAUV_doc(u_x_doc,v_x_doc,r_x_doc,d_x_doc,i);  %和F8中pooldata 作用是一样的
end


theta1 = pinv(X1_doc'*X1_doc)*X1_doc'*Len1;   %Coefficients of the model
theta2 = pinv(X2_doc'*X2_doc)*X2_doc'*Len2;
theta3 = pinv(X3_doc'*X3_doc)*X3_doc'*Len3;   

t1 = toc;  %Calculate training time
