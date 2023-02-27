load HSVACPMCKVLCC2Z1505 HSVACPMCKVLCC2Z1505
load HSVACPMCKVLCC2Z3505 HSVACPMCKVLCC2Z3505
load HSVACPMCKVLCC2Z1010P HSVACPMCKVLCC2Z1010P 

st=1; ov=3200; sample=12;
data = HSVACPMCKVLCC2Z1505(1:sample:ov,:);
set(0,'defaultfigurecolor','w')


t = data(1:end,1);
x = data(st:end,2);
y = data(st:end,3);
psi = data(st:end,4);
u = data(st:end,5);
v = data(st:end,6);
r = data(st:end,7)*pi/180;
phi = data(st:end,8);
d =  data(st:end,9);

load HVSA_z15_05_MTGP HVSA_z15_05_MTGP
load HVSA_z35_05_MTGP HVSA_z35_05_MTGP
load HVSA_z10_10_MTGP HVSA_z10_10_MTGP

load HVSA_z15_05_SMTGP HVSA_z15_05_SMTGP
load HVSA_z35_05_SMTGP HVSA_z35_05_SMTGP

T =HVSA_z15_05_SMTGP(:,1);
u2 = HVSA_z15_05_SMTGP(:,2);
v2 = HVSA_z15_05_SMTGP(:,3);
r2 = HVSA_z15_05_SMTGP(:,4)*pi/180;


RMSE_u_s_N = sqrt(mean((u - u2).^2));

RMSE_v_s_N = sqrt(mean((v - v2).^2));

RMSE_r_s_N = sqrt(mean((r - r2).^2));

plot(T,[r,r2]);



% [R2_u_N,~] = rsquare(u1,u2);
% [R2_u_R,~] = rsquare(u1,u3);
% [R2_u_S,~] = rsquare(u1,u4);
% 
% [R2_v_N,~] = rsquare(v1,v2);
% [R2_v_R,~] = rsquare(v1,v3);
% [R2_v_S,~] = rsquare(v1,v4);
% 
% [R2_r_N,~] = rsquare(r1,r2);
% [R2_r_R,~] = rsquare(r1,r3);
% [R2_r_S,~] = rsquare(r1,r4);