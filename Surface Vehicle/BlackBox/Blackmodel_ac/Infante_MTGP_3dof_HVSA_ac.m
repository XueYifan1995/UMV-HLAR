function [dydt,sm,sPost,uci,vci,rci,U_a] = Infante_MTGP_3dof_HVSA_ac(dt,U,x,spost,sy,dynmodel,rudder,uci,vci,rci) 
%状态变量
tic
dydt = zeros(4,1);
u = U(1); v = U(2); r = U(3);
x_p = x(1); y_p = x(2); 
psi = x(3); d=x(4);
Xs = [u;v;r] ; %trial poionts
% flag_lr=sm;
% %NARX_MT
% [mPost,sPost] = predNIGP(model,Xs,0);
% [mPost,sPost] = predNIGP(model,Xs,2);
%% PILCO
[tileM, tileS] = jointDistribution(Xs,spost,rudder);
[mPost,sPost,~] = gp0(dynmodel, tileM, tileS);
uci = uci+(sum(sPost(1,:))+sum(sPost(:,1))-sPost(1,1))*dt^2;
vci = vci+(sum(sPost(2,:))+sum(sPost(:,2))-sPost(2,2))*dt^2;
rci = rci+(sum(sPost(3,:))+sum(sPost(:,3))-sPost(3,3))*dt^2;
%% pass the values
u_a= mPost(1); 
v_a= mPost(2);
r_a= mPost(3);
U_a = [u_a;v_a;r_a];
%坐标系矩阵转化
dydt(1) = u*cos(psi)-v*sin(psi); %dx   
dydt(2) = u*sin(psi)+v*cos(psi); %dy
dydt(3) = r;                     %dpsi    
% %scheduled maneuver
% dydt(4) = AdjustAngle_t30_r(d,psi,time);
% % dydt(4) = AdjustAngle_t30_r(d,psi,time);
% sm= 1; %psi状态

% % % % %zigzag
[dydt(4),sm] = AdjustAngle_z35_5(d,psi,sy);
t4=toc;
end

