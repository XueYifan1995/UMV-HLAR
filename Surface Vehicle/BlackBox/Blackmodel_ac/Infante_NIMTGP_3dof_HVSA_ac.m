function [dydt,sm,U_a,sPost_u,sPost_v,sPost_r] = Infante_NIMTGP_3dof_HVSA_ac(time,U,x,sy,model,delta_exin)
%状态变量
dydt = zeros(4,1);
u = U(1); v = U(2); r = U(3);
x_p = x(1); y_p = x(2); 
psi = x(3); d=x(4);
% Xs = [u,v,r,d] ; %trial poionts
Xs = [u,v,r,delta_exin] ; %trial poionts
% flag_lr=sm;
% %NARX_MT
[mPost,sPost] = predNIGP(model,Xs,0);
% [mPost,sPost] = predNIGP(model,Xs,2);
%pass the values
u_a= mPost(1); 
v_a= mPost(2);
r_a= mPost(3);
U_a = [u_a;v_a;r_a];
sPost_u = sPost(1);
sPost_v = sPost(2);
sPost_r = sPost(3);
%坐标系矩阵转化
dydt(1) = u*cos(psi)-v*sin(psi); %dx   
dydt(2) = u*sin(psi)+v*cos(psi); %dy
dydt(3) = r;                     %dpsi    
% % % % %zigzag
% [dydt(4),sm ] = AdjustAngle_z15_5(d,psi,sy);
dydt(4)=0;sm=sy;
end

