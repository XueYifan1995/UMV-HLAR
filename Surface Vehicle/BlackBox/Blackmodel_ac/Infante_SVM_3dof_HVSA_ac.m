function [dydt,sm,U_a] = Infante_SVM_3dof_HVSA_ac(U,x,sy,SVM_u,SVM_v,SVM_r)
% tic
%状态变量
dydt = zeros(4,1);
u = U(1); v = U(2); r = U(3);
x_p = x(1); y_p = x(2); 
psi = x(3); d=x(4);
Xs = [u,v,r,d] ; %trial poionts
% flag_lr=sm;
% %sim SVM
u_a= predict(SVM_u,Xs) ; 
v_a= predict(SVM_v,Xs);
r_a= predict(SVM_r,Xs);
U_a = [u_a;v_a;r_a];
%坐标系矩阵转化
dydt(1) = u*cos(psi)-v*sin(psi); %dx   
dydt(2) = u*sin(psi)+v*cos(psi); %dy
dydt(3) = r;                     %dpsi   
% % % %zigzag
[dydt(4),sm] = AdjustAngle_z15_5(d,psi,sy);
% [dydt(4),sm] = AdjustAngle_z30_5(d,psi,sy);
end

