function [dydt,sm,U,sPost_v,sPost_r] = Infante_SONIG_2dof_nouncertain_HVSA(time,U,x,sy,sonig)
%状态变量
tic
dydt = zeros(2,1);
v = U(1); r = U(2);
psi = x(1); d=x(2);
if d<-30
    d=-30;
elseif d>30
    d=30;
end
    
Xs = [d,v,r]' ; %trial poionts
% flag_lr=sm;
% %NARX_MT
[mPost,sPost,~] = makeSonigPrediction(sonig,Xs);
%pass the values
v= mPost(1);
r= mPost(2);
U = [v,r];
sPost_v = sPost(1);
sPost_r = sPost(2);
%坐标系矩阵转化
dydt(1) = r;                     %dpsi    
%scheduled maneuver
% dydt(4) = AdjustAngle_t30_r(d,psi,time);
% % dydt(4) = AdjustAngle_t30_r(d,psi,time);
% sm= 1; %psi状态

% % %zigzag
[dydt(2),sm] = AdjustAngle_z30_5(d,psi,sy);
t4=toc;
end

