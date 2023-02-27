function [dydt,sm,U,sPost_u,sPost_v,sPost_r] = Infante_SONIG_3dof_unCertain(time,U,Variance,x,sy,sonig,i)
%状态变量
tic

dydt = zeros(4,1);
u = U(1); v = U(2); r = U(3);
x_p = x(1); y_p = x(2); 
psi = x(3); d=x(4);
Xs = [u,v,r]' ; %trial poionts

%to distribution
su = sonig.hyp.sx(1);
sx= [sonig.hyp.sx(2);sonig.hyp.sx(3);sonig.hyp.sx(3)];
xDist = createDistribution(Xs, diag(sx.^2));
uDist = createDistribution(d, diag(su.^2));
inputDist = joinDistributions(uDist, xDist);

% flag_lr=sm;
%% Prediction
newxDist = makeSonigStochasticPrediction(sonig, inputDist);
mPost = newxDist.mean;
sPost = sqrt(diag(newxDist.cov));
    
%pass the values
u= mPost(1);
v= mPost(2);
r= mPost(3);
U = [u,v,r];
sPost_u = sPost(1);
sPost_v = sPost(2);
sPost_r = sPost(3);

%% controller
%坐标系矩阵转化
dydt(1) = u*cos(psi)-v*sin(psi); %dx   
dydt(2) = u*sin(psi)+v*cos(psi); %dy
dydt(3) = r;                     %dpsi             
%scheduled maneuver
% dydt(4) = AdjustAngle_t30_r(d,psi,time);
% % dydt(4) = AdjustAngle_t30_r(d,psi,time);
% sm= 1; %psi状态

% % %zigzag
[dydt(4),sm] = AdjustAngle_z30_5(d,psi,sy);

t_pre=toc;
end

