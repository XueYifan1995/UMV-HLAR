function [dydt,sm,U,sPost_r] = Infante_HVSA_1dof_unCertain(time,U,Spost,x,sy,sonig,i)
%状态变量
tic

dydt = zeros(2,1);
r = U(1);
psi = x(1); d=x(2);
Xs = r ; %trial poionts

%to distribution
    su = sonig.hyp.sx(1);
    sx= [sonig.hyp.sx(2)];
    xDist = createDistribution(Xs, diag(sx.^2));
    uDist = createDistribution(d, diag(su.^2));
	inputDist = joinDistributions(uDist, xDist);

%% Prediction
newxDist = makeSonigStochasticPrediction(sonig,inputDist);
mPost = newxDist.mean;
sPost = sqrt(diag(newxDist.cov));
    
%pass the values
r= mPost;
U = r;
sPost_r = sPost;

%% controller
%坐标系矩阵转化
dydt(1) = r;                     %dpsi    
%scheduled maneuver
% dydt(4) = AdjustAngle_t30_r(d,psi,time);
% % dydt(4) = AdjustAngle_t30_r(d,psi,time);
% sm= 1; %psi状态

% % %zigzag
[dydt(2),sm] = AdjustAngle_z30_5(d,psi,sy);

t_pre=toc;
end

