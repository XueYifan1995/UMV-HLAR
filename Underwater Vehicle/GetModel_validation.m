%% Forecasting with Eulerian dispersion - SINDY
xpre = npsAUV_zigzag_2505_005;
x_p=[1,0,0];
xp=x_p;
usesine=0;
polyorder=3;
n=4;
Nvar = 3;
tspan=[0];
u_1 =xpre(4,:)*pi/180;
for k=1:6000   
    t=dt*k;
    tspan = [tspan,t];
    y=[xp(k,:) u_1(k+1)];
    xPool = selfpooldata(y,LibraryType);
    dxPool = xPool*Xi_sphs(:,1:Nvar);
    st = xp(k,:);  dif=dxPool;            
    st_next  = st+ (dt*dif);                 
    xp(k+1,:) = st_next;      
end
xp=xp';%Matrixing the predicted results into a homomorphic matrix

rmse_u_sindy = sqrt(mean((xpre(1,:)-xp(1,:)).^2))  %three dimensions need to be compared with R^2
rmse_v_sindy = sqrt(mean((xpre(2,:)-xp(2,:)).^2))
rmse_r_sindy = sqrt(mean((xpre(3,:)-xp(3,:)).^2))