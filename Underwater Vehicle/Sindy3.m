function [fx,Xi] = Sindy3(xaug_train,xaug_val,dx_train,lamda,LibraryType)
temp=table2array(lamda);
lamda1=temp(1);
lamda2=temp(2);
lamda3=temp(3);


n=4;
Theta = selfpooldata(xaug_train,LibraryType);
Theta_norm = zeros(size(Theta,2),1); 
for i = 1:size(Theta,2)
   Theta_norm(i) = norm(Theta(:,i));
   Theta(:,i) = Theta(:,i)./Theta_norm(i);
end
m = size(Theta,2);

lambda_vec = [lamda1,lamda2,lamda3];


if exist('lambda_vec') == 1
    Xi = sparsifyDynamicsIndependent(Theta,dx_train,lambda_vec,n-1);
else
    Xi = sparsifyDynamics(Theta,dx_train,lambda,n-1);
end


for i = 1:size(Theta,2)
   Xi(i,:) = Xi(i,:)./Theta_norm(i);
end


Nvar = 3;
dxPool = [];

x_p=[1,0,0];
xp=x_p;
dt = 0.05;
u_1 =xaug_val(:,4)*pi/180;
for k=1:length(xaug_val)   %Calculate the predicted acceleration magnitude
    y=[xp(k,:) u_1(k)];
    xPool = selfpooldata(y,LibraryType);
    dxPool = xPool*Xi(:,1:Nvar);
    st = xp(k,:);  dif=dxPool;            
    st_next  = st+ (dt*dif);                
    xp(k+1,:) = st_next;  
end


Rsquare_u = corrcoef(xp(2:end,1),xaug_val(:,1));
Rsquare_v = corrcoef(xp(2:end,2),xaug_val(:,2));
Rsquare_r = corrcoef(xp(2:end,3),xaug_val(:,3));

fx = -(Rsquare_u(1,2)+Rsquare_v(1,2)+Rsquare_r(1,2));   
%The solver finds the minimum value and converts the maximum value into a minimum value problem.

end