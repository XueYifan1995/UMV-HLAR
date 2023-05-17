%% Loading training data
load npsAUV_zigzag_1010_005 npsAUV_zigzag_1010_005
load npsAUV_zigzag_2020_005 npsAUV_zigzag_2020_005
load npsAUV_zigzag_2505_005 npsAUV_zigzag_2505_005
load npsAUV_zigzag_2510_005 npsAUV_zigzag_2510_005

xinput = [npsAUV_zigzag_1010_005,npsAUV_zigzag_2020_005];
x_1_self = xinput(1:3,:);
u = xinput(4,:)*pi/180;

%% Calculate sindy initial parameters
x_1_self = x_1_self';
dt = 0.05;
dx = zeros(length(x_1_self)-5,3);
for i=3:length(x_1_self)-3
        for k=1:size(x_1_self,2)
            dx(i-2,k) = (1/(12*dt))*(-x_1_self(i+2,k)+8*x_1_self(i+1,k)-8*x_1_self(i-1,k)+x_1_self(i-2,k));   
        end
    end

u=u';
xaug = [x_1_self(3:end-3,:) u(3:end-3,:)];
dx(:,size(x_1_self,2)+1) = 0*dx(:,size(x_1_self,2));

n = size(dx,2);


polyorder=3;   %For the actual system, prediction to the third order term is sufficient
usesine = 0;   
LibraryType = 4;

%% Randomly divided data - for long periods
xaug_train = xaug;
xaug_val = npsAUV_zigzag_2510_005';
dx_train = dx;

%% An optimization algorithm is used to determine Xi
lamda1 = optimizableVariable('lamda1',[0.25 0.5]);
lamda2 = optimizableVariable('lamda2',[1e-6 0.2]);
lamda3 = optimizableVariable('lamda3',[1e-6 0.2]);
vars = [lamda1,lamda2,lamda3];
fun = @(lamda) Sindy3(xaug_train,xaug_val,dx_train,lamda,LibraryType);
results = bayesopt(fun,vars,'MaxObjectiveEvaluations',50)

