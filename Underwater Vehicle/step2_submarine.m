load npsAUV_zigzag_2505_005 npsAUV_zigzag_2505_005
load npsAUV_zigzag_3505 npsAUV_zigzag_3505

pre_data= npsAUV_zigzag_2505_005;
order=pre_data(4,:)*pi/180;

%deg angle; rad radian
R2D=180/pi;
D2R=pi/180;
h=0.05;      %Step length: seconds
m= size(pre_data,2);    %Rhythm 
tic
x = zeros(6,1); %Temporary status variables
TEMP_a = zeros(6,1); 
T = zeros(m,1);  %Time
Y = zeros(m,14); %Status Variables
%Variable initialization
u0 = 1; v0 = 0;  r0 = 0;  x0 = 0; y0 = 0; psi0 = 0;d0 = 0;
Initial = [u0; v0; r0;  x0; y0; psi0];




x = Initial;
X1 = zeros(m,2);X2 = zeros(m,7);X3 = zeros(m,7);

for i=1:1:m
    t = h*i;
    T(i,1)=t;
    d_vali=order(i);
    [TEMP_a]=Infante_submarine(x,theta1,theta2,theta3,d_vali)  ; 
    x= x + h.*TEMP_a;
    Y(i,1) = x(1);%u
    Y(i,2) = x(2);%v
    Y(i,3) = x(3);%r
    Y(i,4) = x(4);%x
    Y(i,5) = x(5);%y
    Y(i,6) = x(6);%psi-yaw
end

U_pre  = Y(:,1);
V_pre  = Y(:,2);
R_pre  = Y(:,3);
Xp = Y(:,4);
Yp = Y(:,5);
psi = Y(:,6);
duo = Y(:,7);




duo = duo*180/pi;
psi = psi *180/pi;
t2=toc;

figure
subplot(311),plot(T,U_pre,'linewidth',1.5),xlabel('time (s)'),ylabel('u (m/s)');grid on;hold on
subplot(312),plot(T,V_pre,'linewidth',1.5),xlabel('time (s)'),ylabel('v (m/s)');grid on;hold on
subplot(313),plot(T,R_pre,'linewidth',1.5),xlabel('time (s)'),ylabel('r (deg/s)');grid on;hold on

U_pre_submarine = U_pre;
V_pre_submarine = V_pre;
R_pre_submarine = R_pre;

save U_pre_submarine U_pre_submarine
save V_pre_submarine V_pre_submarine
save R_pre_submarine R_pre_submarine


