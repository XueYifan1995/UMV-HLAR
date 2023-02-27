load HSVACPMCKVLCC2Z1505 HSVACPMCKVLCC2Z1505
load HSVACPMCKVLCC2Z2505 HSVACPMCKVLCC2Z2505
load HSVACPMCKVLCC2Z3005 HSVACPMCKVLCC2Z3005
load HSVACPMCKVLCC2Z3505 HSVACPMCKVLCC2Z3505
load HSVACPMCKVLCC2Z1010P HSVACPMCKVLCC2Z1010P 
load HSVACPMCKVLCC2Z1001 HSVACPMCKVLCC2Z1001
load HSVACPMCKVLCC2Z1501 HSVACPMCKVLCC2Z1501
pre_data= HSVACPMCKVLCC2Z3505;
order=pre_data(:,9)*pi/180;
tic
%deg angle; rad radian
R2D=180/pi;
D2R=pi/180;
h=0.05;      %Step length: seconds
m= 180/0.05;    %Rhythm

x = zeros(6,1); %Temporary status variables
TEMP_a = zeros(6,1); 
T = zeros(m,1);  %Time
Y = zeros(m,14); %Status Variables
%Variable initialization
u0 = 1.179; v0 = 0;  r0 = 0;  x0 = 0; y0 = 0; psi0 = 0;d0 = 0;
Initial = [u0; v0; r0;  x0; y0; psi0];

x = Initial;
X1 = zeros(m,2);X2 = zeros(m,7);X3 = zeros(m,7);

for i=1:1:m
    t = h*i;
    T(i,1)=t;
    d_vali=order(i);
    [TEMP_a]=Infante_3(x,theta1,theta2,theta3,d_vali)  ; 
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
R_pre  = Y(:,3)*180/pi;
Xp = Y(:,4);
Yp = Y(:,5);
psi = Y(:,6);
duo = Y(:,7);

duo = duo*180/pi;
psi = psi *180/pi;

t2 =toc;
figure(3)
subplot(414),plot(T,psi,'linewidth',1.5);hold on
xlabel('time (s)'),ylabel('psi (deg)'),grid on;hold on
subplot(411),plot(T,U_pre,'linewidth',1.5),xlabel('time (s)'),ylabel('u (m/s)');grid on;hold on
subplot(412),plot(T,V_pre,'linewidth',1.5),xlabel('time (s)'),ylabel('v (m/s)');grid on;hold on
subplot(413),plot(T,R_pre,'linewidth',1.5),xlabel('time (s)'),ylabel('r (deg/s)');grid on;hold on

%For drawing comparison charts
U_pre_25 = U_pre;
V_pre_25 = V_pre;
R_pre_25 = R_pre;

save U_pre_25 U_pre_25
save V_pre_25 V_pre_25
save R_pre_25 R_pre_25

