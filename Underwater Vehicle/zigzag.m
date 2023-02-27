function[signal,x]=zigzag(rudder,d_psi,time)
rudder=rudder*pi/180;   %Converts the input angle into radians
d_psi=d_psi*pi/180;
dt=0.05;   %Time for small changes
x0 = [1;0;0;0;0;0;0;0;10;0;0;0];
ui = [rudder;0;0;0;0;1200]; %Set the initial state of the navigator
x(:,1) = x0(1:12);        %Assigning initial values to the robot's state
signal(1)=ui(1,:);
for k=1:time/dt    %Expand the number of simulations
    [xdot(:,k),U(:,k)] = npsauv(x(:,k),ui);
    st = x(:,k);  dif=xdot(:,k);           
    st_next  = st+ (dt*dif);                 
    x(:,k+1) = st_next;      %Update next status
    if x(12,end)>=d_psi   
        ui = [rudder;0;0;0;0;1200];
    elseif -x(12,end)>=d_psi
        ui = [-rudder;0;0;0;0;1200];
    end
    signal(k+1)=ui(1,:);
end
  psi=x(12,2:end)*180/pi;    %Output bow angle
  signal=-signal*180/pi;
  x=x;

end
