function [t,x,y,Vx,Vy,theta] = func_in_air(m,r,alpha,dt,t0,x0,y0,Vx0,Vy0,theta0)
i=1;
%Define initial conditions
g=9.81; %Define acceleration due to gravity (in m/s)
t(i)=t0;
x(i)=x0;
y(i)=y0;
Vx(i)=Vx0;
Vy(i)=Vy0;
theta(i)=theta0;
k=2.16; %Define k value (in kg/m^3)

while y(i)>=0
    V=sqrt((Vx(i))^2+(Vy(i))^2); %calculate velocity
    t(i+1)=(i)*dt; %create time vector
    Vx(i+1)=-k*r^2*abs(sin(alpha-theta(i)))*V^2*cos(theta(i))*(dt/m)+Vx(i); %create Vx vector
    Vy(i+1)=((-k*r^2*abs(sin(alpha-theta(i)))*V^2*sin(theta(i))-m*g))*(dt/m)+Vy(i); %create Vy vector
    x(i+1)=(Vx(i)*dt)+x(i); %create x vector
    y(i+1)=(Vy(i)*dt)+y(i); %create y vector
    theta(i+1)=atan2(Vy(i),Vx(i)); %calculate angle
    i=i+1; %iterate
end
end