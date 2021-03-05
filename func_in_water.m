function [t,x,y,Vx,Vy,theta] = func_in_water(m,r,alpha,dt,t0,x0,y0,Vx0,Vy0,theta0)
i=1;
%Define initial conditions
g=9.81; %Define acceleration due to gravity (in m/s^2)
t(i)=t0;
x(i)=x0;
y(i)=y0;
Vx(i)=Vx0;
Vy(i)=Vy0;
theta(i)=theta0;
rho_water=1000; %Define water density (in kg/m^3)

while y(i)<=0 %iterate while the stone is below the water
    V=sqrt((Vx(i))^2+(Vy(i))^2); %calculate velocity
    S=abs(y(i))/sin(alpha); %calculate S value
    Swet=r^2*(acos(1-(S/r))-(1-(S/r))*sqrt(1-((1-(S/r))^2))); %calculate Swet value
    t(i+1)=(i)*dt; %create time vector
    Vx(i+1)=(-0.5*rho_water*V^2*Swet*sin(alpha-theta(i))*sin(alpha))*(dt/m)+Vx(i); %create Vx vector
    Vy(i+1)=(0.5*rho_water*V^2*Swet*sin(alpha-theta(i))*cos(alpha)-m*g)*(dt/m)+Vy(i); %create Vy vector
    x(i+1)=(Vx(i)*dt)+x(i); %create x vector
    y(i+1)=(Vy(i)*dt)+y(i); %create y vector
    
    if (imag(Vx(i))~=0) && (imag(Vy(i))~=0) %if imaginary velocity, set to zero and delete other entries
        t=[];
        x=[];
        y=[];
        Vx=[];
        Vy=0;
        theta=[];
        return
    end
    
    theta(i+1)=atan2(Vy(i),Vx(i)); %calculate angle
    
    i=i+1; %iterate
    
end

end