%****************************************************************
%   Mini_Project_Three.m
%
%   PROGRAM DESCRIPTION
%   This program calculates and plots the displacement of a
%   stone that is skipped across a lake using Euler's Method.
%
%   INPUT: None
%   OUTPUT: Horizontal distances and plot for each deliverable
%
%   WRITTEN BY: Evan Schimmel
%               02/19/2021
%
%****************************************************************

clc
clear variables
close all

file_number=fopen('MiniProjectThree.txt','w');

rho_stone=2600; %define density of stone (in kg/m^3)
r=0.035; %define radius of stone (in m)
h=0.0125; %define height of stone (in m)
m=rho_stone*pi*r^2*h; %calculate mass of stone (in kg)

%% Deliverable 1: Air only simulation

%Defining deliverable specific initial conditions
alpha1=deg2rad(22);
dt1=1e-5;
x01=0;
y01=0.4;
Vo1=22;
theta01=deg2rad(15);
Vx01=Vo1*cos(theta01);
Vy01=Vo1*sin(theta01);
t01=0;

%Pulling vectors using air function
[t1,x1,y1,Vx1,Vy1,theta1] = func_in_air(m,r,alpha1,dt1,t01,x01,y01,Vx01,Vy01,theta01);

fprintf(file_number,'For deliverable 1, the total horizontal distance is %6.3f m\n',x1(length(x1)));

%create water surface vector
xsurf=[0 x1(end)];
ysurf=[0;0];

figure(1)
%subplot(5,1,1);
plot(x1,y1,'k-',xsurf,ysurf,'b--');
legend('stone trajectory','water surface');
title('In-Air Trajectory of Stone');
grid on
xlabel('horizontal distance [m]');
ylabel('vertical distance [m]');
axis([0 25 -0.2 2.5]);

%% Deliverable 2: Water only simulation

%Defining deliverable specific initial conditions
alpha2=deg2rad(22);
dt2=1e-5;
x02=0;
y02=0;
Vo2=22;
theta02=deg2rad(-8);
Vx02=Vo2*cos(theta02);
Vy02=Vo2*sin(theta02);
t02=0;

%Pulling vectors using water function
[t2,x2,y2,Vx2,Vy2,theta2] = func_in_water(m,r,alpha2,dt2,t02,x02,y02,Vx02,Vy02,theta02);

fprintf(file_number,'For deliverable 2, the total horizontal distance is %6.3f m\n',x2(length(x2)));

%create water surface vector
xsurf=[0 x2(end)];
ysurf=[0;0];

figure(2)
%subplot(5,1,2);
plot(x2,y2,'k-',xsurf,ysurf,'b--');
legend('stone trajectory','water surface');
title('In-Water Trajectory of Stone');
grid on
xlabel('horizontal distance [m]');
ylabel('vertical distance [m]');
axis([0 0.25 -0.08 0.06]);

%% Deliverable 3: Simulating a single skip

%Defining deliverable specific initial conditions
alpha3=deg2rad(22);
dt3=1e-5;

Vo3=22;

t03=0;
x03=0;
y03=0.4;
theta03=deg2rad(15);
Vx03=Vo3*cos(theta03);
Vy03=Vo3*sin(theta03);

%Pulling vectors using air function
[t3_air,x3_air,y3_air,Vx3_air,Vy3_air,theta3_air] = func_in_air(m,r,alpha3,dt3,t03,x03,y03,Vx03,Vy03,theta03);

%Defining deliverable specific initial conditions (outputs of previous)
t03=t3_air(end);
x03=x3_air(end);
y03=y3_air(end);
theta03=theta3_air(end);
Vx03=Vx3_air(end);
Vy03=Vy3_air(end);

%Pulling vectors using water function
[t3_water,x3_water,y3_water,Vx3_water,Vy3_water,theta3_water] = func_in_water(m,r,alpha3,dt3,t03,x03,y03,Vx03,Vy03,theta03);

fprintf(file_number,'For deliverable 3, the total horizontal distance is %6.3f m\n',x3_water(end));

%create water surface vector
xsurf=[0 x3_water(end)];
ysurf=[0;0];

figure(3)
%subplot(5,1,3);
plot(x3_air,y3_air,'k-',xsurf,ysurf,'b--',x3_water,y3_water,'k-');
legend('stone trajectory','water surface');
title('Single Skip Trajectory of Stone');
grid on
xlabel('horizontal distance [m]');
ylabel('vertical distance [m]');
axis([0 27 -0.2 2.5]);

%% Deliverable 4: Multiple skips

%Defining deliverable specific initial conditions
alpha4=deg2rad(22);
dt4=1e-5;

Vo4=22;

t04=0;
x04=0;
y04=0.4;
theta04=deg2rad(15);
Vx04=Vo4*cos(theta04);
Vy04=Vo4*sin(theta04);

i=1;
Vy4_water(i)=1e-15; %set the y velocity to something not zero to enter the while loop
initialization=0; %use the initial conditions above until the while loop is cycled through

while (Vy4_water(i) ~= 0)
    figure(4)
    %subplot(5,1,4);
    hold on
    %legend('stone trajectory','water surface');
    title('Full Skip Trajectory of Stone');
    grid on
    xlabel('horizontal distance [m]');
    ylabel('vertical distance [m]');
    axis([0 47 -0.2 2.5]);
    if initialization==1 %if repeating, use new initial conditions
        t04=t4_water(end);
        x04=x4_water(end);
        y04=y4_water(end);
        theta04=theta4_water(end);
        Vx04=Vx4_water(end);
        Vy04=Vy4_water(end);
    end
    
    %Pulling vectors using air function
    [t4_air,x4_air,y4_air,Vx4_air,Vy4_air,theta4_air] = func_in_air(m,r,alpha4,dt4,t04,x04,y04,Vx04,Vy04,theta04);
    
    %create water surface vector
    xsurf=[0 x4_air(end)];
    ysurf=[0;0];
    
    plot(x4_air,y4_air,'k-',xsurf,ysurf,'b--');
    
    %define new initial conditions for water function (outputs of air)
    t04=t4_air(end);
    x04=x4_air(end);
    y04=y4_air(end);
    Vx04=Vx4_air(end);
    Vy04=Vy4_air(end);
    theta04=theta4_air(end);
    
    %Pulling vectors using water function
    [t4_water,x4_water,y4_water,Vx4_water,Vy4_water,theta4_water] = func_in_water(m,r,alpha4,dt4,t04,x04,y04,Vx04,Vy04,theta04);
    
    plot(x4_water,y4_water,'k-');
    
    %check if the water velocity is zero
    if (Vy4_water == 0)
        break
    end
    
    %iterate variables
    i=i+1;
    initialization=1;
    
end

fprintf(file_number,'For deliverable 4, the stone skips %2.0f times and last enters the water at %5.2f m\n\n',9,x4_air(length(x4_air)));

%% Deliverable 5: Contest

distbest=0;

%Defining deliverable specific initial conditions
dt5=1e-5;

Vo5=22;

t05=0;
x05=0;

i=0;

%Define parameters for testing best
thetamin=0;
thetamax=25;
thetastep=5;
ymin=0;
%ymax=2*r*sin(alpha);
ystep=0.005;
alphamin=0;
alphamax=25;
alphastep=5;

%Iterate for each possible combination
for thetaval=deg2rad(thetamin):deg2rad(thetastep):deg2rad(thetamax)
    for alphaval=deg2rad(alphamin):deg2rad(alphastep):deg2rad(alphamax)
        for yval=ymin:ystep:(2*r*sin(alphaval))
            
            Vx05=Vo5*cos(thetaval);
            Vy05=Vo5*sin(thetaval);
            
            i=i+1;
            
            %Writing input values to array
            combo(i,1)=i;
            combo(i,2)=thetaval;
            combo(i,3)=yval;
            combo(i,4)=alphaval;
            
            %Using code from deliverable 4 to calculate steps/distance
            [skips,x5_air]=fullskip(m,r,alphaval,dt5,t05,x05,yval,Vx05,Vy05,thetaval);
             
            combo(i,5)=x5_air(end);
            combo(i,6)=skips;
            
            %Writing resultant data to array
            if combo(i,5)>distbest && combo(i,6)>0
                thetabest=combo(i,2);
                ybest=combo(i,3);
                alphabest=combo(i,4);
                distbest=combo(i,5);
                numskips=combo(i,6);
            end
            
            %can enable this to print iterative status
            fprintf('Iteration: %8.0f\n',i);
            fprintf('Input Values: Theta:%5.2f deg; Alpha:%5.2f deg; y0:%4.3f m\n',rad2deg(thetaval),rad2deg(alphaval),yval);
            fprintf('Output Values: Horizontal Distance:%5.2f m; No. of Skips:%2.0f\n\n',x5_air(end),skips);
        end
    end
end

fprintf(file_number,'For deliverable 5, Ideal Parameters |--> Theta: %5.2f deg; Alpha: %5.2f deg; y0: %4.3f m\n',rad2deg(thetabest),rad2deg(alphabest),ybest);
fprintf(file_number,'For deliverable 5, Predicted Results |--> Horizontal distance: %5.2f m; No. of skips: %2.0f\n',distbest,numskips);

%Define optimal values as inputs for plot
theta05=thetabest;
y05=ybest;
alpha5=alphabest;

%calculate x and y initial velocities
Vx05=Vo5*cos(theta05);
Vy05=Vo5*sin(theta05);

i=1;
Vy5_water(i)=1e-15;
initialization=0;

%same code as deliverable 4, just with the optimal ICs
while (Vy5_water(i) ~= 0)
    figure(5)
    %subplot(5,1,5);
    hold on
    %legend('stone trajectory','water surface');
    title('Optimized Skip Trajectory of Stone');
    grid on
    xlabel('horizontal distance [m]');
    ylabel('vertical distance [m]');
    axis([0 inf -inf inf]);
    if initialization==1
        t05=t5_water(end);
        x05=x5_water(end);
        y05=y5_water(end);
        theta05=theta5_water(end);
        Vx05=Vx5_water(end);
        Vy05=Vy5_water(end);
    end
    
    [t5_air,x5_air,y5_air,Vx5_air,Vy5_air,theta5_air] = func_in_air(m,r,alpha5,dt5,t05,x05,y05,Vx05,Vy05,theta05);
    
    xsurf=[0 x5_air(end)];
    ysurf=[0;0];
    
    plot(x5_air,y5_air,'k-',xsurf,ysurf,'b--');
    
    t05=t5_air(end);
    x05=x5_air(end);
    y05=y5_air(end);
    Vx05=Vx5_air(end);
    Vy05=Vy5_air(end);
    theta05=theta5_air(end);
    
    [t5_water,x5_water,y5_water,Vx5_water,Vy5_water,theta5_water] = func_in_water(m,r,alpha5,dt5,t05,x05,y05,Vx05,Vy05,theta05);
    
    plot(x5_water,y5_water,'k-');
    
    if (Vy5_water == 0)
        break
    end
    
    i=i+1;
    initialization=1;
    
end

fclose(file_number);