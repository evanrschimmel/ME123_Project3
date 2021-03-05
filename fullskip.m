function [i,x4_air] = fullskip(m,r,alpha4,dt4,t04,x04,y04,Vx04,Vy04,theta04)

%this is all identical code to the deliverable 4 code, just in a function
i=1;
Vy4_water(i)=1e-15;
initialization=0;

while (Vy4_water(i) ~= 0)
    if initialization==1
        t04=t4_water(end);
        x04=x4_water(end);
        y04=y4_water(end);
        theta04=theta4_water(end);
        Vx04=Vx4_water(end);
        Vy04=Vy4_water(end);
    end
    
    [t4_air,x4_air,y4_air,Vx4_air,Vy4_air,theta4_air] = func_in_air(m,r,alpha4,dt4,t04,x04,y04,Vx04,Vy04,theta04);
    
    t04=t4_air(end);
    x04=x4_air(end);
    y04=y4_air(end);
    Vx04=Vx4_air(end);
    Vy04=Vy4_air(end);
    theta04=theta4_air(end);
    
    [t4_water,x4_water,y4_water,Vx4_water,Vy4_water,theta4_water] = func_in_water(m,r,alpha4,dt4,t04,x04,y04,Vx04,Vy04,theta04);
    
    if (Vy4_water == 0)
        break
    end
    
    i=i+1;
    initialization=1;
end
%this accounts for the i+1 above since i is the number of skips sent out
i=i-1;
end