%% 1976 Atmospheric Model

function [temp, p] = atmospheric(elevation)
    if elevation < 11000
        [temp, p] = first(elevation);
    elseif elevation >= 11000 && elevation <20000
        [temp, p] = second(elevation);
    elseif elevation >= 20000 && elevation < 32000
        [temp, p] = third(elevation);
    elseif elevation >= 32000 && elevation < 47000 
        [temp, p] = fourth(elevation);
    else
        [temp, p] = [0, 0];
    end
end

function [temp,p] = first(elevation)
    T_infinity = 288.15; % K
    g = 9.81; % m/sec^2
    R = 287; % J/kg*K
    P_infinity = 101325; % Pa
    k = -6.5/1000;
    temp = T_infinity+k*elevation;
    p = P_infinity*(temp/T_infinity)^(-g/(R*k));
end

function [temp,p] = second(elevation)
    T_infinity = 288.15; % K
    g = 9.81; % m/sec^2
    R = 287; % J/kg*K
    k = -6.5/1000;
    P_infinity = 101325; % Pa
    temp = T_infinity+k*11000;
    p = (P_infinity*(temp/T_infinity)^(-g/(R*k)))*exp(-g*(elevation-11000)/(R*temp));
end

function [temp,p] = third(elevation)
    T_infinity = 288.15; % K
    g = 9.81; % m/sec^2
    R = 287; % J/kg*K
    P_infinity = 101325; % Pa
    k1 = -6.5/1000;
    k = 1/1000;
    temp1 = (T_infinity-(6.5/1000)*11000);
    temp = temp1 +(1/1000)*(elevation-20000);
    p = ((P_infinity*(temp1/T_infinity)^(-g/(R*k1)))*exp(-g*(9000)/(R*temp1)))*(temp/temp1)^(-g/(R*k));
end

function [temp,p] = fourth(elevation)
    T_infinity = 288.15; % K
    g = 9.81; % m/sec^2
    R = 287; % J/kg*K
    P_infinity = 101325; % Pa
    k = 2.8/1000;
    k1 = -6.5/1000;
    k2 = 1/1000;
    temp1 = (T_infinity-(6.5/1000)*11000);
    temp2 = temp1;
    temp3 = temp2 +(1/1000)*(32000-20000);
    temp = temp3+(2.8/1000)*(elevation-32000);
    p = (((P_infinity*(temp1/T_infinity)^(-g/(R*k1)))*exp(-g*(9000)/(R*temp1)))*(temp3/temp2)^(-g/(R*k2)))*(temp/temp3)^(-g/(R*k));
end