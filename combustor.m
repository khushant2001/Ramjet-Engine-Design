%Combustor
clear; clc; close;

% function [mach_out, T_out, P_out. Area_out] = combustor(T_in,P_in,T_0_in,M_in,fuel_mass)
% end

%To-Do
%Change Function Structure to be T_in and T_out
%Added Section to Find Length of Combustor 

%Test
T_in = 300; %k Temperature
T0_in = 273;
P_in = 101325; %Pa Pressure
P_atm = 101325; %Pa Atmospheric Pressure
mach_in = 0.2; %Inlet Mach Numbeer
fuel_mass = 100; %kg Mass of Fuel
gaama = 1.4; %Specific Heat Ratio of Working Fluid

[T_out,P_out,mach_out, combustor_length] = combustor(T_in, P_in, mach_in, T0_in, fuel_mass, gaama, P_atm); %Call Function

function [T_out,P_out,mach_out,combustor_length] = combustor(T_in, P_in, mach_in, T0_in, fuel_mass, gaama, P_atm)
    
    Heat_val = 120*10^3; %kJ/kg Lower Heating Value of H2
    q = Heat_val*fuel_mass; %kJ
    R = 287;

    % Finding the stagnation temperature and pressure
    T_01 = T_in*(1 + (gaama-1)*.5*mach_in^2);
    P_01 = P_in*(1+(gaama-1)*.5*mach_in^2)^(gaama/(gaama-1));
    
    % Find Cp
    cp = gaama*R/(gaama-1);
    T_02 = (q/cp) + T_01;
    % Finding choked properties
    T01_T0star = (gaama+1)*mach_in^2*(2 + (gaama-1)*mach_in^2)/(1 + gaama*mach_in^2)^2;
    T0_star = T_01/T01_T0star;
    P01_P0star = ((1+gaama)/(1+gaama*mach_in^2))*((2+mach_in^2*(gaama-1))/(gaama+1))^(gaama/(gaama-1));
    P0_star = P_01/P01_P0star;
    % Setting up the iterative determination of mach_out
    mach_out = 0;
    tolerance = 0.0001;
    exp1 = T_02/T0_star;
    while 1
        exp2 = (gaama+1)*mach_out^2*(2 + (gaama-1)*mach_out^2)/(1 + gaama*mach_out^2)^2;
        if exp1-exp2 < tolerance || mach_out > 1
            if mach_out > 1
                mach_out = 0;
            end
            break
        else
            mach_out = mach_out + 0.0001;
        end
    end
    
    P_02 = P0_star*((1+gaama)/(1+gaama*mach_out^2))*((2+mach_out^2*(gaama-1))/(gaama+1))^(gaama/(gaama-1));
    T_out = T_02 / (1 + (gaama-1)*.5*mach_out^2);
    P_out = P_02 / (1 + (gaama-1)  *.5*mach_out^2)^(gaama/(gaama-1));

    %Finding Length of Combustor
    tao = 325*(P_out/101325)^(1.6) *exp(-0.8*T0_in/1000)*10^(-6); %Time Fuel took To Burn
    v_in = mach_in*sqrt(gaama*R*T_in); %m/s Inlet Velocity of Fuel
    v_out = mach_out*sqrt(gaama*R*T_out); %m/s Inlet Velocity of Fuel
    combustor_length = 0.5*(v_in+v_out)*tao; %m Length of Combustor
    disp(['Combustor Length = ', num2str(combustor_length)]);

    %Finding Area Out

end
