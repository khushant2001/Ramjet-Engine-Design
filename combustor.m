% This is the code for the combustor. You need to specify the inlet
% conditions: T_in,P_in,Mach_in,phi, area_in

function [T_out,P_out,mach_out,combustor_length,m_dot_fuel,tao] = combustor(T_in, P_in, mach_in,phi,area_in)
    
    %disp("Calculating Properties Across Combustor ...")

    % Declaring constants!
    gaama = 1.4;
    R = 287;
    Heat_val = 120*10^6; %J/kg Lower Heating Value of H2
    
    % Finding the mass flow rates and heat rates
    m_dot = (P_in/(R*T_in))*area_in*mach_in*sqrt(T_in*gaama*R);
    m_dot_fuel = phi*.0291* m_dot;
    q_dot = m_dot_fuel*Heat_val;
    q = q_dot/m_dot;
    %disp(['...Heat Added to the Air Mixture = ', num2str(q)])
    %disp(['...Mass flow rate of fuel = ', num2str(m_dot_fuel)])
    % Finding the stagnation temperature and pressure
    [T_01,P_01] = stagnation_values(T_in,P_in,mach_in);
    
    % Find Cp and T_02
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
    if mach_out == 0
        disp("...ERROR! TOO MUCH HEAT ADDED")
        return
    end
    P_02 = P0_star*((1+gaama)/(1+gaama*mach_out^2))*((2+mach_out^2*(gaama-1))/(gaama+1))^(gaama/(gaama-1));
    T_out = T_02 / (1 + (gaama-1)*.5*mach_out^2);
    P_out = P_02 / (1 + (gaama-1)  *.5*mach_out^2)^(gaama/(gaama-1));

    %Finding Length of Combustor
    tao = 325*(P_in/101325)^(-1.6)*exp(-0.8*T_01/1000)*10^(-6); %Time Fuel took To Burn
    v_in = mach_in*sqrt(gaama*R*T_in); %m/s Inlet Velocity of Fuel
    v_out = mach_out*sqrt(gaama*R*T_out); %m/s Inlet Velocity of Fuel
    combustor_length = 0.5*(v_in+v_out)*tao; %m Length of Combustor
    %disp(['... Combustor Length = ', num2str(combustor_length)]);
end

function [T_02, P_02] = stagnation_values(T,P,M)
    gaama = 1.4;
    T_02 = T*(1 + .5*(gaama-1)*M^2);
    P_02 = P*(1+.5*(gaama-1)*M^2)^(gaama/(gaama-1));
end