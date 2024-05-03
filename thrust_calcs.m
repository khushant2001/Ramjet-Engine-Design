function [thrust,specific_thrust,specfic_fuel_consumption,specfic_impulse] = thrust_calcs(P_in,P_out,T_in,T_out,M_in,M_out,m_dot,A_in,A_out,m_dot_fuel,slope,slope_b,distance,distance_b)
    %disp("Finding the Final Thrust and Other Characteristics ...")
    R = 287;
    gaama = 1.4;
    v_in = M_in*sqrt(R*gaama*T_in);
    v_out = M_out*sqrt(R*gaama*T_out);
    thrust_solo = m_dot*(v_out-v_in) + P_out*A_out - P_in*A_in;

    % Finding the drag caused by the shell
    [~, ~, p_ratio1,~] = shock_relations(M_in,gaama,slope,0,1,0);
    [~, ~, p_ratio2,~] = shock_relations(M_in,gaama,slope_b,0,1,0);
    drag = P_in*p_ratio1*sind(slope)*distance + P_in*p_ratio2*sind(slope_b)*distance_b;
    
    % Determining the total thrust!
    thrust = thrust_solo + drag;
    specific_thrust = thrust/m_dot;
    specfic_fuel_consumption = m_dot_fuel/thrust;
    specfic_impulse = thrust/(9.81*m_dot_fuel);
    %disp(['...Specific Thrust = ', num2str(specific_thrust), ' m/sec'])
    %disp(['...Specific Fuel Consumption = ', num2str(specfic_fuel_consumption), ' sec/m'])
    %disp(['...Specific Impulse = ', num2str(specfic_impulse), ' sec'])
    %disp(['...Total Thrust Produced = ', num2str(thrust), ' N'])
end