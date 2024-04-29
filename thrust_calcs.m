function [thrust] = thrust_calcs(P_in,P2,P_out,T_in,T_out,M_in,M_out,m_dot,A_in,total_A_in,A_out,m_dot_fuel,total_length)
    disp("Finding the Final Thrust and Other Characteristics ...")
    R = 287;
    gaama = 1.4;
    v_in = M_in*sqrt(R*gaama*T_in);
    v_out = M_out*sqrt(R*gaama*T_out);
    thrust_solo = m_dot*(v_out-v_in) - P_out*A_out + P2*A_in + P_in*(total_A_in-A_in);
    drag = 0;
    thrust = thrust_solo - drag;
    specific_thrust = thrust/m_dot;
    specfic_fuel_consumption = m_dot_fuel/thrust;
    specfic_impulse = thrust/(9.81*m_dot_fuel);
    disp(['...Specific Thrust = ', num2str(specific_thrust)])
    disp(['...Specific Fuel Consumption = ', num2str(specfic_fuel_consumption)])
    disp(['...Specific Impulse = ', num2str(specfic_impulse)])
    disp(['...Total Thrust Produced = ', num2str(thrust), ' N'])
end