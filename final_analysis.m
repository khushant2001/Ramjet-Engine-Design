% This is the file that needs to be run in order to get the things solved
% for all stuff!

%% TODO!!
% The converging and diverging nozzle are hardcoded right now. DONT KNOW
% WHATS WRONG!
clc
clear 
close all

%% Inlet conditions
gaama = 1.4;
R = 287;
mach_in = linspace(2.75,3.25,10);

%% Declaring arrays to store values.

thrust = [];
phi = [];
m_dot_fuel = [];
specific_thrust = [];
specfic_fuel_consumption = [];
specfic_impulse = [];
max_temp = [];
max_pressure = [];
combustor_length = [];
m_dot_air = [];
combustor_time = [];
M3_actual = [];

%% Running the analysis
for k = 1:2
    % Importing Area Profiles of diffuser, converging and diverging nozzles!
    if k == 1
        data = load('Area profiles_min.mat');
    else
        data = load('Area profiles_max.mat');
    end
    A_diffuser = data.A_diffuser;
    A_converging = data.A_converging;
    A_nozzle = data.y_wall;
    theta1 = data.theta1;
    theta2 = data.theta2;
    T1 = data.T1;
    P1 = data.P1;
    M4_req = data.M4;
    slope = data.slope;
    slope_b = data.slope_b;
    distance = data.distance;
    distance_b = data.distance_b;

    for j = 1:length(mach_in) 
        disp(j)
    %% Get properties after inlet!
        inlet_area = (A_diffuser(1)-A_diffuser(1+end/2));
        % First step
        [M_12, t_ratio1, p_ratio1,beta1] = shock_relations(mach_in(j),gaama,theta1,0,1,0);
        
        % Second step. Might need move this to find the best optimzation for
        % the given mach number.
        [M_22, t_ratio2, p_ratio2,beta2] = shock_relations(M_12,gaama,theta2,0,1,0);
    
        % Normal shock at the end. 
        [M2, t_ratio, p_ratio,~] = shock_relations(M_22,gaama,0,1,0,0);
    
        T2 = T1*t_ratio1*t_ratio2*t_ratio;
        P2 = P1*p_ratio1*p_ratio2*p_ratio;
        m_dot2 = (P2/(R*T2))*inlet_area*M2*sqrt(T2*gaama*R);
        m_dot_air(k,j) = m_dot2;

    %% Get properties from diffuser
        area_ratio = (A_diffuser(1)-A_diffuser(1+end/2))/(A_diffuser(end/2)-A_diffuser(end));
        M3 = 0;
        tolerance = 0.0001;
        while 1
            exp2 = (M3/M2)*((1+.5*(gaama-1)*M2^2)/(1+.5*(gaama-1)*M3^2))^((gaama+1)/(2*(gaama-1)));
            if area_ratio - exp2 < tolerance
                break
            else
                M3 = M3 + 0.001;
            end
        end
        M3_actual(k,j) = M3;
        % Finding stagnation values at the inlet. 
        T_02 = T2*(1 + .5*(gaama-1)*M2^2);
        P_02 = P2*(1 + .5*(gaama-1)*M2^2)^(gaama/(gaama-1));
    
        T3 = T_02/(1+.5*(gaama-1)*M3^2);
        P3 = P_02/((1+.5*(gaama-1)*M3^2)^(gaama/(gaama-1)));
        
        diffuser_exit_area = A_diffuser(end/2) - A_diffuser(end);
        m_dot3 = (P3/(R*T3))*diffuser_exit_area*M3*sqrt(T3*gaama*R);
    
    %% Get properties from flameholder
        [P3_prime,T3_prime,M3_prime] = flameholder(P3,M3,T3);
        m_dot3_prime = (P3_prime/(R*T3_prime))*diffuser_exit_area*M3_prime*sqrt(gaama*R*T3_prime);
    
    %% Get properties from combustor
        phi_in = 0;
        while 1
            [T4,P4,M4,length_comb,fuel_rate,tao] = combustor(T3_prime,P3_prime,M3_prime,phi_in,diffuser_exit_area);
            if M4 >= M4_req % This number needs to be changed accordingly. 
                phi(k,j) = phi_in;
                m_dot_fuel(k,j) = fuel_rate;
                combustor_length(k,j) = length_comb;
                combustor_time(k,j) = tao;
                break
            else
                phi_in = phi_in+.0001;
            end
        end
        m_dot4 = (P4/(R*T4))*diffuser_exit_area*M4*sqrt(T4*gaama*R);
    
    %% Get properties from converging nozzle
        %disp("Converging section...")
        area_ratio = ((A_converging(1)-A_converging(1+end/2))/(A_converging(end/2)-A_converging(end)));
        M5 = M4;
        tolerance = 0.0001;
        while 1
            exp2 = (M5/M4)*((1+.5*(gaama-1)*M4^2)/(1+.5*(gaama-1)*M5^2))^((gaama+1)/(2*(gaama-1)));
            if area_ratio - exp2 < tolerance || M5 > 1
                break
            else
                M5 = M5 + 0.001;
            end
        end
    
        % Finding stagnation values at the inlet. 
        T_04 = T4*(1 + .5*(gaama-1)*M4^2);
        P_04 = P4*(1 + .5*(gaama-1)*M4^2)^(gaama/(gaama-1));
    
        T5 = T_04/(1+.5*(gaama-1)*M5^2);
        P5 = P_04/((1+.5*(gaama-1)*M5^2)^(gaama/(gaama-1)));
        
        throat_area = A_converging(end/2) - A_converging(end);
        m_dot5 = (P5/(R*T5))*throat_area*M5*sqrt(T5*gaama*R);
    
    %% Get properties from diverging nozzle
        
        %disp("Diverging nozzle...")
        area_ratio = ((A_nozzle(1)-A_nozzle(1+end/2))/(A_nozzle(end/2)-A_diffuser(end)))^2;
        M6 = 1;
        tolerance = 0.0001;
        while 1
            exp2 = (1/M6^2)*((2/(gaama+1))*(1+.5*(gaama-1)*M6^2))^((gaama+1)/(gaama-1));
            if (area_ratio - exp2) < tolerance
                break
            else
                M6 = M6 + 0.001;
            end
        end
        M6 = 3.5;
        % Finding stagnation values at the inlet. 
        T_05 = T5*(1 + .5*(gaama-1)*M5^2);
        P_05 = P5*(1 + .5*(gaama-1)*M5^2)^(gaama/(gaama-1));
    
        T6 = T_05/(1+.5*(gaama-1)*M6^2);
        P6 = P_05/((1+.5*(gaama-1)*M6^2)^(gaama/(gaama-1)));
        exit_area = A_nozzle(end/2)-A_nozzle(end);
        m_dot6 = (P6/(R*T6))*M6*sqrt(gaama*R*T6)*exit_area;
        
    %% Find the overall thurst!!
        [force,one,two,three] = thrust_calcs(P1,P6,T1,T6,mach_in(j),M6,m_dot5,inlet_area,exit_area,m_dot_fuel(k,j),slope,slope_b,distance,distance_b);
        thrust(k,j) = force;
        specific_thrust(k,j) = one;
        specfic_fuel_consumption(k,j) = two;
        specfic_impulse(k,j) = three;
    
        % Finding the max temperature and pressure
        max_pressure(k,j) = max([P1,P2,P3,P4,P5,P6,P3_prime]);
        max_temp(k,j) = max([T1,T2,T3,T4,T5,T6,T3_prime]);
    end
end

figure
set(gcf, 'Color', 'white');
plot(mach_in,phi,LineWidth=5)
legend('M1 = 2.75', 'M2 = 3.25')
set(gca, 'FontSize', 18);
title("Equivalce ratio vs Mach number to Achieve M4 = ", num2str(M4),FontSize=18)
xlabel("Mach number",FontSize=18)
ylabel('$\phi$', 'Interpreter', 'latex',FontSize=18)
grid on
box off

figure
set(gcf, 'Color', 'white');
plot(mach_in,m_dot_fuel,LineWidth=5)
legend('M1 = 2.75', 'M2 = 3.25')
set(gca, 'FontSize', 18);
title("Fuel flow rate required vs Mach number to Achieve M4 = ", num2str(M4),FontSize=18)
xlabel("Mach number",FontSize=18)
ylabel('m_dot [kg/sec]',FontSize=18)
grid on
box off

figure
set(gcf, 'Color', 'white');
plot(mach_in,max_pressure,LineWidth=5)
legend('M1 = 2.75', 'M2 = 3.25')
set(gca, 'FontSize', 18);
title("Max Pressure vs Mach number",FontSize=18)
xlabel("Mach number",FontSize=18)
ylabel("P [Pa]",FontSize=18)
grid on
box off

figure
set(gcf, 'Color', 'white');
plot(mach_in,max_temp,LineWidth=5)
legend('M1 = 2.75', 'M2 = 3.25')
set(gca, 'FontSize', 18);
title("Max Temperature vs Mach number",FontSize=18)
xlabel("Mach number",FontSize=18)
ylabel("T [K]",FontSize=18)
grid on
box off

figure
set(gcf, 'Color', 'white');
plot(mach_in,thrust,LineWidth=5)
legend('M1 = 2.75', 'M2 = 3.25')
set(gca, 'FontSize', 18);
title("Thrust generated vs Mach Number",FontSize=18)
xlabel("Mach number",FontSize=18)
ylabel("Thrust [N]")
grid on
box off

figure
set(gcf, 'Color', 'white');
plot(mach_in,specfic_fuel_consumption,LineWidth=5)
legend('M1 = 2.75', 'M2 = 3.25')
set(gca, 'FontSize', 18);
title("Specific fuel consumption vs Mach Number",FontSize=18)
xlabel("Mach number",FontSize=18)
ylabel("sec/m")
grid on
box off

figure
set(gcf, 'Color', 'white');
plot(mach_in,specific_thrust,LineWidth=5)
legend('M1 = 2.75', 'M2 = 3.25')
set(gca, 'FontSize', 18);
title("Specific Thrust vs Mach Number",FontSize=18)
xlabel("Mach number",FontSize=18)
ylabel("m/sec")
grid on
box off

figure
set(gcf, 'Color', 'white');
plot(mach_in,specfic_impulse,LineWidth=5)
legend('M1 = 2.75', 'M2 = 3.25')
set(gca, 'FontSize', 18);
title("Specific Impulse vs Mach Number",FontSize=18)
xlabel("Mach number",FontSize=18)
ylabel("sec")
grid on
box off

figure
set(gcf, 'Color', 'white');
plot(mach_in,combustor_length,LineWidth=5)
legend('M1 = 2.75', 'M2 = 3.25')
set(gca, 'FontSize', 18);
title("Combustor length vs Mach Number to Achieve M4 = ", num2str(M4),FontSize=18)
xlabel("Mach number",FontSize=18)
ylabel("L [m]")
grid on
box off

figure
set(gcf, 'Color', 'white');
plot(mach_in,m_dot_air,LineWidth=5)
legend('M1 = 2.75', 'M2 = 3.25')
set(gca, 'FontSize', 18);
title("Mass Flow Rate of Air vs Mach number",FontSize=18)
xlabel("Mach number",FontSize=18)
ylabel("M_dot_air [kg/sec]")
grid on
box off

figure
set(gcf, 'Color', 'white');
plot(mach_in,combustor_time,LineWidth=5)
legend('M1 = 2.75', 'M2 = 3.25')
set(gca, 'FontSize', 18);
title("Time to burn fuel vs Mach number",FontSize=18)
xlabel("Mach number",FontSize=18)
ylabel("t [sec]")
grid on
box off

figure
set(gcf, 'Color', 'white');
plot(mach_in,M3_actual,LineWidth=5)
legend('M1 = 2.75', 'M2 = 3.25')
set(gca, 'FontSize', 18);
title("Variation in diffuser's exit mach number designed to be = 0.2",FontSize=18)
xlabel("Mach number",FontSize=18)
ylabel("M3")
grid on
box off