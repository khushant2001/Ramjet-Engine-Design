 %{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ramjet Engine Design %
    % Angel, Christos, Khushant %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%% Things to take care of!
% Get the turning geometry working!
% Make sure energy is same going in and out!

clc
clear 
close all

%% Inlet condtions!

elevation = 16764; %m
mach_in = 3.25; % Having 2 designs and then comparing in the 'final_analysis file'.
[T1, P1] = atmospheric(elevation);  
gaama = 1.4;
R = 287;

%% Imposed conditions!
T_tungsten = 1900; %K. Melting point of tungsten!

% Algorithm implementation!

M6_array = linspace(.25,.5,10);
% Arrays to store data!

thrust = [];
total_length = [];
area_in = [];
air_flow_rate = [];
fuel_flow_rate = [];
phi = [];
T = [];
for i = 1:length(M6_array)
    disp(i)
    %% Inlet
    L1 = .2;
    L2 = 1;
    [T2, P2, M2,x_cowl,y_cowl,x,y,inlet_area,theta1,theta2,inlet_drag] = inlet_design(mach_in,P1,T1,L1,L2);
    m_dot1 = (P1/(R*T1))*y_cowl(2)*mach_in*sqrt(gaama*R*T1);
    m_dot2 = (P2/(R*T2))*inlet_area*M2*sqrt(T2*gaama*R);
    area_in(i) = inlet_area;
 
    %% Diffuser
    
    [x_diffuser,A_diffuser,M3,P3,T3] = diffuser(M2,T2,P2,inlet_area/2,T1,mach_in,.2);
    diffuser_exit_area = 2*A_diffuser(end/2);
    x_diffuser = x_diffuser + x_cowl(end);
    A_diffuser = A_diffuser + (y_cowl(2) - A_diffuser(1));
    m_dot3 = (P3(end)/(R*T3(end)))*diffuser_exit_area*M3(end)*sqrt(T3(end)*gaama*R);
    air_flow_rate(i) = m_dot2;
    %% Flameholder
    
    [P3_prime,T3_prime,M3_prime] = flameholder(P3(end),M3,T3(end));
    m_dot3_prime = (P3_prime/(R*T3_prime))*diffuser_exit_area*M3_prime*sqrt(gaama*R*T3_prime);
    
    %% Combustor
    % TODO! Itarate over the phi to get from the specified mach number to the
    % actual mach number. 
    phi_in = 0;
        while 1
            [T4,P4,M4,length_comb,m_dot_fuel,tao] = combustor(T3_prime,P3_prime,M3_prime,phi_in,diffuser_exit_area);
            if M4 >= M6_array(i) % This number needs to be changed accordingly. 
                fuel_flow_rate(i) = m_dot_fuel;
                phi(i) = phi_in;
                break
            else
                phi_in = phi_in+.0001;
            end
        end
    x_combustor = [x_diffuser(end/2),x_diffuser(end/2)+length_comb];
    y_combustor = [A_diffuser(end/2),A_diffuser(end/2)];
    m_dot4 = (P4/(R*T4))*diffuser_exit_area*M4*sqrt(T4*gaama*R);
    T(i) = T4;
    if T4 > T_tungsten
        disp("... ABOUT TO MELT!!! TOO MUCH HEAT ADDED")
    end
    %% Converging section
    [M5,T5,P5,x_converging,A_converging] = converging_section(M4, T4, P4,diffuser_exit_area/2);
    x_converging = x_converging + x_combustor(end);
    A_converging = A_converging + (y_combustor(end) - A_converging(1));
    throat_area = (A_converging(end/2)-A_converging(end));
    m_dot5 = (P5(end)/(R*T5(end)))*throat_area*M5(end)*sqrt(T5(end)*gaama*R);
    
    %% Nozzle
    n = 20;
    M6_req = 3.5;
    [x_wall,y_wall,P6,T6,M6] = MOC_nozzle(n,M6_req, P5(end), T5(end),throat_area/2);
    exit_area = y_wall(end)*2;
    x_wall = [x_wall,x_wall];
    y_wall = [y_wall,-y_wall];
    x_wall = x_wall + x_converging(end/2);
    y_wall = y_wall + (A_converging(end/2)-y_wall(1));
    
    m_dot6 = (P6/(R*T6))*M6*sqrt(gaama*R*T6)*exit_area;
    
    %% Overall Shell
    [x_shell,y_shell,x_shell_b,y_shell_b,slope,slope_b,distance,distance_b] = outer_shell(x_cowl,y_cowl,x_wall,y_wall,x_diffuser,A_diffuser,x,y);
    
    
    %% Thrust Calcs
    [thrust_got,~,~,~] = thrust_calcs(P1,P6,T1,T6,mach_in,M6,m_dot2,inlet_area,exit_area,m_dot_fuel,slope,abs(slope_b),distance,distance_b,inlet_drag);
    %m_dot_total = (P1/(R*T1))*total_inlet_area*mach_in*sqrt(gaama*R*T1);
    thrust(i) = thrust_got;
    total_length(i) = x_wall(end/2);
end

figure
plot(M6_array,thrust,LineWidth=5)
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 18);
xlabel("M",FontSize=18)
ylabel("N",FontSize=18)
title({'Thrust generated vs Combustor Exit Mach';'(Inlet Mach = 3.25)'},FontSize=18)
grid on
box off

figure
plot(M6_array,total_length,LineWidth=5)
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 18);
xlabel("M",FontSize=18)
ylabel("m",FontSize=18)
title({'Length of the Engine vs Combustor Exit Mach';'(Inlet Mach = 3.25)'},FontSize=18)
grid on
box off

figure
plot(M6_array,air_flow_rate,LineWidth=5)
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 18);
xlabel("M",FontSize=18)
ylabel("kg/sec",FontSize=18)
title({'Air flow rate vs Combustor Exit Mach';'(Inlet Mach = 3.25)'},FontSize=18)
grid on
box off

figure
plot(M6_array,area_in,LineWidth=5)
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 18);
xlabel("M",FontSize=18)
ylabel("m^2",FontSize=18)
title({'Inlet area vs Combustor Exit Mach';'(Inlet Mach = 3.25)'},FontSize=18)
grid on
box off

figure
plot(M6_array,fuel_flow_rate,LineWidth=5)
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 18);
xlabel("M",FontSize=18)
ylabel("kg/sec",FontSize=18)
title({'Fuel flow rate vs Combustor Exit Mach';'(Inlet Mach = 3.25)'},FontSize=18)
grid on
box off

figure
plot(M6_array,phi,LineWidth=5)
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 18);
xlabel("M",FontSize=18)
ylabel('$\phi$', 'Interpreter', 'latex',FontSize=18)
title({'Equivalance Ratio vs Combustor Exit Mach';'(Inlet Mach = 3.25)'},FontSize=18)
grid on
box off
figure
plot(M6_array,T,LineWidth=5)
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 18);
xlabel("M",FontSize=18)
ylabel("K",FontSize=18)
title({'Post-Combustor Temperature vs Combustor Exit Mach';'(Inlet Mach = 3.25)'},FontSize=18)
grid on
box off
%{
%% Plotting the engine!
figure
hold on
plot(x, y, 'ro-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_cowl,y_cowl,'ko-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_diffuser(1:end/2), A_diffuser(1:end/2), 'go-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_diffuser(1+end/2:end), A_diffuser(1+end/2:end), 'go-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_combustor, y_combustor, 'co-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_converging(1:end/2), A_converging(1:end/2), 'bo-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_converging(1+end/2:end), A_converging(1+end/2:end), 'bo-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_wall(1:end/2),y_wall(1:end/2), 'mo-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_wall(1+end/2:end),y_wall(1+end/2:end), 'mo-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_shell,y_shell, 'b*-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_shell_b,y_shell_b, 'b*-', 'MarkerSize', 10, 'LineWidth', 2)
grid on
box off
legend('Inlet','Cowl','Diffuser','','Combustor','Converging Nozzle','','Diverging Nozzle','','Casing','',Location='east')
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 18);
title(["Final Design for Mach in = ",num2str(mach_in)])
xlabel("X [m]")
ylabel("Y [m]")
%}