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
mach_in = 3.25; 
[T1, P1] = atmospheric(elevation);  
gaama = 1.4;
R = 287;

%% Imposed conditions!
T_tungsten = 1900; %K. Melting point of tungsten!

% Algorithm implementation!

%% Inlet
[T2, P2, M2,x_cowl,y_cowl,x,y,inlet_area,theta1,theta2] = inlet_design(mach_in,P1,T1,1,1);
m_dot1 = (P1/(R*T1))*y_cowl(2)*mach_in*sqrt(gaama*R*T1);
m_dot2 = (P2/(R*T2))*inlet_area*M2*sqrt(T2*gaama*R);

%% Diffuser

[x_diffuser,A_diffuser,M3,P3,T3] = diffuser(M2,T2,P2,inlet_area/2,T1,mach_in);
diffuser_exit_area = 2*A_diffuser(end/2);
x_diffuser = x_diffuser + x_cowl(end);
A_diffuser = A_diffuser + (y_cowl(2) - A_diffuser(1));
m_dot3 = (P3(end)/(R*T3(end)))*diffuser_exit_area*M3(end)*sqrt(T3(end)*gaama*R);

%% Flameholder

[P3_prime,T3_prime,M3_prime] = flameholder(P3(end),M3,T3(end));
m_dot3_prime = (P3_prime/(R*T3_prime))*diffuser_exit_area*M3_prime*sqrt(gaama*R*T3_prime);

%% Combustor
% TODO! Itarate over the phi to get from the specified mach number to the
% actual mach number. 
phi = .2;
[T4,P4,M4,combustor_length,m_dot_fuel,tao] = combustor(T3_prime, P3_prime, M3_prime, phi,diffuser_exit_area);
x_combustor = [x_diffuser(end/2),x_diffuser(end/2)+combustor_length];
y_combustor = [A_diffuser(end/2),A_diffuser(end/2)];
m_dot4 = (P4/(R*T4))*diffuser_exit_area*M4*sqrt(T4*gaama*R);

%% Converging section
[M5,T5,P5,x_converging,A_converging] = converging_section(M4, T4, P4,diffuser_exit_area/2);
x_converging = x_converging + x_combustor(end);
A_converging = A_converging + (y_combustor(end) - A_converging(1));
throat_area = (A_converging(end/2)-A_converging(end));
m_dot5 = (P5(end)/(R*T5(end)))*throat_area*M5(end)*sqrt(T5(end)*gaama*R);

%% Nozzle
n = 100;
M6_req = 3.5;
[x_wall,y_wall,P6,T6,M6] = MOC_nozzle(n,M6_req, P5(end), T5(end),throat_area/2);
exit_area = y_wall(end)*2;
x_wall = [x_wall,x_wall];
y_wall = [y_wall,-y_wall];
x_wall = x_wall + x_converging(end/2);
y_wall = y_wall + (A_converging(end/2)-y_wall(1));

m_dot6 = (P6/(R*T6))*M6*sqrt(gaama*R*T6)*exit_area;

%% Saving the Area Profile for the final analysis
filename = 'Area profiles.mat';
save(filename, 'A_diffuser', 'A_converging', 'y_wall','theta1','theta2','T1','P1');

%% Overall Shell

[x_shell,y_shell] = outer_shell(x_cowl,y_cowl,x_wall,y_wall,x_diffuser,A_diffuser);


%% Thrust Calcs
total_inlet_area = 0;%y_shell(end/2)-y_shell(end);
[thrust,~,~,~] = thrust_calcs(P1,P2,P6,T1,T6,mach_in,M6,m_dot2,inlet_area,total_inlet_area,exit_area,m_dot_fuel,x_wall(end));
%m_dot_total = (P1/(R*T1))*total_inlet_area*mach_in*sqrt(gaama*R*T1);

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
%plot(x_shell(1:end/2),y_shell(1:end/2), 'b*-', 'MarkerSize', 10, 'LineWidth', 2)
%plot(x_shell(1+end/2:end),y_shell(1+end/2:end), 'b*-', 'MarkerSize', 10, 'LineWidth', 2)
grid on
box off
legend('Inlet','Cowl','Diffuser','Diffuser','Combustor','Converging Nozzle','Converging Nozzle','Diverging Nozzle','Diverging Nozzle','Casing',Location='east')
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 18);
title("Final Design")
xlabel("X [m]")
ylabel("Y [m]")