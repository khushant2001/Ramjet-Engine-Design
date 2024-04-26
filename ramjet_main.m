%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ramjet Engine Design %
    % Angel, Christos, Khushant %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

% Thing to think about!
% 1. Should have an x and y array passing around the functions to get the
% location of every single component. 
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

%T3 = 1500; % K. To prevent high temperatures! Need to talk about this
T_tungsten = 1900; %K. Melting point of tungsten!

%% Algorithm implementation!

% Inlet
[T2, P2, M2,x_cowl,y_cowl,x,y,inlet_area] = inlet_design(3.25,P1,T1,5,1);

% Diffuser

[x_diffuser,A_diffuser,M3,P3,T3] = diffuser(M2,T2,P2,inlet_area/2,T1,mach_in);
diffuser_exit_area = A_diffuser(end/2);
x_diffuser = x_diffuser + x(end);
A_diffuser = A_diffuser + (y_cowl(2) - A_diffuser(1));

% Flameholder
P_3_prime = flameholder(P3(end),M3);

% Combustor
[T4,P4,M4,combustor_length] = combustor(T3(end), P_3_prime, M3, .007, 1.4);
x_combustor = [x_diffuser(end/2),x_diffuser(end/2)+combustor_length];
y_combustor = [A_diffuser(end/2),A_diffuser(end/2)];

% Converging section
[M5,T5,P5,x_converging,A_converging] = converging_section(M4, T4, P4,diffuser_exit_area);
x_converging = x_converging + x_combustor(end);
A_converging = A_converging + (y_combustor(end) - A_converging(1));

% Nozzle
n = 20;
M6 = 3;
[x_wall,y_wall,P_out,T_out] = MOC_nozzle(n,M6, P5(end), T5(end));
x_wall = [x_wall,x_wall];
y_wall = [y_wall,-y_wall];
x_wall = x_wall + x_converging(end/2);
y_wall = y_wall + (A_converging(end/2)-y_wall(1));
x_wall = [x_converging(end/2),x_wall(1:end/2),x_converging(end),x_wall(1+end/2:end)];
y_wall = [A_converging(end/2),y_wall(1:end/2),A_converging(end),y_wall(1+end/2:end)];

%% Thrust Calcs
exit_area = y_wall(end/2) - y_wall(end);
m_dot = (P1/(R*T1))*inlet_area*mach_in*sqrt(T1*gaama*R);
thrust = thrust_calcs(P1,P_out,T1,T_out,mach_in,M6,m_dot,inlet_area,exit_area,1);

%% Plotting the engine!
figure
hold on
plot(x, y, 'ro-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_cowl,y_cowl,'bo-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_diffuser(1:end/2), A_diffuser(1:end/2), 'go-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_diffuser(1+end/2:end), A_diffuser(1+end/2:end), 'go-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_combustor, y_combustor, 'ro-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_converging(1:end/2), A_converging(1:end/2), 'bo-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_converging(1+end/2:end), A_converging(1+end/2:end), 'bo-', 'MarkerSize', 10, 'LineWidth', 2)
%plot(x_wall(1:end/2),y_wall(1:end/2), 'go-', 'MarkerSize', 10, 'LineWidth', 2)
%plot(x_wall(1+end/2:end),y_wall(1+end/2:end), 'go-', 'MarkerSize', 10, 'LineWidth', 2)