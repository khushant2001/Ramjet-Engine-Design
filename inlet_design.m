%% Function definition for the inlet design:

% Things to think about!
% How to find the angles that increase stag pressure ratio.
% How to find the length of the steps (x and y points) - DONE
% How to find the point where shocks will intersect! - DONE
% How to find the point where the cowl exists! - DONE


function [T_2, P_2, M_2,x_cowl,y_cowl,x,y,area] = inlet_design(mach_in, P_in, T_in,L1,L2)
    gaama = 1.4;
    [~,P_01] = stagnation_values(T_in,P_in,mach_in);

    % First step
    theta1 = 5;
    [M_12, t_ratio1, p_ratio1,beta1] = shock_relations(mach_in,gaama,theta1,0,1,0);
    
    % Second step. Might need move this to find the best optimzation for
    % the given mach number.
    theta2 = 10;
    [M_22, t_ratio2, p_ratio2,beta2] = shock_relations(M_12,gaama,theta2,0,1,0);

    % Normal shock at the end. 
    [M_2, t_ratio, p_ratio,~] = shock_relations(M_22,gaama,0,1,0,0);
    T_2 = T_in*t_ratio1*t_ratio2*t_ratio;
    P_2 = P_in*p_ratio1*p_ratio2*p_ratio;

    % Calling the cowl function!
    [x_cowl,y_cowl] = cowl_design(theta1,theta2,beta1,beta2,L1,L2);

    % Generating points for the inlet
    x = [0,L1*cosd(theta1),L2*cosd(theta2) + L1*cosd(theta1)];
    y = [0,L1*sind(theta1),L2*sind(theta2) + L1*sind(theta1)];

    % Determining the inlet_area
    theta_temp = atand((y_cowl(1) - y(2))/(x_cowl(1) - x(2))) - theta2;
    temp_length = sqrt((x_cowl(1) - x(2))^2 + (y_cowl(1) - y(2))^2);
    area = sind(theta_temp)*temp_length;
    disp(["Intake area = ",num2str(area)]);
end

function [x_cowl,y_cowl] = cowl_design(theta1,theta2,beta1,beta2,L,L2)
    x_cowl = [];
    y_cowl = [];
    tolerance = 0.0001;
    y = 0;
    alpha = sind(beta2-theta1)*L/(sind(180-beta2));
    while 1
        exp = y*tand(beta2)/(alpha*tand(beta2) + y);
        if tand(beta1) - exp < tolerance
            break
        else
            y = y + 0.01;
        end
    end
    x = alpha + y/tand(beta2);
    x_cowl = [x,x+L2*cosd(theta2)];
    y_cowl = [y,y+L2*sind(theta2)];
end

function [T_02, P_02] = stagnation_values(T,P,M)
    gaama = 1.4;
    T_02 = T*(1 + .5*(gaama-1)*M^2);
    P_02 = P*(1+.5*(gaama-1)*M^2)^(gaama/(gaama-1));
end