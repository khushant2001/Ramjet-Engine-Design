%% Function definition for the inlet design:

% Things to think about!
% How to find the angles that increase stag pressure ratio.
% How to find the length of the steps (x and y points) - DONE
% How to find the point where shocks will intersect! - DONE
% How to find the point where the cowl exists! - DONE
% WHAT IS THE MINIMUM amount that the flow needs to be turned!!!

function [T_2, P_2, M_2,x_cowl,y_cowl,x,y,area] = inlet_design(mach_in, P_in, T_in,L1,L2)
    
    disp("Calculating Properties Across Inlet...");
    gaama = 1.4;
    % Finding the optimized angles:
    [theta1,theta2] = find_angles(mach_in,P_in);
    disp(['...Theta1 = ',num2str(theta1)]);
    disp(['...Theta2 = ',num2str(theta2)]);
    % First step
    [M_12, t_ratio1, p_ratio1,beta1] = shock_relations(mach_in,gaama,theta1,0,1,0);
    
    % Second step. Might need move this to find the best optimzation for
    % the given mach number.
    [M_22, t_ratio2, p_ratio2,beta2] = shock_relations(M_12,gaama,theta2,0,1,0);

    % Normal shock at the end. 
    [M_2, t_ratio, p_ratio,~] = shock_relations(M_22,gaama,0,1,0,0);
    
    if imag(M_22) ~= 0 || M_22>M_12
        disp('ERROR!!!! FLOW CANT BE TURNED THIS MUCH.')
        return
    end
    T_2 = T_in*t_ratio1*t_ratio2*t_ratio;
    P_2 = P_in*p_ratio1*p_ratio2*p_ratio;

    % Calling the cowl function!
    [x_cowl,y_cowl] = cowl_design(theta1,theta2,beta1,beta2,L1,L2);

    % Generating points for the inlet
    x = [0,L1*cosd(theta1),x_cowl(end)];
    y = [0,L1*sind(theta1),L1*sind(theta1)+L2*sind(theta2)];

    % Determining the inlet_area
    theta_temp = atand((y_cowl(1) - y(2))/(x_cowl(1) - x(2))) - theta2;
    temp_length = sqrt((x_cowl(1) - x(2))^2 + (y_cowl(1) - y(2))^2);
    
    % TODO: NEEDS TO BE CHECKED
    area = y_cowl(2) - y(3);%sind(theta_temp)*temp_length;
    disp(['...Intake area = ',num2str(area)]);
end

% To find the set of 2 angles that maximize the stagnation pressures. 
function [theta1,theta2] = find_angles(mach_in, P_in)
    gaama = 1.4;
    [~,P_01] = stagnation_values(0,P_in,mach_in);
    theta_1 = linspace(1,30,10);
    theta_2 = linspace(20,20,10);
    ratios = zeros(length(theta_1),length(theta_2));

    % Running a 2d loop!
    for i = 1:length(theta_1)
        for j = 1:length(theta_2)
            [M_12, ~, p_ratio1,~] = shock_relations(mach_in,gaama,theta_1(i),0,1,0);
            [M_22, ~, p_ratio2,~] = shock_relations(M_12,gaama,theta_2(j),0,1,0);
            [~,P_02] = stagnation_values(0,P_in*p_ratio1*p_ratio2,M_22);
            ratios(i,j) = P_02/P_01;
        end
    end

    % Find the set of angles!
    [~, linearIndex] = max(ratios(:));
    [rowIndex, colIndex] = ind2sub(size(ratios), linearIndex);
    theta1 = theta_1(rowIndex);
    theta2 = theta_2(colIndex);
end

% Design the cowl. 
function [x_cowl,y_cowl] = cowl_design(theta1,theta2,beta1,beta2,L,L2)
    disp('Desinging the Cowl...')
    x_cowl = [];
    y_cowl = [];
    tolerance = 0.001;
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

% Find the stagnation values. 
function [T_02, P_02] = stagnation_values(T,P,M)
    gaama = 1.4;
    T_02 = T*(1 + .5*(gaama-1)*M^2);
    P_02 = P*(1+.5*(gaama-1)*M^2)^(gaama/(gaama-1));
end