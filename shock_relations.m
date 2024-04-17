% Function to calculate the flow properties after encountering shocks wave.
function [mach_out, temp_ratio, pressure_ratio] = shock_relations(mach_in, gamma, theta,normal,oblique,expansion)
    % Defining a tolerance metric which will assist in the iterative calc
    tolerance = 0.0001; % Used as a converging metric for the iterative determination of beta

    if oblique == true
        %disp('oblique')
        % Oblique shock wave calculations
        % finding beta from theta, mach_in, and gamma
        beta = 0;
        find_beta = true;

        % Using the beta-theta-mach relationship
        while find_beta
            term1 = 2*cotd(beta);
            term2 = ((mach_in^2)*(sind(beta))^2)-1;
            term3 = (mach_in^2)*(gamma+cosd(beta*2))+2;
            exp = term1*term2/term3; % Right hand side of the beta-theta-mach equation
            
            if tand(theta) - exp < tolerance
                find_beta = false;
            else
                beta = beta + .01;
            end
        end
           
        % Getting all the parameters. 
        mach_in_normal = mach_in*sind(beta);
        pressure_ratio = 1+((2*gamma)/(gamma+1))*(mach_in_normal^2 - 1);
        temp_ratio = (1+((2*gamma)/(gamma+1))*(mach_in_normal^2-1))*(2+mach_in_normal^2*(gamma-1))/(mach_in_normal^2*(gamma+1));
        mach_out_normal = sqrt((1+mach_in_normal^2*.5*(gamma-1))/(gamma*mach_in_normal^2 - .5*(gamma-1)));
        mach_out = mach_out_normal/sind(beta-theta);

    elseif expansion == true
        %disp('expansion')
        % Expansion wave Relations!
        neu_mach_out = theta + calc_neu(gamma, mach_in);
        find_neu_mach_out = true;
        mach_out = 0;
        
        % Using prandel-meyers equation
        while find_neu_mach_out
            term1 = sqrt((gamma+1)/(gamma-1));
            term2 = atand(sqrt((gamma-1)*(mach_out^2-1)/(gamma+1)));
            term3 = atand(sqrt(mach_out^2-1));
            neu = term1*term2-term3;

            if neu_mach_out - neu < tolerance
                find_neu_mach_out = false;
            else
                mach_out = mach_out + .001;
            end
        end

        % Using the isentropic relations to find out other properties.
        p01_p1 = (1 + (gamma-1)*.5*mach_in^2)^(gamma/(gamma-1));
        t01_t1 = 1 + .5*(gamma-1)*mach_in^2;
        p02_p2 = (1 + (gamma-1)*.5*mach_out^2)^(gamma/(gamma-1));
        t02_t2 = 1 + .5*(gamma-1)*mach_out^2;
        temp_ratio = t01_t1/t02_t2;  
        pressure_ratio = p01_p1/p02_p2;

    elseif normal == true
        %disp('normal')
        pressure_ratio = 1+((2*gamma)/(gamma+1))*(mach_in^2 - 1);
        temp_ratio = (1+((2*gamma)/(gamma+1))*(mach_in^2-1))*(2+mach_in^2*(gamma-1))/(mach_in^2*(gamma+1));
        mach_out = sqrt((1+mach_in^2*.5*(gamma-1))/(gamma*mach_in^2 - .5*(gamma-1)));
    end
end

%% Calculate the prandel-myers function value!
function neu = calc_neu(gamma,mach_in)
    term1 = sqrt((gamma+1)/(gamma-1));
    term2 = atand(sqrt((gamma-1)*((mach_in^2)-1)/(gamma+1)));
    term3 = atand(sqrt((mach_in^2)-1));
    neu = term1*term2-term3;
end