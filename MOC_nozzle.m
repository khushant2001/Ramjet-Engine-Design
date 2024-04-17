function [c_plus,wall] = MOC_nozzle(n,gaama,mach_exit, P_star, T_star)

    % Caluclating theta_max for the nozzle!
    theta_max = calc_neu(gaama, mach_exit)*.5;

    % Calculating the stagnation conditions at throat!
    P_0 = P_star/0.528;
    T_0 = T_star/0.833;

    % Generating the characteristic lines: c plus and c minus
    [c_plus, c_minus,wall,centerline] = characteristic_lines(n,gaama,theta_max);
    
    % Looping over the right turning characteristics which are stored as
    % dictionaries in the c_plus cell.
    for i = 1:numel(c_plus)

        % Getting the individual dictionary
        temp = c_plus{i}; 

        % Getting the keys and the values of the entire dictionary!
        k = keys(temp);
        v = values(temp);

        % Looping over the contents of the dictionary. 
        for j = 1:length(k)

            % Accessing the key and the value of the specific element in
            % the dictionary. 
            key = k(j);
            value = v(j);

            % Checking if the point is on the centerline

            if ismember(key, centerline)

                % Getting the values of the left turning characteristics
                % which are required to get the values. 
                c_minus_values = c_minus(j+i-1);

                % Updating the properties of the specific point

                value{1}(2) = 0; % y pos
                value{1}(3) = 0; % theta
                value{1}(4) = c_minus_values{1}(6); % neu
                value{1}(5) = calc_mach(gaama, value{1}(4)); % M
                value{1}(6) = asind(1/value{1}(5)); % nu
                value{1}(7) = value{1}(3) - value{1}(4); % K minus
                value{1}(8) = value{1}(3) + value{1}(4); % K plus
                value{1}(10) = P_0/(1 + .5*(gaama-1)*value{1}(5)^2)^(gaama/(gaama-1)); % P
                value{1}(9) = T_0/(1 + .5*(gaama-1)*value{1}(5)^2); % T

                % Calculating the x and y position for the points on the
                % wall
                if i == 1
                    value{1}(1) = -1/tand(.5*c_minus_values{1}(1) - .5*(c_minus_values{1}(4) + value{1}(6)));
                else
                    previous_wave = c_plus{i-1};
                    temp_keys = keys(previous_wave);
                    index = temp_keys(2);
                    previous_values = previous_wave(index);
                    number = tand(.5*(previous_values{1}(3)) - .5*(value{1}(6) + previous_values{1}(6)));
                    value{1}(1) = -previous_values{1}(2)/number + previous_values{1}(1);
                end

                % Updating the dictionary and the values generated before.
                c_plus{i}(key) = value;
                v(j) = value;

            % Checking if the point is on the wall  

            elseif ismember(key,wall)

                % Accessing the values of the previous point in the
                % dictionary. Needed to get characteristic values.
                v_temp = v(j-1);

                % Updating values in same order as before. 
                value{1}(3) = v_temp{1}(3);
                value{1}(4) = v_temp{1}(4);
                value{1}(5) = calc_mach(gaama,value{1}(4));
                value{1}(6) = asind(1/value{1}(5));
                value{1}(7) = value{1}(3) - value{1}(4);
                value{1}(8) = value{1}(3) + value{1}(4);
                value{1}(10) = P_0/(1 + .5*(gaama-1)*value{1}(5)^2)^(gaama/(gaama-1));
                value{1}(9) = T_0/(1 + .5*(gaama-1)*value{1}(5)^2);

                % Find the x and y coordinates for the wall points
                if i == 1
                    number = tand(value{1}(3) + v_temp{1}(6));
                    number2 = tand(.5*(theta_max + value{1}(3)));
                    A = [1,-number;1,-number2];
                    b = [v_temp{1}(2) - v_temp{1}(1)*number;1];
                else
                    number = tand(value{1}(3) + v_temp{1}(6));
                    previous_wave = c_plus{i-1};
                    
                    index = find(wall == key)-1;
                    index = wall(index);
                    previous_values = previous_wave(index);
                    number2 = tand(.5*(previous_values{1}(3) + v_temp{1}(3)));
                    A = [1,-number;1,-number2];
                    b = [v_temp{1}(2) - v_temp{1}(1)*number;previous_values{1}(2) - previous_values{1}(1)*number2];
                end

                % Solving the 2D matrix used to get x and y values!
                x = linsolve(A,b);
                value{1}(1) = x(2);
                value{1}(2) = x(1);

                % Updating the dictionary and the values. 
                v(j) = value;
                c_plus{i}(key) = value;
            
            % Points other than the centerline and the wall!
            
            else

                % Accessing the previous point characteristics and the left
                % turning characteristic. 
                c_minus_values = c_minus(j+i-1);
                v_temp = v(j-1);

                % Solving the system of equations to get theta and neu
                k_plus = v_temp{1}(7);
                k_minus = c_minus_values{1}(6);
                A = [1,-1;1,1];
                b = [k_plus;k_minus];
                x = linsolve(A,b);
                value{1}(3) = x(1);
                value{1}(4) = x(2);

                % Updating other values as before!
                value{1}(5) = calc_mach(gaama, value{1}(4));
                value{1}(6) = asind(1/value{1}(5));
                value{1}(7) = value{1}(3) - value{1}(4);
                value{1}(8) = value{1}(3) + value{1}(4);
                value{1}(10) = P_0/(1 + .5*(gaama-1)*value{1}(5)^2)^(gaama/(gaama-1));
                value{1}(9) = T_0/(1 + .5*(gaama-1)*value{1}(5)^2);

                % Calculating the x and y position
                number = tand((x(1) + v_temp{1}(3))*.5 + .5*(value{1}(6) + v_temp{1}(6)));
                number2 = tand((x(1) + c_minus_values{1}(1))*.5 - .5*(value{1}(6) + v_temp{1}(6)));
                A = [1,-number;1,-number2];
                b = [v_temp{1}(2) - v_temp{1}(1)*number;1];
                x = linsolve(A,b);
                value{1}(1) = x(2);
                value{1}(2) = x(1);
                c_plus{i}(key) = value;
                v(j) = value;
            end
        end
    end
end

function [c_plus, c_minus,wall,centerline] = characteristic_lines(n,gaama,theta_max)
    c_plus = {}; 
    % Order of c_plus values = [x,y,theta,neu,M,mu,K+,K-,T,P)
    j = 1;
    m = n;

    % Generating arrays to store points at the walls and the centerline. 
    centerline = [];
    wall = [];

    % Populating the cell with dictionaries which contain each wave
    % information. 
    for i = 1:n
        centerline = [centerline,j];
        dict = dictionary;
        for k = 1:m+1
            dict(j) = {NaN(1,10)};
            j = j+1;
        end
        wall = [wall,j-1];
        m = m-1;
        c_plus{i} = dict;
    end
    
    c_minus = dictionary; % All the left turning characteristics information!
    for i = 1:n
        properties = cell(1);
        properties{1}(1) = (theta_max/n)*i; % Theta
        properties{1}(2) = (theta_max/n)*i; % Neu
        properties{1}(3) = calc_mach(gaama,properties{1}(2)); % Mach number
        properties{1}(4) = asind(1/properties{1}(3)); % Wave angle
        properties{1}(5) = properties{1}(1) - properties{1}(2); % K_plus
        properties{1}(6) = properties{1}(1) + properties{1}(2); % K_minus
        c_minus(i) = properties;
    end
end

function neu = calc_neu(gaama,mach_in)
    term1 = sqrt((gaama+1)/(gaama-1));
    term2 = atand(sqrt((gaama-1)*((mach_in^2)-1)/(gaama+1)));
    term3 = atand(sqrt((mach_in^2)-1));
    neu = term1*term2-term3;
end

function M = calc_mach(gaama,neu)
    M = 0;
    tolerance = 0.0001;
    while 1
        if  neu - calc_neu(gaama,M) < tolerance
            break
        else
            M = M + 0.001;
        end
    end
end