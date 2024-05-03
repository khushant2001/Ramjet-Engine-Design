function [x_shell,y_shell] = outer_shell(x_cowl,y_cowl,x_wall,y_wall,x_diffuser,A_diffuser)
    
    disp("Designing the Case of the Engine...")
    theta = zeros(1,length(x_wall)/2-1);
    for i = 1:length(theta)
        theta(i) = atand((y_wall(i+1)-y_wall(i))/(x_wall(i+1)-x_wall(i)));
    end
    theta_req = atand((A_diffuser(end/2)-y_cowl(1))/(x_diffuser(end/2)-x_cowl(1)));
    
    % Calculate absolute differences
    abs_diff = abs(theta - theta_req);
    % Find index of minimum difference
    [minindex, index] = min(abs_diff);
    disp(minindex)
    index = index+1;
    x_shell = [x_cowl(1),x_diffuser(end/2),x_wall(index:end/2)];
    y_shell = [y_cowl(1),A_diffuser(end/2),y_wall(index:end/2)];
end