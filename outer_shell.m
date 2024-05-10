function [x_shell,y_shell,x_shell_b,y_shell_b,slope,slope_b,distance,distance_b] = outer_shell(x_cowl,y_cowl,x_wall,y_wall,x_diffuser,A_diffuser,x,y)
    
    disp("Designing the Case of the Engine...")
    x_shell = [x_cowl(1),x_wall(1),x_wall(end/2)];
    y_shell = [y_cowl(1),y_wall(end/2),y_wall(end/2)];

    slope = atand((y_shell(end)-y_shell(1))/(x_shell(end)-x_shell(1)));
    distance = sqrt((x_shell(2)-x_shell(1))^2 + (y_shell(2)-y_shell(1))^2);
    diffuser_estimate = y_cowl(1) + slope*x_diffuser(end/2);
    if diffuser_estimate < A_diffuser(end/2)
        disp("... ERROR NEED MORE SPACE")
        %return
    end
    x_shell_b = [x(1),x_wall(round(end/8)),x_wall(end/2)];
    y_shell_b = [y(1),y_wall(end),y_wall(end)];
    distance_b = sqrt((x_shell_b(2)-x_shell_b(1))^2 + (y_shell_b(2)-y_shell_b(1))^2);
    slope_b = atand((y_shell_b(end)-y_shell_b(1))/(x_shell_b(end)-x_shell_b(1)));
end