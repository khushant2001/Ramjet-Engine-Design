
function [mach_out, T_out, P_out, Area_out] = combustor(T_in,P_in,T_0_in,M_in,fuel_mass)

end

function [P2,T2,mach_out] = rayleigh(mach_in,T1,P1,gaama,q)
    
    R = 287;
    % Finding the stagnation temperature and pressure
    T_01 = T1*(1 + (gaama-1)*.5*mach_in^2);
    P_01 = P1*(1+(gaama-1)*.5*mach_in^2)^(gaama/(gaama-1));

    % Find Cp
    cp = gaama*R/(gaama-1);

    T_02 = (q/cp) + T_01;

    % Finding choked properties
    T01_T0star = (gaama+1)*mach_in^2*(2 + (gaama-1)*mach_in^2)/(1 + gaama*mach_in^2)^2;
    T0_star = T_01/T01_T0star;
    P01_P0star = ((1+gaama)/(1+gaama*mach_in^2))*((2+mach_in^2*(gaama-1))/(gaama+1))^(gaama/(gaama-1));
    P0_star = P_01/P01_P0star;

    % Setting up the iterative determination of mach_out
    mach_out = 0;
    tolerance = 0.0001;
    exp1 = T_02/T0_star;
    while 1
        exp2 = (gaama+1)*mach_out^2*(2 + (gaama-1)*mach_out^2)/(1 + gaama*mach_out^2)^2;
        if exp1-exp2 < tolerance || mach_out > 1
            if mach_out > 1
                mach_out = 0;
            end
            break
        else
            mach_out = mach_out + 0.0001;
        end
    end
    
    P_02 = P0_star*((1+gaama)/(1+gaama*mach_out^2))*((2+mach_out^2*(gaama-1))/(gaama+1))^(gaama/(gaama-1));
    T2 = T_02 / (1 + (gaama-1)*.5*mach_out^2);
    P2 = P_02 / (1 + (gaama-1)  *.5*mach_out^2)^(gaama/(gaama-1));
end