function [theta] = BTMeq_FindT(beta,M1)
    
    % Function used to solve for theta given beta and M1
    
    a = (M1^2 * sin(beta)^2) - 1;
    b = M1^2 * (1.4 + cos(2*beta)) + 2;
    
    tantheta = 2 * cot(beta) * a/b;
    
    theta = atan(tantheta);


end

