function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
    
    %% Constants
    c_avg = (c_t + c_r)/2;
    S = b * c_avg;
    AR = b^2/S;

    %% Structs for the win
    TIP.a0 = a0_t; % cross-sec lift slope at tip (a_0)
    TIP.c = c_t; % chord length at tip
    TIP.aero = aero_t*(pi/180); % zero lift aoa at tip (alpha_L=0)
    TIP.geo = geo_t*(pi/180); % geometric aoa at tip (alpha)
    
    ROOT.a0 = a0_r; % cross-sec lift slope at tip
    ROOT.c = c_r; % chord length at root
    ROOT.aero = aero_r*(pi/180); % zero lift aoa at root
    ROOT.geo = geo_r*(pi/180); % geometric aoa at root
    
    %% The loop
    
    % The matrices
    A_mat = zeros(N);
    b_mat = zeros(N,1);
    
    for i = 1:N
        
        % Theta and y at panel i
        theta_i = i*pi/(2*N);
        y_theta = b/2 * cos(theta_i);
        
        % lc-slope, Chord, aero and geo aoa at panel i
        a0_y = ROOT.a0 - (ROOT.a0 - TIP.a0)*(y_theta/(b/2));
        c_y = ROOT.c - (ROOT.c - TIP.c)*(y_theta/(b/2));
        aero_y = ROOT.aero - (ROOT.aero - TIP.aero)*(y_theta/(b/2));
        geo_y = ROOT.geo - (ROOT.geo - TIP.geo)*(y_theta/(b/2));
        
        % Updating the b matrix
        b_mat(i) = geo_y - aero_y;
        
        for j = 1:N
            
            % n is odd
            n = 2*j - 1;
            
            % Updating the A matrix
            A_mat(i,j) = 4*b*sin(n*theta_i)/(a0_y*c_y) + n*sin(n*theta_i)/sin(theta_i);
            
        end
        
    end
    
    x_mat = A_mat^(-1) * b_mat; % x_mat has all the A_(2n-1) coefficients
    
    %% Calculating e,c_L,c_Di
    
    % c_L
    c_L = x_mat(1) * pi * AR;
    
    % delta
    delta = 0;
    for j = 2:N
        coeff = 2*j-1;
        temp = (coeff) * (x_mat(j)/x_mat(1))^2;
        delta = delta + temp;
    end
    
    % e
    e = 1 / (1+delta);
    
    % c_Di
    c_Di = c_L^2 / (pi*e*AR);
    
end

