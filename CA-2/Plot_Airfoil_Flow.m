function Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N)
    
    
    %% Define Domain
    
    b = c + 0.5;
    xmin = -b;
    xmax = b;
    ymin = -b;
    ymax = b;

    %% Define Number of Grid Points

    nx = 100; % steps in the x direction
    ny = 100; % steps in the y direction

    %% Create mesh over domain using number of grid points specified
    [x,y] = meshgrid(linspace(xmin, xmax, nx), linspace(ymin, ymax, ny));

    %% Making vectors for the x-coordinates
    x_vec = linspace(0, c, N+1);
    x_c_vec = x_vec/c;
    
    % Step size between consecutive points
    delta_x = c/N;
    
    
    %% Vortex sheet strength
    
    % Converting alpha to radians
    alpha = alpha * pi/180;
    
    % Finding the control points - the center of each panel
    ctrl_xc = x_c_vec(2:end) - delta_x/2;
    x_ctrl = ctrl_xc*c;
    
    % Finding gamma at each control point
    gam_x = alpha * V_inf * sqrt((1 - ctrl_xc) ./ ctrl_xc);
    
    
    %% Numerically calculating the vortex sheet
    
    % Pre-allocating for vector
    Gamma = zeros(1,N);
    
    for i = [1:N]
        
        % Finding the strength Gamma of Vortex i
        Gamma(i) = gam_x(i) * delta_x;
    end
    
    
    %% Stream function for uniform flow
    
    
    Psi_uni = V_inf * (x*sin(-alpha) + y*cos(-alpha));
    
    %% Stream function as affected by each Vortex
    
    % Pre-allocation
    Psi_Gamma = zeros(nx, ny);
    
    for i = 1:nx
        for j = 1:ny
            
            % This will be the value of phi at (i,j)
            psi_ij = 0;
            
            for k = 1:length(Gamma)
                % Current x coordinate
                x_curr = x_ctrl(k);
                
                % Distance from panel k to (x,y)
                radius = sqrt( (y(i,j))^2 + (x_curr-x(i,j))^2 );
                
                % psi increment
                psi_ij = psi_ij + Gamma(k)/(2*pi)*log(radius);
                
            end
            
            % Adding it to the array
            Psi_Gamma(i,j) = psi_ij;
            
        end
    end
    
    %% Superimpose stream functions and plot
    
    % Summing them up
    StreamFunction = Psi_uni + Psi_Gamma;
    
    % Determine color levels for contours
    levmin = StreamFunction(1,nx); 
    levmax = StreamFunction(ny,nx/2);
    levels = linspace(levmin,levmax,50)';

    % Plot streamfunction at levels
    contour(x,y,StreamFunction,levels, 'LineWidth', 1)
    
    % Plotting airfoil for reference
    hold on;
    plot(x_vec, zeros(1,N+1), 'k', 'LineWidth', 1);
    
    xlabel("x");
    ylabel("y");
    title("Problem 1: Streamlines");
    
    %% Plotting the equipotential lines

    % Using the same contours for stream function
    levmin = StreamFunction(1, nx); 
    levmax = StreamFunction(ny, nx/2);
    levels = linspace(levmin, levmax, 50)';

    % Relating Stream Function to Velocity Potential
    % d_phi/dx = v = -d_psi/d_y 
    % d_phi/dy = u = d_psi/dx
    
    % Plotting: Since stream function is related to velocity potential as
    % shown in the above comments, we just switch the contour coordinates
    % and negate the x coordinates
    figure;
    contour(y, -x, StreamFunction, levels, 'LineWidth', 1)

    % Plotting airfoil for reference
    hold on;
    plot(x_vec, zeros(1,N+1), 'k', 'LineWidth', 1);
    
    xlabel("x");
    ylabel("y");
    title("Problem 1: Velocity Potential");
    
    % Setting the limits for the plot
    xlim([xmin, xmax]);
    ylim([ymin,ymax]);
    
    %% Pressure contours
    
    % Pre-allocation
    V_Psi = zeros(20);
    
    for i = 1:nx
        for j = 1:ny
            
            % Variables that will be used
            Vx_ij = 0;
            Vy_ij = 0;
            x_ij = x(i,j);
            y_ij = y(i,j);
            
            
            for k = 1:length(Gamma)
                
                % Current control point x-coordinate
                x_curr = x_ctrl(k);
                
                % Distance from panel k to (x,y)
                radius = sqrt( (y_ij)^2 + (x_curr-x_ij)^2 );
                
                % Angle between chord line and radius
                theta_k = atan((y_ij) / (x_ij-x_curr));
                if x_ij-x_curr < 0
                    theta_k = theta_k + pi/2;
                end
                
                % Contribution of panel k to each velocity component
                Vx_k = Gamma(k) / (2*pi*radius) * cos(theta_k - pi/2);
                Vy_k = Gamma(k) / (2*pi*radius) * sin(theta_k - pi/2);
                
                % Incrementing each velocity component with the new contribution
                Vx_ij = Vx_ij + Vx_k;
                Vy_ij = Vy_ij + Vy_k;
                
                
            end
            
            % Adding uniform flow contributions
            Vx_ij = Vx_ij + V_inf * cos(alpha);
            Vy_ij = Vy_ij + V_inf * sin(alpha);
            
            % The final speed as obtained from the Stream function
            V_Psi(i,j) = sqrt(Vx_ij^2 + Vy_ij^2);
            
        end
    end
    
    % C_p from bernoulli's equation
    C_p = 1 - (V_Psi/V_inf).^2;
    
    % Not optimal color-allocation but enough to see the important contours
    levmin = min(C_p, [], 'all');
    levmax = max(C_p, [], 'all');
    levels = linspace(levmin, levmax, 200)';

    % Plotting the Pressure contours
    figure;
    contour(x,y,C_p,levels, 'LineWidth', 1);
    
    hold on;
    
    % Airfoil for reference
    plot(x_vec, zeros(1,N+1), 'k', 'LineWidth', 1);
    
    hold off;
    xlabel("x");
    ylabel("y");
    title("Problem 1: Pressure contours");
    
end

