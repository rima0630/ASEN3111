%% ASEN 3111 - Computational Assignment 1 - Main
%
% Author: Rishi Mayekar
% Date Created: August 31, 2021   |   Date Last Edited: September 13, 2021
% Purpose: To understand the use of integration methods, which can be used
% to find sectional coefficients of lift and drag, and use the results to 
% compare the two methods. 

%% Housekeeping

clc; clear all; close all;

%% ------------------ PROBLEM #1 ------------------

fprintf("PROBLEM 1 \n");

    %% Analytically determining sectional coefficients of lift and drag
    syms x
    expr_Cp = -4*(sin(x)^2 +sin(x));
    expr_cl = -0.5 * expr_Cp * sin(x);
    expr_cd = -0.5 * expr_Cp * cos(x);

    c_l = int(expr_cl, 0, 2*pi);
    c_d = int(expr_cd, 0, 2*pi);
    
    panels = 1000;

    %% Trapezoidal rule
    
    % Initializing vectors
    trap_cl = zeros(1,panels);
    trap_cd = zeros(1,panels);
    
    fprintf("Computing Trapezoidal Rule...\n");

    for N = [1:panels]
        
        % Vector of N divisions from 0 to 2*pi
        theta = linspace(0, 2*pi, N);
        
        % "f(x)" for cd and cl both
        values_cl = 2 * (sin(theta).^3 + sin(theta).^2);
        values_cd = 2 * (sin(theta).^2 + sin(theta)) .* cos(theta);
        
        % temp values needed for calculations
        temp_cl = 0;
        temp_cd = 0;
        step = 2*pi/N;

        for i = 1:N-1
            % calculating the i'th panel
            incr_cl = 0.5 * step * (values_cl(i) + values_cl(i+1));
            incr_cd = 0.5 * step * (values_cd(i) + values_cd(i+1));
            
            % Adding the calculated panel to the answer
            temp_cl = temp_cl + incr_cl;
            temp_cd = temp_cd + incr_cd;
        end
        
        % updating the vector with all comp trap values for each N
        trap_cl(N) = temp_cl;
        trap_cd(N) = temp_cd;

    end
    
    fprintf("Trapezoidal Rule Complete! \n\n")
    
    % Ploting
    plot([1:panels], trap_cl, 'r');
    hold on;
    plot([1:panels], trap_cd, 'b');
    title('Effect of number of panels on c_d and c_l with Trapezoidal Rule');
    xlabel('Number of panels (N)');
    ylabel('Coefficient');
    legend('Coefficient of lift', 'Coefficient of drag')
    
    
    %% Simpson's rule
    
    % Initializing vectors
    simp_cd = zeros(1,panels);
    simp_cl = zeros(1,panels);
    
    fprintf("Computing Simpson's Rule...\n");
    
    for N = [1:panels]
        
        % Vector of N divisions from 0 to 2*pi
        theta = linspace(0, 2*pi, N);
        
        % "f(x)" for cd and cl both
        values_cl = 2 * (sin(theta).^3 + sin(theta).^2);
        values_cd = 2 * (sin(theta).^2 + sin(theta)) .* cos(theta);
        
        % temp values needed for calculations
        temp_cl = 0;
        temp_cd = 0;
        step = pi/N;

        for i = 2:N-1
            % calculating the i'th panel
            incr_cl = values_cl(i-1) + 4*values_cl(i) + values_cl(i+1);
            incr_cd = values_cd(i-1) + 4*values_cd(i) + values_cd(i+1);
            
            % Adding the calculated panel to the answer
            temp_cl = temp_cl + (step/3)*incr_cl;
            temp_cd = temp_cd + (step/3)*incr_cd;
        end
        
        % updating the vector with all comp trap values for each N
        simp_cl(N) = temp_cl;
        simp_cd(N) = temp_cd;

    end
    
    fprintf("Simpson's Rule Complete! \n\n")
    
    % Ploting
    figure;
    plot([1:panels], simp_cl, 'r');
    hold on;
    plot([1:panels], simp_cd, 'b');
    title("Effect of number of panels on c_d and c_l with Simpson's Rule");
    xlabel('Number of panels (N)');
    ylabel('Coefficient');
    legend('Coefficient of lift', 'Coefficient of drag')
    
    
    %% Errors for both methods
    
    threshold = 2*pi - (2*pi * (0.2/100));
    simp_ind = find(simp_cl >= threshold, 1 );
    trap_ind = find(trap_cl >= threshold, 1 );
    
    fprintf("PROBLEM 1 RESULTS: \n");
    fprintf("It takes %d panels for the Complex Trapezoidal rule and %d panels for Simpson's Rule \nto achieve 2/10 percent error in the predicted sectional lift.\n\n\n", trap_ind, simp_ind);
    
    

%% ------------------ PROBLEM #2 ------------------

fprintf("PROBLEM 2 \n");

    load Data_CA1_Mayekar_Rishi;

    xc_low = Cp_lower.knots(7:107);
    xc_up = Cp_upper.knots(7:107);

    %% constants

    c = 1.5; % m
    aoa = 9; % degree
    v_inf = 40; % m/s
    rho_inf = 1.225; % kg/m^3
    P_inf = 101300; % Pa
    q_inf = 0.5 * rho_inf * v_inf^2; % kg/ms^2
    t = 12/100; % Specific to NACA 0012

    %% Trapezoidal
    
    % Converting from percent to x values
    xc_low = xc_low * c;
    xc_up = xc_up * c;
    
    % Setting max panels to test, and initializing D and L vectors
    N_max = 1000;
    D = zeros(1, N_max);
    L = zeros(1, N_max);
    
    fprintf("Computing Drag and Lift using Trapezoidal rule... \n");

    for i = [1:N_max]

        N = i;
        
        % Setting up the x values
        x_trap = linspace(0,c,N+1);
        x_c_trap = x_trap / c;
        delta_x = diff(x_trap);
        
        % Step size and temp variables
        step = c/N;
        F_incr = 0;
        F_temp = 0;
        
        % Initializing Normal and Axial variables
        N_low = 0;
        N_up = 0;
        A_low = 0;
        A_up = 0;
        
        % Airfoil shape equation with x_c vector plugged in
        y_t = t/0.2 * c * (0.2969*sqrt(x_c_trap) - 0.1260*x_c_trap ...
              - 0.3516*x_c_trap.^2 + 0.2843*x_c_trap.^3 - 0.1036*x_c_trap.^4);
        delta_y = diff(y_t);
        
        % Extracting Cp values for each x/c value
        Cp_up = fnval(Cp_upper, x_c_trap);
        Cp_low = fnval(Cp_lower, x_c_trap);
        
        % Finding Pressure using C_p
        P_up = (Cp_up * q_inf) + P_inf;
        P_low = (Cp_low * q_inf) + P_inf;

        for j = 1:N-1
            
            % Delta values for the panel
            dx = delta_x(j);
            dy = delta_y(j);
            
            % The increments for the upper and lower surface trapezoids
            trap_up = 0.5 * (P_up(j) + P_up(j+1));
            trap_low = 0.5 * (P_low(j) + P_low(j+1));
            
            % Adding the increments to the Normal values with dx
            N_up = N_up + trap_up * dx;
            N_low = N_low + trap_low * dx;
            
            % Adding the increments to the Axial values with dy
            A_up = A_up + trap_up * dy;
            A_low = A_low + trap_low * dy;

        end

        % Calculating N' and A'
        N_pr = N_low - N_up;
        A_pr = A_low + A_up;
        
        % Finding the new Lift and Drag, adding them to the L and D vectors
        L(i) = N_pr*cosd(aoa) - A_pr*sind(aoa);
        D(i) = N_pr*sind(aoa) + A_pr*cosd(aoa);


    end
    
    fprintf("Lift and Drag successfully computed \n\n");
    
    % Lift and Drag plot
    figure;
    N_vec = linspace(1,N_max,N_max);
    plot(N_vec, L);
    hold on;
    plot (N_vec, D);
    legend("Lift", "Drag");
    xlabel("Number of Panels")
    ylabel("Lift and Drag per unit span (N/m)")
    title("Lift and Drag estimates per unit span using different number of panels")


    %% Error
    
    % Using the last values to find the error thresholds
    E_1 = (1 - 5/100) * L(end);
    E_2 = (1 - 1/100) * L(end);
    E_3 = (1 - 0.2/100) * L(end);
    
    % Looking for what N value the Lift crosses each threshold
    N_1 = find(L >= E_1, 1);
    N_2 = find(L >= E_2, 1);
    N_3 = find(L >= E_3, 1);


    %% Printing the results
    
    fprintf(" PROBLEM 2 RESULTS: \n");
    fprintf("The lift per unit span of the airfoil: %2d N/m \n", L(end));
    fprintf("The drag per unit span of the airfoil: %2d N/m \n \n", D(end));
    fprintf("The Trapezoidal rule was tested for up to %2d panels... \n", N_max);
    fprintf("Number of panels needed to achieve 5%% relative error: %2d \n", N_1);
    fprintf("Number of panels needed to achieve 1%% relative error: %2d \n", N_2);
    fprintf("Number of panels needed to achieve 0.2%% relative error: %2d \n", N_3);
    
