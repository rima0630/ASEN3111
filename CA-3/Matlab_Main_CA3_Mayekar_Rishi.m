%% ASEN 3111 - Computational Assignment 3 - Main
%
% Author: Rishi Mayekar
% Date Created: October 12, 2021   |   Date Last Edited: November 1, 2021
%
% Purpose: In this CA, we first construct a function to model Prandtl's
% Lifting Line Theory with numerical summations in order to find the
% coefficient of lift, induced drag and the span efficiency factor for a
% finite wing with specific geometry. We then use this function to
% determine the lift and drag experienced by a specific form of wing with
% geometric twist and variable airfoil along the wingspan. We then
% calculate how many odd coefficients were needed to lower the error down
% to certain percents. Lastly, we apply the function to calculate how span
% efficiency varies with the root-to-tip chord ratio with different Aspect
% Ratios.

clc; clear; close all;

%% Problem 2

v_inf = 130 * 5280 / 3600; % [ft/sec]
rho_inf = 0.002377; % [lbs/ft^3]

b = 80; % Span in ft

c_t = 4; % Tip chord
c_r = 12; % Root chord

S = (c_t + c_r)/2 * b;

geo_t = 1; % Tip geometric aoa
geo_r = 6; % Root geometric aoa
N_max = 200; % Number of panels

% Tip is 0012, root is 2412
% Values taken from results of CA-2
a0_t = 0.11733 * (180/pi); % Cross sec lift slope at tip [rad^-1]
a0_r = 0.11708 * (180/pi); % Cross sec lift slope at root [rad^-1]

aero_t = 0; % 0-lift aoa at tip
aero_r = -2; % 0-lift aoa at root

% Pre-allocating e, CL and CDi vectors
e = zeros(1, N_max-1);
c_L = zeros(1, N_max-1);
c_Di = zeros(1, N_max-1);

% Doing PLLT for N number of coefficients from 2 to N_max
for N = 2:N_max
    [e(N-1), c_L(N-1), c_Di(N-1)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);
end

N = [2:N_max];

% Calculating Lift and Drag using coefficients
q = 0.5 * rho_inf * v_inf^2;
Di = c_Di * q * S; % [lbs]
L = c_L * q * S; % [lbs]

% Error Matrices
Err_e = abs(e-e(end))/e(end) * 100;
Err_CDi = abs(c_Di-c_Di(end))/c_Di(end) * 100;
Err_CL = abs(c_L-c_L(end))/c_L(end) * 100;

N_CDi_5p = 0;
N_CDi_2p = 0;
N_CDi_02p = 0;

N_CL_5p = 0;
N_CL_2p = 0;
N_CL_02p = 0;

% To find the number of panels needed for each percent.
for i = 1:N_max-1
    
    if N_CDi_5p == 0 && Err_CDi(i) < 5
        N_CDi_5p = i+1;
    end
    
    if N_CDi_2p == 0 && Err_CDi(i) < 2
        N_CDi_2p = i+1;
    end
    
    if N_CDi_02p == 0 && Err_CDi(i) < 0.2
        N_CDi_02p = i+1;
    end
    
    if N_CL_5p == 0 && Err_CL(i) < 5
        N_CL_5p = i+1;
    end
    
    if N_CL_2p == 0 && Err_CL(i) < 2
        N_CL_2p = i+1;
    end
    
    if N_CL_02p == 0 && Err_CL(i) < 0.2
        N_CL_02p = i+1;
    end
    
end

% Printing the minimum no. of coefficients for each error percent
disp("Number of odd coefficients needed to calculate C_L to... ")
disp("...5% accuracy: " + N_CL_5p);
disp("...2% accuracy: " + N_CL_2p);
disp("...0.2% accuracy: " + N_CL_02p);

disp("Number of odd coefficients needed to calculate C_Di to... ")
disp("...5% accuracy: " + N_CDi_5p);
disp("...2% accuracy: " + N_CDi_2p);
disp("...0.2% accuracy: " + N_CDi_02p);

%% Problem 3

cr = 12; % Span in ft, keeping constant

AR_vec = [4 6 8 10]; % ARs provided to plot

% Define linear vector of c_t
ct_min = 0;
ct_max = 12;
ct_n = 100;
c_t = linspace(ct_min, ct_max, ct_n);

geo_t = 6; % Tip geometric aoa
geo_r = 6; % Root geometric aoa
N = 50; % Number of panels

% Tip is 0012, root is 2412
% Values taken from results of CA-2 sent in Canvas Announcement
a0_t = 6.283185; % Cross sec lift slope at tip [rad^-1]
a0_r = 6.283185; % Cross sec lift slope at root [rad^-1]
aero_t = 0; % 0-lift aoa at tip
aero_r = -2.077239; % 0-lift aoa at root

% Pre-allocating e, CL and CDi vectors
e = zeros(1, ct_n);
c_L = zeros(1, ct_n);
c_Di = zeros(1, ct_n);

figure;
hold on;

for AR = AR_vec
    
    % Iterating through the tip chord array
    for i = 1:ct_n
        
        ct_curr = c_t(i);
        
        % Calculate Corresponding Wingspan
        % Area (S) = (cr + ct) * b/2
        % AR = b^2 / S
        % AR = b / ((cr + ct) * 0.5)
        b_curr = AR * 0.5 * (cr + ct_curr);
        
        % Function call for each element
        [e(i), c_L(i), c_Di(i)] = PLLT(b_curr,a0_t,a0_r,ct_curr,cr,aero_t,aero_r,geo_t,geo_r,N);
        
    end
    
    % Plotting
    plot(c_t./c_r, e);
    
end

title("Span efficienfy vs. c_t/c_r for various Aspect Ratios");
xlabel("c_t/c_r")
ylabel("Span efficiency factor")
legend("AR = 4", "AR = 6", "AR = 8", "AR = 10")
