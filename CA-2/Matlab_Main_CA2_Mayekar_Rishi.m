%% ASEN 3111 - Computational Assignment 2 - Main
%
% Author: Rishi Mayekar
% Contributors: Elijah Stubbs, Ketan Kamut
% Date Created: September 14, 2021   |   Date Last Edited: October 12, 2021
%
% Purpose: In this CA, we first explore the flat plate model and plotting
% the streamlines, potential and pressure contours for one. We then look
% into NACA airfoils, starting with the classic symmetric NACA 0012
% airfoil. We find a nominal number of panels to use before using the
% vortex-panel method to estimate the contribution of vortex flow by adding
% up each vortext sheet's contribution at any point in space to find the
% velocity induced at the point. We use this to find sectional lift
% coefficient vs angle of attack slope, and then do the same with cambered
% NACA airfoils too, before eventually looking into the effect of flaps!

%%
clc; clear all; close all

%% ------------------ PROBLEM #1 ------------------

disp("Working on Problem 1...");

% Thin Airfoil Theory
c = 1.5;
alpha = 6;
V_inf = 30;
p_inf = 1.013e5;
rho_inf = 0.225;
N = 3;

% Function call
Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N);

% For lower N-Values, the streamlines and the equipotential lines seemed to
% be a little sharp especially at the regions where the lines curved the
% most. The pressure contours showed this effect even more so, and it was
% actually more evident to see the vortext created by each panel when the 
% number of panels was less than 50.

disp("Problem 1 complete!");
disp(" ");

%% ------------------ PROBLEM #2 ------------------

disp("Commencing Problem 2...");

% Parameters for the Vortex Panel function
C = 1.5;
t = 12;
ALPHA = 0;
VINF = 30;
FLAG = 0;
N_MAX = 200;
N_MIN = 20;

disp("Calculating nominal panel count");

% Pre-allocating error matrix
Error = zeros(1,N_MAX - N_MIN + 1);

for N = N_MIN:N_MAX
    
    % Obtaining Boundary points for NACA 0012
    [XB, YB] = NACA_Airfoils(0, 0, t, C, N);

    % Function call
    [CL,CP,XC] = Vortex_Panel(XB,YB,VINF,ALPHA,FLAG);
    
    % We are comparing the mean CP to 0
    Error(N - N_MIN + 1) = mean(abs(CP));
end

% Plotting the error
figure;
plot([N_MIN:N_MAX], Error);
[E, I] = min(Error);
xlabel("Number of panels (N)");
ylabel("Mean coefficient of pressure");
title("Problem 2: Mean C_p vs Number of Panels");

% It seems like the CL hovers around 0 after a certain number of panels
N_OPT = (I-1) + N_MIN; 
disp("Optimal number of panels needed to minimize the error is: " + N_OPT);

% Using the optimal number of panels to construct another NACA 0012
[XB, YB] = NACA_Airfoils(0, 0, t, C, N_OPT);

% Vectors that shall be plotted
CLvec = zeros(1,4);
ALPHA = [-6, 0, 6, 9];

disp("Plotting C_L vs Alpha for NACA 0012");

% Setting up for the cp vs xc plots
FLAG = 1;

for i = 1:4
    
    % Plotting c[-xc for each of the four alphas
    [CL,CP,XC] = Vortex_Panel(XB,YB,VINF,ALPHA(i),FLAG);
    
    % Unique title for each plot
    title("Problem 2: C_p over x/c for " + ALPHA(i) + " degree AOA");
    
    % Saving the CL for each in the vector
    CLvec(i) = CL;
end

% Plotting the Coefficient of lift vs angle of attack for the 4 angles
figure;
plot(ALPHA, CLvec);
grid on;
xlabel('Angle of Attack (Degrees)');
ylabel('Coefficient of Lift');
title('Problem 2: Coefficient of lift vs Angle of Attack for NACA 0012');

disp("Problem 2 Complete!");
disp(" ");

%% ------------------ PROBLEM #3 ------------------

disp("Working on Problem 3...");

% The four airfoils in vector format:
NACA = [0 0 12 ; 2 4 12 ; 4 4 12 ; 2 4 24];
C = 1.5;
N = 200;

% Setting which angles of attack to plot
A_MIN = -10;
A_MAX = 20;

% Pre-allocating, as always
CLVEC = zeros(1,A_MAX - A_MIN + 1);

figure;
grid on;
hold on;

disp("Plotting CL vs Alpha for the 4 NACA Airfoils");

for i = 1:4
    
    % Obtaining the boundary points of
    [XB, YB] = NACA_Airfoils(NACA(i,1), NACA(i,2), NACA(i,3), C, N);
    
    FLAG = 0;
    VINF = 30;
    
    % Alpha vector
    ALPHA = [A_MIN:A_MAX];
    
    for A_CURR = ALPHA
        
        % Calling function and saving each CL into the vector
        [CL,CP,XC] = Vortex_Panel(XB,YB,VINF,A_CURR,FLAG);
        CLVEC(A_CURR - A_MIN + 1) = CL;
    end
    
    % The Lift vs Alpha plot
    plot([A_MIN:A_MAX], CLVEC);
    
    % Finding the closest angle of attack of C_L = 0
    [~, I_CL0] = min(abs(CLVEC));
    ALPHA_CL0 = ALPHA(I_CL0);
    
    % Lift curve slope
    dCl_da = (CLVEC(end) - CLVEC(1)) / (A_MAX - A_MIN);
    
    % The results for each airfoil
    disp("Alpha for CL=0 for NACA " + NACA(i,1) + NACA(i,2) + NACA(i,3) + " Airfoil: " + ALPHA_CL0 );
    disp("Lift-curve slope for NACA " + NACA(i,1) + NACA(i,2) + NACA(i,3) + " Airfoil: " + dCl_da );
    
end

title("Problem 3: NACA Airfoils CL");
xlabel("Angle of Attack (Degrees)");
ylabel("Coefficient of Lift");
legend("0012", "2412", "4412", "2424");

disp("Problem 3 complete!")
disp(" ");

%% ---------------- BONUS PROBLEM ------------------

disp("Working on bonus problem...");

% Constants
C = 1.5;
ALPHA = 6;
N = 200;
DEL = 20;
x_rot = 0.8*C;

% Creating the two airfoils
[xb1, yb1] = NACA_Airfoils(0, 0, 12, C*0.75, N);
[xb2, yb2] = NACA_Airfoils(0, 0, 12, C*0.2, N);
xb2 = xb2 + 0.8*C;

% Creating NACA 0012 with no flap
[XB, YB] = NACA_Airfoils(0, 0, 12, C, N);

figure;
FLAG = 0;
VINF = 30;

% Pre-allocation
CLTVEC = zeros(1,A_MAX - A_MIN + 1);
CL1VEC = zeros(1,A_MAX - A_MIN + 1);

disp("Plotting CL vs Alpha for with and w/o flap")

for ALPHA = [A_MIN:A_MAX]

    % Calling function and saving CL at ALPHA for main airfoil
    [CL1,~,~] = Vortex_Panel(xb1,yb1,VINF,ALPHA,FLAG);
    
    % Doing the same for flap, but adjusting ALPHA by adding DEL
    [CL2,~,~] = Vortex_Panel(xb2,yb2,VINF,ALPHA+DEL,FLAG);
    
    % Original no-flap airfoil
    [CL,~,~] = Vortex_Panel(XB,YB,VINF,ALPHA,FLAG);
    
    % Saving the CL values
    CLTVEC(ALPHA - A_MIN + 1) = CL1 + CL2;
    CLVEC(ALPHA - A_MIN + 1) = CL;
end

% Plotting
plot([A_MIN:A_MAX], CLTVEC);
hold on;
plot([A_MIN:A_MAX], CLVEC);
legend("With flap", "Without flap");
grid on;
title("Bonus: C_L vs Angle of Attack for NACA 0012 with and without flap");
xlabel("Angle of Attack (degrees)");
ylabel("Coefficient of Lift");

disp("Computational Assignment 2 complete!")
