%% ASEN 3111 - Computational Assignment 3 - Main
%
% Author: Rishi Mayekar
% Date Created: November 2, 2021   |   Date Last Edited: November 29, 2021
%
% Purpose: In this CA, we first use matlab functions to verify the numbers
% referenced from the Appendix tables and B-T-M diagram in the textbook
% examples. This gives good practice for using the function. Problem 2
% presents us with a diamond airfoil for which we must calculate the lift
% and wave drag coefficients using the pressure ratios and shock/expansion
% fan designation with the help of the functions given. Lastly, we then
% produce a plot of the relation between beta and theta for various mach
% numbers when it comes to oblique shocks.

clc; clear; close all;

%% Problem 1 

fprintf("<strong>----------PROBLEM #1----------</strong>\n\n")

% EX 8.8
    % Givens
    M = 3.5;
    P_stat = 0.3; % atm
    T_stat = 180; % K

    % Find Properties at given M
    [p0op, t0ot, ~] = isentropic(M);
    
    % Printing
    fprintf("<strong>Example 8.8</strong>:\n")
    fprintf('P0/P = %2.3f\nT0/T = %2.3f\n\n', p0op, t0ot);

% EX 8.9
    [p0op, ~, ~] = isentropic(0.6);
    M1 = isentropicFindM(1.4, 1.69);
    
    % Printing
    fprintf("<strong>Example 8.9</strong>:\n")
    fprintf('P0/P = %2.3f\nM = %2.3f\n\n', p0op, M1);
    
% EX 8.11
    M1 = 2;
    [M4n, p2op1, ~, t2ot1, ~, ~] = shock_calc(M1);
    
    % Printing
    fprintf("<strong>Example 8.11</strong>:\n")
    fprintf('P2/P1 = %2.3f\nT2/T1 = %2.3f\nM2 = %2.4f\n\n', p2op1, t2ot1, M4n);

% EX 8.12
    M1 = 2;
    [p01op1, ~, ~] = isentropic(M1);
    [~, ~, ~, ~, ~, p02op01] = shock_calc(M1);
    
    % Printing
    fprintf("<strong>Example 8.12 (a)</strong>:\n")
    fprintf('P01/P1 = %1.3f\nP02/P01 = %1.4f\n\n', p01op1, p02op01);
    
    M1 = 4;
    [p01op1, ~, ~] = isentropic(M1);
    [~, ~, ~, ~, ~, p02op01] = shock_calc(M1);
    
    % Printing
    fprintf("<strong>Example 8.12 (b)</strong>:\n")
    fprintf('P01/P1 = %1.3f\nP02/P01 = %1.4f\n\n', p01op1, p02op01);

% EX 8.13
    Minf = 2;
    [p0infopinf, t0infotinf, ~] =isentropic(Minf);
    [~, ~, ~, ~, ~, p01op0inf] = shock_calc(Minf);
    
    % Printing
    fprintf("<strong>Example 8.13</strong>:\n")
    fprintf('P0inf/Pinf = %1.3f\nT0inf/Tinf = %1.1f\nP01/P0inf = %1.4f\n'...
        , p0infopinf, t0infotinf, p01op0inf);
    
    M2 = 0.2;
    [p02op2, t02ot2, ~] =isentropic(M2);
    fprintf('P02/P2 = %1.3f\nT02/T2 = %1.3f\n\n', p02op2, t02ot2);
    
% EX 8.14
    M = 10;
    [p0infopinf, t0infotinf, ~] =isentropic(M);
    [~, ~, ~, ~, ~, p01op0inf] = shock_calc(M);
    
    % Printing
    fprintf("<strong>Example 8.14</strong>:\n")
    fprintf('P0inf/Pinf = %5.0f\nT0inf/Tinf = %1.0f\nP01/P0inf = %.3d\n'...
        , p0infopinf, t0infotinf, p01op0inf);
    
    M2 = 0.2;
    [p02op2, t02ot2, ~] =isentropic(M2);
    fprintf('P02/P2 = %1.3f\nT02/T2 = %1.3f\n\n', p02op2, t02ot2);
    
% EX 8.15
    M1 = 2;
    [M4n,p2op1,rho2orho1,t2ot1,~,~] = shock_calc(M1);
    
    % Printing
    fprintf("<strong>Example 8.15</strong>:\n")
    fprintf('M1 = %d\nP2/P1 = %1.1f\nM2 = %1.4f\nRho2/Rho1 = %1.3f\nT2/T1 = %1.3f\n\n'...
        , M1, p2op1, M4n, rho2orho1, t2ot1);
    
% EX 8.18
    M1 = 3.5;
    [M4n,~,~,t2ot1,~,~] = shock_calc(M1);
    
    % Printing
    fprintf("<strong>Example 8.18</strong>:\n")
    fprintf('M2 = %1.4f\nT2/T1 = %1.3f\n\n', M4n, t2ot1);
    
% EX 8.21
    M1 = 3.53;
    [M4n,~,~,~,~,~] = shock_calc(M1);
    
     % Printing
    fprintf("<strong>Example 8.21</strong>:\n")
    fprintf('M2 = %1.2f\n\n', M4n);
    
% EX 8.23
    M1 = 8;
    [M4n, p2op1, ~, ~, ~, ~] = shock_calc(M1);
    [p02op2, ~, ~] =isentropic(M4n);
    
    p02op1 = p02op2 * p2op1;
    
     % Printing
    fprintf("<strong>Example 8.23</strong>:\n")
    fprintf('P02/P1 = %1.2f\n\n', p02op1);
    
% EX 9.2
    Mn1 = 1.6;
    [Mn2,p2op1,~,t2ot1,~,p02op01] = shock_calc(Mn1);
    
     % Printing
    fprintf("<strong>Example 9.2</strong>:\n")
    fprintf('M_n,2 = %1.4f\np2/p1 = %1.2f\nT2/T1 = %1.3f\nP02/P01 = %1.4f\n',... 
        Mn2, p2op1, t2ot1, p02op01);
    
    M1 = 2;
    [p01op1, t01ot1, ~] = isentropic(M1);
    fprintf('P01/P1 = %1.4f\nT01/T1 = %1.1f\n\n', p01op1, t01ot1);
    
% EX 9.3
    Mn1 = 1.2;
    [Mn2, p2op1, ~, t2ot1, ~, ~] = shock_calc(Mn1);
    
     % Printing
    fprintf("<strong>Example 9.3</strong>:\n")
    fprintf('M_n,2 = %1.4f\np2/p1 = %1.3f\nT2/T1 = %1.3f\n\n',... 
        Mn2, p2op1, t2ot1);
    
% EX 9.4
    p2op1 = 3;
    [Mn1] = shock_calcGetM1(p2op1);
    
    % Printing
    fprintf("<strong>Example 9.4</strong>:\n")
    fprintf('M_n,1 = %1.2f\n\n', Mn1);
    
% EX 9.5
    M = 3;
    [~,~,~,~,~,p02op01] = shock_calc(M);
    
    % Printing
    fprintf("<strong>Example 9.5</strong>:\n")
    fprintf('Case 1:\nP02/P01 = %1.4f\n', p02op01);
    
    Mn1 = 1.93;
    [Mn2,~,~,~,~,p02op01] = shock_calc(Mn1);
    fprintf('Case 2:\nP02/P01 = %1.4f\nM_n,2 = %1.3f\n', p02op01, Mn2);
    
    theta = 22 * pi/180;
    M1 = 3;
    beta = BTMeq(theta, M1)*180/pi;
    fprintf('beta = %2.0f\n', beta);
    
    M2 = 1.9;
    [~,~,~,~,~,p03op02] = shock_calc(M2);
    fprintf('P03/P02 = %1.4f\n\n', p03op02);
    
% EX 9.6
    M1 = 5;
    theta = 15 * pi/180;
    beta = BTMeq(theta, M1)*180/pi;
    
    % Printing
    fprintf("<strong>Example 9.6</strong>:\n")
    fprintf('beta = %2.1f\n', beta);
    
    % Mn1 = M1 * sind(beta);
    Mn1 = 2.05;
    [~,p2op1,~,~,~,~] = shock_calc(Mn1);
    fprintf('P2/P1 = %1.3f\n\n', p2op1);
    
% EX 9.7
    M1 = 3.6;
    theta = 10 * pi/180;
    beta = BTMeq(theta, M1);
    Mn1 = M1*sin(beta);
    
    [Mn2,p2op1,~,t2ot1,~,~] = shock_calc(Mn1);
    
    % Printing
    fprintf("<strong>Example 9.7</strong>:\n")
    fprintf('beta = %2.0f\nMn2 = %0.4f\nP2/P1 = %1.2f\nT2/T1 = %1.3f\n',...
        beta*180/pi, Mn2, p2op1, t2ot1);
    
    Mn2 = 2.96*sind(27.3);
    [Mn3,p3op2,~,t3ot2,~,~] = shock_calc(Mn2);
    
    % Printing
    fprintf('\nMn2 = %0.3f\nP3/P2 = %1.3f\nT3/T2 = %1.3f\nMn3 = %1.4f\n\n',...
        Mn2, p3op2, t3ot2, Mn3);
    
% EX 9.8
    Mn1 = 6.9;
    [~, p2op1, ~, t2ot1, ~, ~] = shock_calc(Mn1);
    
     % Printing
    fprintf("<strong>Example 9.8</strong>:\n")
    fprintf('P2/P1 = %2.2f\nT2/T1 = %2.1f\n\n',...
        p2op1, t2ot1);
   
% EX 9.9
    M1 = 1.5; theta = 15*pi/180;
    v1 = PMeq_findv(M1);
    v2 = v1 + theta;
    M2 = PMeq_findM(v2, 2);
    
    % Printing
    fprintf("<strong>Example 9.9</strong>:\n")
    fprintf('v1 = %2.2f\nM2 = %1.1f\n\n',...
        v1*180/pi, M2);
    
    [p01op1, t01ot1, ~] = isentropic(M1);
    [p02op2, t02ot2, ~] = isentropic(M2);
    fprintf('P01/P1 = %1.3f\nT01/T1 = %1.2f\nP02/P2 = %1.2f\nT02oT2 = %1.1f\n\n',...
        p01op1, t01ot1, p02op2, t02ot2);
    
% EX 9.10
    M1 = 10; theta = 15*pi/180;
    v1 = PMeq_findv(M1);
    v2 = v1 - theta;
    M2 = PMeq_findM(v2, 2);
    
    [p01op1, ~, ~] = isentropic(M1);
    [p02op2, ~, ~] = isentropic(M2);

    % Printing
    fprintf("<strong>Example 9.10</strong>:\n")
    fprintf('v1 = %2.1f\nM2 = %1.1f\n', v1*180/pi, M2);
    fprintf('P01/P1 = %1.3f\nP02/P2 = %1.2f\n\n', p01op1, p02op2);

% EX 9.11
    M1 = 10; theta = 15*pi/180;
    beta = BTMeq(theta, M1);
    Mn1 = M1*sin(beta);
    
    [Mn2, p2op1, ~, ~, ~, p02op01] = shock_calc(Mn1);
    M2 = Mn2/sin(beta - theta);
    
    % Printing
    fprintf("<strong>Example 9.11</strong>:\n")
    fprintf('beta = %1.0f\nP2/P1 = %2.2f\nP02/P01 = %1.4f\nMn2 = %1.4f\n\n'...
        , beta*180/pi, p2op1, p02op01, Mn2);
    
    [p01op1, ~, ~] = isentropic(M1);
    [p02op2, ~, ~] = isentropic(M2);
    
    fprintf('P01/P1 = %1.0f\nP02/P2 = %2.2f\n\n', p01op1, p02op2);
    
% EX 9.12
    M1 = 3; alpha = 5*pi/180;
    v1 = PMeq_findv(M1);
    v2 = v1 + alpha;
    M2 = PMeq_findM(v2, 2);
    
    [p01op1, ~, ~] = isentropic(M1);
    [p02op2, ~, ~] = isentropic(M2);
    
     % Printing
    fprintf("<strong>Example 9.12</strong>:\n")
    fprintf('v1 = %1.2f\nM2 = %1.2f\nP01/P1 = %2.2f\nP02/P2 = %2.0f\n\n'...
        , v1*180/pi, M2, p01op1, p02op2);
    
    beta = BTMeq(alpha, M1);
    Mn1 = M1*sin(beta);
    [~, p3op1, ~, ~, ~, ~] = shock_calc(Mn1);
    fprintf('beta = %2.1f\nP3/P1 = %1.3f\n\n', beta*180/pi, p3op1);
    
 % EX 9.13
    M1 = 7; alpha = 10*pi/180;
    v1 = PMeq_findv(M1);
    
    v2 = v1 + alpha;
    M2 = PMeq_findM(v2, 2);
    
    [p01op1, ~, ~] = isentropic(M1);
    [p02op2, ~, ~] = isentropic(M2);
    
     % Printing
    fprintf("<strong>Example 9.13</strong>:\n")
    fprintf('v1 = %2.2f\nM2 = %1.2f\nP01/P1 = %1.2f\nP02/P2 = %1.2f\n\n',...
        v1*180/pi, M2, p01op1, p02op2);
    
    
    beta = BTMeq(alpha, M1)*180/pi;
    Mn1 = 1.99;
    [~, p3op1, ~, ~, ~, ~] = shock_calc(Mn1);
    
    v2 = 95.97*pi/180;
    M2 = PMeq_findM(v2, 2);
    [p02op2, ~, ~] = isentropic(M2);
    
    fprintf('beta = %2.1f\nP3/P1 = %1.2f\nM2 = %1.1f\nP02/P2 = %1.4f\n\n',...
        beta, p3op1, M2, p02op2);
    
    [p01op1, ~, ~] = isentropic(M1); 
    theta = 15*pi/180;
    beta = BTMeq(theta, M1)*180/pi;
    
    Mn1 = M1*sin(beta);
    [~, p3op1, ~, ~, ~, ~] = shock_calc(Mn1);
    
    fprintf('P01/P1 = %1.4f\nbeta = %2.1f\nP3/P1 = %1.2f\n\n',...
        p01op1, beta, p3op1);
    

% EX 10.1
    AeoAstar = 10.25;
    Me = nozzle_area(AeoAstar);
    [p0ope, t0ote, ~] = isentropic(Me); 
    
     % Printing
    fprintf("<strong>Example 10.1</strong>:\n")
    fprintf('Me = %2.2f\nP0/Pe = %1.0f\nT0/Te = %1.2f\n\n',...
        Me, p0ope, t0ote);
    
% EX 10.2
    AeoAstar = 2;
    [~, Me] = nozzle(AeoAstar);
    [p0ope, t0ote, ~] = isentropic(Me); 
    
    % Printing
    fprintf("<strong>Example 10.2</strong>:\n")
    fprintf('Me = %2.1f\nP0/Pe = %1.2f\nT0/Te = %1.3f\n\n',...
        Me, p0ope, t0ote);
    
    AeoAstar = 2;
    [Me, ~] = nozzle(AeoAstar);
    [p0ope, t0ote, ~] = isentropic(Me);
    
    fprintf('Me = %2.1f\nP0/Pe = %1.3f\nT0/Te = %1.3f\n\n',...
        Me, p0ope, t0ote);

% EX 10.3
    p0ope = 1.028;
    Me = isentropicFindM(1.4, p0ope);
    AeoAstar = AratFindA(1.4, Me);
    
    AtoAstar = 1.482;
    [Mt, ~] = nozzle(AtoAstar);
    
    % Printing
    fprintf("<strong>Example 10.3</strong>:\n")
    fprintf('Me = %2.1f\nAe/A* = %1.3f\nMt = %1.2f\n\n',...
        Me, AeoAstar, Mt);

% EX 10.6
    M = 2;
    [~, ~, ~, ~, ~, p02op01] = shock_calc(M);
    
     % Printing
    fprintf("<strong>Example 10.6</strong>:\n")
    fprintf('P02/P01 = %1.4f\n\n', p02op01);
    
    

%% Problem 2
% Consider a diamond-wedge airfoil such as shown in Figure 9.36, with a
% half-angle ε = 10◦. The airfoil is at an angle of attack α = 15◦ to a Mach 3
% freestream. Calculate the lift and wave-drag coefficients for the airfoil

clear;
fprintf("\n\n<strong>----------PROBLEM #2----------</strong>\n")

% Given constants
epsilon = 10*pi/180;
alpha = 15*pi/180;
M1 = 3;
gamma = 1.4;

% Region numbering is 1 for freestream, 2 for top left, 3 for top right, 
% 4 for bottom left and 5 for bottom right

% Region 1 to 2 - PM Expansion fan
theta1 = alpha - epsilon;
v_M1 = PMeq_findv(M1);
v_M2 = theta1 + v_M1;
M2 = PMeq_findM(v_M2,2);

[p01op1, ~, ~] = isentropic(M1);
[p02op2, ~, ~] = isentropic(M2);
p2op1 = p01op1/p02op2;

% Region 2 to 3 - PM Expansion fan
theta2 = 2*epsilon;
v_M3 = theta2 + v_M2;
M3 = PMeq_findM(v_M3,2);

[p03op3, ~, ~] = isentropic(M3);
p3op1 = p01op1/p03op3;

% Region 1 to 4 - Oblique shock
theta3 = alpha + epsilon;
beta3 = BTMeq(theta3,M1);

Mn1 = M1*sin(beta3);
[M4n, p4op1, ~, ~, ~, ~] = shock_calc(Mn1);
M4 = M4n/sin(beta3-theta3);
p4op1 = 1*p4op1;

% Region 4 to 5 - PM Expansion fan
theta4 = 2*epsilon;
v_M4 = PMeq_findv(M4);
v_M5 = theta4 + v_M4;
M5 = PMeq_findM(v_M5,2);

[p05op5, ~, ~] = isentropic(M5);
p5op1 = p01op1/p05op5;

% Normal and axial coefficients
cn = (p5op1 + p4op1 - p3op1 - p2op1) / (gamma * M1^2);
ca = (p2op1 + p4op1 - p3op1 - p5op1) * tan(epsilon) / (gamma * M1^2);

% Lift and drag coefficients
cl = cn*cos(alpha) - ca*sin(alpha);
cd = cn*sin(alpha) + ca*cos(alpha);

fprintf("\nLift coefficient: " + cl + "\nWave Drag coefficient: " + cd + "\n");


%% Problem 3
% Using given MATLAB functions to reproduce Anderson Figure 9.9

clear;
fprintf("\n\n<strong>----------PROBLEM #3----------</strong>\n")

% Defining all the mach functions to be plotted
M_array = [1:0.05:1.5 1.6:0.1:2 2.2:0.2:4 4.5 5 6 8 10];
beta_array = (0:1:90)*pi/180;

figure;
hold on;
xlim([0 60]);

fprintf("Plotting Theta v. Beta curves...\n")

% A for loop to plot each curve
for i = 1:length(M_array)    
    theta_array = zeros(1,91);
    
    % Another for loop to create a row of thetas for the given betas
    for j = 1:91
        theta_array(j) = BTMeq_FindT(beta_array(j), M_array(i));
    end
    
    % Plotting and creating the legend entry
    plot(theta_array*180/pi, beta_array*180/pi);
    legend_array(i) = "M = " + M_array(i);
    
    
end

% Formatting
title("The \beta-\theta-M Diagram")
xlabel("Deflection Angle - \theta degrees")
ylabel("Shockwave Angle - \beta degrees")
legend(legend_array);
    
