function [XB,YB] = NACA_Airfoils(M, P, T, C, N)
    
    % Initializing an un-cambered mean chordline
    XB = linspace(0, C, (N+2)/2);
    XB_C = XB/C;
    T = T/100; % Converting thickness
    
    % Formula for thickness of any airfoil
    A = [0.2969, -0.1260, -0.3516, 0.2843, -0.1036];
    YT = T/0.2 * C * (A(1)*sqrt(XB_C) + A(2)*XB_C + A(3)*(XB_C).^2 + A(4)*XB_C.^3 + A(5)*XB_C.^4);

    if M == 0 && P == 0 % Symmetrical
        
        % For symmetrical airfoils it is just the airfoil thickness
        % Panels go from trailing edge to leading edge back to trailing edge
        XB = [flip(XB), XB(2:end)];
        YB = [-flip(YT), YT(2:end)];
       
    else 
        
        % Location and magnitude of maximum camber from NACA number
        M = M/100;
        P = P/10;
        
        % Pre-allocating size of arrays for camber line 
        YC = zeros(1, length(XB_C));
        dYC_dx = zeros(1, length(XB_C));
        
        % A counter for indexing
        i = 1;
        
        for xc = XB_C
            
            % Evaluating the camber line functions piecewise
            if xc <= P
                YC(i) = (M/(P^2)) * (2*P*xc - xc^2);
                dYC_dx(i) = (2*M/P^2) * (P - xc);
            elseif xc > P && xc <= 1
                YC(i) = (M/(1-P)^2) * ((1 - 2*P) + (2*P*xc) - xc^2);
                dYC_dx(i) = (2*M/(1-P)^2) * (P - xc);
            end
            
            % Counter update
            i = i+1;
        end
        
        % Calculating theta, needed for boundary layer calculation
        THETA = atan(dYC_dx);
        
        % X-Coords of BL, Panels go from TE to LE back to TE
        XU = XB - YT.*sin(THETA);
        XL = XB + YT.*sin(THETA);
        XB = [flip(XL), XU(2:end)];
        
        % Y-Coords of BL
        YU = YC + YT.*cos(THETA);
        YL = YC - YT.*cos(THETA);
        YB = [flip(YL), YU(2:end)];
        
        
    end
    
     
end

