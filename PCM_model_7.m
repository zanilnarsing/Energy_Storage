% PCM Model 7.0
% Fixed Boundary layer width
% Attempt to include Charging, Discharge, and Storage
% This is the current latest version as of 8/9/17

%...Define constants and initial conditions

% Info obtained from the Accumulator... Coded in here for proof of concept,
% this section should be replaced when implementing into main code



time_charge = 50;                      % [s]
time_store = 25;                       % [s]
time_discharge = 50;                   % [s]
num_cycle = 3;
temp_i = 300;                           % [K]
temp_max = 500;                         % [K]
charging_first = 1;                     % 1 = charging first, 0 = discharging first

T_steam = Temp_vector(time_charge, time_store, time_discharge, num_cycle,temp_i,temp_max,charging_first);

% User Input conditions
Ri = 0.4046;                        % [m] inner diameter of PCM, outer diameter of pipe
Ro = 3;                             % [m] outer diameter of PCM , just a guess for now
nr = 99;                            % number of nodes in the radial direction
dr = (Ro - Ri)/(nr-1);              % [m] Radial step
r = Ri:dr:Ro;                       % setting up vector for radius steps
l = 300000;                         % [m] length of pipe


% Time Conditions
nt = 100000;
t = 1:(time_charge + time_store + time_discharge)*num_cycle;
dt = 5;

% Initial Temp Variables *** Look at again***                                    
T_init = 290;                       % [K] temperature of the PCM initially
T_inf = 290;                        % [K] temperature of the surroundings
T_melt = 315;                       % [K} *assumed* melting point                                    

% Material Properties (hard code) (needs updating to an usable PCM)
rho_s = 814;                        % [kg/m^3] Density of Solid PCM
rho_l = 774;                        % Density of Liquid PCM                                    
k_s = 0.38;                         % [W/m*K] Thermal Conductivity of Solid PCM
k_l = 0.35;                         % Thermal Conductivity of Liquid PCM
cp_s = 1900;                        % [J/kg*K] Specific heat capacity of solid PCM
cp_l = 2200;                        % Specific heat capacity of liquid PCM
alpha_s = k_s / (rho_s*cp_s);       % [m2/s] thermal diffusivity of solid PCM
alpha_l = k_l / (rho_l*cp_l);       % [m2/s] thermal diffusivity of liquid PCM
L = 0.242e6;                        % [J/kg] latent heat of fusion

%... meta time loop, following the accumulator - this outer loop included
%in the code for proof of concept. This outer loop should NOT be included
%when implementing the code in the main accu. model

accu_time = size(T_steam);
for a = 1:accu_time
    T_wall = T_steam(a);                

    %... Initializing vectors
    T = zeros(nr,nt);                   % Initializing Temperature Profile
    T(1,:) = T_wall;                      
    T(:,1) = T_init;
    q_in = zeros(nr,nt);                % Initializing Vector for heat energy [J] in
    q_in(1) = 0;                        % First entry is zero, no temperature gradient initially
    R_melt = zeros(1,nt);
    R_melt(1) = 0;
    n_melt = zeros(1,nt);
    n_melt(1) = 0;
    R_pc = zeros(1,nt);
    n_pc = zeros(1,nt);
    Qdot_in = zeros(1,nt);
    Qdot_in(1) = 0;
    Qdot_out = zeros(1,nt);
    Qdot_out(1) = 0;
    Qdot_store = zeros(1,nt);
    Qdot_store(1) = 0;
    R_pc = zeros(1,nt);
    R_pc(1) = 0;

    ifmelted = 0;                       % used to determine where to draw previous time step Temp values from. Initializing at 0. 


    %... Beginning Time Loop(m)
    for m = 2:nt %beginning at nt=2 b/c nt=1 is already established in inital conditions
        if ifmelted == 0 %case to determine where to draw previous Temp values from

            %... Beginning Space Loop (n)
            for n = 2:(nr-1)            
                r_plus = (r(n)+ r(n+1))/2; %used for calc during diffusion
                r_minus = (r(n)+ r(n-1))/2;
                %... Solid Heat Diffusion
                T(n,m) = T(n,m-1) + ((alpha_s*dt)/(r(n)*dr))*(((r_plus/dr)*(T(n+1,m-1) - T(n,m-1))) - ((r_minus/dr)*(T(n,m-1) - T(n-1,m-1))));
                if n == nr-1 %case for outer edge of pipe
                    T(n,m) = T(n,m-1) + ((alpha_s*dt)/(r(n)*dr))*(((r_plus/dr)*(T(n,m-1) - T(n-1,m-1))) - ((r_minus/dr)*(T(n-1,m-1) - T(n-2,m-1))));
                end
                %... Calculating Energy [J] into PCM
                %q_in(n,m) = pi()*rho_s*cp_s*l*(Ro^2 - Ri^2)*(T(n-1,m) - T(n,m)); %Finding q [J] between each node
                %q_in_total = sum(q_in); %[J], Finding total q for each time step, summing across radial nodes
            end
            Temp_check = T(2:nr,m);%Storing all temp values for time step to check melting case (excluding wall temperature)
            [Temp_max, n_max] = max(Temp_check);%finding the max temp in Temp_check and the index where it occurs. *note, the index in T is n_max+1*
            if Temp_max > T_melt %case if melting

                ifmelted = 1; %changing if melting case

                if n_pc(m-1) == 0 %if n_pc doesn't exist
                    n_pc(m) = n_max+1; %note; n_max index is 1 less than n_pc index
                    R_pc(m) = r(n_pc(m)); %adding location of Phase Change Radius
                    n_melt(m) = n_pc(m) + 1; %Fixed width of interface, 1 index length (ie. interface width = nr)            
                    R_melt(m) = r(n_melt(m)); %adding location of Melting Radius 

                else %case where n_pc > 0
                    Qdot_in(m) = 2*pi()*k_l*l*(T((n_pc(m-1)-1),(m-1)) - T_melt)/(log(R_pc(m-1)/r(n_pc(m-1)-1))); %finding heat in
                    Qdot_out(m) = 2*pi()*k_s*l*(T_melt - T((n_melt(m-1)+1),(m-1)))/(log(r(((n_melt(m-1)+1)))/R_melt(m-1)));%finding heat out
                    Qdot_store(m) = Qdot_in(m) - Qdot_out(m); %finding heat stored (this heat moves interface)
                    R_pc(m) = sqrt(((R_pc(m-1))^2) - ((Qdot_store(m) * dt)/(L*rho_s*l*pi()))); %finding new value of phase change radius
                    n_pc_exact(m) = 1 + ((R_pc(m) - Ri)/(dr)); %%exact location of new phase change interface index
                    n_pc(m) = floor(n_pc_exact(m)); %integer value of new phase change interface index (for simplicity)
                    R_pc(m) = r(n_pc(m)); %finding location of new phase change radius
                    n_melt(m) = n_pc(m) + 1; %finding new index of  melting radius (fixed interface width)
                    R_melt(m) = r(n_melt(m));%finding new location of melting radius


                end

                % Setting up Temperatures in sections ***(Creating T array)***

                %... Interface Temperatures (Temp at R_pc and R_melt)
                T(n_pc(m), m) = T_melt; 
                T(n_melt(m), m) = T_melt;

                %... Solid Temperatures (below melting point)
                for y = (n_melt(m)+1):(nr-1) % Solid Diffusion (locations under melting temperature)
                    r_plus = (r(y)+ r(y+1))/2; %used for diffusion calculations
                    r_minus = (r(y)+ r(y-1))/2;
                    T(y,m) = T(y,m-1) + ((alpha_s*dt)/(r(y)*dr))*(((r_plus/dr)*(T(y+1,m-1) - T(y,m-1))) - ((r_minus/dr)*(T(y,m-1) - T(y-1,m-1))));
                    if y == nr-1 %Boundary condition
                        T(y,m) = T(y,m-1) + ((alpha_s*dt)/(r(y)*dr))*(((r_plus/dr)*(T(y,m-1) - T(y-1,m-1))) - ((r_minus/dr)*(T(y-1,m-1) - T(y-2,m-1))));
                    end
                end

                %... Liquid Temperatures (above melting point)
                if n_pc(m)>= 3 %Liquid diffusion only neccessary if n_pc is greater than 3, otherwise all values known already
                    for u = 2:(n_pc(m)-1) % Liquid Diffusion (locations above melting temperature)
                        r_plus = (r(u)+ r(u+1))/2;
                        r_minus = (r(u)+ r(u-1))/2;
                        T(u,m) = T(u,m-1) + ((alpha_l*dt)/(r(u)*dr))*(((r_plus/dr)*(T(u+1,m-1) - T(u,m-1))) - ((r_minus/dr)*(T(u,m-1) - T(u-1,m-1))));
                        if u == nr-1
                            T(u,m) = T(u,m-1) + ((alpha_l*dt)/(r(u)*dr))*(((r_plus/dr)*(T(u,m-1) - T(u-1,m-1))) - ((r_minus/dr)*(T(u-1,m-1) - T(u-2,m-1))));
                        end

                    end
                end
            else
                ifmelted = 0; %case for not melting
            end
        else
            %... Case in which the PCM already began to melt
            Qdot_out(m) = 2*pi()*k_s*l*(T_melt - T((nr-1),(m-1)))/(log(r(nr)/R_melt(m-1)));
            Qdot_in(m) = 2*pi()*k_l*l*(T(1,(m-1)) - T_melt)/(log(R_pc(m-1)/r(1))); %finding heat in
            Qdot_store(m) = Qdot_in(m) - Qdot_out(m); %finding heat stored (this heat moves interface)
            R_pc_exact(m) = sqrt(((R_pc(m-1))^2) + ((Qdot_store(m) * dt)/(L*rho_s*l*pi()))); %finding new value of phase change radius
            n_pc_exact(m) = 1 + ((R_pc_exact(m) - Ri)/(dr)); %%exact location of new phase change interface index
            n_pc(m) = ceil(n_pc_exact(m)); %integer value of new phase change interface index (for simplicity)
            R_pc(m) = r(n_pc(m)); %finding location of new phase change radius
            n_melt(m) = n_pc(m) + 1; %finding new index of  melting radius (fixed interface width)
            R_melt(m) = r(n_melt(m));%finding new location of melting radius    

            check1 = ((R_pc(m-1))^2) + ((Qdot_store(m) * dt)/(L*rho_s*l*pi()));

            % Setting up Temperatures in sections ***(Creating T array)***

            %... Interface Temperatures (Temp at R_pc and R_melt)
            T(n_pc(m), m) = T_melt; 
            T(n_melt(m), m) = T_melt;

            %... Solid Temperatures (below melting point)
            for y = (n_melt(m)+1):(nr-1) % Solid Diffusion (locations under melting temperature)
                r_plus = (r(y)+ r(y+1))/2; %used for diffusion calculations
                r_minus = (r(y)+ r(y-1))/2;
                T(y,m) = T(y,m-1) + ((alpha_s*dt)/(r(y)*dr))*(((r_plus/dr)*(T(y+1,m-1) - T(y,m-1))) - ((r_minus/dr)*(T(y,m-1) - T(y-1,m-1))));
                if y == nr-1 %Boundary condition
                    T(y,m) = T(y,m-1) + ((alpha_s*dt)/(r(y)*dr))*(((r_plus/dr)*(T(y,m-1) - T(y-1,m-1))) - ((r_minus/dr)*(T(y-1,m-1) - T(y-2,m-1))));
                end
            end

            %... Liquid Temperatures (above melting point)
            if n_pc(m)>= 3 %Liquid diffusion only neccessary if n_pc is greater than 3, otherwise all values known already
                for u = 2:(n_pc(m)-1) % Liquid Diffusion (locations above melting temperature)
                    r_plus = (r(u)+ r(u+1))/2;
                    r_minus = (r(u)+ r(u-1))/2;
                    T(u,m) = T(u,m-1) + ((alpha_l*dt)/(r(u)*dr))*(((r_plus/dr)*(T(u+1,m-1) - T(u,m-1))) - ((r_minus/dr)*(T(u,m-1) - T(u-1,m-1))));
                    if u == nr-1
                        T(u,m) = T(u,m-1) + ((alpha_l*dt)/(r(u)*dr))*(((r_plus/dr)*(T(u,m-1) - T(u-1,m-1))) - ((r_minus/dr)*(T(u-1,m-1) - T(u-2,m-1))));
                    end
                end
            end
        end
    end
end



%... Creating Plots
figure(2)
N_curves = 5;
ir = round(linspace(1,nt,N_curves));
for p=1:N_curves
    plot(r,T(:,ir(p)))
    hold on;
end
title('Temperature Profiles During Conduction in PCM');
xlabel('Radial Location (m)');
ylabel('Temperature (K)');