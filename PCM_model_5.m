%new test, using Temp_check
% PCM Model 4.0

%...Define constants and initial conditions
Ri = 0.4046; %[m] inner diameter of PCM, outer diameter of pipe
Ro = 3; % [m] outer diameter of PCM , just a guess for now
nr = 99; %number of nodes in the radial direction
dr = (Ro - Ri)/(nr-1); % [m] Radial step
r = Ri:dr:Ro; %setting up vector for radius steps
l = 300000; %[m] length of pipe

tmax = 3000; % [s] max time
dt = 5;
nt = (tmax/dt) + 1; % nodes in time variable
t = 1:nt; % creating a time vector **note the error in using time steps, some inaccuracy**


%z = zeros(nt,1);%Mass Fraction of Liquid
rho_s = 2.0256e3; % [kg/m^3] Density of Solid PCM
rho_l = 2.0256e3; %Density of Liquid PCM **assume same as Solid Density for now**
%rho_m = %Average Density of Mixture PCM

k_s = 5; % [W/m*K] ***this is just a guessed value*** Thermal Conductivity of Solid PCM
k_l = .1; %Thermal Conductivity of Liquid PCM
%k_m = zeros(nt,1);
%k_m(i) = z(i)*(k_l) + (1-z)*k_s;

cp_s = 141.9; % [J/kg*K] Specific heat capacity of solid PCM
cp_l = 85; %Specific heat capacity of liquid PCM
%cp_m = %Avg specific heat capacity of mixture PCM

alpha_s = k_s / (rho_s*cp_s); % [m2/s] thermal diffusivity of solid PCM
alpha_l = k_l / (rho_l*cp_l); %[m2/s] thermal diffusivity of liquid PCM

L = 1.775e6; % [J/kg] latent heat of fusion

h = 2.79; %[W/m2*K] Free convection coefficient


T_init = 300; %[K] temperature of the PCM initially
T_inf = 300; %[K] temperature of the surroundings
T_melt = 400; %[K} *assumed* melting point

%... variables to calculate wall temperature
t_charge = 60; %[s] Time to reach full temp
T_max = 500; %[K] Max temperature reached during charging cycle
t_store = tmax - t_charge; %time max temperature is held

%Finding the inner wall temperature at each time step
for a = 1:nt
    if a <= t_charge
        T_Ri(a) = T_init + (a-1)*((T_max-T_init)/t_charge);
    else
        T_Ri(a) = T_max;
    end
end


T = zeros(nr,nt); %Initializing Temperature Profile
T(1,:) = T_Ri; % Initial Conditions
T(:,1) = T_init;

q_in = zeros(nr,nt); % Initializing Vector for heat energy [J] in
q_in(1) = 0; % First entry is zero, no temperature gradient initially


R_melt = zeros(1,nt);
R_melt(1) = 0;
Qdot_in = zeros(1,nt);
Qdot_in(1) = 0;
Qdot_out = zeros(1,nt);
Qdot_out(1) = 0;
Qdot_store = zeros(1,nt);
Qdot_store(1) = 0;
R_pc = zeros(1,nt);
R_pc(1) = 0;

count1 = 1;
count2 = 1;
count3 = 1;
count4 = 1;


%... Beginning Main Loop (time loop)
for m = 2:nt
    for n = 2:(nr-1)
        r_plus = (r(n)+ r(n+1))/2;
        r_minus = (r(n)+ r(n-1))/2;
        %... Solid Heat Diffusion
        T(n,m) = T(n,m-1) + ((alpha_s*dt)/(r(n)*dr))*(((r_plus/dr)*(T(n+1,m-1) - T(n,m-1))) - ((r_minus/dr)*(T(n,m-1) - T(n-1,m-1))));
        if n == nr-1
            T(n,m) = T(n,m-1) + ((alpha_s*dt)/(r(n)*dr))*(((r_plus/dr)*(T(n,m-1) - T(n-1,m-1))) - ((r_minus/dr)*(T(n-1,m-1) - T(n-2,m-1))));
        end
        %... Calculating Energy [J] into PCM
        %q_in(n,m) = pi()*rho_s*cp_s*l*(Ro^2 - Ri^2)*(T(n-1,m) - T(n,m)); %Finding q [J] between each node
        %q_in_total = sum(q_in); %[J], Finding total q for each time step, summing across radial nodes
    end
    Temp_check = T(2:nr,m); %Storing all temp values for time step to check melting case
    [Temp_max, n_max] = max(Temp_check); %finding the max temp in Temp_check and the index where it occurs. *note, the index in T is n_max+1*
    if Temp_max > T_melt
        R_melt(m) = r(n_max+1) - (((T_melt - Temp_max)/(T(n_max+2) - Temp_max))*(r(n_max+2) - r(n_max+1)));
        n_melt_exact(m) = 1 + ((R_melt(m) - Ri)/(dr));
        n_melt(m) = ceil(n_melt_exact(m));
        Qdot_in(m) = 2*pi()*k_s*l*(T(n_max,m) - T_melt)/(log(R_melt(m)/r(n_max)));
        Qdot_out(m) = 2*pi()*k_s*l*(T_melt - T(n_melt(m),m))/(log(r((n_melt(m)))/R_melt(m)));
        Qdot_store(m) = Qdot_in(m) - Qdot_out(m);
        R_pc(m) = sqrt(((R_melt(m)).^2) - ((Qdot_store(m) .* dt)/(L*rho_s*l*pi())));
        n_pc_exact(m) = 1 + ((R_pc(m) - Ri)/(dr));
        n_pc(m) = floor(n_pc_exact(m));
        if n_pc(m) == 1
            n_pc(m) = 2;
        end
        
        
        if R_pc(m) < Ri
            for i = 2:(n_melt(m))
                T(i,m) = T_melt;
                count3 = count3+1;
            end
            for y = (n_melt(m)+1):(nr-1)
                r_plus = (r(y)+ r(y+1))/2;
                r_minus = (r(y)+ r(y-1))/2;
                T(y,m) = T(y,m-1) + ((alpha_s*dt)/(r(y)*dr))*(((r_plus/dr)*(T(y+1,m-1) - T(y,m-1))) - ((r_minus/dr)*(T(y,m-1) - T(y-1,m-1))));
                if y == nr-1
                    T(y,m) = T(y,m-1) + ((alpha_s*dt)/(r(y)*dr))*(((r_plus/dr)*(T(y,m-1) - T(y-1,m-1))) - ((r_minus/dr)*(T(y-1,m-1) - T(y-2,m-1))));
                end
                count4 = count4+1;
            end
        else
            if n_pc(m) <= 3;
                T(n_pc(m),m) = T_melt;
                for y = n_melt(m):(nr-1)
                    r_plus = (r(y)+ r(y+1))/2;
                    r_minus = (r(y)+ r(y-1))/2;
                    T(y,m) = T(y,m-1) + ((alpha_s*dt)/(r(y)*dr))*(((r_plus/dr)*(T(y+1,m-1) - T(y,m-1))) - ((r_minus/dr)*(T(y,m-1) - T(y-1,m-1))));
                    if y == nr-1
                        T(y,m) = T(y,m-1) + ((alpha_s*dt)/(r(y)*dr))*(((r_plus/dr)*(T(y,m-1) - T(y-1,m-1))) - ((r_minus/dr)*(T(y-1,m-1) - T(y-2,m-1))));
                    end
                end
            else
                for u = 2:(n_pc(m)-1)
                    r_plus = (r(u)+ r(u+1))/2;
                    r_minus = (r(u)+ r(u-1))/2;
                    T(u,m) = T(u,m-1) + ((alpha_l*dt)/(r(u)*dr))*(((r_plus/dr)*(T(u+1,m-1) - T(u,m-1))) - ((r_minus/dr)*(T(u,m-1) - T(u-1,m-1))));
                    if u == nr-1
                        T(u,m) = T(u,m-1) + ((alpha_l*dt)/(r(u)*dr))*(((r_plus/dr)*(T(u,m-1) - T(u-1,m-1))) - ((r_minus/dr)*(T(u-1,m-1) - T(u-2,m-1))));
                    end
                end
                for i = n_pc(m):n_melt
                    T(i,m) = T_melt;
                end
                for y = (n_melt(m)+1):(nr-1)
                    r_plus = (r(y)+ r(y+1))/2;
                    r_minus = (r(y)+ r(y-1))/2;
                    T(y,m) = T(y,m-1) + ((alpha_s*dt)/(r(y)*dr))*(((r_plus/dr)*(T(y+1,m-1) - T(y,m-1))) - ((r_minus/dr)*(T(y,m-1) - T(y-1,m-1))));
                    if y == nr-1
                        T(y,m) = T(y,m-1) + ((alpha_s*dt)/(r(y)*dr))*(((r_plus/dr)*(T(y,m-1) - T(y-1,m-1))) - ((r_minus/dr)*(T(y-1,m-1) - T(y-2,m-1))));
                    end
                end
            end
        end
    end
end