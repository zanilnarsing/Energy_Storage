%% Function to create a test set of temperature vectors for given cycle characteristics

function [temp] = Temp_vector(time_charge, time_store, time_discharge, num_cycle,temp_min,temp_max,charging_first)

time_onecyc = time_charge+time_store+time_discharge;


time_tot = (time_charge + time_store + time_discharge)*num_cycle;


temp = zeros(time_tot+1,1);
time = [0:time_tot]';


i = 1; 


if charging_first == 1
    while i <= num_cycle
                
        start = time_onecyc*(i-1) + 1;
        
        for c = start:(time_charge+(time_onecyc*(i-1)))
            temp(c,1) = temp_min + (c-start)*((temp_max-temp_min)/time_charge);
        end
        
        for s = c:(c+time_store)
            temp(s,1) = temp_max;
        end
        
        for d = s:(s+time_discharge)
            temp(d,1) = temp_max - (d-s)*((temp_max-temp_min)/time_discharge);
        end
        
        i = i+1;
    end
else
    while i <= num_cycle
               
        start = time_onecyc*(i-1) + 1;
        
        for d = start:(time_discharge+(time_onecyc*(i-1)))
            temp(d,1) = temp_max - (d-start)*((temp_max-temp_min)/time_discharge);
        end
        
        for c = d:(d+time_charge)
            temp(c,1) = temp_min + (c-d)*((temp_max-temp_min)/time_charge);
        end
        
        for s = c:(c+time_store)
            temp(s,1) = temp_max;
        end
        
        i = i+1;
        
    end
end
figure(1)
plot(time,temp);
title('Temperature Profile of Steam During Cycles');
xlabel('Time(s)');
ylabel('Temperature (K)');