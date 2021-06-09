function[q,init_mode] = ADS_init(SS_avail,MAG_avail,MAG_reading,mag_eci,SS_reading,sun_eci)

% Function to initialize ADS

% INPUTS: SS_avail- flag to denote SS availability (0 -> unavailable)
%       MAG_avail - flag to denote MAG availability (0 -> unavailable)
%       MAG_reading - measured mag body vector
%       mag_eci - mag ECI vector
%       SS_reading - measured sun body vector
%       sun_eci - sun ECI vector
% OUTPUTS: q -> initialization quaternion
%          init_mode -> ADS initialization mode (0-> coarse, 1-> normal)

    if SS_avail && MAG_avail == 1
            [q,init_mode] = Normal_init(MAG_reading,mag_eci,SS_reading,sun_eci);
    else 
            [q,init_mode] = Coarse_init();
    end
    
    function[q,init_mode]= Normal_init(MAG_reading,mag_eci,SS_reading,sun_eci)
        q = QUEST(MAG_reading,mag_eci,SS_reading,sun_eci,0.4,0.6); %row vector
        init_mode = 1;
    end
    function[q,init_mode]= Coarse_init()
        q = [0.5 0.5 0.5 0.5];
        init_mode = 0;
    end

end


    




