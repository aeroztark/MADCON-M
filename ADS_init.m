function[q,init_mode] = ADS_init(SS_avail,MAG_avail,MAG_reading,mag_eci,SS_reading,sun_eci)

% Function to initialize ADS

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


    




