close all;
clear all;
clc;

%defining dimensions
for T_air_in = 10:10:50
        num_fins= 300; 
        fin_length = 1.5;
        fin_width = 0.2;
        tube_length = 1.5;
        fin_thickness = 0.001;
        num_tubes= 620;
        num_rows = 5;
        num_columns = num_tubes/num_rows;
        tube_diam = 0.00795; %meters (m)
        
        tube_spacing_lengthwise = fin_length/(num_columns+1);
        tube_spacing_widthwise = fin_width/(num_rows+1);
        transverse_pitch = tube_diam + tube_spacing_lengthwise;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %defining fluid properties for ethylene glycol
        mass_flow_ethyl = 1.2; %kg/s
        mass_flow_ethyl_per_tube = mass_flow_ethyl/num_tubes;
        
        T_ethyl_in = 90;
        T_ethyl_out = 60;
        
        cp_ethyl = 3641.5; %J/KgC;
        dyn_visc_ethyl = 0.000906;
        K_ethyl = 0.3947; %thermal conductivity of fluid
        Pr_ethyl = 9.82;
        density_ethyl = 1045;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        total_heat_to_be_removed = mass_flow_ethyl*cp_ethyl*(T_ethyl_in-T_ethyl_out);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %defining fluid properties for air
        %T_air_in = 30;
        Temp_diff = 20;
        T_air_out = T_air_in + Temp_diff;
        cp_air = 1008;
        dyn_visc_air = 1.918E-05;
        K_air = 0.02662;
        Pr_air = 0.7255;
        density_air = 1.127;
        
        %mass_flow_air = total_heat_to_be_removed/(cp_air*(T_air_out - T_air_in)); %kg/s
        mass_flow_air = 53.5;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Inside the tube: convective heat transfer coefficient
        v_inside = ((4*mass_flow_ethyl_per_tube)/(density_ethyl*pi*(tube_diam)^2));
        Re_inside = (density_ethyl*tube_diam*v_inside)/dyn_visc_ethyl;
        
        if Re_inside < 2300
            Nu_inside = 4.36;
        elseif Re_inside>2300 && Re_inside<10000
            f= ((0.79*log(Re_inside))-1.64)^-2;
            Nu_numerator = (f/8)*(Re_inside-1000)*Pr_ethyl;
            Nu_denominator= 1+(12.7*((f/8)^0.5)*((Pr_ethyl^(2/3))-1));
            Nu_inside = Nu_numerator/Nu_denominator;
        %     Nu_inside = ((f/8)*(Re-1000)*Pr_ethyl)/(1+(12.7*((f/8)^0.5)*((Pr_ethyl)^(2/3))-1));
        
        else
            Nu_inside = 0.023*((Re_inside)^0.8)*(9.82^(1/3));
        end
        
        h_inside = (K_ethyl*Nu_inside)/tube_diam;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Outside the tube: convective heat transfer coefficient
        
        v_outside = mass_flow_air/(density_air*tube_length*fin_length);
        v_outside_max = (transverse_pitch*v_outside)/(transverse_pitch-tube_diam);
        Re_outside = (density_air*tube_diam*v_outside_max)/dyn_visc_air;
        Nu_outside = 0.27*(Re_outside^0.63)*(Pr_air^0.36);
        Nu_corrected = 0.93*Nu_outside ;
        h_outside = (Nu_corrected*K_air)/tube_diam;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %calculating relevant areas
        
        A_inside = num_tubes*pi*tube_diam*tube_length;
        A_unfin = num_tubes*((pi*tube_diam*tube_length)-(num_fins*pi*tube_diam*fin_thickness));
        A_fin= num_fins*((2*((fin_length*fin_width)+(fin_width*fin_thickness)+(fin_length*fin_thickness))-(2*num_fins*pi*(tube_diam^2)*0.25)));
        A_outside = A_unfin + A_fin;
        
        A_nofin= A_inside;
        
        %Calculating efficiency of the fin
        K_aluminum = 401;
        a=(2*h_outside)/(K_aluminum*fin_thickness);
        m = sqrt(a);
        Lc= fin_length+(fin_thickness/2);
        fin_efficiency = tanh(m*Lc)/(m*Lc);
        % fin_efficiency = 0.7;
        
        A_s = A_unfin + (fin_efficiency*A_fin);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Calculating R
        R_f=0.00035;
        
        R_thermal= (1/(h_inside*A_inside))+(1/(h_outside*A_outside))+(R_f/A_inside);
        
        U = 1/(A_s*R_thermal);
        U_inside = (1/(A_inside*R_thermal));
        U_outside = (1/(A_outside*R_thermal));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculating NTU
        
        C_ethyl= mass_flow_ethyl*cp_ethyl;
        C_air = mass_flow_air*cp_air;
        
        if C_ethyl<C_air
            C_min= C_ethyl;
            C_max=C_air;
        else 
            C_min = C_air;
            C_max = C_ethyl;
        end
        
        NTU = (U*A_s)/C_min;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculating LMTD
        
        LMTD = ((T_ethyl_in-T_air_out)-(T_ethyl_out-T_air_in))/(log((T_ethyl_in-T_air_out)/(T_ethyl_out-T_air_in)));
        
        P= (T_ethyl_out-T_ethyl_in)/(T_air_in-T_ethyl_in);
        R = (T_air_in-T_air_out)/(T_ethyl_out-T_ethyl_in);
    
        F=0.95;
        
        LMTD_corrected = F*LMTD;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculating effectiveness
        
        delta_T_max = T_ethyl_in-T_air_in;
        Q_dot_max = C_min * delta_T_max;
        Effectiveness_1 = total_heat_to_be_removed/Q_dot_max;
        
        c = C_min/C_max;
        exp_term = ((NTU^0.22)/c)*((exp(-c*(NTU^0.78)))-1);
        Effectiveness_2 = 1 - exp(exp_term);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Finding U actual
       
        U_actual = total_heat_to_be_removed/(A_s*LMTD_corrected);
        U_percent_diff = ((U-U_actual)/U_actual)*100;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clc;
fprintf('Values for Ambient air temperature: %i \n', T_air_in);
fprintf('Physics U= %f \n', U);
fprintf('Theoretical/Math U= %f \n', U_actual);
fprintf('U percentage difference = %f \n', U_percent_diff);
fprintf('Effectiveness= %f \n', Effectiveness_2);
disp("__________________________________________________________")
%fprintf('NTU= %f \n', NTU);
%fprintf('Fin Efficiency= %f \n', fin_efficiency);
end
disp("Done")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Off-design performance




