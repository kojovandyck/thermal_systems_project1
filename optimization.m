close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREVIOUS OPTIMAL SPECIFICATIONS
% num_fins= 650; 
% fin_length = 1;
% fin_width = 0.2;
% tube_length = 1.5411;
% fin_thickness = 0.001;
% num_tubes= 600;
% num_rows = 20;
% num_columns = num_tubes/num_rows;
% tube_diam = 0.00635; %meters (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


num_fins= 650; 
tube_length = 1.5411;
fin_thickness = 0.001;
num_tubes= 600;
num_rows = 20;
num_columns = num_tubes/num_rows;
tube_diam = 0.00635; %meters (m)
transverse_pitch = 4*tube_diam;
longitudinal_pitch=  3*tube_diam;

tube_spacing_lengthwise = transverse_pitch - tube_diam ;
tube_spacing_widthwise = longitudinal_pitch - tube_diam;
fin_length = tube_spacing_lengthwise*(num_columns+1);
fin_width = tube_spacing_widthwise*(num_rows+1);

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

math_Q_rate = mass_flow_ethyl*cp_ethyl*(T_ethyl_in-T_ethyl_out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%defining fluid properties for air
T_air_in = 30;
Temp_diff = 20;
T_air_out = T_air_in + Temp_diff;
cp_air = 1008;
dyn_visc_air = 1.918E-05;
K_air = 0.02662;
Pr_air = 0.7255;
density_air = 1.127;
mass_flow_air = math_Q_rate/(cp_air*(T_air_out - T_air_in)); %kg/s
%mass_flow_air = 53.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Constants to Define
%LTMD Correction Factor
P= (T_ethyl_out-T_ethyl_in)/(T_air_in-T_ethyl_in);
R = (T_air_in-T_air_out)/(T_ethyl_out-T_ethyl_in);
F=0.95;
% Graph is on page 664

%No of Rows Correction Factor
row_number_correction_factor = 1;
% 3 rows = 0.86, 4 rows = 0.9, 5 rows = 0.93
% 7 rows = 0.96, 10 rows =0.98, 13 rows = 0.99

%Pressure drop and Pumping Power Constants
friction_factor = 0.17;
correction_factor = 0.8;
% Data can be found on pg 448

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
Nu_corrected = row_number_correction_factor*Nu_outside ;
h_outside = (Nu_corrected*K_air)/tube_diam;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculating relevant areas
A_inside = num_tubes*pi*tube_diam*tube_length;
A_unfin = num_tubes*((pi*tube_diam*tube_length)-(num_fins*pi*tube_diam...
    *fin_thickness));
A_fin= num_fins*((2*((fin_length*fin_width)+(fin_width*fin_thickness)...
    +(fin_length*fin_thickness))-(2*num_fins*pi*(tube_diam^2)*0.25)));
A_nofin= A_inside;
A_outside = A_unfin + A_fin;
    
%Calculating efficiency of the fin
K_copper = 401;
a=(2*h_outside)/(K_copper*fin_thickness);
m = sqrt(a);
area_of_unit_cell = (fin_length*fin_width)/num_tubes;
side_of_unit_cell = sqrt(area_of_unit_cell);

Lc= ((side_of_unit_cell-tube_diam)/2)+(fin_thickness/2);


fin_efficiency = tanh(m*Lc)/(m*Lc);

%Effective External Surface Area
A_s = A_unfin + (fin_efficiency*A_fin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating R
R_f=0.00035;
R_thermal= (1/(h_inside*A_inside))+(1/(h_outside*A_s))+(R_f/A_inside);

U = 1/(A_s*R_thermal);
%U_inside = (1/(A_inside*R_thermal));
%U_outside = (1/(A_s*R_thermal));
    
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

LMTD = ((T_ethyl_in-T_air_out)-(T_ethyl_out-T_air_in))/...
    (log((T_ethyl_in-T_air_out)/(T_ethyl_out-T_air_in)));



LMTD_corrected = F*LMTD;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating effectiveness

delta_T_max = T_ethyl_in-T_air_in;
Q_dot_max = C_min * delta_T_max;
physics_Q_rate = U*A_s*LMTD_corrected;
Effectiveness_1 = physics_Q_rate/Q_dot_max;

c = C_min/C_max;
exp_term = ((NTU^0.22)/c)*((exp(-c*(NTU^0.78)))-1);
Effectiveness_2 = 1 - exp(exp_term);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finding U math

U_math = math_Q_rate/(A_s*LMTD_corrected);
%U_perc_diff = abs(((U-U_math)/U_math)*100);
U_perc_diff = ((U-U_math)/U_math)*100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pressure drop and Pumping Power
P_L = longitudinal_pitch/tube_diam;
P_T = transverse_pitch/tube_diam;
correction_value_finder=(P_T-1)/(P_L-1);

pressure_drop = 0.5*num_rows*friction_factor*correction_factor*...
density_air*((v_outside_max)^2);
pumping_power = ((mass_flow_air*pressure_drop)/density_air);

volume = (fin_length*fin_width*fin_thickness*num_fins);
%(0.25*pi*(tube_diam^2)*tube_length)+...
SA_V_Ratio = A_outside/volume;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Changing Parameter (Pick a parameter to change)
% fprintf('Tube Length = %f, U percentage diff = %f \n', tube_length,U_perc_diff);

% values(count,1) = tube_length;
% values(count,2) =  U_perc_diff;
% count = count+1;
% %end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Plot Graph
% plot(values(:,1),values(:,2), 'red','LineWidth',3);
% xlabel('Tube Length')
% ylabel('U Percentage Difference')
% 
% hold on
% yline(0,'LineWidth',3);
% 
% hold off
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp(values);


fprintf('Values for Calculations (Not results)\n');
fprintf('_________________________________________________\n');
fprintf('P_L= %f \n', P_L);
fprintf('Correction Value Finder = %f \n', correction_value_finder);
fprintf('Reynolds Number Outside = %f \n', Re_outside);
fprintf('P = %f \n', P);
fprintf('R= %f \n', R);
fprintf('\n');
fprintf('Results under Normal Ambient Air \n');
disp('___________________________________________________');
fprintf('Dimensions of Heat Exchanger \n');
fprintf('Length of Fin (m) = %f \n', fin_length);
fprintf('Width of Fin (m)= %f \n', fin_width);
fprintf('Thickness of Fin (m)= %f \n', fin_thickness);
fprintf('Tube Length (m)= %f \n', tube_length);
fprintf('Tube Diameter (m)= %f \n', tube_diam);
fprintf('Longitudinal Pitch (m)= %f \n', longitudinal_pitch);
fprintf('Transverse Pitch (m)= %f \n', transverse_pitch);
fprintf('Volume = %f \n', volume);
fprintf('SA = %f \n', A_outside);
fprintf('SA:V = %f \n', SA_V_Ratio);
fprintf('\n');

fprintf('Performance Results \n');
fprintf('Mass flow rate of air = %f \n', mass_flow_air);
fprintf('Required Air Velocity = %f \n', v_outside);
fprintf('Ambient Air Temperature In: %f C \n', T_air_in);
fprintf('Exit Air Temperature: %f C \n', T_air_out);
fprintf('Q_Physics= %f W \n', physics_Q_rate);
fprintf('Q_Expected %f W \n', math_Q_rate);
fprintf('Physics U= %f W/(m2°C) \n', U);
fprintf('Expected U= %f W/(m2°C) \n', U_math);
fprintf('Effective Heat Transfer Area= %f \n', A_s);
fprintf('U percentage diff = %f \n', U_perc_diff);
fprintf('LMTD (inclusive of correction factor) = %f \n', LMTD_corrected);
fprintf('Fin Efficiency= %f \n', fin_efficiency);
fprintf('Effectiveness1= %f \n', Effectiveness_1);
fprintf('Effectiveness2= %f \n', Effectiveness_2);
fprintf('NTU= %f \n', NTU);
fprintf('Pressure Drop= %f Pa \n', pressure_drop);
fprintf('Pumping Power= %f kW \n', pumping_power/1000);


fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Off-design performance - different code
fprintf('Off-design Performance Evaluation \n');
fprintf('_________________________________________________\n');

for i = 10:10:50
    %Initialization 
    T_air_in = i;
    Q_dot_max2 = C_min*(T_ethyl_in-T_air_in);
    Q_dot2 = Effectiveness_1*Q_dot_max2;
    T_air_out = (Q_dot2/(C_air))+T_air_in;
    T_ethyl_out = T_ethyl_in-(Q_dot2/(C_ethyl));

    %Printing Exit Air Temp
    fprintf('Ambient Air Temperature In: %f C \n', T_air_in);
    fprintf('Exit Air Temperature: %f C \n', T_air_out);
    fprintf('Exit Ethyl Glycol Temperature: %f C \n', T_ethyl_out);
    disp('___________________________________________________');

end

disp("Done")





