function [T1 T2 TS TF]=HX_performance(UA_hx, mdot_hx, HXfluid, q, Thx_in, mdot_boiler, UA_primary, mdot_primary)

T_boiler_guess = 40;	%guess for average boiler temperature
T_primary_guess = 50;	%guess for average primary circuit temperature  

%%%%%%%%%% Cp Tables taken from NIST %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mol_CO2 = 44.0095;		%g/mol - molecular weight of CO2
% Temperature (K) vs specific heat capacity (J/mol/K) table for CO2 taken from NIST
T_cp_CO2 = [...
			298	37.12;
			300	37.22;
			400	41.34];

mol_H2O = 18.0153;		%g/mol - molecular weight of water
% Temperature (K) vs specific heat capacity (J/mol/K) table for water taken from NIST
T_cp_H20 = [...
			298	75.38;
			300	75.35;
			400	76.74];	

mol_N2 = 28.0134;			%g/mol - molecular weight of nitrogen
% Temperature (K) vs specific heat capacity (J/mol/K) table for nitrogen taken from NIST
T_cp_N2 = [...
			200	29.11;
			300	29.12;
			400	29.25];	
			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			
% HX heat balance
if HXfluid == "water"
	cp_hx = 1000*(interp1(T_cp_H20(:,1),T_cp_H20(:,2),Thx_in+273.15))/mol_H2O;
elseif HXfluid == "air"
	cp_hx = 1000*(interp1(T_cp_N2(:,1),T_cp_N2(:,2),Thx_in+273.15))/mol_N2;
end
dT_hx = q/mdot_hx/cp_hx;
Thx_out = Thx_in + dT_hx;
Thx_out_old = Thx_out;
Thx_err = 10;

while Thx_err > 0.1
	Thx_ave = mean([Thx_in, Thx_out_old]);
	if HXfluid == "water"
		cp_hx = 1000*(interp1(T_cp_H20(:,1),T_cp_H20(:,2),Thx_ave+273.15))/mol_H2O;
	elseif HXfluid == "air"
		cp_hx = 1000*(interp1(T_cp_N2(:,1),T_cp_N2(:,2),Thx_ave+273.15))/mol_N2;
	end
	dT_hx = q/mdot_hx/cp_hx;
	Thx_out = Thx_in + dT_hx;
	Thx_err = abs(Thx_out-Thx_out_old);
	Thx_out_old = Thx_out;
end

% Boiler heat balance
cp_boiler = 1000*(interp1(T_cp_H20(:,1),T_cp_H20(:,2),T_boiler_guess+273.15))/mol_H2O;
dT_boiler = q/mdot_boiler/cp_boiler;

a_boiler = 1/mdot_boiler/cp_boiler - 1/mdot_hx/cp_hx;

eta_boiler = exp(-1*UA_hx*a_boiler);

TS = (Thx_in + dT_boiler - eta_boiler*Thx_out)/(1 - eta_boiler);
TF = TS - dT_boiler;

T_boiler_ave_old = mean([TS, TF]);
T_boiler_err = abs(T_boiler_ave_old - T_boiler_guess);

while T_boiler_err > 0.1
	cp_boiler = 1000*(interp1(T_cp_H20(:,1),T_cp_H20(:,2),T_boiler_ave_old+273.15))/mol_H2O;
	dT_boiler = q/mdot_boiler/cp_boiler;
	a_boiler = 1/mdot_boiler/cp_boiler - 1/mdot_hx/cp_hx;
	eta_boiler = exp(-1*UA_hx*a_boiler);
	TS = (Thx_in + dT_boiler - eta_boiler*Thx_out)/(1 - eta_boiler);
	TF = TS - dT_boiler;
	T_boiler_ave = mean([TS, TF]);
	T_boiler_err = abs(T_boiler_ave-T_boiler_ave_old);
	T_boiler_ave_old = T_boiler_ave;
end

% Primary side heat balance

cp_primary = 1000*(interp1(T_cp_CO2(:,1),T_cp_CO2(:,2),T_primary_guess+273.15))/mol_CO2;
dT_primary = q/mdot_primary/cp_primary;

a_primary = 1/mdot_primary/cp_primary - 1/mdot_boiler/cp_boiler;

eta_primary = exp(-1*UA_primary*a_primary);

T2 = (TF + dT_primary - eta_primary*TS)/(1 - eta_primary);
T1 = T2 - dT_primary;

T_primary_ave_old = mean([T1, T2]);
T_primary_err = abs(T_primary_ave_old - T_primary_guess);

while T_primary_err > 0.1
	cp_primary = 1000*(interp1(T_cp_CO2(:,1),T_cp_CO2(:,2),T_primary_ave_old+273.15))/mol_CO2;
	dT_primary = q/mdot_primary/cp_primary;
	a_primary = 1/mdot_primary/cp_primary - 1/mdot_boiler/cp_boiler;
	eta_primary = exp(-1*UA_primary*a_primary);
	T2 = (TF + dT_primary - eta_primary*TS)/(1 - eta_primary);
	T1 = T2 - dT_primary;
	T_primary_ave = mean([T1, T2]);
	T_primary_err = abs(T_primary_ave-T_primary_ave_old);
	T_primary_ave_old = T_primary_ave;
end

end