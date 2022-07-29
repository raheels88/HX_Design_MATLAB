function [Th2 Tc1 Tc2 UA_hx mdot_hx]=HX_sizer(P_r,Th1,Thx_in,Thx_out,hx_fluid,F)
% function [T_outH T_outC]=primary_heatex_out(T_inH,T_inC,mdot_H,mdot_C,Cp_h,Cp_c,h)
%  - Summary: 
%    Calculates primary circuit outlet temperatures for coolant gas and boiler
%
% Outputs: Th2 		= cooled reactor gas tempereature i.e. T1 (degC)
%	       Tc1		= warmed boiler water out temperature (degC)
%          Tc2 		= cool boiler water in temperature (degC)
%          UA_hx	= UA for the new heat exchanger (W/degC)
%          mdot_hx	= mass flow rate for the new heat exchanger (kg/s)
%%%%%%%%%%%
% Inputs:  P_r		= reactor gas pressure (barg) - should be between 0.1 and 0.5
%          Th1		= hot reactor gas temperature i.e. T2 (degC)
%          Thx_in	= heat exchanger fluid inlet temperature (degC)
%          Thx_out	= heat exchanger fluid allowable outlet temperature (degC)
%          hx_fluid	= heat exchanger fluid (either "air" or "water")
%          F		= LMTD correction factor - depends on type of HX (typically 0.9)
%

q = 3e6;			%W    - overall decay heat
UA_primary = 9e5; 	%W/K  - UA for primary circuit heat removal
mdot_c = 24;		%kg/s - boiler mass flow rate from 2 quads
Tc_guess = 40;		%degC - inital guess for average boiler circuit temperature


% Pressure (barg) vs mdot_h (kg/s) table taken from /tasks/TSK-4961/3.computer_runs/item_075
% results files in run066 (0barg) and run022 (0.5barg)
P_mdot_h = [... 
				0	119.05; 
				0.5	191.15];

mdot_h = interp1(P_mdot_h(:,1),P_mdot_h(:,2),P_r);	%kg/s - linear interpolation of reactor mass flow rate

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

cp_h = 1000*(interp1(T_cp_CO2(:,1),T_cp_CO2(:,2),Th1+273.15))/mol_CO2;	%interpolation of specific heat capacity and conversion to J/kg/degC
dTh = q/mdot_h/cp_h;
Th2 = Th1 - dTh;
Th2_old = Th2;
Th_err = 10;

while Th_err > 0.1
	Tha = mean([Th1, Th2_old]);
	cp_h = 1000*(interp1(T_cp_CO2(:,1),T_cp_CO2(:,2),Tha+273.15))/mol_CO2;
	dTh = q/mdot_h/cp_h;
	Th2 = Th1 - dTh;
	Th_err = abs(Th2-Th2_old);
	Th2_old = Th2;
end

cp_c = 1000*(interp1(T_cp_H20(:,1),T_cp_H20(:,2),Tc_guess+273.15))/mol_H2O;
dTc = q/mdot_c/cp_c;

a = 1/mdot_h/cp_h - 1/mdot_c/cp_c;

eta = exp(-1*UA_primary*a);

Tc1 = (Th2 + dTc - eta*Th1)/(1 - eta);
Tc2 = Tc1 - dTc;

Tca_old = mean([Tc1, Tc2]);
Tc_err = abs(Tca_old - Tc_guess);

while Tc_err > 0.1
	cp_c = 1000*(interp1(T_cp_H20(:,1),T_cp_H20(:,2),Tca_old+273.15))/mol_H2O;
	dTc = q/mdot_c/cp_c;
	a = 1/mdot_h/cp_h - 1/mdot_c/cp_c;
	eta = exp(-1*UA_primary*a);
	Tc1 = (Th2 + dTc - eta*Th1)/(1 - eta);
	Tc2 = Tc1 - dTc;
	Tca = mean([Tc1, Tc2]);
	Tc_err = abs(Tca-Tca_old);
	Tca_old = Tca;
end

%varargout{1}=[Th1, Th2, Tc1, Tc2];

[UA_hx mdot_hx] = hx_params(Tc1,Tc2,Thx_out,Thx_in,hx_fluid,q,F);

end