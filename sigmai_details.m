%%%%%% About this code
% sigmai_details.m is a MATLAB function which calculates  two dimensional charge density (sigmai), spontaneus polarization (P_SP_tot), piezo-electric polarization (P_PE_tot),
% lattice constant (aAlInGaN), bandgap (Eg_AlInGaN),and conduction band discontinuty (del_E_C) based on AL and In mole fractions.
% This code also gives different options for choosing the bowing parameters used for calculation of AlInGaN bandgap as discussed in [1].

% [1] H. R. Mojaver, F. Manouchehri, and P. Valizadeh, “Theoretical evaluation of two dimensional electron gas characteristics of quaternary AlxInyGa1-x-yN/GaN hetero-junctions,” J. Appl. Phys., vol. 119, no. 15, pp. 154502-1–154502-7, Apr. 2016.
%%%%%%%

function [sigmai,P_SP_tot,P_PE_tot,aAlInGaN,Eg_AlInGaN,del_E_C] = sigmai_details(X,Y)

%Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 1.38e-23;                                                              %J/K
T = 300;
h = 6.625e-34;                                                             %Js
h_bar = h/(2*pi);                                                     
e = 1.6e-19;                                                               %C

%Varaiables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = 1-X-Y;

%Lattice Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aGaN = 3.189e-10;                                                    
aAlN = 3.112e-10;    
aInN = 3.538e-10;
aAlInGaN = X*aAlN + Y*aInN + Z*aGaN;

%BandGap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ref 21 of [1]
 Eg_GaN = 3.45;
 Eg_AlN = 6.21;
 Eg_InN = 0.68;
 b_InGaN=1.72;
 b_AlGaN=0.9;
 b_AlInN=6.43/(1+1.21*X^2);
 Eg_AlInGaN = X*Eg_AlN + Y*Eg_InN +(1-X-Y)*Eg_GaN -b_AlGaN*X*(1-X)-b_InGaN*Y*(1-Y)-(b_AlInN-b_AlGaN-b_InGaN)*X*Y;
 
% %Ref 22 of [1]
% b_AlGaN=1.4;
% b_InGaN=1.37;
% b_InAlN=2.914;
% C=2.59;
% Eg_AlInGaN= X*Eg_AlN+(1-X-Y)*Eg_GaN+Y*Eg_InN-b_InGaN*(1-X-Y)*(X+Y)-b_InAlN*X*(1-X)-b_AlGaN*X*(1-X-Y)+(b_InGaN+b_InAlN)*X*(1-X-Y)-C*X*(1-X-Y)*Y;

% % 
% %Ref 23 of [1]
% b_InAlN=3.7;
% b_InGaN=1.4;
% b_GaAlN=0.71;
% Eg_AlInGaN=(1-X-Y)*Eg_GaN + X*Eg_AlN + Y*Eg_InN-(1-X)*Y*Z*b_InGaN-(1-Z)*X*Y*b_InAlN-(1-Y)*X*Z*b_GaAlN;

%Del_E_C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VBO_AlN_GaN = 0.34;
% VBO_InN_GaN = -0.62;
% 
% VBO_AlInGaN = X*VBO_AlN_GaN+Y*VBO_InN_GaN;
% CBO_AlInGaN = (Eg_AlInGaN - Eg_GaN)-VBO_AlInGaN;
% del_E_C = CBO_AlInGaN;
del_E_C=0.7*(Eg_AlInGaN-Eg_GaN);

%Polarization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a0XY = (aGaN - aAlInGaN)/(aAlInGaN);            %basal strain between AlInGaN and GaN

%%%%Spontaneous%%%%
P_SP_GaN = -0.0339;  %[Piprek, Bernardini 2007][Bernardini 2002]                                                         %C/m2
P_SP_AlN = -0.0898;                                                     
P_SP_InN = -0.0413;  
b_AlGaN = 0.0191;   %[Piprek, Bernardini 2007]
b_InGaN = 0.0378;
b_AlInN = 0.0709;


P_SP_AlInGaN = X*P_SP_AlN + Y*P_SP_InN + Z*P_SP_GaN+ b_AlGaN*X*Z + b_InGaN*Y*Z + b_AlInN*X*Y;

%%%%Piezo-Electric%%%%
if a0XY>0
    P_PE_AlN = -1.808 * a0XY - 7.888 * a0XY^2;
else
    P_PE_AlN = -1.808 * a0XY + 5.624 * a0XY^2;
end
P_PE_GaN = -0.918 * a0XY + 9.541 * a0XY^2;
P_PE_InN = -1.373 * a0XY + 7.559 * a0XY^2;

P_PE_AlInGaN=(X*P_PE_AlN+Y*P_PE_InN+Z*P_PE_GaN);

%%%%Total%%%%
P_PE_tot=P_PE_GaN-P_PE_AlInGaN;
P_SP_tot = P_SP_GaN - P_SP_AlInGaN;
sigmai = P_PE_tot + P_SP_tot;


end


