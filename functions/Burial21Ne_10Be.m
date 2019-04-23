function [ Tb, s_Tb, Ero, s_Ero ] = Burial21Ne_10Be( X, D_X, Y, D_Y, Lat, Z, Parametres )


%This function computes burial ages from 21Ne-10Be data, assuming steady state erosion

%Inputs: 
% X: 21Ne concentration (at.g-1)
% sX: Error on 21Ne concentration (at.g-1)
% Y: 10Be concentration (at.g-1)
% sY: Error on 10Be concentration (at.g-1)
% Lat: latitude
% Z: altitude (m)
% Parameters 

%Outputs:
% Tb: burial age (Ma)
% s_Tb: uncertainty (Ma)
% Ero: erosion before burial (mm.kyr-1)
% s_ero: uncertainty (mm.kyr-1)

%10Be radioactive constant
Lambda_Y = Parametres(3);

% P_X: 21Ne production rate (at.g-1.yr-1)
P_X = Parametres(2); 

% P_Y: 10Be production rate (at.g-1.yr-1)
P_Y = Parametres(4);

Density = 2.7; Attenuation_length = 160; 
Mu = Density / Attenuation_length;


% f is computed using the function: "StoneFactL.m" no VDM assumed for this figure as it doesn'make much sense
gmr = -0.03417;
dtdz = 0.0065;
SLP = 1013.25;
%Pressure at elevation z
Pk = SLP .* exp((gmr./dtdz) .* (log(288.15) - log(288.15 - (Z.*dtdz))));

%Calculation of Stone factor
f = StoneFactorL(Lat,Pk,SLP);

%Calculation of the burial age Tb
Tb = (-1/Lambda_Y)*log( (Y/P_Y)* ( (P_X/X) + (Lambda_Y/f)) );


%Calculation of uncertainty
s_Tb= (1/Lambda_Y) * ( ( D_Y / Y )^2 + ( (D_X * P_X) / (X^2*(P_X/X + (Lambda_Y / f)))   )^2  )^(1/2);

%Calculation of preburial erosion (mm/kyr)
Ero = f*P_X / (Mu*X)*10000;

% Uncertainty on Erosion
s_Ero = Ero * D_X/X;