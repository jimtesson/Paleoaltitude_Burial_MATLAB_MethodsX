function [ Tbm, s_Tb, Ero, s_Ero, Tburial] = Burial26Al_10Be( X, D_X, Y, D_Y, Lat, Z, Parametres )


%This function computes burial ages from 26Al-10Be data, assuming steady
%state erosion before burial

%Inputs:
% X: 10Be concentration (at.g-1)
% sX: Error on 10Be concentration (at.g-1)
% Y: 26Al concentration (at.g-1)
% sY: Error on 26Al concentration (at.g-1)
% Lat: latitude
% Z: altitude (m)
% Parameters


%Outputs:
% Tburial_dis (Ma): vector with the probabilty distribution of burial ages 
% Tbm (Ma): mean burial age 
% s_Tb (Ma): uncertainty
% Ero (cm.yr-1): erosion before burial 
% s_ero (cm.yr-1): uncertainty

%10Be radioactive constant
Lambda_X = Parametres(1); 


% P_X: 10Be production rate (at.g-1.yr-1)
P_X = Parametres(2); 

%26Al radioactive constant
Lambda_Y = Parametres(3);

% P_Y: 26Al production rate (at.g-1.yr-1)
P_Y = Parametres(4);

Density = 2.7; Attenuation_length = 160; 
Mu = Density / Attenuation_length;


% f is computed using the function: "StoneFactL.m" no VDM assumed for this figure as it doesn'make much sense
gmr = -0.03417;
dtdz = 0.0065;
SLP = 1013.25;

%Pressure at elevation Z
Pk = SLP .* exp((gmr./dtdz) .* (log(288.15) - log(288.15 - (Z.*dtdz))));

%Calculation of Stone factor
f = StoneFactorL(Lat,Pk,SLP);

%Calculation of the burial age Tb with a Monte Carlo method

%Number of draws
n = 5e3;

    Tburial_dis = zeros(n,1);

    for i = 1:n
    
    Xrand = normrnd(X,D_X);
    Yrand = normrnd(Y,D_Y);
    
    if Xrand > 0 && Yrand >0
    
    %Definition of the burial age equation
    Burial = @(t) f.* ((P_X/Xrand)*exp(-Lambda_X*t) - (P_Y/Yrand).*exp(-Lambda_Y*t)) - Lambda_X + Lambda_Y;
    option = optimset('MaxFunEvals',10000, 'MaxIter', 10000); 
 
    %Calculation of a first guess value using the simplified formula of the
    %burial age
     
    Tini = 1/(Lambda_X-Lambda_Y)*log((Y/X)/(P_Y/P_X));
    
    Solution = fzero(@(t) Burial(t),Tini, option);
 
    Tburial_dis(i)=Solution;
    
    end
    
    end
    
Tbm = nanmean(Tburial_dis);
s_Tb = nanstd(Tburial_dis);

% Calculation of the solution without uncertainties for plot
    Burial = @(t) f.* ((P_X/X)*exp(-Lambda_X*t) - (P_Y/Y).*exp(-Lambda_Y*t)) - Lambda_X + Lambda_Y;
    option = optimset('MaxFunEvals',10000, 'MaxIter', 10000); 
    Solution = fzero(@(t) Burial(t),Tini,option);
    Tburial = Solution;

% Calculation of preburial erosion (mm/kyr)

Xcor = X*exp(Lambda_X*Tbm);
Ero = (1/Mu)* (f*P_X / Xcor - Lambda_X)*10000;

% Uncertainty on Erosion
s_Ero = Ero * D_X/X;