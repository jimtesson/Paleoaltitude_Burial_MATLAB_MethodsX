function [Tint, s_Tint ] = Integration_time(C, s_C, P_SLHL, f, Lamb)

% Input:
% C (at/g): Measured cosmogenic nuclide concentration
% s_C (at/g): Uncertainty of the measured cosmogenic nuclide concentration
% P_SLHL (at/g/yr): 
% f: Stone Factor
% Lamb (yr-1): radioactive decay constant

% Output:
% Tint (yr): integration time of the computed elevation - average residence time of the
% analyzed cosmogneic nuclide in the mineral 
% s_Tint (yr): uncertainty attached to Tint


%clear all
%close all

%f = 4;
%P_SLHL = 1; %(at/g/yr)
%This production rate is artificially set to 1, but its value has actually
%no impact on the computed integration time

%Lamb = log(2)/(1.39e6);
%Lamb = log(2)/(0.7e6);
%Lamb = log(2)/(5300);

%C_mes = f.*P_SLHL/Lamb/1.0000001;

% Definition of the equation linking time and nuclide concentration C
%t = @(x) (-1/Lamb).*log(1-(Lamb.*x/(f.*P_SLHL)));

if C < f*P_SLHL/Lamb
    
%Tint = (1/C)* integral(t, 0, C)

Tint = (1/Lamb) + (1/Lamb) * log(1 - (Lamb*C/(f*P_SLHL))) * (f*P_SLHL/(Lamb*C) - 1);

s_Tint = abs((s_C/(Lamb*C)) * (log(1 - (Lamb*C/(f*P_SLHL)))*f*P_SLHL/(Lamb*C) + 1));
 
 

else
    
Tint = 1/Lamb;
s_Tint = NaN;

end
    
%figure
%fplot(t, [0 C_mes], 100);
%xlabel('[C(at/g)]'); ylabel('t(yr)');

% Exposure time
%t_exp = (-1/Lamb).*log(1-(Lamb.*C/(f.*P_SLHL)));
%r = Tint/t_exp;
