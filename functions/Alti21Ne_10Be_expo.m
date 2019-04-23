function [ Z, Z_min, Z_max, Tint, s_Tint,  StoneFactor ] = Alti21Ne_10Be_expo( X, D_X, Y, D_Y, Lat, Parametres )

%This function computes elevations from 21Ne-10Be data, assuming no
%erosion - continuous exposure

%Inputs: 
% X: 21Ne concentration (at.g-1)
% sX: Error on 21Ne concentration (at.g-1)
% Y: 10Be concentration (at.g-1)
% sY: Error on 10Be concentration (at.g-1)
% Lat: latitude
% Parameters

%Outputs:
% Z: mean altitude (m)
% Zmin: min altitude (m)
% Zmax: max altitude (m)
% Tint: average residence time of the shortest half life nuclide 10Be (yr)
% s_Tint: uncertainty on Tint (yr)
% StoneFactor: facteur de Stone

%%% Stone factor computation %%%

P_X = Parametres(2);
Lambda_Y = Parametres(3) ;
P_Y = Parametres(4);

A = Lambda_Y * X / P_X;
B = Lambda_Y * Y / P_Y;
D_A = Lambda_Y * D_X / P_X;
D_B = Lambda_Y * D_Y / P_Y;


if A > B + sqrt(D_A^2 + D_B^2) %Case 1 : min, max and mean altitude can be computed
    StoneFactor_max = fzero(@(f) f*(1-exp(-A/f))-B - sqrt(exp(-2*A/f) * D_A^2 + D_B^2), [0, 20000]);
    StoneFactor_min = fzero(@(f) f*(1-exp(-A/f))-B + sqrt(exp(-2*A/f) * D_A^2 + D_B^2), [0, 20000]);
    StoneFactor = fzero(@(f) f*(1-exp(-A/f))-B, [0,20000]);

elseif A < B - sqrt(D_A^2 + D_B^2) %Case 4, not enough Ne, no possible computation
    Z_min = NaN; Z_max = NaN; Z = NaN; t_max = NaN; D_t_max = NaN; StoneFactor = NaN;
    return

elseif A < B %Case 3: computation of minimal elevation
    StoneFactor_min = fzero(@(f) f*(1-exp(-A/f))-B + sqrt(exp(-2*A/f) * D_A^2 + D_B^2), [0, 20000]);
    StoneFactor = NaN; StoneFactor_max = NaN;

else %Case 2, minimal and mean elevation computed
    StoneFactor_min = fzero(@(f) f*(1-exp(-A/f))-B + sqrt(exp(-2*A/f) * D_A^2 + D_B^2), [0, 20000]);
    StoneFactor = fzero(@(f) f*(1-exp(-A/f))-B, [StoneFactor_min, 20000]);
    StoneFactor_max = NaN;
end

%%% Computation of integration time = average residence time of 10Be
%%% (shortest half-life nuclide)

if not(isnan(StoneFactor))
    
    [Tint, s_Tint ] = Integration_time(Y, D_Y, P_Y, StoneFactor, Lambda_Y); 
    % Tint (yr): integration time - average residence time of the analyzed cosmogenic nuclide in the mineral 
else
    Tint = NaN;
    s_Tint = NaN;
end


%%% Computation of elevation %%%

if isnan(Tint)
    Tint > 9000000;
    t_c = 9000000;
else
    t_c = Tint;
end

% Definition of a function that accounts fo the sample latitude and for 
% time variations in Earth magnetic field

fun = @(x) StoneFactCT(t_c/1000, Lat, x);

    if isnan(StoneFactor_min)
    Z_min = NaN;
    else
    Z_min = int16(fzero(@(x) fun(x)- StoneFactor_min, [-5000, 20000]));
    end

    if isnan(StoneFactor)
    Z = NaN;
    else
    Z = int16(fzero(@(x) fun(x) - StoneFactor, [-5000, 20000]));
    end

    if isnan(StoneFactor_max)
    Z_max = NaN;
    else
    Z_max =  int16(fzero(@(x) fun(x) - StoneFactor_max, [-5000, 20000]));
    end

end