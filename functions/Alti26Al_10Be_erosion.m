function [ Z, Z_min, Z_max, Erosion, s_Erosion, Tint, s_Tint, StoneFactor ] = Alti26Al_10Be_erosion( X, D_X, Y, D_Y, Lat, Parametres )

%This function computes elevations from 26Al-10Be data, assuming steady state erosion

%Inputs: 
% X: 10Be concentration (at.g-1)
% sX: Error on 10Be concentration (at.g-1)
% Y: 26Al concentration (at.g-1)
% sY: Error on 26Al concentration (at.g-1)
% Lat: latitude
% Parameters

%Outputs:
% Z: mean altitude (m)
% Zmin: min altitude (m)
% Zmax: max altitude (m)
% Erosion (m/Ma)
% s_Erosion: uncertainty of Erosion (m/Ma)
% Tint: average residence time of the shortest half life nuclide 10Be (yr)
% s_Tint: uncertainty on Tint (yr)
% StoneFactor: facteur de Stone


Lambda_X = Parametres(1);
P_X = Parametres(2);
Lambda_Y = Parametres(3); 
P_Y = Parametres(4);
Mu = Parametres(5);


%%% Stone factor computation %%%

StoneFactor = ((Lambda_X - Lambda_Y) * X * Y) / (P_X * Y - P_Y * X);
D_StoneFactor = (Lambda_Y - Lambda_X) / (P_Y*X - P_X*Y)^2 * sqrt((D_Y*X^2*P_Y)^2 + (P_X*Y^2*D_X)^2);

if StoneFactor < 0
    if StoneFactor + D_StoneFactor > 0
        StoneFactor_min = max(X*Lambda_X/P_X, Y*Lambda_Y/P_Y);
    else
        StoneFactor_min = NaN;
    end
    StoneFactor_max = NaN;
else
    StoneFactor_min = max(X*Lambda_X/P_X, Y*Lambda_Y/P_Y);
    StoneFactor_min = max(StoneFactor_min, StoneFactor - D_StoneFactor);
    StoneFactor_max = StoneFactor + D_StoneFactor;
end

if Y/X < P_Y*Lambda_X / (P_X * Lambda_Y)
    StoneFactor_min = NaN;
    StoneFactor = NaN;
    if StoneFactor_max < max(X*Lambda_X/P_X, Y*Lambda_Y/P_Y)
        StoneFactor_max = NaN;
    end
end



%%% Computation of erosion from 10Be (m/Ma) %%
Erosion = 10000 * (1/Mu)*(StoneFactor*P_Y/Y - Lambda_Y );
s_Erosion = Erosion * D_Y/Y;



%%% Computation of integration time = average residence time of 26Al
%%% (shortest half-life nuclide)

if not(isnan(StoneFactor))
    
    [Tint, s_Tint ] = Integration_time(Y, D_Y, P_Y, StoneFactor, Lambda_Y); 
    % Tint (yr): integration time - average residence time of the analyzed cosmogenic nuclide in the mineral 
else
    Tint = NaN;
    s_Tint = NaN;
end


if StoneFactor < 0
    StoneFactor = NaN;
    Erosion = NaN; 
    Tint = NaN;
    s_Tint = NaN;
end


%%% Computation of elevation %%%

% Geomagnetic database time set
if isnan(Tint)
    Tint > 9000000;
    t_c = 9000000;
else
    t_c = Tint;
end

% Solving Elevation Z using fzero function

fun = @(x) StoneFactCTMoise(t_c/1000, Lat, x);

if isnan(StoneFactor_min)
    Z_min = NaN;
else
    Z_min = int16(fzero(@(x) fun(x) - StoneFactor_min, [-5000, 20000]));
end

if isnan(StoneFactor)
    Z = NaN;
else
    Z = int16(fzero(@(x) fun(x) - StoneFactor, [-5000, 20000]));
end

if isnan(StoneFactor_max)
    Z_max = NaN;
elseif StoneFactor_max > 250
    Z_max = Inf;
else
    Z_max = int16(fzero(@(x) fun(x) - StoneFactor_max, [-5000, 20000]));
end

end
