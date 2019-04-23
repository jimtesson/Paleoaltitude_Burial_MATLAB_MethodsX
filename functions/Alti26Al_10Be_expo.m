function [ Z, Z_min, Z_max, Tint, s_Tint,  StoneFactor ] = Alti26Al_10Be_expo( X, D_X, Y, D_Y, Lat, Parametres )

%This function computes elevations from 10Be-21Ne data, assuming no
%erosion - continuous exposure

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
% Tint: average residence time of the shortest half life nuclide 26Al (yr)
% s_Tint: uncertainty on Tint (yr)
% StoneFactor: facteur de Stone

%%% Stone factor computation %%%

Lambda_X = Parametres(1);
P_X = Parametres(2);
Lambda_Y = Parametres(3);
P_Y = Parametres(4);


r = Lambda_X/Lambda_Y;

A = Lambda_Y * Y / P_Y;
B = Lambda_X * X / P_X;
D_A = Lambda_Y * D_Y / P_Y;
D_B = Lambda_X * D_X / P_X;


g = @(x) -B - x .* (-A ./ x + 1) .^ r + x;
D_g = @(x) sqrt (D_A .^ 2 .* r .^ 2 .* (-A ./ x + 1) .^ (2 * r) ./ (-A ./ x + 1) .^ 2 + D_B .^ 2);


% Computation of max value for  Stone factor

if r*A-B+sqrt((r*D_A)^2+D_B^2)<0
    StoneFactor_max = fzero(@(x) g(x)+D_g(x), [A+0.0001, 20000]);
   if StoneFactor_max < max(A,B)
      StoneFactor_max = NaN;
    end
else
    StoneFactor_max = NaN;
end



% Computation of Stone factor

if A > B && r*A < B
    StoneFactor = fzero(@(x) g(x), [A, 20000]);
else
    StoneFactor = NaN ;
end

% Computation of min value for Stone factor

if not(isnan(StoneFactor))
    if r*A-B-sqrt((r*D_A)^2+D_B^2)<0
        D_h_min = @(f) -A .* r .* (-A ./ f + 1) .^ r ./ (f .* (-A ./ f + 1)) - (-A ./ f + 1) .^ r + 1 - (A.* D_A .^ 2 .* r .^ 3 .* (-A ./ f + 1) .^ (2 * r) ./ (f .^ 2 .* (-A ./ f + 1) .^ 3) - A .* D_A .^ 2 .* r.^ 2 .* (-A ./ f + 1) .^ (2 * r) ./ (f .^ 2 .* (-A ./ f + 1) .^ 3)) ./ sqrt (D_A .^ 2 .* r .^ 2 .* (-A ./ f + 1) .^ (2 * r) ./ (-A ./ f + 1) .^ 2 + D_B .^ 2);
        x_inv = fzero(D_h_min, [A+0.0001, 20000]);
        if g(x_inv)-D_g(x_inv) > 0
            StoneFactor_min = fzero(@(x) g(x)-D_g(x), [x_inv, 2000]);
        else
            StoneFactor_min = max(A,B);
        end
    else
        StoneFactor_min = max(A,B);
    end
else
    %Sample in forbidden zone
    StoneFactor_min = NaN;
end


%Computation of the Stone factor uncertainty
if isnan(StoneFactor_max) && not(isnan(StoneFactor)) && not(isnan(StoneFactor_min))
    D_f = (StoneFactor - StoneFactor_min);
elseif isnan(StoneFactor_min) && not(isnan(StoneFactor)) && not(isnan(StoneFactor_max))
    D_f = (StoneFactor_max - StoneFactor);
elseif not(isnan(StoneFactor))
    D_f = (StoneFactor_max - StoneFactor_min)/2;
else
    D_f = NaN;
end

%%% Computation of integration time = average residence time of 26Al
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


