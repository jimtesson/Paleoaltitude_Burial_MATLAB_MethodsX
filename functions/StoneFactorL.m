function [ StoneFact, StFact_neut, StFact_muon ] = StoneFactorL( Latitude,P,SLP)
%--------Description------------------------------------------------------
% This function computes the scaling factors of Lal-Stone (Lal, 1991;
% Stone, 2000) for a given latitude, atmosheric pressure and with respect
% to a sea level pressure.
%
% Input : 
%          Latitude  : Sample latitude (positive decimal degree)
%          P         : Air pressure
%          SLP       : Sea level pressure (in hPa) 
%
% Output  :
%          StoneFact : Lal-Stone scaling factor
% 
% Code written by LCP Martin, PH Blard and J Lave
% Centre de Recherches Petrographiques et Geochimiques (CRPG-CNRS), France
% blard@crpg.cnrs-nancy.fr
% Program desciprtion provided in Martin et al. (2017)
% 
% Copyright 2015, CNRS-Universite de Lorraine
% All rights reserved
%
% This file is part of the CREp program.
% CREp is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. See <http://www.gnu.org/licenses/>
%-------------------------------------------------------------------------

% Parameters

% The Spallation vs Muonic production ratio is the global value taken from
% Braucher et al.,  (2013)
neutron = 0.9937;
muon = 1 - neutron;

% Scaling factor calculation
St =    [0,	31.8518,	250.3193,	-0.083393,	7.4260E-05,	-2.2397E-08,	0.5870;
        10,	34.3699,	258.4759,	-0.089807,	7.9457E-05,	-2.3697E-08,	0.6000;
        20,	40.3153,	308.9894,	-0.106248,	9.4508E-05,	-2.8234E-08,	0.6780;
        30,	42.0983,	512.6857,	-0.120551,	1.1752E-04,	-3.8809E-08,	0.8330;
        40,	56.7333,	649.1343,	-0.160859,	1.5463E-04,	-5.0330E-08,	0.9330;
        50,	69.0720,	832.4566,	-0.199252,	1.9391E-04,	-6.3653E-08,	1.0000;
        60,	71.8733,	863.1927,	-0.207069,	2.0127E-04,	-6.6043E-08,	1.0000];

if Latitude<60;
    infer=floor(Latitude/10)+1;
    super=infer+1;
    dec=10*infer-Latitude;
else
    infer=6;
    super=7;
    dec=0;
end
a = (dec*St(infer,2)+(10-dec)*St(super,2))/10;
b = (dec*St(infer,3)+(10-dec)*St(super,3))/10;
c = (dec*St(infer,4)+(10-dec)*St(super,4))/10;
d = (dec*St(infer,5)+(10-dec)*St(super,5))/10;
e = (dec*St(infer,6)+(10-dec)*St(super,6))/10;
M = (dec*St(infer,7)+(10-dec)*St(super,7))/10;

StoneFact=((neutron*(a+b.*exp(-P/150)+c.*P+d.*P.^2+e.*P.^3))+muon.*M.*exp((SLP-P)/242));

StFact_neut = a+b.*exp(-P/150)+c.*P+d.*P.^2+e.*P.^3;

StFact_muon = M.*exp((SLP-P)/242);

end
