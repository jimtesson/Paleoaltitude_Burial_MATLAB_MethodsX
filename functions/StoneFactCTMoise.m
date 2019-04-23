function [FactStoneCorrige] = StoneFactCT(AgeFact,Latitude,Altitude)

%--------Description------------------------------------------------------
% This functions computes Lal-Stone scaling factors (Lal, 1991; Stone 2000)
% with the geomagnetic time correction proposed by Nishiizumi et al. (1989)
% as described in Balco et al. (2008). It is implemented for two atmosphere
% models. Different Virtual Dipolar Moment reconstruction can be input if
% they fit the required format.
%
% Input : 
%        AgeFact   : Non corrected age of the calibration object (in ka)                
%        Latitude  : Latitude of the calibration object (in decimal degrees, <0 if S)
%        Altitude  : Elevation (in masl)
%      
% Output  : 
%        FactStoneCorrige : Lal-Stone scaling factor accounting for geomagnetic correction
%
% IMPORTANT : Requires Matlab 2009 or any more recent versions.

% Code written by LCP Martin, PH Blard and J Lave
% Centre de Recherches Petrographiques et Geochimiques (CRPG-CNRS), France
% email: blard@crpg.cnrs-nancy.fr
% Program description in Martin et al., (In Prep)
% 
% Copyright 2015, CNRS-Universite de Lorraine
% All rights reserved
%
% This file is part of the CREp program (Martin et al., QG, 2017).
% CREp is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. See <http://www.gnu.org/licenses/>
%-------------------------------------------------------------------------


%------------------IMPORTATION VDM----------------------------------------

load GMDB

PaleoMag=GMDB.Musch;
VecAgeReel1=PaleoMag(1,:);
PaleoVDM=PaleoMag(2,:);

% Block if the age is too old
if AgeFact>VecAgeReel1(end);
    FactStoneCorrige=NaN;
    return
end

[~,Dates]=size(VecAgeReel1);
if VecAgeReel1(end)>(1.5*AgeFact); % Chope if too long
    VecIndice=find(VecAgeReel1>(1.5*AgeFact));
    VecAgeReel1=VecAgeReel1(1:VecIndice(1));
    PaleoVDM=PaleoVDM(1:VecIndice(1));
    Dates=length(VecAgeReel1);
end

%------------------PALEOMAGNETIC CORRECTIONS-------------------------------

% Paleomag correction is applied followong the Nishiizumi et al (1989)
% formula

CosLambdaM=(PaleoVDM/PaleoVDM(1)).^0.25*cos(pi*Latitude/180);

% For cos(LambdaM) > 1, we assume that  LambdaM = 0.01.

LambdaM=zeros(1,Dates);
for k=1:length(CosLambdaM);
    if CosLambdaM(k)<=1;
        LambdaM(1,k)=acos(CosLambdaM(k))*180/pi;
    else
        LambdaM(1,k)=0.01;
    end
end

%--------------------SPATIAL CORRECTION-------------------------------------
% Computation of the Stone factor for each paleolatitude LambdaM

StoneFact=zeros(1,Dates);
gmr = -0.03417;
dtdz = 0.0065;
P=1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (Altitude.*dtdz)) ) );
for k=1:length(LambdaM);
    StoneFact(k)=StoneFactorL(LambdaM(k),P,1013.25);
end


%------------------TIME INTEGRATION OF STONE FACTORS-------------------
StoneFactMoy=MoyenneIntegrV2(VecAgeReel1,StoneFact);

%---Interpolation
if AgeFact<=VecAgeReel1(end-1);
    NbPoints=100;
    Mil=floor(NbPoints/2);
    VectT=linspace(AgeFact*0.90,AgeFact*1.10,NbPoints);
    VectSt=interp1(VecAgeReel1,StoneFactMoy,VectT,'spline');
    %Picking of correct value
    VecTFin=[VectT(Mil-1) VectT(Mil) VectT(Mil+1)];
    VecTFin=abs(VecTFin-AgeFact);
    [~,Indice]=min(VecTFin);
    Indice=Mil+Indice-2;
    %Output sorting
    FactStoneCorrige=VectSt(Indice);
else
    FactStoneCorrige=(AgeFact-VecAgeReel1(end-1))*((StoneFactMoy(end)-StoneFactMoy(end-1))/(VecAgeReel1(end)-VecAgeReel1(end-1)))+StoneFactMoy(end-1);
end

end











