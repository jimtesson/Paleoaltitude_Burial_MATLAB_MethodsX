function [ YMoyVect] = MoyenneIntegrV2( XVect,YVect )
%--------Description------------------------------------------------------
% This function computes time average values of the YVect=f(XVect)
% function for each time step of the XVect vector.
%
% Input : 
%           XVect    : time vector
%           YVect    : corresponding signal
%
% Output  : 
%           YMoyVect : time average values of the input data (for each time step)
%
% Code written by LCP Martin, PH Blard and J Lavé
% Centre de Recherches Pétrographiques et Géochimiques (CRPG-CNRS), France
% blard@crpg.cnrs-nancy.fr
% Program desciprtion provided in Martin et al., (In Prep)
% 
% Copyright 2015, CNRS-Université de Lorraine
% All rights reserved
%
% This file is part of the CREp program.
% CREp is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. See <http://www.gnu.org/licenses/>
%-------------------------------------------------------------------------

YIntegr=zeros(1,length(YVect));
YMoyVect=zeros(1,length(YVect));

for k=1:(length(YVect)-1);
    DeltaX=XVect(k+1)-XVect(k);
    Trapeze=0.5*(YVect(k+1)+YVect(k));
    YIntegr(1,k+1)=YIntegr(1,k)+DeltaX*Trapeze;
    YMoyVect(k+1)=YIntegr(k+1)/(XVect(k+1)-XVect(1));
end

end