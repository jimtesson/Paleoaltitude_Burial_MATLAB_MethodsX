function [ u, v ] = Ellipse_incertitude_lin( Num_x, Den_x, Num_y, Den_y, cor_mat, confidence )
%UNTITLED Summary of this function goes here
%   Draw and return coordinates of incertitude ellipse of correlated ratio of variable Num_x/Den_x in x and Num_y/Den_y in y with given confidence
%   Correlation matrix defined symetric and like this :
%   
%           Num_x  Den_x  Num_y  Den_y      where r_ij is the correlation between variable i and j
%   Num_x    r11    r12    r13    r14       r_ij = r_ji
%   Den_x    r21    r22    r23    r24       r_ij = 1 if i and j are the same variable
%   Num_y    r31    r32    r33    r34       r_ij = 0 if variable i and j are uncorrelated or i or j is a constant
%   Den_y    r41    r42    r43    r44       for example r23 is correlation between Den_x and Num_y
%   
%   Confidence in range ]0,1[
%   If no confidence given or confidence equal 0, ellipse is first order (k=1 ; P = 0.39)

% Ameliorer : prendre en compte correlation dans le calcul de D_X et D_Y

if nargin < 6
   confidence = 0;
end

% Data loading
Value_Num_x = Num_x(1); Value_Den_x = Den_x(1); Value_Num_y = Num_y(1); Value_Den_y = Den_y(1);
D_Num_x = Num_x(2); D_Den_x = Den_x(2); D_Num_y = Num_y(2); D_Den_y = Den_y(2);

% coefficient of variation
v_Num_x = Num_x(2)/Num_x(1);
v_Num_y = Num_y(2)/Num_y(1);
v_Den_x = Den_x(2)/Den_x(1);
v_Den_y = Den_y(2)/Den_y(1);

% correlation calculation
% https://en.wikipedia.org/wiki/Spurious_correlation_of_ratios - Pearson 1897
% Warning ! Num_x = x1 ; Den_x = x3 ; Num_y = x2 ; Den_y = x4 : r_i2 <-> r_i3 and r_2j <-> r_3j
% This equation is an approximation of the real correlation
cor = (cor_mat(1,3) * v_Num_x * v_Num_y - cor_mat(1,4) * v_Num_x * v_Den_y - cor_mat(3,2) * v_Num_y * v_Den_x + cor_mat(2,4) * v_Den_x * v_Den_y) / sqrt(v_Num_x^2 + v_Den_x^2 - 2*cor_mat(1,2)*v_Num_x*v_Den_x) / sqrt(v_Num_y^2 + v_Den_y^2 - 2*cor_mat(3,4)*v_Num_y*v_Den_y);

% ratio error
% Classis propagation error formula (sum of squared partial derivative error power 1/2)
% D_X = 1 / Value_Den_x^2 * sqrt(Value_Den_x^2 * D_Num_x^2 + D_Den_x^2 * Value_Num_x^2)
% D_Y = 1 / Value_Den_y^2 * sqrt(Value_Den_y^2 * D_Num_y^2 + D_Den_y^2 * Value_Num_y^2)
% Formula with correlation
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty
D_X = 1 / Value_Den_x * sqrt(D_Num_x^2 + Value_Num_x^2/Value_Den_x^2*D_Den_x^2 - 2*cor_mat(1,2)*Value_Num_x/Value_Den_x*D_Den_x*D_Num_x);
D_Y = 1 / Value_Den_y * sqrt(D_Num_y^2 + Value_Num_y^2/Value_Den_y^2*D_Den_y^2 - 2*cor_mat(3,4)*Value_Num_y/Value_Den_y*D_Den_y*D_Num_y);

% ellipse parameter
% http://web.cortial.net/optim/annexe2.html
cen = [Value_Num_x / Value_Den_x ; Value_Num_y / Value_Den_y];
alpha = 0.5 * atan(2 * cor * D_X * D_Y / (D_X^2 - D_Y^2));
half_axis_x = sqrt(D_X^2*cos(alpha)^2 + cor*D_X*D_Y*sin(2*alpha) + D_Y^2*sin(alpha)^2);
half_axis_y = sqrt(D_X^2*sin(alpha)^2 - cor*D_X*D_Y*sin(2*alpha) + D_Y^2*cos(alpha)^2);

% Ellipse order calculation
% http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/#Source_Code
% and http://web.cortial.net/optim/annexe2.html
if confidence == 0
   k = 1;
else
   k = sqrt(2 * log(1/(1-confidence)));
end

% draw ellipse in a linear plot
% calculated with ellipse parametric equation 
t = 0 : pi/100 : 2*pi;
cost = cos(t);
sint = sin(t);

u = cen(1) + k*half_axis_x * cos(alpha) * cost - k*half_axis_y * sin(alpha) * sint;
v = cen(2) + k*half_axis_x * cost * sin(alpha) + k*half_axis_y * sint * cos(alpha);

v(v < 0) = 1E-9; % Trunquate elipse for negative values.

end
