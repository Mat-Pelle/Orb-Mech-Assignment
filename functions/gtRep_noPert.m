function [a] = gtRep_noPert(k,m,gt_we,Earth_mu)

% rep_groundtrack computes the value of the semi-major axis to obtain repeating groundtrack

%PROTOTYPE: 
%     [a,T] = rep_groundtrack(k,m,we,mu)
%
% INPUT:
%     k [1]          Revolution of the satelite
%     m [1]          Rotations of the planet
%     gt_we [1]      Angular velocity of the Earth [rad/s]
%     Earth_mu [1]   Gravitational constant of the Earth [km^3/s^2]
%     
% OUTPUT:
%    gtRep_a [1]        Semimajor axis [km]
%
% CONTRIBUTORS
%       Bertolini Edoardo
%       Busi Silvia
%       Muylle Julia
%       Pellegrini Matias
%
% VERSIONS
%
% 30/11/2020: First Version

%Earth’s rotational period 
Te = 2*pi/gt_we;

%s/c’s orbital period
T  = Te*m/k;

%Semi-major axis for repeating groundtrack
a  = ((T/(2*pi))^2*Earth_mu)^(1/3);

