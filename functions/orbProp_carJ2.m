function dy = orbProp_carJ2 (t,y,Earth_mu,Re,J2)

%odefun_J2 gives the derivative of the state vector using the dynamic equation of 2BP with J2
%perturbation
%
%PROTOTYPE: 
%     dy = odefun_J2 (t,y,mu,J2)
% 
% INPUT:
%     t           Time [s]
%     y [6x1]     State Vector containing position and velocity vectors concatanate
%     mu [1]      Gravitational constant of the Earth [km^3/s^2]
%     J2 [1]      J2 effect
%     
% OUTPUT:
%     dy [6x1]     Derivative fo the state vector
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

r      = y(1:3);
r_norm = norm(r);
%Re     = UT_astroConstants(23);
k      = 3/2*J2*Earth_mu*Re^2/r_norm^4;

aj(1) = k*(r(1)./r_norm.*(5*r(3).^2./r_norm^2-1));
aj(2) = k*(r(2)./r_norm.*(5*r(3).^2./r_norm^2-1));
aj(3) = k*(r(3)./r_norm.*(5*r(3).^2./r_norm^2-3));

dy = [y(4:6); -Earth_mu/(norm(r)^3).*y(1:3)+ aj'];

end