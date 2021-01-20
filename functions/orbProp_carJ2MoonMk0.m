function dy = orbProp_carJ2Moon(t,y,mu1,mu2,j2,Re,date)

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

r = [y(1) y(2) y(3)];
rmod  = norm(r);

%% perturbations

%J2
aj2 = (3/2)*((j2*mu1*Re^2)/rmod^4)*[((y(1)/rmod)*(5*(y(3)^2/rmod^2)-1))
                                    ((y(2)/rmod)*(5*(y(3)^2/rmod^2)-1))
                                    ((y(3)/rmod)*(5*(y(3)^2/rmod^2)-3))];
%Moon
date = date+(t/(3600*24));
[Moon_r, ~] = ephMoon(date); %Initial Moon position and velocity vectors, in [km] & [km/s], respectively

aMoon = mu2*(((y-Moon_r)/((norm(r-Moon_r))^3))-(Moon_r/((norm(Moon_r))^3)));
%  aMoon = mu2*[(((y(1)-Moon_r(1))/((norm(r-Moon_r))^3))-(Moon_r(1)/((norm(Moon_r))^3)))
%               (((y(2)-Moon_r(2))/((norm(r-Moon_r))^3))-(Moon_r(2)/((norm(Moon_r))^3)))
%               (((y(3)-Moon_r(3))/((norm(r-Moon_r))^3))-(Moon_r(3)/((norm(Moon_r))^3)))];

%% Main function
dy = [y(4)
      y(5)
      y(6)
      (-mu1/rmod^3)*y(1)+aj2(1)+aMoon(1)
      (-mu1/rmod^3)*y(2)+aj2(2)+aMoon(2)
      (-mu1/rmod^3)*y(3)+aj2(3)+aMoon(3)];
end