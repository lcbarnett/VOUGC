function [t,y] = lorenz(T,y0,parms,dt)

% DESCRIPTION:
%
% Integrate the classical Lorenz system (1) numerically
%
% PARAMETERS:
%
% T      - maximum time (integrates from tt = 0 to tt = T
% y0     - initial values (3-vector)
% parms  - Lorenz parameters (3-vector, containing sigma, rho and beta parameters - see (1))
% dt     - sampling interval
%
% RETURN VALUE:
%
% t     - sample time values
% y     - variable values
%
% REFERENCES:
%
% (1) https://en.wikipedia.org/wiki/Lorenz_system
%
% (C) Lionel Barnett, 2025

% parms = [10.00, 28.00, 8.0/3.0]; - reasonable default parms
% y0 =  = [1,1,1];                 - reasonable default initial values

s = parms(1); % sigma
r = parms(2); % rho
b = parms(3); % beta

[t,y] = ode45(@lorenz_fun,0:dt:T,y0);

function dydt = lorenz_fun(tt,yy)
	dydt = [s*(yy(2)-yy(1)); yy(1)*(r-yy(3))-yy(2); yy(1)*yy(2)-b*yy(3)];
end

end
