function [t,y] = lorenz(T,y0,p,dt)

% DESCRIPTION:
%
% Integrate the classical Lorenz system (1) numerically
%
% PARAMETERS:
%
% T      - maximum time (integrates from t = 0 to t = T)
% y0     - initial values (3-vector)
% p      - Lorenz system parameters (3-vector, containing sigma, rho and beta parameters - see (1))
% dt     - sampling interval of returned values (may not be same as integration step size!)
%
% RETURN VALUE:
%
% t      - sample time values
% y      - Lorenz variable values
%
% REFERENCES:
%
% (1) https://en.wikipedia.org/wiki/Lorenz_system
%
% (C) Lionel Barnett, 2025

% Some nice defaults:

if nargin < 2 || isempty(y0), y0 = [1; 1; 1];               end
if nargin < 3 || isempty(p),  p  = [10.00; 28.00; 8.0/3.0]; end
if nargin < 4 || isempty(dt), dt = 0.01;                    end

s = p(1); % sigma
r = p(2); % rho
b = p(3); % beta

[t,y] = ode45(@lorenz_fun,0:dt:T,y0);

function dydt = lorenz_fun(tt,yy)
	dydt = [s*(yy(2)-yy(1)); yy(1)*(r-yy(3))-yy(2); yy(1)*yy(2)-b*yy(3)];
end

end
