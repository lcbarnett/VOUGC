function J = lorenz_jacobian(y,p)

% DESCRIPTION:
%
% Integrate the classical Lorenz system (1) numerically
%
% PARAMETERS:
%
% y  - Lorenz variable values - point in phase space (3-vector)
% p  - Lorenz syatem parameters (3-vector, containing sigma, rho and beta parameters - see (1))
%
% RETURN VALUE:
%
% J  - the Jacobian matrix at y (3 x 3 matrix)
%
% REFERENCES:
%
% (1) https://en.wikipedia.org/wiki/Lorenz_system
%
% (C) Lionel Barnett, 2025

J = [ ...
	-p(1)       p(1)   0    ; ...
	 p(2)-y(3) -1     -y(1) ; ...
	 y(2)       y(1)  -p(3)   ...
];
