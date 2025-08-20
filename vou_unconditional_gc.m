function [F,err1,err2] = vou_unconditional_gc(A,V,x,y)

% DESCRIPTION:
%
% Calculate time-domain unconditional Granger causality rate for a vector
% Ornstein-Uhlenbeck (VOU) process X(t):
%
%	dX(t) = AX(t)dt + dW(t)
%
% with autoregression matrix A, and W(t) a Wiener process with dW(t) ~ N(0,Vdt).
% NOTE: A need not be stable (i.e., may have eigenvalues with nonnegative real
% part).
%
% Source/target variables may be multivariate. Uses a state-space method which
% involves solving an associated continuous-time algebraic Riccati equation (CARE).
%
% PARAMETERS:
%
% A     - VOU coefficients matrix
% V     - VOU Wiener process covariance matrix
% x     - multi-index of target variable
% y     - multi-index of source variable
%
% RETURN VALUES:
%
% F     - Granger causality rate from y to x, conditional on other variables
% err   - CARE error report number (zero if no error)
%
% Possible CARE errors are (see 'icare' in the Matlab Control System Toolbox):
%
%	0 - No errors (unique solution is accurate)
%	1 - Solution accuracy is poor
%	2 - Solution not finite
%	3 - No solution found (Hamiltonian spectrum has imaginery eigenvalues)
%	4 - No solution found ("pencil" is singular)
%
% REFERENCES:
%
% (1) L. Barnett and A. K. Seth (2015): Granger causality for state-space models, Phys. Rev. E 91(4) Rapid Communication.
% (2) L. Barnett and A. K. Seth (2016): Detectability of Granger causality for subsampled continuous-time neurophysiological processes, J. Neurosci. Methods 275.
% (3) L. Barnett (2017): Granger causality rate for a vector Ornstein-Uhlenbeck process (working notes).
%
% (C) Lionel Barnett, 2024

n = size(A,1);
x = x(:)'; % vectorise
y = y(:)'; % vectorise
z = 1:n; z([x y]) = []; % indices of remaining (conditioning) variables

err1 = 0;
err2 = 0;

[F1,err1] = vou_conditional_gc(A,V,x,z); % F(z -> x | y)
if err1 ~= 0
	F = NaN;
	return
end

[F2,err2] = vou_conditional_gc(A,V,x,[y z]); % F(yz -> x) - actually unconditional
if err2 ~= 0
	F = NaN;
	return
end

F = F2 - F1; % F = F(yz -> x) - F(z -> x | y) - Geweke identity
