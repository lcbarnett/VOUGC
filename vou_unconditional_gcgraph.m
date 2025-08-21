function [F,err1,err2] = vou_unconditional_gcgraph(A,V)

% DESCRIPTION:
%
% Calculate time-domain unconditional Granger-causal graph for a vector
% Ornstein-Uhlenbeck (VOU) process X(t):
%
%	dX(t) = AX(t)dt + dW(t)
%
% with autoregression matrix A, and W(t) a Wiener process with dW(t) ~ N(0,Vdt).
% NOTE: A need not be stable (i.e., may have eigenvalues with nonnegative real
% part).
%
% PARAMETERS:
%
% A     - VOU coefficients matrix
% V     - VOU Wiener process covariance matrix (or empty for identity matrix)
%
% RETURN VALUE:
%
% F     - pairwise-unconditional Granger causality rates (unconditional Granger-causal graph)
%
% REFERENCES:
%
% (1) L. Barnett and A. K. Seth (2015): Granger causality for state-space models, Phys. Rev. E 91(4) Rapid Communication.
% (2) L. Barnett and A. K. Seth (2016): Detectability of Granger causality for subsampled continuous-time neurophysiological processes, J. Neurosci. Methods 275.
% (3) L. Barnett (2017): Granger causality rate for a vector Ornstein-Uhlenbeck process (working notes).
%
% (C) Lionel Barnett, 2024

if nargin < 2
	V = [];
end

n = size(A,1);

err1 = 0;
err2 = 0;

F = nan(n);
for i = 1:n
	oi = 1:n; oi(i) = []; % omit i
	[F1,err1] = vou_conditional_gc(A,V,i,oi); % F([i] -> i) - actually unconditional
	if err1 ~= 0
		F = NaN;
		return
	end
	for j = 1:n
		if j == i, continue; end
		oij = 1:n; oij([i j]) = []; % omit i,j
		[F2,err2] = vou_conditional_gc(A,V,i,oij); % F([ij] -> i | j)
		if err2 ~= 0
			F = NaN;
			return
		end
		F(i,j) = F1 - F2; % F(i,j) = F([i] -> i) - F([ij] -> i | j) - Geweke identity
	end
end
