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

[n, n1]  = size(A); assert(n1 == n, 'VOU coefficients matrix must be square');
[n1,n2]  = size(V); assert(n1 == n2,'VOU covariance matrix must be square');
                    assert(n1 == n, 'VOU covariance matrix must be same size as coefficients matrix');

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),     'Some target indices out of range');
assert(all(y >=1 & y <= n),     'Some source indices out of range');
assert(isempty(intersect(x,y)), 'Source/target multi-indices overlap');

z = 1:n; z([x y]) = []; % indices of remaining (conditioning) variables
r = [x z];              % indices for the reduced system

err1 = 0;
err2 = 0;

% First we calculate the conditional GC rate F(y -> x | z)

if length(y) == 1 % P scalar, so CARE is a quadratic equation.

	L = chol(V(r,r));
	AOL = A(r,y)'/L;
	VOL = V(y,r)/L;
	a = AOL*AOL';
	b = AOL*VOL'-A(y,y);
	c = VOL*VOL'-V(y,y);
	P = (sqrt(b^2-a*c)-b)/a;
	F1 = P*trace(V(x,x)\(A(x,y)*A(x,y)'));

else

	[P,~,~,rep] = icare(A(y,y)',A(r,y)',V(y,y),V(r,r),V(r,y)');
	err1 = rep.Report;
	if err ~= 0
		F = NaN;
		return
	end
	F1 = trace(V(x,x)\(A(x,y)*P*A(x,y)'));

end

% Now we calculate the unconditional GC rate F(yz -> x)

o = [z y]; % indices for the "other" system

[P,~,~,rep] = icare(A(o,o)',A(x,o)',V(o,o),V(x,x),V(x,o)');
err2 = rep.Report;
if err ~= 0
	F = NaN;
	return
end
F1 = trace(V(x,x)\(A(x,o)*P*A(x,o)'));

% F = F(yz -> x) - F(y -> x | z)

F = f2 - F1;
