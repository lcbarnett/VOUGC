function stabilizable = is_stabilizable(A,B,tol)

% PBH test for stabilzability
%
% DUALITY: (A,B) stabilizable iff (A',B') detectable

if nargin < 3 || isempty(tol)
	tol = 1e-10;
end

n = size(A,1);

aeigs = eig(A); % all eigenvalues

% Filter for unstable eigenvalues in continuous-time (Real part >= 0)
% A small tolerance accounts for numerical precision issues
ueigs = aeigs(real(aeigs) > -tol);

% Remove duplicate eigenvalues to avoid redundant rank checks
ueigs = unique(ueigs);

% Evaluate the PBH rank condition for each unstable eigenvalue
stabilizable = true;
for i = 1:length(ueigs)
	% Check if the PBH matrix [λ*I-A,B] has full row rank
	if rank([ueigs(i)*eye(n)-A,B]) < n
		stabilizable = false;
		break;
	end
end
