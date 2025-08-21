% Demonstration script for Granger causality map of a (noiseless) Lorenz system
%
% Demo parameters (override on command line)

defvar('Te',     10                      ); % equilibriation time
defvar('T',      200                     ); % sample time
defvar('dt',     0.01                    ); % sampling time increment
defvar('y0',     [1; 1; 1]               ); % initial values
defvar('p',      [10.00; 28.00; 8.0/3.0] ); % Lorenz system parameters

% Simulate the Lorenz system

[t,y] = lorenz(Te+T,y0,p,dt);

n = round(T/dt);

t = t(end-n+1:end);
y = y(end-n+1:end,:);

% Calculate Jacobian matrices along simulation trajectory (attractor)

J = zeros(3,3,n);
for k = 1:n
	J(:,:,k) = lorenz_jacobian(y(k,:),p);
end

% Calculate maximum of real part of eigenvalues of the Jacobians (< 0 at locally-stable points)

lam = zeros(n,1);
for k = 1:n
	lam(k) = max(real(eig(J(:,:,k))));
end

% Calculate pairwise-conditional local Granger causality rates (assume identity noise covariance matrix)

R = zeros(3,3,n);
for k = 1:n
	R(:,:,k) = vou_conditional_gcgraph(J(:,:,k));
end

% Calculate average of local GC rates across the attractor (global causality rate)

RAVG = mean(R,3);
fprintf('\nAverage GC rate =\n');
disp(RAVG);

% Plot eigenvalues on attractor

z1 = [y(:,1) nan(n,1)]; % for Matlab 'patch' function (don't ask)
z2 = [y(:,2) nan(n,1)];
z3 = [y(:,3) nan(n,1)];

figure(1); clf;
patch(z1,z2,z3,[lam nan(n,1)],'EdgeColor','interp','FaceColor','none');
colorbar;
title('Stability (maximum of real part of Jacobian eigenvalues)');

% Plot GC maps on the attractor

figure(2); clf;
maxr = ceil(nanmax(R(:)));
sgtitle(sprintf('Granger causality maps for Lorenz system\n'));
for i = 1:3
	for j = 1:3
		if j == i, continue; end
		subplot(3,3,3*(i-1)+j)
		patch(z1,z2,z3,[squeeze(R(i,j,:)) nan(n,1)],'EdgeColor','interp','FaceColor','none');
		clim([0,maxr]);
		colorbar;
		xlabel(sprintf('R(%d -> %d)',j,i));
	end
end
