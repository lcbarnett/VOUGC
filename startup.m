%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VOUGC Toolbox start-up script
%
% Initialise VOUGC Toolbox. This script is run automatically if Matlab is started
% in the toolbox root (installation) directory.
%
% (C) Lionel Barnett, 2025. See the README file in the installation directory for
% an overview, and the LICENSE file for licensing terms.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('[VOUGC startup] Initialising VOUGC toolbox.\n');

global vougc_root
vougc_root = fileparts(mfilename('fullpath')); % directory containing this file
addpath(vougc_root);
addpath(fullfile(vougc_root,'utils'));

% Check if we have the icare function from the Control System Toolbox.

if exist('icare') ~= 2
	fprintf(2,'[VOUGC startup]\n');
    fprintf(2,'[VOUGC startup] WARNING: The ''icare'' function from the Matlab Control System Toolbox does\n');
    fprintf(2,'[VOUGC startup]          not appear to be present. Some functionality will be unavailable.\n');
	fprintf(2,'[VOUGC startup]\n');
end

fprintf('[VOUGC startup] Initialisation complete (you may re-run ''startup'' at any time).\n');
