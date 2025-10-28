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

% Skip initialisation if already initialised (unless vougc_forceinit is set)

global vougc_initialised
if vougc_initialised
	if exist('vougc_forceinit') && vougc_forceinit
		fprintf('[VOUGC startup] Already initialised - forcing re-initialisation\n');
		vougc_initialised = false;
		clear vougc_forceinit
	else
		fprintf('[VOUGC startup] Already initialised - skipping\n');
		return
	end
else
	fprintf('[VOUGC startup] Initialising\n');
end

global vougc_root
vougc_root = fileparts(mfilename('fullpath')); % directory containing this file
addpath(vougc_root);
addpath(fullfile(vougc_root,'demo'));
fprintf('[VOUGC startup] Added appropriate paths\n');

% Check if we have the icare function from the Control System Toolbox.

if exist('icare') ~= 2
	fprintf(2,'[VOUGC startup]\n');
    fprintf(2,'[VOUGC startup] WARNING: The ''icare'' function from the Matlab Control System Toolbox does\n');
    fprintf(2,'[VOUGC startup]          not appear to be present. Some functionality will be unavailable.\n');
	fprintf(2,'[VOUGC startup]\n');
end

vougc_initialised = true;
fprintf('[VOUGC startup] Initialisation complete (you may re-run ''startup'' at any time)\n');
