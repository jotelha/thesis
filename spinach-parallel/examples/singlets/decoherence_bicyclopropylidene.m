% Long-lived spin states in the bicyclopropylidene molecule
% (8 protons, 65536-dimensional Liouville space). The relaxation
% superoperator accounts for every dipolar coupling and every
% CSA tensor in the system.
%
% The calculation takes about an hour and requires 24GB of RAM.
%
% hannah.hogben@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function decoherence_bicyclopropylidene()

% Read the spin system (coordinates, chemical shifts,
% J-couplings and CSAs) from a vacuum DFT calculation
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\bicyclopropylidene.log'),{{'H','1H'}},31.8);

% Set magnet field to 1.0 Tesla
sys.magnet=1.0;

% Tighten up the tolerances
sys.tols.prox_cutoff=Inf;
sys.tols.rlx_integration=1e-5;

% Set relaxation theory parameters
inter.relaxation='redfield';
inter.rlx_keep='full';
inter.tau_c=100e-12;
spin_system=create(sys,inter);

% Use complete basis set
bas.mode='complete';
spin_system=basis(spin_system,bas);

% Build the relaxation superoperator
R=r_superop(spin_system);

% List twenty smallest magnitude eigenvalues
disp(eigs(R-speye(size(R)),20,'SM')+1);

end