% A function that generates the equilibrium density matrix at a given
% temperature.
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function rho=hs_equilibrium(spin_system,contributions)

%% Preamble

% Click forward for output
spin_system=click(spin_system,'forward');

% Inform the user
report(spin_system,'hs_equilibrium: computing the thermal equilibrium density matrix...');

% Set the defaults
if ~exist('contributions','var'), contributions='isotropic'; end

% Generate the Hamiltonian
switch contributions
    
    case 'isotropic'
        report(spin_system,'hs_equilibrium: using the isotropic part of the Hamiltonian...');
        H=hs_hamiltonian(secularity(spin_system,'keep_all'));
    
    case 'full'
        report(spin_system,'hs_equilibrium: using the full Hamiltonian at the input orientation...');
        [H_iso,Q]=hs_hamiltonian(secularity(spin_system,'keep_all'));
        H=H_iso+orientation(Q,[0 0 0]);
    
    otherwise
        error('hs_equilibrium: unknown contribution specification.');

end

% Decide the approximation
if spin_system.rlx.temperature==0
    
    % If zero is supplied, use the high temperature approximation
    report(spin_system,'hs_equilibrium: WARNING - high temperature approximation, trace ignored.')
    rho=-H./max(diag(H));
    
else
    
    % Temperature factor
    beta=spin_system.tols.hbar/(spin_system.tols.kbol*spin_system.rlx.temperature);
    
    % Compute density matrix as a propagator in imaginary time
    rho=propagator(spin_system,H,1i*beta);
    
    % Normalize with respect to the trace
    rho=rho./trace(rho);
    
end

end

% Blessed is the man that walketh not in the counsel of the wicked.
%
% Psalm 1:1-3

