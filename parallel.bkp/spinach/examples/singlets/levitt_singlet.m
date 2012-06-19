% A demonstration that the two-spin singet state is immune 
% to dipolar relaxation.
%
% ilya.kuprov@oerc.ox.ac.uk

function levitt_singlet()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};

inter.coordinates={[0.0 0.0 0.0]
                   [0.5 0.6 0.7]};

inter.relaxation='redfield';
inter.rlx_keep='full';
inter.tau_c=5e-9;
spin_system=create(sys,inter);

% Basis specification
bas.mode='complete';
spin_system=basis(spin_system,bas);

% Relaxation superoperator
R=r_superop(spin_system);

% Singlet relaxation rate
S=singlet(spin_system,1,2);
disp(['Singlet state relaxation rate: ' num2str(S'*R*S)]);

end