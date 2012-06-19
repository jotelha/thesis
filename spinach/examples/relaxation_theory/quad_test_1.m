% Quadrupolar relaxation.
%
% ilya.kuprov@oerc.ox.ac.uk

function quad_test_1()

%% System specification
sys.magnet=14.1;
sys.isotopes={'14N'};
inter.coupling.eigs={[1e4 1e4 -2e4]};
inter.coupling.euler={[0 pi/3 0]};
inter.relaxation='redfield';
inter.rlx_keep='full';
inter.tau_c=2e-9;

spin_system=create(sys,inter);

%% Basis specification
bas.mode='complete';
spin_system=basis(spin_system,bas);

%% Relaxation superoperator
disp(full(r_superop(spin_system)));

end