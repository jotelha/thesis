% Relaxation superoperator test - DD.
%
% ilya.kuprov@oerc.ox.ac.uk

function dd_test_1()

%% System specification
sys.magnet=14.1;
sys.isotopes={'1H','13C'};

inter.coordinates={[0.0 0.0 0.0]
                   [0.0 0.0 1.02]};

inter.relaxation='redfield';
inter.rlx_keep='full';
inter.tau_c=5e-9;
spin_system=create(sys,inter);

%% Basis specification
bas.mode='complete';
spin_system=basis(spin_system,bas);

%% Relaxation superoperator
disp(full(r_superop(spin_system)));

end