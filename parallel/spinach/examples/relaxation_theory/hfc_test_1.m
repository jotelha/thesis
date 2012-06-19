% Relaxation superoperator test - hyperfine coupling.
%
% ilya.kuprov@oerc.ox.ac.uk

function hfc_test_1()

%% System specification
sys.magnet=14.1;
sys.isotopes={'1H','E'};

inter.coordinates={[0.0 0.0 0.0 ]
                   [0.0 0.0 1.5]};

inter.relaxation='redfield';
inter.rlx_keep='full';
inter.tau_c=10e-12;
spin_system=create(sys,inter);

%% Basis specification
bas.mode='complete';
spin_system=basis(spin_system,bas);

%% Relaxation superoperator
disp(full(r_superop(spin_system)));

end