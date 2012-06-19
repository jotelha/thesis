% CSA-CSA cross-correlation.
%
% ilya.kuprov@oerc.ox.ac.uk

function csa_test_1()

%% System specification
sys.magnet=14.1;
sys.isotopes={'1H','13C'};
inter.zeeman.eigs={[7 15 -22]
                   [11 18 -29]};
inter.zeeman.euler={[pi/5 pi/3 pi/11]
                    [pi/6 pi/7 pi/15]};
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