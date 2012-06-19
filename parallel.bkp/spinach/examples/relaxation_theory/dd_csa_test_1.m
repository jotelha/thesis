% DD-CSA cross-correlation in a two-spin system.
%
% ilya.kuprov@oerc.ox.ac.uk

function dd_csa_test_1()

% Spin system
sys.isotopes={'1H','13C'};

% Basis set
bas.mode='complete';

% Interactions
sys.magnet=14.1;
inter.zeeman.eigs={[7  15 -22]
                   [11 18 -29]};
inter.zeeman.euler={[pi/3 pi/4 pi/5]
                    [pi/6 pi/7 pi/8]};
inter.coordinates={[0 0 0]
                   [0 0 1.02]};
               
% Relaxation theory
inter.relaxation='redfield';
inter.rlx_keep='full';
inter.tau_c=1e-9;

% Spinach code
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
disp(full(r_superop(spin_system)));

end