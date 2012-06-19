% Self-consistency test for dipolar relaxation. The relaxation
% superoperator should not depend on the system orientation.
%
% ilya.kuprov@oerc.ox.ac.uk

function dd_test_2()

%% System specification
sys.magnet=14.1;
sys.tols.rlx_integration=1e-5;
sys.isotopes={'1H','1H','1H'};

R=euler2dcm([pi/3 pi/4 pi/5]);

inter.coordinates={[-0.2230    0.7893   -0.5721]*R
                   [-0.6752   -0.3561    0.6460]*R
                   [ 0.8982   -0.4333   -0.0739]*R};
               
inter.relaxation='redfield';
inter.rlx_keep='full';
inter.tau_c=1e-9;
spin_system=create(sys,inter);

%% Basis specification
bas.mode='complete';
spin_system=basis(spin_system,bas);

%% Relaxation superoperator
R=r_superop(spin_system);
disp(R);

end

