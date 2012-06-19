% Commutation test 1: all output must be zero.
%
% ilya.kuprov@oerc.ox.ac.uk

function commutation_test_1()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={2.5};
spin_system=create(sys,inter);

% Basis specification
bas.mode='complete';
spin_system=basis(spin_system,bas);

% Generate single-spin commutation superoperators
Lp=operator(spin_system,'L+',1);
Lm=operator(spin_system,'L-',1);
Lz=operator(spin_system,'Lz',1);
Lx=(Lp+Lm)/2;
Ly=(Lp-Lm)/2i;

% Check the commutation relations
disp(norm(Lz*Lp-Lp*Lz-Lp,'fro'));
disp(norm(Lz*Lm-Lm*Lz+Lm,'fro'));
disp(norm(Lx*Ly-Ly*Lx-1i*Lz,'fro'));

end