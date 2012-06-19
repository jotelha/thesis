% Normaization test 1: all output must be zero.
%
% ilya.kuprov@oerc.ox.ac.uk

function normalization_test_1()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','235U'};
inter.zeeman.scalar={2.5 1.0};
inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{2,1}=10;
spin_system=create(sys,inter);

% Basis specification
bas.mode='complete';
spin_system=basis(spin_system,bas);

% Single-spin state norm test
Up=state(spin_system,'L+',2);
Um=state(spin_system,'L-',2);
Uz=state(spin_system,'Lz',2);
Ux=(Up+Um)/2; Uy=(Up-Um)/2i;

Hp=state(spin_system,'L+',1);
Hm=state(spin_system,'L-',1);
Hz=state(spin_system,'Lz',1);
Hx=(Hp+Hm)/2; Hy=(Hp-Hm)/2i;

disp(norm(Ux)-norm(Uy));
disp(norm(Uy)-norm(Uz));
disp(norm(Uz)-norm(Ux));
disp(norm(Hx)-norm(Hy));
disp(norm(Hy)-norm(Hz));
disp(norm(Hz)-norm(Hx));

end