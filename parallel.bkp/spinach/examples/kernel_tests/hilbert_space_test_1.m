% Hilbert space commutation test 1: all output must be zero.
%
% ilya.kuprov@oerc.ox.ac.uk

function hilbert_space_test_1()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','235U'};
inter.zeeman.scalar={2.5 1.0};
inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{2,1}=10;
spin_system=create(sys,inter);

% Single-spin operator commutation test
Up=hs_operator(spin_system,'L+',2);
Um=hs_operator(spin_system,'L-',2);
Uz=hs_operator(spin_system,'Lz',2);
Ux=(Up+Um)/2;
Uy=(Up-Um)/2i;

disp(norm(full(Uz*Up-Up*Uz-Up)));
disp(norm(full(Uz*Um-Um*Uz+Um)));
disp(norm(full(Ux*Uy-Uy*Ux-1i*Uz)));

% Two-spin operator commutation test
HpUp=hs_operator(spin_system,{'L+','L+'},{1,2});
HmUm=hs_operator(spin_system,{'L-','L-'},{1,2});
Hz=hs_operator(spin_system,'Lz',1);

disp(norm(full(Uz*HpUp-HpUp*Uz-HpUp)));
disp(norm(full(Hz*HpUp-HpUp*Hz-HpUp)));
disp(norm(full(Uz*HmUm-HmUm*Uz+HmUm)));
disp(norm(full(Hz*HmUm-HmUm*Hz+HmUm)));

end