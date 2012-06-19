% Commutation test 3: all output must be zero.
%
% ilya.kuprov@oerc.ox.ac.uk

function commutation_test_3()

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

% Single-spin operator commutation test
Up=operator(spin_system,'L+',2);
Um=operator(spin_system,'L-',2);
Uz=operator(spin_system,'Lz',2);
Ux=(Up+Um)/2;
Uy=(Up-Um)/2i;

disp(norm(full(Uz*Up-Up*Uz-Up)));
disp(norm(full(Uz*Um-Um*Uz+Um)));
disp(norm(full(Ux*Uy-Uy*Ux-1i*Uz)));

% Two-spin operator commutation test
HpUp=operator(spin_system,{'L+','L+'},{1,2});
HmUm=operator(spin_system,{'L-','L-'},{1,2});
Hz=operator(spin_system,'Lz',1);

disp(norm(full(Uz*HpUp-HpUp*Uz-HpUp)));
disp(norm(full(Hz*HpUp-HpUp*Hz-HpUp)));
disp(norm(full(Uz*HmUm-HmUm*Uz+HmUm)));
disp(norm(full(Hz*HmUm-HmUm*Hz+HmUm)));

% Product superoperator commutation test
HpUm_left=-2*p_superop(spin_system,[3 1],'left');
HpUm_right=-2*p_superop(spin_system,[3 1],'right');

disp(norm(full(HpUm_left*HpUm_right-HpUm_right*HpUm_left)));

end