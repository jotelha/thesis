% Symmetry function test 1.
%
% ilya.kuprov@oerc.ox.ac.uk

function symmetry_test_1()

% Spin system
sys.isotopes ={'E','E','1H', '1H', '1H', '1H'};

% Symmetry
sys.sym_spins={[3 4 5 6]};
sys.sym_group={'D2h'};
sys.sym_a1g_only=0;

% Basis set
bas.mode='complete';
bas.projections=0;

% Interactions
inter.zeeman.scalar={2.002 2.002 0 0 0 0};
inter.coupling.scalar=num2cell(mt2hz([0      0   0.295 0.295 0.295 0.295
                                      0      0     0     0     0     0
                                      0.295  0     0     0     0     0
                                      0.295  0     0     0     0     0
                                      0.295  0     0     0     0     0
                                      0.295  0     0     0     0     0]));

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Hamiltonian superoperator
H=h_superop(secularity(spin_system,'lowfield'));

% Symmetry factorization
S=horzcat(spin_system.bas.irrep.projector);

% Plotting
subplot(1,2,1); spy(abs(H)>1e-6);
subplot(1,2,2); spy(abs(S'*H*S)>1e-6);

end