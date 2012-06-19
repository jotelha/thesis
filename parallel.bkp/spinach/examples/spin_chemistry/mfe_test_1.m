% Recombination kinetics test 1.
%
% hannah.hogben@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function mfe_test_1()

%% System specification

sys.isotopes ={'E','E','1H','1H','1H','1H'};
sys.sym_spins={[3 4]};
sys.sym_group={'S2'};

inter.chem.rp_theory='hore-jones';
inter.chem.rp_spins=[1 2];
inter.chem.rp_rates=[1e6 0];

inter.zeeman.scalar={2.002 2.002 0 0 0 0};
inter.coupling.scalar = num2cell(mt2hz([0       0     0.195  0.195  0    0
                                        0       0     0      0     -1.3  0.2
                                        0.195   0     0      0      0    0  
                                        0.195   0     0      0      0    0  
                                        0      -1.3   0      0      0    0     
                                        0       0.2   0      0      0    0  ]/2));

bas.mode='complete';

spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

%% Simulation

S=mfe(spin_system,1e-3,2e-9,1000);
plot(real(S),'r-');

end