% MARY function test 1.
%
% hannah.hogben@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function mary_test_1()

% Set the common parameters
sys.isotopes ={'E','E','1H','1H','1H','1H'};
sys.sym_spins={[3 4]};
sys.sym_group={'S2'};
bas.mode='complete';
bas.projections=0;

% Set the couplings
inter.zeeman.scalar={2.002 2.002 0 0 0 0};
inter.coupling.scalar = num2cell(mt2hz([0       0     0.195  0.195  0    0
                                        0       0     0      0     -1.3  0.2
                                        0.195   0     0      0      0    0  
                                        0.195   0     0      0      0    0  
                                        0      -1.3   0      0      0    0     
                                        0       0.2   0      0      0    0  ]/2));

% Set the fields and kinetics parameters
kinetics=[0.176 0.880 1.76 3.52 8.8 17.6 35.2 52.8]*1e6;
fields=1e-3*(0:0.01:5);

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
M=mary(spin_system,fields,kinetics);

% Plot the answer
plot(fields,M,'r-');

end