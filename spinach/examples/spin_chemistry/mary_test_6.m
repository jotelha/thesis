% MARY function test 6. Phosphorus radical, including relaxation on 31P.
%
% hannah.hogben@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function mary_test_6()

% Set the common parameters
sys.isotopes ={'E','E','31P','1H','1H','1H','1H','1H','1H'};
sys.sym_spins={[4 5],[6 7],[8 9]};
sys.sym_group={'S2','S2','S2'};
bas.mode='complete';
bas.projections=0;

% Set the couplings
inter.zeeman.scalar={2.002 2.002 0 0 0 0 0 0 0};
inter.coupling.scalar = num2cell(mt2hz([0       0     69.5   0      0      0     0     0      0
                                        0       0     0      1.526  1.526  0.490 0.490 0.175  0.175
                                        69.5    0     0      0      0      0     0     0      0
                                        0       1.526 0      0      0      0     0     0      0 
                                        0       1.526 0      0      0      0     0     0      0 
                                        0       0.490 0      0      0      0     0     0      0 
                                        0       0.490 0      0      0      0     0     0      0 
                                        0       0.175 0      0      0      0     0     0      0 
                                        0       0.175 0      0      0      0     0     0      0 ]/2));
inter.coupling.matrix{1,3}=mt2hz([-11.0 0 0; 0 -11.0 0; 0 0 22.0]);
inter.coupling.matrix{9,9}=[];

% Relaxation theory
inter.tau_c=2.4e-12;     
inter.relaxation='redfield';
inter.rlx_keep='full';

% Set the fields and kinetics parameters
kinetics=1e6;
fields=1e-3*(0:0.1:5);

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
M=mary(spin_system,fields,kinetics);

% Plot the answer
plot(fields*1000,M,'g-');

end