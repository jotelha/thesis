% Figure 3 from Chris Timmel's paper (Chemical Physics Letters 298 (1-3), 7-14).
%
% ilya.kuprov@oerc.ox.ac.uk

function mary_test_3()

% Spin system and basis
sys.isotopes={'E','E','1H','1H'};
bas.mode='complete';
bas.projections=0;

% Couplings
inter.zeeman.scalar={2.0023 2.0044 0 0};
inter.coupling.scalar{1,3}=gauss2mhz(35)*1e6;
inter.coupling.scalar{1,4}=gauss2mhz(30)*1e6;
inter.coupling.scalar{4,4}=0;

% Magnetic field and kinetics
log_mT=-5:0.1:3; fields=1e-3*10.^log_mT; 
kinetics=[0.1 1.0 10.0 100.0 1000.0]*1e6;

% Spinach run
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
M=mary(spin_system,fields,kinetics);

% Plotting
plot(log_mT,M)

end