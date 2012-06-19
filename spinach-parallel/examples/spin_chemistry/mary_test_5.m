% Figure 5 from Chris Timmel's paper (Molecular Physics 95 (1), 71-89).
%
% ilya.kuprov@oerc.ox.ac.uk

function mary_test_5()

% Spin system and basis
sys.isotopes={'E','E','1H','1H'};
bas.mode='complete';
bas.projections=0;

% Couplings
inter.zeeman.scalar={2.0023 2.0023 0 0};
inter.coupling.scalar{1,3}=gauss2mhz(20)*1e6;
inter.coupling.scalar{2,4}=gauss2mhz(20)*1e6;
inter.coupling.scalar{4,4}=0; 

% Magnetic field and kinetics
fields=linspace(0,0.5*20/1e4,200);
kinetics=2*pi*[0.0001 0.001 0.01 0.02 0.05 0.1 0.2]*gauss2mhz(20)*1e6;

% Spinach run
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
M=mary(spin_system,fields,kinetics);

% Plotting
plot(linspace(0,0.5,200),M)

end