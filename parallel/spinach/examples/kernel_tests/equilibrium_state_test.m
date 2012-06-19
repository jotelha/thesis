% Thermal equilibrium state at 4.2 Kelvin, in both Hilbert and Liouville
% space.
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function equilibrium_state_test()

%% Set up the spin system
sys.magnet=14.1;
sys.isotopes={'E','1H','13C','15N'};
inter.zeeman.scalar={2.0023 1.0 2.0 3.0};
inter.coupling.scalar=cell(length(sys.isotopes));
inter.coupling.scalar{1,2}=1e6;
inter.coupling.scalar{2,3}=1e6;
inter.coupling.scalar{1,3}=1e3;
inter.coupling.scalar{3,4}=1e3;
inter.coupling.scalar{1,4}=1e2;
inter.temperature=4.2;
bas.mode='complete';

%% Run spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

%% Get the equilibrium state in Hilbert space
rho=hs_equilibrium(spin_system);
gens=hs_generators(spin_system);
disp(trace(gens(1).Lz*rho));
disp(trace(gens(2).Lz*rho));
disp(trace(gens(3).Lz*rho));
disp(trace(gens(4).Lz*rho));

%% Get the equilibrium state in Liouville space
rho=equilibrium(spin_system);
disp(state(spin_system,'Lz',1)'*rho);
disp(state(spin_system,'Lz',2)'*rho);
disp(state(spin_system,'Lz',3)'*rho);
disp(state(spin_system,'Lz',4)'*rho);

end