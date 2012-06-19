% A 1H-15N cross-relaxation experiment with lots of protons.
%
% ilya.kuprov@oerc.ox.ac.uk

function cross_relaxation()

% Read the spin system parameters (vacuum DFT calculation)
options.min_j=1.0;
[sys,inter]=g03_to_spinach(g03_parse('../molecules/tma.log'),{{'H','1H'},{'N','15N'}},[31.8 0],options);

% Magnetic field
sys.magnet=14.1;

% Basis set
bas.mode='IK-1';
bas.level=3;
bas.space_level=3;

% Relaxation theory
inter.relaxation='redfield';
inter.equilibrium='thermal';
inter.temperature=298;
inter.tau_c=100e-12;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Initial state - thermal equilibrium
rho=equilibrium(spin_system);

% Detection state - Lz on the nitrogen
coil=state(spin_system,'Lz','15N');

% Static Liouvillian
L=h_superop(secularity(spin_system,'nmr'));

% Proton pulse operator
Lp=operator(spin_system,'L+','1H');
Lx=(Lp+Lp')/2;

% Relaxation superoperator
R=r_superop(spin_system);

% An inversion pulse on protons
rho=step(spin_system,Lx,rho,pi);

% Evolution
answer=evolution(spin_system,L+1i*R,coil,rho,1e-3,10000,'observable');

% Plotting
plot(real(answer)');

end