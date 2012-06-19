% Cross-polarization experiment in the doubly rotating frame. A single
% nitrogen-15 in a bath of 12 protons scattered on a 2 Angstrom radius
% sphere around it.
%
% WARNING - the applicability of state space restriction to this system has
% not been established.
%
% WARNING - absorptive state space boundary is used
%
% ilya.kuprov@oerc.ox.ac.uk

function cross_polarization_test_1()

% System specification
sys.magnet=9.394;
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','15N'};
sys.regime='crystal';

% Case-specific tolerance settings
sys.tols.zte_tol=1e-10;
sys.tols.prox_cutoff=3.0;

% Interactions
inter.zeeman.scalar={0 0 0 0 0 0 0 0 0 0 0 0 0};
inter.coordinates={[-2.51887974   -0.99807637   -0.87365101]
                   [-1.09220609   -0.49369054   -0.00000060]
                   [-2.51887819   -0.99807636    0.87365165]
                   [-2.54376425    2.65945619   -1.20122894]
                   [-1.11551509    1.65289357   -1.19927242]
                   [-2.54058796    1.14568573   -2.07390040]
                   [-4.50219437    0.45366792    0.87353250]
                   [-4.50217572    1.96708252    0.00023756]
                   [-4.50219487    0.45407982   -0.87377011]
                   [-2.54233875    1.14692134    2.07390117]
                   [-1.11551334    1.65103779    1.20034313]
                   [-2.54201521    2.66007643    1.20015742]
                   [-2.67552180    0.95825426    0.00000000]};

% Basis - currently set to five-spin orders
bas.mode='IK-1';
bas.level=5;
bas.space_level=5;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up an absorptive outer boundary
spin_system.rlx.boundary='absorptive';
spin_system.rlx.boundary_damping_rate=2e4;

% Get the dipolar Liouvillian at a random orientation
[L,Q]=h_superop(secularity(spin_system,'nmr'));
L=L+orientation(Q,[0.1 0.2 0.3])+r_superop(spin_system);

% Get the irradiation operators
H_plus=operator(spin_system,'L+','1H');
Hx=(H_plus+H_plus')/2; Hy=(H_plus-H_plus')/2i;
N_plus=operator(spin_system,'L+','15N');
Nx=(N_plus+N_plus')/2; Ny=(N_plus-N_plus')/2i;

% Add the irradiation operators with matched power levels
power_1h=2*pi*5e4;
power_15n=2*pi*5e4;
L=L+power_1h*Hy+power_15n*Nx;

% Set the initial condition
rho=equilibrium(spin_system);

% Set the detection states
coil_x=(state(spin_system,'L+','15N')+...
        state(spin_system,'L-','15N'))/2;

% Appy the excitation pulses
rho=step(spin_system,Hx,rho,pi/2);
rho=step(spin_system,Ny,rho,pi/2);

% Run the evolution for 0.01 seconds and watch the nitrogen
duration=0.0005; [timestep,nsteps]=stepsize(L,duration);
nitrogen_x=evolution(spin_system,L,coil_x,rho,timestep,nsteps,'observable');

% Plot the answer
plot(linspace(0,duration,nsteps+1),real(nitrogen_x));

end