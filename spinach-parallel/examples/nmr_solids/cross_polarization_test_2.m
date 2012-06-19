% Cross-polarization experiment in the doubly rotating frame. A single
% nitrogen-15 in a bath of 17 protons taken from a polyalanine chain
% geometry.
%
% WARNING - the applicability of state space restriction to this system has
% not been fully established.
%
% WARNING - absorptive state space boundary is used
%
% ilya.kuprov@oerc.ox.ac.uk
% aanevzor@ncsu.edu

function cross_polarization_test_2()

% System specification
sys.magnet=9.394;
sys.isotopes={'15N','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'};
sys.regime='crystal';

% Case-specific tolerance settings
sys.tols.zte_tol=1e-9;
sys.tols.prox_cutoff=3.0;

% Interactions
inter.zeeman.scalar={0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0};
inter.coordinates={[ 3.371   1.462   0.000] 
                   [-0.524   0.894   0.000] 
                   [-0.543  -0.938   0.000]
                   [ 1.847  -0.534   0.928]
                   [ 3.058  -0.939  -1.274]
                   [ 1.571  -1.903  -1.181]
                   [ 1.610  -0.425  -2.172]
                   [ 3.917   0.530   0.000]
                   [ 5.139   2.584  -0.049]
                   [ 4.080   4.509  -1.336]
                   [ 3.838   2.969  -2.184]
                   [ 2.498   3.695  -1.263]
                   [ 3.953   1.736   2.265]
                   [ 3.600   2.602   4.464]
                   [ 1.325   3.326   3.311]
                   [ 4.464   4.954   4.915]
                   [ 5.601   4.043   3.904]
                   [ 4.542   5.267   3.165]};

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

% Get the dipolar Liouvillian at the input orientation
[L,Q]=h_superop(secularity(spin_system,'nmr'));
L=L+orientation(Q,[0 0 0])+r_superop(spin_system);

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

% Run the evolution for 50 microseconds and watch the nitrogen
duration=0.0005; [timestep,nsteps]=stepsize(L,duration);
nitrogen_x=evolution(spin_system,L,coil_x,rho,timestep,nsteps,'observable');

% Plot the answer
plot(linspace(0,duration,nsteps+1),real(nitrogen_x));
drawnow;

end

