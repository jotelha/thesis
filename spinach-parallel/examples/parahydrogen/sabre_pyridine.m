% SABRE experiment simulation for Eibe Duecker and Lars Kuhn. Set to
% reproduce Figure 3b from 
%
%        http://dx.doi.org/10.1021/ja903601p
%
% ilya.kuprov@oerc.ox.ac.uk

function sabre_pyridine()

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H'};

% Chemical shifts
inter.zeeman.scalar={8.54 7.44 7.86 7.44 8.54 -23.5 -23.5};

% Couplings inside pyridine
inter.coupling.scalar=cell(7);
inter.coupling.scalar{1,2}= 4.88; 
inter.coupling.scalar{4,5}= 4.88;
inter.coupling.scalar{1,4}= 1.00; 
inter.coupling.scalar{2,5}= 1.00;
inter.coupling.scalar{1,3}= 1.84; 
inter.coupling.scalar{3,5}= 1.84;
inter.coupling.scalar{1,5}=-0.13;
inter.coupling.scalar{2,3}= 7.67; 
inter.coupling.scalar{3,4}= 7.67;
inter.coupling.scalar{2,4}= 1.37;

% Couplings of the hydride group - please check
inter.coupling.scalar{1,6}=1.12;
inter.coupling.scalar{1,7}=1.02;
inter.coupling.scalar{6,7}=7.00;

% Basis set
bas.mode='complete';

% Set the magnet to 25 mT
sys.magnet=25e-3;

% Do the housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set the secularity assumptions
spin_system=secularity(spin_system,'nmr');

% Get the Hamiltonian commutation superoperator
H=h_superop(spin_system);

% Start in a singlet state
rho=singlet(spin_system,6,7);

% Determine the optimal time step
[timestep,nsteps]=stepsize(H,2.5);

% Evolve the system for 2.5 seconds
rho=evolution(spin_system,H,[],rho,timestep,nsteps,'final');

% Disconnect the parahydrogen
[H,rho]=decouple(spin_system,H,rho,[6 7]);

% Determine the optimal time step
[timestep,nsteps]=stepsize(H,2.5);

% Evolve the system for a further 2.5 seconds
rho=evolution(spin_system,H,[],rho,timestep,nsteps,'final');

% Get the Zeeman Hamiltonian (by killing all couplings)
spin_system.inter.coupling.matrix=cell(spin_system.comp.nspins);
H_zeeman=decouple(spin_system,h_superop(spin_system),[],[6 7]);

% Lift the field exponentially in 1024 steps 
% from 25 mT to 7.05 T over 5 seconds
 for n=2:1024
     H_zeeman=1.005524881*H_zeeman; H_total=H+H_zeeman;
     rho=step(spin_system,H_total,rho,5/1024);
 end

% Determine the optimal time step
[timestep,nsteps]=stepsize(H_total,1.0);

% Evolve the system for a further second in high field
rho=evolution(spin_system,H_total,[],rho,timestep,nsteps,'final');

% Pulse-acquire parameters
parameters.sweep=2000;
parameters.npoints=8192;
parameters.zerofill=65536;
parameters.spins='1H';
parameters.axis_units='Hz';
parameters.invert_axis=1;

% Pulse-acquire sequence
fid=pulse_acquire(spin_system,parameters,H_total,rho);

% Apodization
fid=apodization(fid,'exponential-1d',20);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
plot_1d(spin_system,real(spectrum),parameters);

