% DD-CSA cross-correlation -- a reproduction of Fig 5a from
% Grace and Kumar, JMR, 2005.
%
% ilya.kuprov@oerc.ox.ac.uk

function dd_csa_test_2()

% Read the spin system parameters (vacuum DFT calculation)
[sys,inter]=g03_to_spinach(g03_parse('../molecules/fdnb.log'),{{'H','1H'},{'F','19F'}},[32.0 270.0]);

% Set up the calculation
sys.magnet=9.4;                  % Magnet induction
sys.tols.prox_cutoff=5;          % Increase proximity cutoff to 5 Angstrom
inter.relaxation='redfield';     % Redfield relaxation theory
inter.rlx_keep='secular';        % Rotating-frame version of Redfield theory
inter.equilibrium='thermal';     % Relax to thermal equilibrium
inter.tau_c=9.6e-12;             % Correlation time
bas.mode='complete';             % Complete basis set

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set simulation parameters
parameters.offset=-530;          % Axis offset
parameters.sweep_width=70;       % Sweep width
parameters.f1.npoints=1024;      % Number of simulation points
parameters.f1.zerofill=2048;     % Zerofill fid to twice the length 
parameters.nuclei='19F';         % Run the experiment on 19F

% Get the Hamiltonian superoperator, add Redfield superoperator, subtract
% the offset Zeeman terms
L=   h_superop(secularity(spin_system,'nmr'))+...
  1i*r_superop(spin_system)-...
     offset(spin_system,parameters.nuclei,parameters.offset);

% Get the L+ commutation superoperator on fluorine
Lp=operator(spin_system,'L+',parameters.nuclei);

% Set detection state to L+
coil=state(spin_system,'L+',parameters.nuclei);

% Calculate the time step of the simulation
timestep=1/parameters.sweep_width;

% Do not overwrite the plots
hold on;

% Loop over mixing times
for t_mix=[0.1 1.4 1.6 1.8 2.0 2.2 2.4 10]
    
    % Set the state vector to thermal equilibrium
    rho=equilibrium(spin_system);
    
    % Apply the inversion pulse
    rho=step(spin_system,(Lp+Lp')/2,rho,pi);
    
    % Run the mixing time
    rho=evolution(spin_system,L,coil,rho,t_mix,1,'final');
    
    % Apply the detection pulse
    rho=step(spin_system,(Lp-Lp')/2i,rho,pi/2);
    
    % Run the detection period
    fid=evolution(spin_system,L,coil,rho,timestep,parameters.f1.npoints-1,'observable');
    
    % Perform the Fourier transform
    spectrum=fftshift(fft(apodization(fid,'exponential-1d',40),parameters.f1.zerofill));
    
    % Plot the spectrum
    plot(real(spectrum)); drawnow;
    
end

% Invert the X axis
set(gca,'XDir','reverse'); axis tight;

