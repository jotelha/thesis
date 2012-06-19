% W-band pulsed ESR spectrum of parafluorotoluene radical.
%
% ilya.kuprov@oerc.ox.ac.uk

function relaxation_parafluorotoluene()

% Read the spin system properties (vacuum DFT calculation)
options.no_xyz=1;
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\parafluorotoluene.log'),{{'E','E'},{'H','1H'},{'F','19F'}},[],options);

% Ignore small HFC anisotropies
sys.tols.inter_cutoff=1e5;

% Set the simulation parameters
sys.magnet=0.33;
bas.mode='ESR-2';
inter.relaxation='redfield';
inter.rlx_keep='secular';
inter.tau_c=50e-12;

% Set the sequence parameters
parameters.offset=0;
parameters.sweep=3e8;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.spins='E';
parameters.axis_units='Gauss';
parameters.derivative=1;

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=pulse_acquire(spin_system,parameters);

% Apodization
fid=apodization(fid,'none-1d');

% Perform Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plot the spectrum
plot_1d(spin_system,real(spectrum),parameters);

end

