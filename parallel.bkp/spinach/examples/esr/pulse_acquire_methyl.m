% X-band pulsed ESR spectrum of methyl radical.
%
% ilya.kuprov@oerc.ox.ac.uk

function pulse_acquire_methyl()

% Read the spin system properties (vacuum DFT calculation)
options.no_xyz=1;
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\methyl.log'),{{'E','E'},{'H','1H'}},[],options);

% Set the simulation parameters
sys.magnet=0.33;
bas.mode='ESR-1';
inter.relaxation='damp';
inter.damp_rate=5e7;

% Set the sequence parameters
parameters.offset=0;
parameters.sweep=5e8;
parameters.npoints=128;
parameters.zerofill=512;
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

