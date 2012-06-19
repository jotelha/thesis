% X-band pulsed ESR spectrum of perfluoropyrene cation radical, computed using
% brute force operator algebra in the full 4,194,304 dimensional Liouville space.
% 
% This calculation requires at least 12GB of RAM and tests the implementation of
% trajectory-level state space restriction tools supplied with Spinach.
%
% ilya.kuprov@oerc.ox.ac.uk

function pulse_acquire_perfluoropyrene()

% Read the spin system properties (vacuum DFT calculation)
options.no_xyz=1;
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\pyrene_cation.log'),{{'E','E'},{'H','1H'}},[],options);

% Set the simulation parameters
sys.magnet=0.33;
bas.mode='complete';

% Set the sequence parameters
parameters.offset=0;
parameters.sweep=3e8;
parameters.npoints=2048;
parameters.zerofill=4096;
parameters.spins='E';
parameters.axis_units='Gauss';
parameters.derivative=1;

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=pulse_acquire(spin_system,parameters);

% Apodization
fid=apodization(fid,'exponential-1d',10);

% Perform Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plot the spectrum
plot_1d(spin_system,real(spectrum),parameters);

end

