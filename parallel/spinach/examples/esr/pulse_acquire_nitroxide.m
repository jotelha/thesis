% W-band pulsed powder ESR spectrum of nitroxide radical.
%
% ilya.kuprov@oerc.ox.ac.uk

function pulse_acquire_nitroxide()

% Read the spin system properties (vacuum DFT calculation)
options.no_xyz=1;
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\nitroxide.log'),{{'E','E'},{'N','14N'}},[],options);

% Set the simulation parameters
sys.magnet=3.5;
sys.regime='powder';
sys.disable={'zte','pt','krylov'};
bas.mode='ESR-1';

% Set the sequence parameters
parameters.offset=-2e8;
parameters.sweep=1e9;
parameters.npoints=64;
parameters.zerofill=512;
parameters.spins='E';
parameters.axis_units='Gauss';
parameters.derivative=1;

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=pulse_acquire(spin_system,parameters);

% Apodization
fid=apodization(fid,'crisp-1d');

% Perform Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plot the spectrum
plot_1d(spin_system,real(spectrum),parameters);

end

