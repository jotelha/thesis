% Quadrupole powder pattern.
%
% ilya.kuprov@oerc.ox.ac.uk

function solids_test_4()

% System specification
sys.magnet=14.1;
sys.isotopes={'235U'};
sys.regime='powder';
sys.disable={'zte','pt','krylov'};
inter.coupling.matrix={[400 0 0; 0 400 0; 0 0 -800]};

% Basis set
bas.mode='complete';

% Pulse-acquire parameters
parameters.sweep=20000;
parameters.npoints=128;
parameters.zerofill=512;
parameters.spins='235U';
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Spinach code
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=pulse_acquire(spin_system,parameters);

% Apodization
fid=apodization(fid,'crisp-1d');

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
plot_1d(spin_system,-real(spectrum),parameters);

end