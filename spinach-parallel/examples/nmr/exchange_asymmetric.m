% Two-spin asymmetric chemical exchange pattern.
%
% ilya.kuprov@oerc.ox.ac.uk

function exchange_asymmetric()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,3};
inter.chem.exchange=[0 2000; 500 0];
spin_system=create(sys,inter);

% Basis specification
bas.mode='complete';
spin_system=basis(spin_system,bas);

% Pulse-acquire parameters
parameters.offset=900;
parameters.sweep=5000;
parameters.npoints=512;
parameters.zerofill=1024;
parameters.spins='1H';
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Running the sequence
fid=pulse_acquire(spin_system,parameters);

% Apodization
fid=apodization(fid,'crisp-1d');

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
plot_1d(spin_system,real(spectrum),parameters);

end