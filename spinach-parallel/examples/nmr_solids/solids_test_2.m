% Two-spin CSA powder pattern.
%
% ilya.kuprov@oerc.ox.ac.uk

function solids_test_2()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
sys.regime='powder';
sys.disable={'pt','zte','krylov'};
inter.zeeman.eigs={[-2 -2 4]-5,[-1 -3 4]+5};
inter.zeeman.euler={[0 0 0],[0 0 0]};

% Basis set
bas.mode='complete';

% Pulse-acquire parameters
parameters.sweep=20000; %Hz
parameters.npoints=128;
parameters.zerofill=512;
parameters.spins='1H';
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
plot_1d(spin_system,real(spectrum),parameters);

end