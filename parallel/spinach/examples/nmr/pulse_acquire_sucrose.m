% 1H NMR spectrum of sucrose (magnetic parameters read in from a DFT calculation),
% including Redfield relaxation superoperator.
%
% ilya.kuprov@oerc.ox.ac.uk

function pulse_acquire_sucrose()

% Read the spin system properties (vacuum DFT calculation)
options.min_j=1.0;
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\sucrose.log'),{{'H','1H'}},31.8,options);

% Set the simulation parameters
sys.magnet=14.1;
bas.mode='IK-2';
inter.relaxation='redfield';
inter.rlx_keep='kite';
inter.tau_c=200e-12;

% Sequence parameters
parameters.offset=2800;
parameters.sweep=6500;
parameters.npoints=8192;
parameters.zerofill=65536;
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

