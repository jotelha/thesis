% HMQC spectrum of sucrose (magnetic parameters computed with DFT).
%
% matthew.krzystyniak@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function hmqc_sucrose()

% Read the spin system properties (vacuum DFT calculation)
options.min_j=3.0; options.no_xyz=1; options.dilute={'13C'};
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\sucrose.log'),{{'H','1H'},{'C','13C'}},[31.8 182.1],options);

% Set the simulation parameters
sys.magnet=5.9;
bas.mode='IK-2';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.sweep_f1=3000;
parameters.sweep_f2=1000;
parameters.offset_f1=4250;
parameters.offset_f2=1200;
parameters.npoints_f1=256;
parameters.npoints_f2=256;
parameters.zerofill_f1=512;
parameters.zerofill_f2=512;
parameters.spins_f1='13C';
parameters.spins_f2='1H';
parameters.J=140;
parameters.axis_units='ppm';

% HMQC simulation
fid=hmqc(spin_system,parameters);

% Apodization
fid=apodization(fid,'sinebell-2d');

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill_f2,parameters.zerofill_f1));

% Plotting
contour_plot(spin_system,abs(spectrum),parameters);

end

