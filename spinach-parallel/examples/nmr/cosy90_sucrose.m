% COSY spectrum of sucrose (magnetic parameters computed with DFT).
%
% ilya.kuprov@oerc.ox.ac.uk

function cosy90_sucrose()

% Read the spin system properties (vacuum DFT calculation)
options.min_j=2.0; options.no_xyz=1;
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\sucrose.log'),{{'H','1H'}},31.8,options);

% Set the simulation parameters
sys.magnet=5.9;
bas.mode='IK-2';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.offset=800;
parameters.sweep=1700;
parameters.npoints_f1=512;
parameters.npoints_f2=512;
parameters.zerofill_f1=2048;
parameters.zerofill_f2=2048;
parameters.spins='1H';
parameters.axis_units='ppm';

% COSY simulation
fid=cosy(spin_system,parameters);

% Apodization
fid=apodization(fid,'sinebell-2d');

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill_f2,parameters.zerofill_f1));

% Plotting
contour_plot(spin_system,real(spectrum),parameters,20,[0.02 0.2],2,256,6);

end

