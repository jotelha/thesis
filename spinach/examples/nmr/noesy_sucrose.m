% NOESY spectrum of sucrose (magnetic parameters computed with DFT).
%
% ilya.kuprov@oerc.ox.ac.uk

function noesy_sucrose()

% Read the spin system properties (vacuum DFT calculation)
options.min_j=1.0;
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\sucrose.log'),{{'H','1H'}},31.8,options);

% Set the simulation parameters
sys.magnet=5.9;
inter.relaxation='redfield';
inter.tau_c=200e-12;
bas.mode='IK-2';
bas.space_level=3;
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
parameters.tmix=0.5;
parameters.axis_units='ppm';

% Hypercomplex NOESY simulation
[cos_term,sin_term]=noesy(spin_system,parameters);

% Apodization
cos_term=apodization(cos_term,'gaussian-2d',5);
sin_term=apodization(sin_term,'gaussian-2d',5);

% F2 Fourier transform
f1_cos=real(fftshift(fft(cos_term,parameters.zerofill_f2,1),1));
f1_sin=real(fftshift(fft(sin_term,parameters.zerofill_f2,1),1));

% F1 Fourier transform
spectrum=fftshift(fft(f1_cos-1i*f1_sin,parameters.zerofill_f1,2),2);

% Plotting
contour_plot(spin_system,-real(spectrum),parameters,20,[0.05 0.2],2,256,6);

end

