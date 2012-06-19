% HSQC spectrum of sucrose (magnetic parameters computed with DFT).
%
% matthew.krzystyniak@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function hsqc_sucrose()

% Read the spin system properties (vacuum DFT calculation)
options.min_j=3.0; options.no_xyz=1; options.dilute={'13C'};
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\sucrose.log'),{{'H','1H'},{'C','13C'}},[31.8 182.1],options);

% Set the magnet field
sys.magnet=5.9;

% Set the basis
bas.mode='IK-2';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.sweep_f1=2500;
parameters.sweep_f2=950;
parameters.offset_f1=4300;
parameters.offset_f2=1100;
parameters.npoints_f1=128;
parameters.npoints_f2=128;
parameters.zerofill_f1=512;
parameters.zerofill_f2=512;
parameters.spins_f1='13C';
parameters.spins_f2='1H';
parameters.J=140;
parameters.axis_units='ppm';

% HSQC simulation
[P_term,N_term]=hsqc(spin_system,parameters);

% Apodization
P_term=apodization(P_term,'sinebell-2d');
N_term=apodization(N_term,'sinebell-2d');

% F2 Fourier transform (directly detected dimension)
f1_P=fftshift(fft(P_term,parameters.zerofill_f2,1),1);
f1_N=fftshift(fft(N_term,parameters.zerofill_f2,1),1);

% Form States signal
fid=f1_P+conj(f1_N);

% F1 Fourier transform (indirectly detected dimension)
spectrum=fftshift(fft(fid,parameters.zerofill_f1,2),2);

% Plotting
contour_plot(spin_system,real(spectrum),parameters,20,[0.05 1.0],2,256,6,'positive');

end

