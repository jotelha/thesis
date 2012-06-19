% HSQC spectrum of a C-H system in the presence of RDCs.
%
% ilya.kuprov@oerc.ox.ac.uk

function rdc_test_1()

% Spin system parameters
sys.magnet=5.9;
sys.isotopes={'1H','13C'};
inter.zeeman.scalar={5.0 65.0};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=140;
inter.coordinates={[0.0 0.0 0.0];
                   [0.6 0.7 0.8]};
inter.order_matrix=diag([1e-3 2e-3 -3e-3]);

% Basis set
bas.mode='complete';

% Spinach housekeeping
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
parameters.residual='on';

% HSQC simulation
fid=hsqc(spin_system,parameters);

% Apodization
fid=apodization(fid,'sinebell-2d');

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill_f2,parameters.zerofill_f1));

% Plotting
contour_plot(spin_system,abs(spectrum),parameters);

end

