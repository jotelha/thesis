% Figure 8.26 from Andrew Derome's "Modern NMR Techniques for Chemistry Research".
%
% ilya.kuprov@oerc.ox.ac.uk

function cosy90_derome()

% Spin system and interactions
sys.magnet=16.1;
sys.isotopes={'1H','1H','1H'};
inter.zeeman.scalar={3.70 3.92 4.50};
inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{2,3}=12;
inter.coupling.scalar{1,3}=4;
inter.coupling.scalar{3,3}=0;

% Basis set
bas.mode='complete';

% Sequence parameters
parameters.offset=2800;
parameters.sweep=700;
parameters.npoints_f1=1024;
parameters.npoints_f2=1024;
parameters.zerofill_f1=2048;
parameters.zerofill_f2=2048;
parameters.spins='1H';
parameters.axis_units='ppm';

% Spinach code
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=cosy(spin_system,parameters);

% Apodization
fid=apodization(fid,'sinebell-2d');

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill_f2,parameters.zerofill_f1));

% Plotting
contour_plot(spin_system,real(spectrum),parameters,20,[0.05 1.0],2,256,6);

end

