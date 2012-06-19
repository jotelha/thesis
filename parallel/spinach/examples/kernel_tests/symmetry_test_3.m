% 1H NMR spectrum of valine with six equivalent spins.
%
% ilya.kuprov@oerc.ox.ac.uk

function symmetry_test_3()

sys.magnet=11.7;
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H'};
inter.zeeman.scalar={3.5950 2.2580 1.0270 1.0270 1.0270 0.9760 0.9760 0.9760};
inter.coupling.scalar{1,2}=4.34;
inter.coupling.scalar{2,3}=7.00;
inter.coupling.scalar{2,4}=7.00;
inter.coupling.scalar{2,5}=7.00;
inter.coupling.scalar{2,6}=7.00;
inter.coupling.scalar{2,7}=7.00;
inter.coupling.scalar{2,8}=7.00;
inter.coupling.scalar{8,8}=0.00;

sys.sym_spins={[3 4 5 6 7 8]};
sys.sym_group={'S6'};

% Set the simulation parameters
bas.mode='IK-2';    

% Sequence parameters
parameters.offset=1000;
parameters.sweep=3000;
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
fid=apodization(fid,'exponential-1d',5);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
plot_1d(spin_system,real(spectrum),parameters);

end

