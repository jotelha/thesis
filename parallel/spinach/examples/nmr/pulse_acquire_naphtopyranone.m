% NMR spectrum of 3-phenylmethylene-1H,3H-naphtho-[1,8-c,d]-pyran-1-one.
%
% Source: Penchav, et al., Spec. Acta Part A, 78 (2011) 559-565.
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function pulse_acquire_naphtopyranone()

% Magnetic induction
sys.magnet=14.095;

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'};

% Chemical shifts
inter.zeeman.scalar={8.345,7.741,8.097,8.354,7.784,8.330,7.059,7.941,7.466,7.326,7.466,7.941};

% Scalar couplings
inter.coupling.scalar=cell(12,12);
inter.coupling.scalar{1,2}=7.8;
inter.coupling.scalar{1,3}=0.9;
inter.coupling.scalar{2,3}=7.8;
inter.coupling.scalar{4,5}=8.4;
inter.coupling.scalar{4,6}=1.2;
inter.coupling.scalar{5,6}=7.2;
inter.coupling.scalar{8,9}=7.8;
inter.coupling.scalar{8,10}=1.2;
inter.coupling.scalar{9,10}=7.8;
inter.coupling.scalar{10,11}=7.8;
inter.coupling.scalar{10,12}=1.2;
inter.coupling.scalar{11,12}=7.8;

% Basis set
bas.mode='IK-2';

% Sequence parameters
parameters.offset=4600;
parameters.sweep=1000;
parameters.npoints=8192/2;
parameters.zerofill=65536/2;
parameters.spins='1H';
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Spinach code
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=pulse_acquire(spin_system,parameters);

% Apodization
fid=apodization(fid,'gaussian-1d',10);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
plot_1d(spin_system,real(spectrum),parameters);

end
