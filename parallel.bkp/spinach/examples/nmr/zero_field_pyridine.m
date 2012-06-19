% Dima Budker's zero-field spectroscopy -- 15N pyridine.
%
% ilya.kuprov@oerc.ox.ac.uk

function zero_field_pyridine()

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','15N'};

% Interactions
inter.zeeman.scalar={0 0 0 0 0 0};
inter.coupling.scalar{1,2}=4.9; 
inter.coupling.scalar{4,5}=4.9;
inter.coupling.scalar{1,4}=1.0; 
inter.coupling.scalar{2,5}=1.0;
inter.coupling.scalar{1,3}=1.8; 
inter.coupling.scalar{3,5}=1.8;
inter.coupling.scalar{1,5}=-0.1;
inter.coupling.scalar{2,3}=7.7; 
inter.coupling.scalar{3,4}=7.7;
inter.coupling.scalar{2,4}=1.4;
inter.coupling.scalar{1,6}=-10.8; 
inter.coupling.scalar{5,6}=-10.8;
inter.coupling.scalar{2,6}=-1.5; 
inter.coupling.scalar{4,6}=-1.5;
inter.coupling.scalar{6,6}=0;

% Basis set
bas.mode='complete';

% Sequence parameters
parameters.sweep=150;
parameters.npoints=512;
parameters.zerofill=1024;

% Spinach code
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=zerofield(spin_system,parameters);

% Apodization
fid=apodization(fid-mean(fid),'exponential-1d',10);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
plot_1d(spin_system,real(spectrum),parameters);

end

