% Dima Budker's zero-field spectroscopy -- 13C methanol.
%
% ilya.kuprov@oerc.ox.ac.uk

function zero_field_methanol()

% Spin system
sys.isotopes={'1H','1H','1H','13C'};
sys.sym_group={'C3v'};
sys.sym_spins={[1 2 3]};

% Interactions
sys.magnet=0;
inter.zeeman.scalar={0 0 0 0};
inter.coupling.scalar{1,4}=140;
inter.coupling.scalar{2,4}=140;
inter.coupling.scalar{3,4}=140;
inter.coupling.scalar{4,4}=0;

% Basis set
bas.mode='complete';

% Sequence parameters
parameters.sweep=1000;
parameters.npoints=4096;
parameters.zerofill=4096;

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

