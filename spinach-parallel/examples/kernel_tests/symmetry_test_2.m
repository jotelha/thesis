% A highly symmetric spin system from Andres Castillo.
%
% ilya.kuprov@oerc.ox.ac.uk

function symmetry_test_2()

% Spin system specification
sys.magnet=9.4;
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'};

inter.zeeman.scalar={0.89 0.89 0.89 0.895 0.895 0.895 1.16 1.16 1.16 1.2 1.39 1.71 3.85};
inter.coupling.scalar=num2cell(...
      [  0         0         0         0         0         0         0         0         0         0         0    6.7200         0
         0         0         0         0         0         0         0         0         0         0         0    6.7200         0
         0         0         0         0         0         0         0         0         0         0         0    6.7200         0
         0         0         0         0         0         0         0         0         0         0         0    6.6400         0
         0         0         0         0         0         0         0         0         0         0         0    6.6400         0
         0         0         0         0         0         0         0         0         0         0         0    6.6400         0
         0         0         0         0         0         0         0         0         0         0         0         0    6.0800
         0         0         0         0         0         0         0         0         0         0         0         0    6.0800
         0         0         0         0         0         0         0         0         0         0         0         0    6.0800
         0         0         0         0         0         0         0         0         0         0   14.0850    5.1000    8.3000
         0         0         0         0         0         0         0         0         0   14.0850         0    8.3050    6.4000
    6.7200    6.7200    6.7200    6.6400    6.6400    6.6400         0         0         0    5.1000    8.3050         0         0
         0         0         0         0         0         0    6.0800    6.0800    6.0800    8.3000    6.4000         0         0]/2);

% Symmetry
sys.sym_group={'C3v','C3v','C3v'};
sys.sym_spins={[1 2 3],[4 5 6],[7 8 9]};

% Basis set
bas.mode='IK-2';
     
% Sequence parameters
parameters.offset=800;
parameters.sweep=2000;
parameters.npoints=2048;
parameters.zerofill=8196;
parameters.spins='1H';
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Spinach code
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=pulse_acquire(spin_system,parameters);

% Apodization
fid=apodization(fid,'crisp-1d');

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
plot_1d(spin_system,real(spectrum),parameters);     

end
