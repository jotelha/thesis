% COSY-45 spectrum of rotanone.
% 
% R. Nunlist and J.Ralph, J. Heterocyclic Chem. 25, 351 (1988)
%
% matthew.krzystyniak@oerc.ox.ac.uk

function cosy45_rotanone()

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','1H'};
sys.sym_group={'C3v','C3v','C3v'};
sys.sym_spins={[14 15 16],[17 18 19],[20 21 22]};

% Interactions
sys.magnet=5.9;
inter.zeeman.scalar={6.72 6.40 4.13 4.56 4.89 6.46 7.79 3.79 2.91...
                     3.27 5.19 4.89 5.03 1.72 1.72 1.72 3.72 3.72...
                     3.72 3.76 3.76 3.76};
inter.coupling.scalar{3,4}=12.1; 
inter.coupling.scalar{4,5}=3.1; 
inter.coupling.scalar{3,5}=1.0; 
inter.coupling.scalar{3,8}=1.0; 
inter.coupling.scalar{1,8}=1.0;
inter.coupling.scalar{6,7}=8.6; 
inter.coupling.scalar{5,8}=4.1; 
inter.coupling.scalar{7,9}=0.7; 
inter.coupling.scalar{7,10}=0.7; 
inter.coupling.scalar{9,10}=15.8;
inter.coupling.scalar{10,11}=9.8; 
inter.coupling.scalar{9,11}=8.1; 
inter.coupling.scalar{13,14}=1.5; 
inter.coupling.scalar{12,14}=0.9; 
inter.coupling.scalar{22,22}=0;

% Relaxation theory
inter.relaxation='t1_t2';
inter.r1_rates=1.0*ones(22,1);
inter.r2_rates=1.5*ones(22,1);

% Basis set
bas.mode='IK-2';

% Sequence parameters
parameters.offset=1200;
parameters.sweep=2000;
parameters.npoints_f1=512;
parameters.npoints_f2=512;
parameters.zerofill_f1=2048;
parameters.zerofill_f2=2048;
parameters.spins='1H';
parameters.axis_units='ppm';
parameters.angle=pi/4;

% Spinach code
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=cosy(spin_system,parameters);

% Apodization
fid=apodization(fid,'sinebell-2d');

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill_f2,parameters.zerofill_f1));

% Plotting
contour_plot(spin_system,abs(spectrum),parameters,20,[0.01 0.2],2,256,6,'both');

end

