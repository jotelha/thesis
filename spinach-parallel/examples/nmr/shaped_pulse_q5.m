% 90-degree Q5 pulse on a chain of 31 strongly coupled protons.
%
% ilya.kuprov@oerc.ox.ac.uk

function shaped_pulse_q5()

% Magnetic field
sys.magnet=14.1;

% Isotopes
sys.isotopes={ '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H',...
               '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H',...
               '1H','1H','1H','1H','1H','1H','1H','1H','1H'};

% Zeeman interactions
inter.zeeman.scalar=num2cell(linspace(-4,4,31));

% Couplings
inter.coupling.scalar=cell(31);
for n=1:30
    inter.coupling.scalar{n,n+1}=10;
end

% Basis set
bas.mode='IK-2';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=secularity(spin_system,'nmr');
L=h_superop(spin_system);

% Initial state
rho=equilibrium(spin_system);

% Detection state
coil=state(spin_system,'L+','1H');

% Q5 pulse setup
npoints=200;
duration=0.012;
offset=960;
[amplitudes,phases,integral]=read_wave([spin_system.sys.root_dir '\etc\wave\q5_1000.pk'],npoints);

% Pulse calibration
calibration_factor=8*npoints*(pi/2)/(integral*duration);
amplitudes=calibration_factor*amplitudes;

% Pulse execution
rho=shaped_pulse(spin_system,L,rho,'1H',offset,phases,amplitudes,duration);

% Detection
sweep=7000; npoints=2048; zerofill=16384;
fid=evolution(spin_system,L,coil,rho,1/sweep,npoints,'observable');

% Apodization
fid=apodization(fid,'crisp-1d');

% Fourier transform
spectrum=fftshift(fft(fid,zerofill));

% Plotting
plot(imag(spectrum));

end

