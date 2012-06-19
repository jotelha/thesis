% 1H NMR spectrum of a protein fragment. The calculation includes
% the full Redfield relaxation superoperator, taking into account 
% all DD, CSA and DD-CSA terms.
%
% kpervushin@ntu.edu.sg
% ilya.kuprov@oerc.ox.ac.uk

function pulse_acquire_peptide()

% Magnetic induction
sys.magnet=16.4;

% Isotopes
sys.isotopes={'15N','1H', '13C','1H','13C','15N','1H', '13C',...
              '1H', '13C','1H', '1H','13C','1H', '13C','1H',...
              '1H', '1H', '13C','1H','1H', '1H', '13C','15N',...
              '1H', '13C','1H', '1H','1H', '1H', '1H', '1H'}; 

% Chemical shielding tensors
inter.zeeman.eigs=  {[19.2   91.6   -110.8] + 119.8
                     [ 1.0    6.4     -7.5] +   8.2
                     [ 0.4   14.3    -14.7] +  56.6
                     [ 5.0    5.0    -10.0] +   4.3
                     [97.6   -4.8    -92.8] + 176.3
                     [43.3   82.9   -124.1] + 122.2
                     [ 5.0    6.4     -6.9] +   8.2
                     [ 2.9   12.7    -15.4] +  55.6
                     [ 5.0    5.0    -10.0] +   4.3
                     [-2.0   11.7     -9.8] +  42.3
                     [ 5.0    5.0    -10.0] +   1.6
                     [ 5.0    5.0    -10.0] +   1.5
                     [-2.0   11.7     -9.8] +  26.8
                     [ 5.0    5.0    -10.0] +   1.5
                     [-2.0   11.7     -9.8] +  24.7
                     [ 5.0    5.0    -10.0] +   0.7
                     [ 5.0    5.0    -10.0] +   0.7
                     [ 5.0    5.0    -10.0] +   0.7
                     [-2.0   11.7     -9.8] +  24.1
                     [ 5.0    5.0    -10.0] +   0.7
                     [ 5.0    5.0    -10.0] +   0.7
                     [ 5.0    5.0    -10.0] +   0.7
                     [-5.6   95.0    -90.1] + 176.9
                     [-136.0 76.7     59.1] + 120.7
                     [-4.6    4.1      0.5] +   8.3
                     [-12.5   11.7     0.7] +  57.4
                     [ 5.0    5.0    -10.0] +   4.3
                     [ 5.0    5.0    -10.0] +   2.0
                     [ 5.0    5.0    -10.0] +   0.8
                     [ 5.0    5.0    -10.0] +   7.1
                     [ 5.0    5.0    -10.0] +   0.8
                     [ 5.0    5.0    -10.0] +   7.0};

% Euler angles for shielding tensors (set to zero for lack of data)
inter.zeeman.euler=cell(32,1);

% Coupling constants
inter.coupling.scalar{1,2}=-92; 
inter.coupling.scalar{1,3}=-11; 
inter.coupling.scalar{3,4}=140; 
inter.coupling.scalar{2,4}=8; 
inter.coupling.scalar{3,5}=55; 
inter.coupling.scalar{5,6}=15; 
inter.coupling.scalar{6,7}=-92; 
inter.coupling.scalar{6,8}=-11; 
inter.coupling.scalar{8,9}=140; 
inter.coupling.scalar{7,9}=8; 
inter.coupling.scalar{8,10}=35; 
inter.coupling.scalar{10,11}=135; 
inter.coupling.scalar{25,27}=8;
inter.coupling.scalar{10,12}=135; 
inter.coupling.scalar{9,12}=-3; 
inter.coupling.scalar{9,11}=12; 
inter.coupling.scalar{10,13}=30; 
inter.coupling.scalar{13,14}=130; 
inter.coupling.scalar{11,14}=-3;
inter.coupling.scalar{12,14}=12; 
inter.coupling.scalar{13,15}=20; 
inter.coupling.scalar{13,19}=20; 
inter.coupling.scalar{15,16}=120; 
inter.coupling.scalar{15,17}=120; 
inter.coupling.scalar{15,18}=120;
inter.coupling.scalar{19,20}=120; 
inter.coupling.scalar{19,21}=120; 
inter.coupling.scalar{19,22}=120; 
inter.coupling.scalar{8,23}=55; 
inter.coupling.scalar{23,24}=15; 
inter.coupling.scalar{24,25}=-92;
inter.coupling.scalar{24,26}=-11; 
inter.coupling.scalar{26,27}=140;
inter.coupling.scalar{32,32}=0;

% Coordinates
inter.coordinates= {[-3.138 -16.334   5.186]
                    [-2.456 -17.063   5.031]
                    [-3.030 -15.529   6.404]
                    [-3.925 -14.912   6.487]
                    [-1.814 -14.631   6.357]
                    [-1.984 -13.367   6.738] 
                    [-2.896 -13.043   7.027]
                    [-0.868 -12.424   6.749]
                    [ 0.029 -13.000   6.522]
                    [-1.029 -11.366   5.656]
                    [-1.945 -10.820   5.881]
                    [-0.175 -10.695   5.742]
                    [-1.115 -11.782   4.187]
                    [-2.023 -12.370   4.050]
                    [-1.201 -10.547   3.315]
                    [-1.291 -10.842   2.269]
                    [-2.074  -9.961   3.602]
                    [-0.301  -9.947   3.446]
                    [ 0.061 -12.628   3.769]
                    [-0.039 -12.919   2.723]
                    [ 0.980 -12.058   3.897]
                    [ 0.100 -13.524   4.388]
                    [-0.644 -11.752   8.105]
                    [ 0.626 -11.697   8.508]
                    [ 1.333 -12.157   7.952] 
                    [ 1.064 -11.002   9.722]
                    [ 0.318 -10.251   9.982]
                    [-6.635 -14.491   7.235]
                    [ 3.330 -12.701   8.884]
                    [ 2.550 -13.892   6.500]
                    [-5.600 -15.99    6.600]
                    [-4.900 -13.199   7.550]};

% Basis set
bas.mode='IK-2';
bas.space_level=2;

% Relaxation theory
inter.relaxation='redfield';
inter.tau_c=1e-9;

% Pulse sequence parameters
parameters.offset=3000;
parameters.sweep=10000;
parameters.npoints=8192;
parameters.zerofill=32768;
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



