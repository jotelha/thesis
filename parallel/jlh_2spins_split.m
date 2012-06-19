% W-band pulsed powder ESR spectrum of nitroxide radical.
%
% ilya.kuprov@oerc.ox.ac.uk

function jlh_2spins_split()

%jlh - set simulation data
title = 'jlh_2spins_split';
nodes = 3;
ppn = 1;

% Set the simulation parameters
sys.magnet=3.356;
sys.regime='powder';
bas.mode='ESR-2';
sys.tols.grid_rank=131;


% Interactions
sys.isotopes={'E','14N'};
inter.zeeman.matrix=cell(2,1);
inter.zeeman.matrix{1,1}=[2.0064953    0.0000748 0.0000187; 0.0001123    2.0057462    0.0000510; 0.0000024    0.0000522    2.0023775];
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{1,2}=1e6*[1.9000       0.0730 -0.4967; 0.0730       1.9759      -0.8768; -0.4967      -0.8768      77.1529];
%N-S distance 1.80 A, from N-S bond scan using orca

%inter.coupling.scalar=cell(2,2);
%inter.coupling.scalar{1,2}=1e6*27.4502;


% Set the sequence parameters
parameters.offset=10;
parameters.sweep=1e9;
parameters.npoints=1024;
parameters.zerofill=2048;
parameters.spins='E';
parameters.axis_units='Gauss';
parameters.derivative=0;

% jlh - parallel extension
parameters.nodes = nodes;
parameters.ppn = ppn;
right_now = datestr(now,'yyyymmdd_HHMMSS');
realtitle = [title '_' right_now];
parameters.title = realtitle;

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

%jlh - keep track of execution time

jlh_master_pulse_acquire(spin_system,parameters);
end