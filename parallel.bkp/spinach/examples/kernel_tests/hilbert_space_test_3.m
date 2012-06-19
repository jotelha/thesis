% Parallelization test: multi-threaded evaluation of observables in Hilbert
% space time propagation. This test will only run on a Matlab cluster with
% at least 128 workers. Change 'amazon' to your own parallel configuration
% name in line 56.
% 
% Spin system of 3-phenylmethylene-1H,3H-naphtho-[1,8-c,d]-pyran-1-one.
% Source: Penchav, et al., Spec. Acta Part A, 78 (2011) 559-565.
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function hilbert_space_test_3()

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

% Spinach housekeeping
spin_system=create(sys,inter);

% Set conditions to high-field NMR
spin_system=secularity(spin_system,'nmr');

% Get the Hamiltonian
H=hs_hamiltonian(spin_system);

% Get the initial state
rho=(hs_operator(spin_system,'L+','all')+hs_operator(spin_system,'L-','all'))/2;

% Get the time stepping parameters
trajectory_duration=0.2;
[timestep,nsteps]=stepsize(H,trajectory_duration);

% Benchmark parallel propagation
for n=[128 64 32 16 8 4 2 1]
    matlabpool('amazon',n); pause(10);
    tic
        hs_evolution(spin_system,H,rho,rho,timestep,nsteps,'observable');
    toc
    matlabpool('close');
end

end
