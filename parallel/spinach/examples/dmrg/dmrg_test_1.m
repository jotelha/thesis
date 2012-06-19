% Liquid state pulse-acquire simulation of wtphene spin system. Relaxation
% switched on to keep the high spin orders down.
%
% A standard NMR pulse-acquire experiment is propagated forward using Spinach
% CI representation. At each step, the state vector is compressed into the MPS
% representation and MPS dimension statistics are printed to the console.
%
% ilya.kuprov@oerc.ox.ac.uk

function dmrg_test_1()

% Read the spin system properties (vacuum DFT calculation)
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\wtphene.log'),{{'H','1H'},{'F','19F'}},[31.8 400]);

% Set the simulation parameters
sys.magnet=14.1;
bas.mode='complete';

% Relaxation settings
inter.relaxation='t1_t2';
inter.r1_rates=50*ones(1,8);
inter.r2_rates=50*ones(1,8);

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Operators, states and timing
L=h_superop(secularity(spin_system,'nmr'))+1i*r_superop(spin_system);
rho=state(spin_system,'L+',1); [timestep,nsteps]=stepsize(L,0.05);

% MPS dimension statistics at each step
for n=1:nsteps
    
    % Take a good step forward
    rho=step(spin_system,L,rho,10*timestep);
    
    % Get the state vector into a tensor form
    cits=spinach2cits(spin_system,full(rho));
    
    % Get the state tensor into an MPS
    cits2mps(spin_system,cits,1e-6); disp(' ');
    
end

end

