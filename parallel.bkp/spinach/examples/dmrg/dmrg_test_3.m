% Solid state pulse-acquire simulation of wtphene spin system in the absence
% of relaxation.
%
% A standard NMR pulse-acquire experiment is propagated forward using Spinach
% CI representation. At each step, the state vector is compressed into the MPS
% representation and MPS dimension statistics are printed to the console.
%
% ilya.kuprov@oerc.ox.ac.uk

function dmrg_test_3()

% Read the spin system properties (vacuum DFT calculation)
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\wtphene.log'),{{'H','1H'},{'F','19F'}},[31.8 400]);

% Set the simulation parameters
sys.magnet=14.1;
sys.tols.prox_cutoff=Inf;
bas.mode='complete';

% Spinach code
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Operators, states and timing
[H_iso,Q]=h_superop(secularity(spin_system,'nmr'));
H=H_iso+orientation(Q,[pi/3 pi/4 pi/5]);
rho=state(spin_system,'L+','1H');
[timestep,nsteps]=stepsize(H,0.05);

% MPS dimension statistics at each step
for n=1:nsteps
    
    % Take a good step forward
    rho=step(spin_system,H,rho,10*timestep);
    
    % Get the state vector into a tensor form
    cits=spinach2cits(spin_system,full(rho));
    
    % Get the state tensor into an MPS
    cits2mps(spin_system,cits,1e-6); disp(' ');
    
end

end

