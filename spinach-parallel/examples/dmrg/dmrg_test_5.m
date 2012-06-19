% Time evolution of a simple linear spin chain with nearest-neighbor
% coupling at zero field. Relaxation is switched on to keep the high
% spin orders down.
%
% The L+ state on the first spin propagated forward using Spinach CI 
% representation. At each step, the state vector is compressed into 
% the MPS representation and MPS dimension statistics are printed to
% the console.
%
% ilya.kuprov@oerc.ox.ac.uk

function dmrg_test_5()

% Magnetic field
sys.magnet=14.1;

% Isotopes
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H'};

% Couplings
inter.coupling.scalar=cell(8);
for n=1:7
    inter.coupling.scalar{n,n+1}=10+n;
end

% Relaxation settings
inter.relaxation='t1_t2';
inter.r1_rates=10*ones(1,8);
inter.r2_rates=10*ones(1,8);

% Basis set
bas.mode='complete';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Trajectory generation and analysis
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

