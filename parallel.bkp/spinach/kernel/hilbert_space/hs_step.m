% A propagation step function in hilbert space. Arguments
%
%      H   - the Hamiltonian to be used for propagation
%      rho - density matrix to be propagated
%      timestep - the length of the time step to take
%
% luke.edwards@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function rho=hs_step(spin_system,H,rho,timestep)

% Decide the propagator calculation method
if (length(H)>spin_system.tols.small_matrix)
    
    % For large matrices use sparse exponentiation
    P=propagator(spin_system,H,timestep);
    
else
    
    % For small matrices use Matlab's expm
    P=expm(full(-1i*H*timestep));
    
end

% Perform the propagation step
rho=P*rho*P';

end

% "Faith in the supernatural begins as faith in the superiority of others."
%
% Ayn Rand

