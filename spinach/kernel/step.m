% A propagation step function using Krylov propagation for large matrices
% and Matlab's expm for small matrices. Arguments
%
%      L   - the Liouvillian to be used for propagation
%      rho - the state vector (or a horizontal stack thereof) to be propagated
%      time_step - the length of the time step to take
%
% ilya.kuprov@oerc.ox.ac.uk

function rho=step(spin_system,L,rho,time_step)

% Patch in thermal relaxation
if strcmp(spin_system.rlx.equilibrium,'thermal')
    rho(1,:)=1;
end

% Perform the propagation step
if ~any(strcmp(spin_system.sys.disable,'expv'))
    
    % Use Krylov propagation
    for n=1:size(rho,2)
        rho(:,n)=expv(-1i*L*time_step,rho(:,n));
    end
    
else
    
    % On special request, use expm
    rho=expm(full(-1i*L*time_step))*rho;
    
end

% Patch out thermal relaxation
if strcmp(spin_system.rlx.equilibrium,'thermal')
    rho(1,:)=0;
end

end

% Chemistry is physics without thought; mathematics is physics without purpose.
%
% Physics folklore
