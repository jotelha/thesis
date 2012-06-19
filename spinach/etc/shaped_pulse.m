% General shaped pulse function. Input variables:
%
% L:            the drift Liouvillian
% offset:       pulse offset frequency       (Hz)
% phases:       [phi1 phi2 ... phiN]         (degrees)
% amplitudes:   [a1 a2 ... aN]               (rad/s)
% duration:     pulse duration               (seconds)
%
% spins:       '1H', '13C' etc
%
% ilya.kuprov@oerc.ox.ac.uk

function answer=shaped_pulse(spin_system,L,rho,spins,offset,phases,amplitudes,duration)

% Get the basic operators
Lm=operator(spin_system,'L-',spins);
Lx=(Lm'+Lm)/2; Ly=(Lm'-Lm)/2i; Lz=(Lm'*Lm-Lm*Lm')/2;

% Decide the time step
nsteps=length(amplitudes);
dt=duration/nsteps;

% Patch in thermal relaxation
if strcmp(spin_system.rlx.equilibrium,'thermal')
    rho(1,:)=1;
end

% Preallocate the answer
answer=zeros(size(rho));

% Ignore symmetry if lacking or forbidden
if (~isfield(spin_system.bas,'irrep'))||any(strcmp(spin_system.sys.disable,'symmetry'))
    spin_system.bas.irrep.projector=speye(size(L));
end

% Loop over the symmetry subspaces
for subs=1:numel(spin_system.bas.irrep)
    
    % Project the operators into the current subspace
    L_subs=spin_system.bas.irrep(subs).projector'*L*spin_system.bas.irrep(subs).projector;
    Lx_subs=spin_system.bas.irrep(subs).projector'*Lx*spin_system.bas.irrep(subs).projector;
    Ly_subs=spin_system.bas.irrep(subs).projector'*Ly*spin_system.bas.irrep(subs).projector;
    Lz_subs=spin_system.bas.irrep(subs).projector'*Lz*spin_system.bas.irrep(subs).projector;
    
    % Project the state vector into the current subspace
    rho_subs=spin_system.bas.irrep(subs).projector'*rho;
    
    % Run the propagation
    for n=1:nsteps
        
        % Generate the evolution step operator
        pulse_op=amplitudes(n)*(Lx_subs*cosd(phases(n))+Ly_subs*sind(phases(n)));
        total_op=pulse_op+L_subs+2*pi*offset*Lz_subs;
        
        % Apply the pulse chunk
        rho_subs=step(spin_system,total_op,rho_subs,dt);
        
    end
    
    % Add the subspace to the total
    answer=answer+spin_system.bas.irrep(subs).projector*rho_subs;
    
end

% Patch out thermal relaxation
if strcmp(spin_system.rlx.equilibrium,'thermal')
    answer(1,:)=0;
end

end

% When IK proposed the fibre etching technique featured in the 2004 JMR
% paper (http://dx.doi.org/10.1016/j.jmr.2004.08.017), he could see terror
% in his supervisors's eyes - Peter Hore reasonably thought that the
% notoriously eccentric Russian student could not possibly be trusted with
% boiling hydrofluoric acid in a high-power laser lab. The Chemistry 
% Department safety officer held a similar view. It is not entirely clear
% how IK got hold of several milliliters of concentrated HF and a heat gun
% on a Saturday night in the PTCL Teaching Lab, but the photographs of the
% resulting fibre tip were left on Peter's table on Monday. The paper was
% accepted by JMR without revisions.

