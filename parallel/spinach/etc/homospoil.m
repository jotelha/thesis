% Emulates a strong homospoil pulse -- only the zero-frequency states
% survive the process.
%
% The function supports horizontal stacks of state vectors.
%
% ilya.kuprov@oerc.ox.ac.uk

function rho=homospoil(spin_system,rho,zqc_flag)

% Pull the projection information from the basis
[~,M]=lin2lm(spin_system.bas.basis);

% Set the defaults
if nargin==2
    zqc_flag='destroy';
end

% Filter the state vector
switch zqc_flag
    case 'keep'
        % Find the states that have zero carrier frequency and kill everything else
        rho(abs(sum(repmat(spin_system.inter.basefrq,spin_system.bas.nstates,1).*M,2))>1e-6,:)=0;
    case 'destroy'
        % Find the longitudinal states and kill everything else
        rho(sum(abs(M),2)>0,:)=0;
    otherwise
        error('homospoil: unknown ZQC flag.');
end

end

% Rocket science has been mythologized all out of proportion to its true difficulty.
%
% John Carmack 