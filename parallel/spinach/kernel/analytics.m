% Norm dynamics analysis function. Analyzes the system trajectory and
% plots the time dependence of the density matrix norm, partitioned
% into subspaces according to the correlation order of the correspon-
% ding states. This is useful in deciding whether a particular basis
% set is large enough to correctly reproduce the system dynamics.
%
% gareth.charnock@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function analytics(spin_system,trajectory)

% Determine the spin order figures for each state
state_orders=sum(logical(spin_system.bas.basis),2);

% Preallocate the norm trajectory array
norm_trajectory=zeros(max(state_orders),size(trajectory,2));

% Loop over the state orders that are present
for n=1:max(state_orders)
    
    % Get the subspace mask
    subspace_mask=(state_orders==n);
    
    % Get the part of the trajectory belonging to the subspace
    subspace_trajectory=trajectory(subspace_mask,:);
    
    % Get the norm of the trajectory
    norm_trajectory(n,:)=sqrt(sum(subspace_trajectory.*conj(subspace_trajectory),1));
    
end

% Do the plotting
plot(norm_trajectory');

end

% And this is the whole shabby secret: to some men, the sight of an
% achievement is a reproach, a reminder that their own lives are
% irrational, and that there is no loophole - no escape from reason
% and reality.
%
% Ayn Rand

