% Zero track elimination function. Inspects the first few steps in the
% system trajectory and drops the states that did not get populated to a
% user-specified tolerance. Arguments:
%
%      L -   the Liouvillian to be used for time propagation
%      rho - the initial state to be used for time propagation
%
%      projector - projector matrix into the reduced space, to be used as
%                  follows: L=P'*L*P,   rho=P'*rho;
%
% Further information is available in IK's JMR paper on the subject:
%
%      http://dx.doi.org/10.1016/j.jmr.2008.08.008
%
% ilya.kuprov@oerc.ox.ac.uk

function projector=zte(spin_system,L,rho)

% Click forward for output
spin_system=click(spin_system,'forward');

% Run Zero Track Elimination
if any(strcmp(spin_system.sys.disable,'zte'))
    
    % Skip if instructed to do so by the user
    report(spin_system,'zte: WARNING - zero track elimination disabled, basis set left unchanged.');
    
    % Return a unit matrix
    projector=speye(size(L));

elseif nnz(rho)/numel(rho)>spin_system.tols.zte_maxden
    
    % Skip if the benefit is likely to be minor
    report(spin_system,'zte: WARNING - too few zeros in the state vector, basis set left unchanged.');
    
    % Return a unit matrix
    projector=speye(size(L));
    
elseif norm(rho,'inf')<spin_system.tols.zte_tol
    
    % Skip if the state vector norm is too small for Krylov procedure
    report(spin_system,'zte: WARNING - state vector norm below drop tolerance, basis set left unchanged.');
    
    % Return a unit matrix
    projector=speye(size(L));
    
else
    
    % Estimate the Larmor time step
    timestep=stepsize(L);
    
    % Report to the user
    report(spin_system,['zte: dropping states with amplitudes below ' num2str(spin_system.tols.zte_tol)...
                        ' within the first ' num2str(timestep*spin_system.tols.zte_nsteps) ' seconds of the trajectory.']);
    report(spin_system,['zte: ' num2str(spin_system.tols.zte_nsteps) ' steps taken, ' num2str(timestep) ' seconds each.']);
    
    % Compute trajectory steps with the Krylov technique
    trajectory=zeros(length(rho),spin_system.tols.zte_nsteps);
    for n=1:spin_system.tols.zte_nsteps
        trajectory(:,n)=rho;
        rho=step(spin_system,L,rho,timestep);
    end
    
    % Determine which tracks are zero
    zero_track_mask=(max(abs(trajectory),[],2)<spin_system.tols.zte_tol);
    
    % Always keep the unit state
    zero_track_mask(1)=false();
    
    % Take a unit matrix and delete the columns corresponding to zero tracks
    projector=speye(size(L));
    projector(:,zero_track_mask)=[];
     
    % Report back to the user
    report(spin_system,['zte: old active space dimension ' num2str(size(L,2))]);
    report(spin_system,['zte: new active space dimension ' num2str(size(projector,2))]);
    
end

end

% Every great scientific truth goes through three stages. First, people say
% it conflicts with the Bible. Next they say it had been discovered before.
% Lastly they say they always believed it. 
%
% Louis Agassiz

