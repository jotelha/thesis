% Time evolution function. Performs all types of time propagation with
% automatic trajectory level state space restriction. Arguments:
%
%      L      - the Liouvillian to be used during evolution
%
%      rho    - the initial state vector or a horizontal stack thereof
% 
%      output - is a string giving the type of evolution that is required
%
%                'final' - returns the final state vector or a stack thereof
%
%                'trajectory' - returns the stack of state vectors giving
%                               the trajectory of the system starting from
%                               rho with the user-specified number of steps
%                               and step length.
%
%                'refocus' - evolves the first vector for zero steps, 
%                            second vector for one step, third vector for
%                            two steps, etc., consistent with the second 
%                            stage of evolution in the indirect dimension
%                            after a refocusing pulse.
%
%                'observable' - returns the time dynamics of an observable
%                               as a vector (if starting from a single initial
%                               state) or a matrix (of starting from a stack
%                               of initial states).
%
%                'multichannel' - returns the time dynamics of several
%                                 observables as rows of a matrix. Note
%                                 that destination state screening may be
%                                 less efficient when there are multiple
%                                 destinations to screen against.
%
%      coil   - the detection state, used when 'observable' is specified as
%               the output option. If 'multichannel' is selected, the coil
%               should contain multiple columns corresponding to individual
%               observable vectors.
%
%      destination - the state to be used for destination state screening.
%
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk
% ohad.levinkron@weizmann.ac.il

function answer=evolution(spin_system,L,coil,rho,timestep,nsteps,output,destination)

% Click forward for output
spin_system=click(spin_system,'forward');

% Normalize the coils
if ~isempty(coil)
    for n=1:size(coil,2)
        coil(:,n)=coil(:,n)/norm(coil(:,n));
    end
    %jlh - reducing 
    %report(spin_system,'evolution: coils have been normalized.');
end

% Decide the screening procedure
if any(strcmp(spin_system.sys.disable,'dss'))
    report(spin_system,'evolution: WARNING - destination state screening disabled by the user.');
    screen_state=rho;
elseif exist('destination','var')
    report(spin_system,'evolution: using the user-supplied state as destination filter...');
    screen_state=destination;
elseif any(strcmp(output,{'observable','dec'}))
    %jlh - reducing output
    %report(spin_system,'evolution: using coil state as destination filter...');
    screen_state=coil;
else
    screen_state=rho;
end

% Patch in thermal relaxation
if strcmp(spin_system.rlx.equilibrium,'thermal')
    rho(1,:)=1; screen_state(1,:)=1;
end

% Reduce the problem dimension
%jtlh - reduce output
% report(spin_system,'evolution: trying to reduce the problem dimension...');
projectors=reduce(spin_system,L,screen_state);
        
% Run the evolution
switch output
    
    case 'final'
        
        % Preallocate the answer
        answer=zeros(size(rho));
               
        % Loop over independent subspaces
        report(spin_system,'evolution: propagating the system...');
        for subs=1:length(projectors)
            
            % Project into the current subspace
            L_subs=projectors{subs}'*L*projectors{subs};
            rho_subs=projectors{subs}'*rho;
            
            % Propagate the system
            if (nnz(L_subs)<spin_system.tols.krylov_switchover)||any(strcmp(spin_system.sys.disable,'krylov'))
                
                % Use exponential propagation
                P=propagator(spin_system,L_subs,timestep);
                for n=1:nsteps
                    rho_subs=P*rho_subs;
                end
                
            else
                
                % For very large subspaces use Krylov propagation
                report(spin_system,['evolution: WARNING - Krylov propagation used in subspace #' num2str(subs)]);
                for n=1:nsteps
                    rho_subs=step(spin_system,L_subs,rho_subs,timestep);
                end
                
            end
            
            % Add the subspace to the total
            answer=answer+projectors{subs}*rho_subs;
            
        end
        
    case 'trajectory'
        
        % Preallocate the answer and set the starting point
        answer=zeros(size(rho,1),nsteps+1); answer(:,1)=rho;
               
        % Loop over independent subspaces
        report(spin_system,'evolution: propagating the system...');
        for subs=1:length(projectors)
            
            % Project into the current subspace
            L_subs=projectors{subs}'*L*projectors{subs};
            rho_subs=projectors{subs}'*rho;
            
            % Propagate the system
            if (nnz(L_subs)<spin_system.tols.krylov_switchover)||any(strcmp(spin_system.sys.disable,'krylov'))
                
                % Use exponential propagation
                P=propagator(spin_system,L_subs,timestep);
                for n=1:nsteps
                    rho_subs=P*rho_subs;
                    answer(:,n+1)=answer(:,n+1)+projectors{subs}*rho_subs;
                end
                
            else
                
                % For very large subspaces use Krylov propagation
                report(spin_system,['evolution: WARNING - Krylov propagation used in subspace #' num2str(subs)]);
                for n=1:nsteps
                    rho_subs=step(spin_system,L_subs,rho_subs,timestep);
                    answer(:,n+1)=answer(:,n+1)+projectors{subs}*rho_subs;
                end
                
            end
            
        end
        
    case 'refocus'
        
        % Preallocate the answer
        answer=zeros(size(rho));
        
        % The first point is ready
        answer(:,1)=rho(:,1);
               
        % Loop over independent subspaces
        report(spin_system,'evolution: propagating the system...');
        for subs=1:length(projectors)
            
            % Project into the current subspace
            L_subs=projectors{subs}'*L*projectors{subs};
            rho_subs=projectors{subs}'*rho;
            
            % Propagate the system
            if (nnz(L_subs)<spin_system.tols.krylov_switchover)||any(strcmp(spin_system.sys.disable,'krylov'))
                
                % Use exponential propagation
                P=propagator(spin_system,L_subs,timestep);
                for n=2:size(rho,2)
                    rho_subs(:,n:end)=P*rho_subs(:,n:end);                  
                end
                answer(:,2:end)=answer(:,2:end)+projectors{subs}*rho_subs(:,2:end);
                
            else
                
                % For very large subspaces use Krylov propagation
                report(spin_system,['evolution: WARNING - Krylov propagation used in subspace #' num2str(subs)]);
                for n=2:size(rho,2)
                    rho_subs(:,n:end)=step(spin_system,L_subs,rho_subs(:,n:end),timestep);                   
                end
                answer(:,2:end)=answer(:,2:end)+projectors{subs}*rho_subs(:,2:end);
                
            end
        
        end
        
    case 'observable'
        
        % Preallocate the answer
        answer=zeros(nsteps+1,size(rho,2));       
        
        % Loop over independent subspaces
        %jlh -reduce output
        %report(spin_system,'evolution: propagating the system...');
        for subs=1:length(projectors)
            
            % Project into the current subspace
            L_subs=projectors{subs}'*L*projectors{subs};
            rho_subs=projectors{subs}'*rho;
            coil_subs=projectors{subs}'*coil;
            
            % The first point does not require propagation
            answer(1,:)=answer(1,:)+coil_subs'*rho_subs;
            
            % Propagate the system
            if (nnz(L_subs)<spin_system.tols.krylov_switchover)||any(strcmp(spin_system.sys.disable,'krylov'))

                % Use exponential propagation
                P=propagator(spin_system,L_subs,timestep);
                for n=1:nsteps
                    rho_subs=P*rho_subs;
                    answer(n+1,:)=answer(n+1,:)+coil_subs'*rho_subs;
                end
                
            else
                
                % For very large subspaces use Krylov propagation
                report(spin_system,['evolution: WARNING - Krylov propagation used in subspace #' num2str(subs)]);
                for n=1:nsteps
                    rho_subs=step(spin_system,L_subs,rho_subs,timestep);
                    answer(n+1,:)=answer(n+1,:)+coil_subs'*rho_subs;
                end
                
            end
            
        end
        
    case 'dec'
        
        % Preallocate the answer
        answer=zeros(nsteps+1,size(rho,2));
        
        report(spin_system,'evolution: using DEC algorithm to calculate the observable...')
        for subs=1:length(projectors)
            
            L_subs=projectors{subs}'*L*projectors{subs};
            rho_subs=projectors{subs}'*rho;
            coil_subs=projectors{subs}'*coil;
            
            % Generate fid using DEC algorithm
            answer=answer+dec(spin_system,L_subs,coil_subs,rho_subs,timestep,nsteps).';
            
        end
        
    case 'multichannel'
        
        % Preallocate the answer
        answer=zeros(size(coil,2),nsteps+1);
        
        % Reduce the problem dimension
        report(spin_system,'evolution: trying to reduce the problem dimension...');
        projectors=reduce(spin_system,L,rho);
        
        % Loop over independent subspaces
        report(spin_system,'evolution: propagating the system...');
        for subs=1:length(projectors)
            
            % Project into the current subspace
            L_subs=projectors{subs}'*L*projectors{subs};
            rho_subs=projectors{subs}'*rho;
            coil_subs=projectors{subs}'*coil;
            
            % The first point does not require propagation
            answer(:,1)=answer(:,1)+coil_subs'*rho_subs;
            
            % Propagate the system
            if (nnz(L_subs)<spin_system.tols.krylov_switchover)||any(strcmp(spin_system.sys.disable,'krylov'))
                
                % Use exponential propagation
                P=propagator(spin_system,L_subs,timestep);
                for n=1:nsteps
                    rho_subs=P*rho_subs;
                    answer(:,n+1)=answer(:,n+1)+coil_subs'*rho_subs;
                end
                
            else
                
                % For very large subspaces use Krylov propagation
                report(spin_system,['evolution: WARNING - Krylov propagation used in subspace #' num2str(subs)]);
                for n=1:nsteps
                    rho_subs=step(spin_system,L_subs,rho_subs,timestep);
                    answer(:,n+1)=answer(:,n+1)+coil_subs'*rho_subs;
                end
                
            end
            
        end
        
    otherwise
        
        error('evolution: invalid output option.');

end

% Patch out thermal relaxation
if strcmp(spin_system.rlx.equilibrium,'thermal')
    switch output
        case {'final','trajectory','refocused_trajectory','refocus'}
            answer(1,:)=0;
        case {'observable','multichannel','dec'}
            % Do nothing
    end
end

end

% Degrees of ability vary, but the basic principle remains the same: the
% degree of a man's independence, initiative and personal love for his work
% determines his talent as a worker and his worth as a man. Independence is
% the only gauge of human virtue and value. What a man is and makes of
% himself; not what he has or hasn't done for others. There is no
% substitute for personal dignity. There is no standard of personal dignity
% except independence.
%
% Ayn Rand, "The Fountainhead"

