% Symmetry and trajectory-level state space restriction for a
% user-specified Liouvillian L and initial state rho. Applies all
% reduction methods (unless disabled during the call to create.m) and
% returns a cell array of projectors into a set of independent reduced
% subspaces. The projectors are to be used as follows:
%
%     L=P'*L*P;     rho=P'*rho;
%
% ilya.kuprov@oerc.ox.ac.uk

function projectors=reduce(spin_system,L,rho)

%% Preliminaries

% Click forward for output
spin_system=click(spin_system,'forward');

% If a stack is supplied, choose a representative state vector
if size(rho,2)>1
    rho=mean(abs(rho),2);
end

%% Symmetry factorization
if any(strcmp(spin_system.sys.disable,'symmetry'))
    
    % Issue a reminder to the user
    report(spin_system,'reduce: WARNING - symmetry factorization disabled by the user.');
    
    % Pass the Liouvillian and the density matrix on unchanged
    L_symm={L}; rho_symm={rho};
    
    % Projector is a unit matrix
    projectors_symm={speye(size(L))};
    
elseif ~isfield(spin_system.bas,'irrep')
    
    % Inform the user
    %jlh - reduce output
    % report(spin_system,'reduce: no symmetry information available.');
    
    % Pass the Liouvillian and the density matrix on unchanged
    L_symm={L}; rho_symm={rho};
    
    % Projector is a unit matrix
    projectors_symm={speye(size(L))};
    
else
    
    % Start with empty cell arrays (dimensions are hard to predict)
    L_symm={}; rho_symm={}; projectors_symm={};
    
    % Loop over irreducible representations
    for n=1:length(spin_system.bas.irrep)
        
        if spin_system.bas.irrep(n).dimension==0
            
            % Ignore empty irreducible representations
            report(spin_system,['reduce: irreducible representation #' num2str(n) ' has no states in it - dropped.']);
            
        elseif negligible(spin_system.bas.irrep(n).projector'*rho,'vector',spin_system.tols.irrep_drop)
            
            % Ignore unpopulated irreducible representations
            report(spin_system,['reduce: irreducible representation #' num2str(n) ', dimension ' num2str(spin_system.bas.irrep(n).dimension) ' is unpopulated - dropped.']);
            
        else
            
            % Keep populated irreducible representations
            report(spin_system,['reduce: irreducible representation #' num2str(n) ', dimension ' num2str(spin_system.bas.irrep(n).dimension) ' is active - kept.']);
            
            % Create the irrep projector array
            projectors_symm=[projectors_symm {spin_system.bas.irrep(n).projector}]; %#ok<AGROW>
            
            % Create the irrep Liouvillian array
            L_symm=[L_symm {spin_system.bas.irrep(n).projector'*L*spin_system.bas.irrep(n).projector}]; %#ok<AGROW>
            
            % Create the irrep state vector array
            rho_symm=[rho_symm {spin_system.bas.irrep(n).projector'*rho}]; %#ok<AGROW>
            
        end
        
    end
 
end

%% Zero Track Elimination

% Preallocate the cell arays
L_zte=cell(size(L_symm)); rho_zte=cell(size(rho_symm)); projectors_zte=cell(size(projectors_symm));

% Loop over the independent subspaces received from above
for n=1:length(L_symm)
    
    % Report to the user
    %jlh - reduce output
    %report(spin_system,['reduce: subspace #' num2str(n) ', attempting zero track elimination...']);
    
    % Run Zero Track Elimination
    zte_projector=zte(spin_system,L_symm{n},rho_symm{n});
    
    % Update the subspace Liouvillian
    L_zte{n}=zte_projector'*L_symm{n}*zte_projector;
    
    % Update the subspace state vector
    rho_zte{n}=zte_projector'*rho_symm{n};
    
    % Update the subspace projector
    projectors_zte{n}=projectors_symm{n}*zte_projector;
    
end

%% Path Tracing
if any(strcmp(spin_system.sys.disable,'pt'))
    
    % Issue a reminder to the user
    report(spin_system,'reduce: WARNING - path tracing disabled by the user.');
    
    % Pass the projectors on unchanged
    projectors=projectors_zte;
    
else
    
    % Start with an empty cell array (dimensions are hard to predict)
    projectors={};
    
    % Loop over independent subspaces
    for n=1:length(projectors_zte)
        
        % Inform the user
        %jlh - reduce output
        %report(spin_system,['reduce: path-tracing subspace #' num2str(n) '...']);
        
        % Run the path tracing
        subspace_projectors=path_trace(spin_system,L_zte{n});
        
        % Set the counters for dropped subspaces
        n_dropped_subspaces=0; total_dropped_dimensions=0;
        
        % Loop over the subspaces
        for k=1:length(subspace_projectors)
            
            % Inspect the current subspace for triviality
            if negligible(subspace_projectors{k}'*rho_zte{n},'vector',spin_system.tols.path_drop)
            
                % If the subspace is negligible, ignore it and update the counters
                total_dropped_dimensions=total_dropped_dimensions+size(subspace_projectors{k},2);
                n_dropped_subspaces=n_dropped_subspaces+1;
            
            else
                
                % If the subspace is useful, update the projector array
                projectors=[projectors {projectors_zte{n}*subspace_projectors{k}}]; %#ok<AGROW>
            
            end
            
        end
        
        % Update the user on the dropped subspaces, if any
        if n_dropped_subspaces>0
            report(spin_system,['reduce: dropped ' num2str(n_dropped_subspaces) ' empty subspaces with a total dimension of ' num2str(total_dropped_dimensions)]);
        end
        
    end

end

end

% For centuries, the battle of morality was fought between those who
% claimed that your life belongs to God and those who claimed that it
% belongs to your neighbours � between those who preached that the good is
% self-sacrifice for the sake of ghosts in heaven and those who preached
% that the good is self-sacrifice for the sake of incompetents on earth.
% And no one came to say that your life belongs to you and the good is to
% live it.
%
% Ayn Rand, "Atlas Shrugged"

