% Chemical kinetics superoperator. All adjustable parameters are specified
% in the call to create.m function.
%
% ilya.kuprov@oerc.ox.ac.uk

function K=k_superop(spin_system)

% Preallocate the kinetics superoperator
K=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,2*spin_system.bas.nstates);

% Process chemical kinetics
if ~isempty(spin_system.chem.rate)
    
    % Inform the user
    if strcmp(spin_system.chem.type,'intramolecular')
        report(spin_system,'k_superop: intermolecular exchange, multi-spin orders involving exchanging spins will be kept.');
    else
        report(spin_system,'k_superop: intramolecular exchange, multi-spin orders involving exchanging spins will be damped.');
    end
    
    % Separate the basis into single- and multi-spin orders
    if any(spin_system.chem.rate)
        sso_index=find(sum(logical(spin_system.bas.basis),2)==1);
        sso_basis=spin_system.bas.basis(sso_index,:);
        mso_index=find(sum(logical(spin_system.bas.basis),2)>1);
        mso_basis=spin_system.bas.basis(mso_index,:);
    end
    
    % Loop over the rates
    for n=1:length(spin_system.chem.rate)
        
        % Classify the states into those to be added to and those to be subtracted from
        source_index=find(sso_basis(:,spin_system.chem.from(n))~=0);
        destin_index=find(sso_basis(:,spin_system.chem.to(n))~=0);
        deplet_index=find(mso_basis(:,spin_system.chem.from(n))~=0);
        if length(source_index)~=length(destin_index)
            error('k_superop: restricted state space breach in the kinetics superoperator.');
        end
        
        % Deplete the source single-spin orders
        K=K-sparse(sso_index(source_index),sso_index(source_index),spin_system.chem.rate(n)*ones(size(source_index)),spin_system.bas.nstates,spin_system.bas.nstates);
        
        % Populate the destination single-spin orders
        K=K+sparse(sso_index(destin_index),sso_index(source_index),spin_system.chem.rate(n)*ones(size(source_index)),spin_system.bas.nstates,spin_system.bas.nstates);
        
        % Deplete the source multi-spin orders
        K=K-sparse(mso_index(deplet_index),mso_index(deplet_index),spin_system.chem.rate(n)*ones(size(deplet_index)),spin_system.bas.nstates,spin_system.bas.nstates);
        
        % Populate the destination multi-spin orders
        if strcmp(spin_system.chem.type,'intramolecular')
            K=K+sparse(mso_index(destin_index),mso_index(source_index),spin_system.chem.rate(n)*ones(size(source_index)),spin_system.bas.nstates,spin_system.bas.nstates);
        end
        
    end
    
end

% Process radical pair recombination
if any(strcmp(spin_system.chem.rp_theory,{'haberkorn','hore-jones'}))
    
    % Get the singlet and triplet projector product superoperators
    unit=speye(spin_system.bas.nstates);
    LzSz_l=operator(spin_system,{'Lz','Lz'},num2cell(spin_system.chem.rp_spins),'left');
    LzSz_r=operator(spin_system,{'Lz','Lz'},num2cell(spin_system.chem.rp_spins),'right');
    LpSm_l=operator(spin_system,{'L+','L-'},num2cell(spin_system.chem.rp_spins),'left');
    LpSm_r=operator(spin_system,{'L+','L-'},num2cell(spin_system.chem.rp_spins),'right');
    LmSp_l=operator(spin_system,{'L-','L+'},num2cell(spin_system.chem.rp_spins),'left');
    LmSp_r=operator(spin_system,{'L-','L+'},num2cell(spin_system.chem.rp_spins),'right');
    singlet_l=unit/4-(LzSz_l+0.5*(LpSm_l+LmSp_l)); triplet_l=unit-singlet_l;
    singlet_r=unit/4-(LzSz_r+0.5*(LpSm_r+LmSp_r)); triplet_r=unit-singlet_r;
    
    % Assemble the recombination kinetics superoperator
    switch spin_system.chem.rp_theory
        
        case 'haberkorn'
            
            K=K-0.5*(spin_system.chem.rp_rates(1)*(singlet_l+singlet_r)+spin_system.chem.rp_rates(2)*(triplet_l+triplet_r));
            
        case 'hore-jones'
            
            K=K-sum(spin_system.chem.rp_rates)*unit+spin_system.chem.rp_rates(1)*triplet_l*triplet_r+spin_system.chem.rp_rates(2)*singlet_l*singlet_r;
            
    end
    
end

end

% The egoist is the absolute sense is not the man who sacrifices others.
% He is the man who stands above the need of using others in any manner. He
% does not function through them. He is not concerned with them in any
% primary matter. Not in his aim, not in his motive, not in his thinking,
% not in his desires, not in the source of his energy. He does not exist
% for any other man - and he asks no other man to exist for him. This is the
% only form of brotherhood and mutual respect possible between men.
%
% Ayn Rand, "The Fountainhead"


