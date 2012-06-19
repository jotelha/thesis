% Basis set control. See the documentation for the input parameters.
%
% ilya.kuprov@oerc.ox.ac.uk

function spin_system=basis(spin_system,bas)

%% Preliminaries

% Show the banner
banner(spin_system,'basis_banner');

% Validate the input
validate(spin_system,bas,'basis_input')

% Set the defaults
if ~isfield(bas,'mode'), bas.mode='IK-2'; end
if ~isfield(bas,'level'), bas.level=4; end
if ~isfield(bas,'space_level'), bas.space_level=2; end
if ~isfield(bas,'cnes'), bas.cnes=0; end
if ~isfield(bas,'manual'), bas.manual=[]; end

%% Connectivity analysis

if strcmp(bas.mode,'IK-1')
    
    % Update the user
    report(spin_system,'basis: IK-1 connectivity-adaptive basis set.');
    
    % Run the connectivity analysis
    report(spin_system,['basis: direct product states up to (and including) order ' num2str(bas.level) ' between directly coupled spins.']);
    report(spin_system,['basis: ACS(' num2str(bas.level) ') partitioning of coherent coupling graph...']);
    coherent_subgraphs=dfpt(spin_system.inter.conmatrix,bas.level);
    report(spin_system,['basis: ' num2str(size(coherent_subgraphs,1)) ' subgraphs generated.']);

    % Run the proximity analysis
    report(spin_system,['basis: direct product states up to (and including) order ' num2str(bas.space_level) ' between spins within ' num2str(spin_system.tols.prox_cutoff) ' Angstrom of each other.']);
    report(spin_system,['basis: ACS(' num2str(bas.space_level) ') partitioning of proximity graph...']);
    proximity_subgraphs=dfpt(spin_system.inter.proxmatrix,bas.space_level);
    report(spin_system,['basis: ' num2str(size(proximity_subgraphs,1)) ' subgraphs generated.']);
    
end
        
if strcmp(bas.mode,'IK-2')
    
    % Update the user
    report(spin_system,'basis: IK-2 connectivity-adaptive basis set.');
    
    % Create the cluster list based on the connectivity matrix
    report(spin_system,'basis: direct product states involving all nearest neighbours for each spin.');
    coherent_subgraphs=spin_system.inter.conmatrix|eye(size(spin_system.inter.conmatrix));
    report(spin_system,['basis: ' num2str(size(coherent_subgraphs,1)) ' subgraphs generated.']);
    
    % Run the proximity analysis
    report(spin_system,['basis: direct product states up to (and including) order ' num2str(bas.space_level) ' between spins within ' num2str(spin_system.tols.prox_cutoff) ' Angstrom of each other.']);
    report(spin_system,['basis: ACS(' num2str(bas.space_level) ') partitioning of proximity graph...']);
    proximity_subgraphs=dfpt(spin_system.inter.proxmatrix,bas.space_level);
    report(spin_system,['basis: ' num2str(size(proximity_subgraphs,1)) ' subgraphs generated.']);
    
end
        
if strcmp(bas.mode,'IK-0')
    
    % Update the user
    report(spin_system,['basis: IK-0 basis set - all direct product states up to (and including) order ' num2str(bas.level)]);
    
    % Generate the cluster list using all possible picks
    picks=combnk(1:spin_system.comp.nspins,bas.level); nsubgraphs=size(picks,1);
    coherent_subgraphs=zeros(nsubgraphs,spin_system.comp.nspins);
    for n=1:nsubgraphs
        coherent_subgraphs(n,picks(n,:))=1;
    end
    coherent_subgraphs=logical(coherent_subgraphs);
    
    % Do not run proximity analysis
    proximity_subgraphs=[];
    
end

if any(strcmp(bas.mode,{'IK-1','IK-2','IK-0'}))
    
    % Include user-specified subgraphs
    if isempty(bas.manual)
        report(spin_system,'basis: no subgraphs specified manually by the user.');
        manual_subgraphs=[];
    else
        manual_subgraphs=bas.manual;
        report(spin_system,['basis: added ' num2str(size(manual_subgraphs,1)) ' subgraphs specified by the user.']);
    end
    
    % Run the CNES procedure
    if bas.cnes
        report(spin_system,'basis: CNES partitioning of coherent coupling graph...');
        cnes_subgraphs=spin_system.inter.conmatrix|eye(size(spin_system.inter.conmatrix));
        report(spin_system,['basis: ' num2str(size(cnes_subgraphs,1)) ' subgraphs generated.']);
    else
        report(spin_system,'basis: CNES augmentation disabled.');
        cnes_subgraphs=[];
    end
    
    % Remove identical subgraphs and those completely contained within other subgraphs
    report(spin_system,'basis: removing identical subgraphs...');
    subgraphs=unique([coherent_subgraphs; proximity_subgraphs; cnes_subgraphs; manual_subgraphs],'rows');

end

if any(strcmp(bas.mode,{'IK-1','IK-2'}))||bas.cnes||(~isempty(bas.manual))
    
    % Remove the subraphs completely enclosed by other subgraphs
    report(spin_system,'basis: removing subgraphs completely enclosed by other subgraphs...');
    positions_to_remove=zeros(size(subgraphs,1),1);
    for n=1:size(subgraphs,1)
        for k=1:size(subgraphs,1)
            if (nnz(subgraphs(k,:))<nnz(subgraphs(n,:)))&&(nnz(or(subgraphs(k,:),subgraphs(n,:)))==nnz(subgraphs(n,:)))
                positions_to_remove(k)=1;
            end
        end
    end
    subgraphs(logical(positions_to_remove),:)=[];
    
end
    
if any(strcmp(bas.mode,{'IK-1','IK-2','IK-0'}))
    
    % Write to the data structure and tell the user
    spin_system.bas.subgraphs=subgraphs;
    spin_system.bas.nsubgraphs=size(spin_system.bas.subgraphs,1);
    subgraph_sizes=sum(spin_system.bas.subgraphs,2);
    for n=min(subgraph_sizes):max(subgraph_sizes)
        if nnz(subgraph_sizes==n)>0
            report(spin_system,['basis: generated ' num2str(nnz(subgraph_sizes==n)) ' subgraphs of size ' num2str(n)]);
        end
    end
    report(spin_system,['basis: a total of ' num2str(size(spin_system.bas.subgraphs,1)) ' subgraphs generated.']);
    
end

%% Basis set generation

% For connectivity based basis sets, merge the subgraph basis sets
if any(strcmp(bas.mode,{'IK-1','IK-2','IK-0'}))
    
    % Preallocate the basis descriptor array
    nstates=0;
    for n=1:size(subgraphs,1)
        nstates=nstates+prod(spin_system.comp.mults(logical(subgraphs(n,:))))^2;
    end
    basis_spec=zeros(nstates,spin_system.comp.nspins,'single');
    
    % Populate the basis descriptor array
    list_position=1;
    for n=1:size(subgraphs,1)
        
        % Determine the total number of states in the current subgraph
        nstates=prod(spin_system.comp.mults(logical(subgraphs(n,:))))^2;
        
        % Determine which spins belong to the current cluster
        spins_involved=find(subgraphs(n,:));
        
        % Build the basis table
        for k=1:length(spins_involved)
            basis_spec(list_position:(list_position+nstates-1),spins_involved(k))=...
            kron(kron(ones(prod(spin_system.comp.mults(spins_involved(1:(k-1))))^2,1),...
                               (0:(spin_system.comp.mults(spins_involved(k))^2-1))'),...
                      ones(prod(spin_system.comp.mults(spins_involved((k+1):end)))^2,1));
        end
        list_position=list_position+nstates;
        
    end
    
end

% For the exact basis sets use explicit generation
if strcmp(bas.mode,'complete')
    
    % Determine the total number of states and preallocate the basis table
    basis_spec=zeros(prod(spin_system.comp.mults)^2,spin_system.comp.nspins);
    
    % Build the basis table
    for k=1:spin_system.comp.nspins
        basis_spec(:,k)=kron(kron(ones(prod(spin_system.comp.mults(1:(k-1)))^2,1),(0:(spin_system.comp.mults(k)^2-1))'),ones(prod(spin_system.comp.mults((k+1):end))^2,1)); 
    end
    
end

if strcmp(bas.mode,'ESR-1')
    
    report(spin_system,'basis: ESR-1 basis set - complete on electrons, {E,Lz} on nuclei.');
    
    % Prepare the direct product components
    components=cell(1,spin_system.comp.nspins);
    for n=1:spin_system.comp.nspins
        if strcmp(spin_system.comp.isotopes{n}(1),'E')
            components{n}=(0:(spin_system.comp.mults(n)^2-1))';
        else
            components{n}=((1:spin_system.comp.mults(n)).*((1:spin_system.comp.mults(n))-1))';
        end
    end

    % Preallocate the basis table
    basis_spec=zeros(prod(cellfun(@numel,components)),spin_system.comp.nspins);
    
    % Include {E,Lz} for nuclei and a complete basis set for electrons
    for n=1:spin_system.comp.nspins
        dimension_before=prod(cellfun(@numel,components(1:(n-1))));
        dimension_after=prod(cellfun(@numel,components((n+1):end)));
        basis_spec(:,n)=kron(kron(ones(dimension_before,1),components{n}),ones(dimension_after,1));
    end

end

if strcmp(bas.mode,'ESR-2')
    
    % Analyze the coupling structure
    spins_with_full_basis=zeros(1,spin_system.comp.nspins);
    for n=1:spin_system.comp.nspins
        
        % User complete basis set on anything that has an anisotropic coupling
        for k=1:spin_system.comp.nspins
            if (norm(spin_system.inter.coupling.matrix{n,k})>spin_system.tols.inter_cutoff)&&...
                norm(spin_system.inter.coupling.matrix{n,k}-eye(3)*trace(spin_system.inter.coupling.matrix{n,k})/3)>spin_system.tols.inter_cutoff
                spins_with_full_basis(n)=1; spins_with_full_basis(k)=1;
            end
        end
        
        % User complete basis set on all electrons
        if strcmp(spin_system.comp.isotopes{n}(1),'E')
            spins_with_full_basis(n)=1;
        end
        
        % Inform the user
        if spins_with_full_basis(n)
            report(spin_system,['basis: [ESR-2] complete basis on spin ' num2str(n)]);
        else
            report(spin_system,['basis: [ESR-2] zero-quantum basis on spin ' num2str(n)]);
        end
        
    end
    
    % Prepare the direct product components
    components=cell(1,spin_system.comp.nspins);
    for n=1:spin_system.comp.nspins
        if spins_with_full_basis(n)
            components{n}=(0:(spin_system.comp.mults(n)^2-1))';
        else
            components{n}=((1:spin_system.comp.mults(n)).*((1:spin_system.comp.mults(n))-1))';
        end
    end

    % Preallocate the basis table
    basis_spec=zeros(prod(cellfun(@numel,components)),spin_system.comp.nspins);
    
    % Include {E,Lz} for nuclei and a complete basis set for electrons
    for n=1:spin_system.comp.nspins
        dimension_before=prod(cellfun(@numel,components(1:(n-1))));
        dimension_after=prod(cellfun(@numel,components((n+1):end)));
        basis_spec(:,n)=kron(kron(ones(dimension_before,1),components{n}),ones(dimension_after,1));
    end

end

% Apply the coherence level filter
if isfield(bas,'projections')
    
    % Compute the projection quantum numbers for each basis element
    [~,M]=lin2lm(basis_spec);
    
    % Keep the projection quantum numbers indicated
    state_mask=zeros(size(basis_spec,1),1);
    projection_numbers=sum(M,2);
    for n=bas.projections
        state_mask=state_mask|(projection_numbers==n);
    end
    basis_spec(~state_mask,:)=[];
    
    % Inform the user
    report(spin_system,['basis: projection filter applied for M=[' num2str(bas.projections) ']']);
    
end

% Write the data structure
spin_system.bas.basis=unique(basis_spec,'rows');
spin_system.bas.nstates=size(spin_system.bas.basis,1);

% Print the summary
summary(spin_system,'basis');

%% Symmetry treatment

if any(strcmp(spin_system.sys.disable,'symmetry'))
    
    % Issue a reminder to the user
    report(spin_system,'basis: WARNING - symmetry factorization disabled by the user.');
    
else
    
    % Get the group data
    if length(spin_system.comp.sym.group)>1
        
        % Lift the constituent groups from the database
        ngroups=length(spin_system.comp.sym.group);
        groups=cell(1,ngroups);
        for n=1:ngroups
            groups{n}=symmetry_group(spin_system.comp.sym.group{n});
        end
        
        % Compute the direct product character table
        group.characters=1;
        for n=1:ngroups
            group.characters=kron(group.characters,groups{n}.characters);
        end
        group.irrep_dims=group.characters(:,1)';
        group.n_irreps=size(group.characters,1);
        report(spin_system,['basis: ' num2str(group.n_irreps) ' irreps in the group direct product.']);
        report(spin_system,['basis: dimensions of the irreps ' num2str(group.irrep_dims)]);
        
        % Compute the direct product element list
        group.elements=groups{1}.elements; group.order=groups{1}.order;
        for n=2:ngroups
            group.elements=[kron(group.elements,ones(groups{n}.order,1)) kron(ones(group.order,1),groups{n}.elements+size(group.elements,2))];
            group.order=group.order*groups{n}.order;
        end
        report(spin_system,['basis: ' num2str(size(group.elements,1)) ' symmetry operations in the direct product.']);
        
        % Concatenate the spin lists
        spins=horzcat(spin_system.comp.sym.spins{:});
        
    elseif length(spin_system.comp.sym.group)==1
        
        % Lift the group from the database
        spins=spin_system.comp.sym.spins{1};
        group=symmetry_group(spin_system.comp.sym.group{1});
        report(spin_system,['basis: ' num2str(group.n_irreps) ' irreps in the symmetry group.']);
        report(spin_system,['basis: dimensions of the irreps ' num2str(group.irrep_dims)]);
        
    else
        
        % Remind the user that symmetry is not operational
        report(spin_system,'basis: no symmetry information available.');
        
    end
    
    % Make sure the basis is sorted
    if ~issorted(spin_system.bas.basis,'rows')
        error('basis: the basis is not correctly sorted.');
    end
    
    % Run the SALC procedure
    if exist('group','var')
        
        % Compute the permutation table
        permutation_table=zeros(spin_system.bas.nstates,group.order);
        for n=1:group.order
            group_element=1:spin_system.comp.nspins;
            group_element(spins)=group_element(spins(group.elements(n,:)));
            permuted_basis=spin_system.bas.basis(:,group_element);
            [~,index]=sortrows(permuted_basis);
            permutation_table(:,n)=index;
        end
        
        % Compute the irreducible representations
        if spin_system.comp.sym.a1g_only
            
            % Inform the user
            report(spin_system,'basis: Liouville space symmetry mode - fully symmetric irrep only.');
            
            % Prune the permutation table
            symmetry_related_states=unique(sort(permutation_table,2,'ascend'),'rows');
            dimension=size(symmetry_related_states,1);
            
            % Populate the coefficient matrix
            index=unique([kron(ones(group.order,1),(1:dimension)') symmetry_related_states(:) ones(dimension*group.order,1)],'rows');
            coeff_matrix=sparse(index(:,1),index(:,2),index(:,3),dimension,spin_system.bas.nstates);
            
            % Normalize the coefficient matrix
            norms=sqrt(sum(coeff_matrix.^2,2));
            coeff_matrix=spdiags(norms.^(-1),0,dimension,dimension)*coeff_matrix;
            
            % Report back to the user
            report(spin_system,['basis: A1g irrep, ' num2str(dimension) ' states.']);
            
            % Return the projector and dimension
            spin_system.bas.irrep.projector=coeff_matrix';
            spin_system.bas.irrep.dimension=dimension;
            
        else
            
            % Inform the user
            report(spin_system,'basis: WARNING - full symmetry treatment, all irreps will be included.');
            
            for n=1:group.n_irreps
                
                % Compute the coefficient matrix
                coeff_matrix=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,spin_system.bas.nstates*group.order);
                for k=1:spin_system.bas.nstates
                    
                    % Populate the coefficient matrix
                    for m=1:group.order
                        coeff_matrix(k,permutation_table(k,m))=coeff_matrix(k,permutation_table(k,m))+group.characters(n,m); %#ok<SPRIX>
                    end
                    
                    % Unify the signs
                    non_zero_elements=nonzeros(coeff_matrix(k,:));
                    if any(non_zero_elements)
                        coeff_matrix(k,:)=coeff_matrix(k,:)*sign(non_zero_elements(1)); %#ok<SPRIX>
                    end
                    
                end
                
                % Remove identical rows and all-zero rows
                coeff_matrix=unique(coeff_matrix,'rows');
                coeff_matrix(sum(abs(coeff_matrix),2)<spin_system.tols.liouv_zero,:)=[];
                
                % If the irrep is 2+ dimensional, orthonormalize the SALCs. Otherwise,
                % just normalize the SALCs.
                if group.irrep_dims(n)>1
                    report(spin_system,['basis: WARNING - ' num2str(group.irrep_dims(n)) ' dimensional irrep encountered. Will have to orthogonalize (slow).']);
                    coeff_matrix=orth(full(coeff_matrix'))';
                else
                    for k=1:size(coeff_matrix,1)
                        coeff_matrix(k,:)=coeff_matrix(k,:)/norm(coeff_matrix(k,:));
                    end
                end
                
                % Inform the user and write the irrep into the data structure
                report(spin_system,['basis: Irrep #' num2str(n) ', ' num2str(size(coeff_matrix,1)) ' states.']);
                spin_system.bas.irrep(n).projector=coeff_matrix';
                spin_system.bas.irrep(n).dimension=size(coeff_matrix,1);
                
            end
            
        end
        
    end
       
end

end

% In 1969, Robert Rathbun Wilson, the US physicist who headed Fermilab, the
% world's highest-energy particle accelerator laboratory, addressed the
% Congressional Joint Committee on Atomic Energy. Rhode Island Senator John
% Pastore asked Wilson to spell out what research into high-energy particle
% physics would do to improve the defence of the United States. Wilson gave
% a reply that went down in scientific history. Fermilab, he said, had
% "nothing to do directly with defending our country, except to make it
% worth defending".
%
% http://www.theregister.co.uk/2009/02/09/woudhuysen_energise_1/

