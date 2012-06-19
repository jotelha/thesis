% Sided product superoperator. Returns superoperators corresponding to
% right or left multiplication of a density matrix by a user-specified
% operator. Arguments:
%
%     opspec - Spinach operator specification (see the manual)
%     side   - 'left' or 'right'
%
% If multiple lines are present in opspec, returns the sum of the cor-
% responding superoperators. 
%
% ilya.kuprov@oerc.ox.ac.uk
% hannah.hogben@chem.ox.ac.uk

function A=p_superop(spin_system,opspec,side)

% If operator caching is enabled, consult the cache first
if isfield(spin_system.sys,'identifier')
    
    % Generate a unique operator identifier
    identifier=[spin_system.sys.identifier sprintf('_%i',opspec(:)) side];
    filename=[spin_system.sys.root_dir spin_system.sys.slash 'scratch' spin_system.sys.slash identifier '.mat'];
    
    % If the file is present, load the operator from disk
    if exist(filename,'file')
        
        % Load the operator
        load(filename);
        
        if exist('A','var')
            
            % Run a cursory check on dimensions
            if ~all(size(A)==[spin_system.bas.nstates spin_system.bas.nstates]) %#ok<NODEF>
                error('p_superop: dimension mismatch between the current basis and the cached superoperators.');
            end
            
            % Return the operator
            return
            
        end
        
    end
    
end

% Make sure the basis is sorted
if ~issorted(spin_system.bas.basis,'rows')
    error('p_superop: the basis set is not sorted, please sort your basis.');
end

% Preallocate the answer
A=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,size(opspec,1)*spin_system.bas.nstates);

% Loop over the lines in opspec
for k=1:size(opspec,1)
    
    % Determine the active spins
    active_spins=find(opspec(k,:));
    
    % Decide how to proceed
    if isempty(active_spins)
        
        % For unit state, just return the unit matrix
        L=speye(spin_system.bas.nstates);
        
    else
        
        % Lift the structure coefficients for the individual spins
        source=cell(1,length(active_spins));
        destin=cell(1,length(active_spins));
        struct=cell(1,length(active_spins));
        for n=1:length(active_spins)
            [product_table_left,product_table_right]=ist_product_table(spin_system,spin_system.comp.mults(active_spins(n)));
            switch side
                case 'left'
                    product_table=sparse(squeeze(product_table_left(opspec(k,active_spins(n))+1,:,:)));
                case 'right'
                    product_table=sparse(squeeze(product_table_right(opspec(k,active_spins(n))+1,:,:)));
                otherwise
                    error('p_superop: invalid side specification.');
            end
            product_table=product_table.*(abs(product_table)>spin_system.tols.liouv_zero);
            [destin{n},source{n},struct{n}]=find(product_table); source{n}=source{n}-1; destin{n}=destin{n}-1;
        end
        
        % Compute the structure coefficients for the active sub-algebra
        from=source{1}; to=destin{1}; coeff=struct{1};
        for n=2:length(active_spins)
            from=[kron(from,ones(size(source{n},1),1)) kron(ones(size(from,1),1),source{n})];
            to=[kron(to,ones(size(destin{n},1),1)) kron(ones(size(to,1),1),destin{n})];
            coeff=kron(coeff,struct{n});
        end
        
        % Preallocate the superoperator
        L=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,spin_system.bas.nstates);
        
        % Lift the basis columns corresponding to the active spins
        basis_cols=spin_system.bas.basis(:,active_spins);
        
        % Compute the superoperator
        for n=1:size(from,1)
            
            % Retrieve the source subspace
            source_subspace_index=true(size(basis_cols,1),1);
            for m=1:size(from,2)
                source_subspace_index=and(source_subspace_index,(basis_cols(:,m)==from(n,m)));
            end
            source_subspace=spin_system.bas.basis(source_subspace_index,:);
            source_subspace_index=find(source_subspace_index);
            source_subspace(:,active_spins)=[];
            
            % Retrieve the destination subspace
            destin_subspace_index=true(size(basis_cols,1),1);
            for m=1:size(to,2)
                destin_subspace_index=and(destin_subspace_index,(basis_cols(:,m)==to(n,m)));
            end
            destin_subspace=spin_system.bas.basis(destin_subspace_index,:);
            destin_subspace_index=find(destin_subspace_index);
            destin_subspace(:,active_spins)=[];
            
            % Decide on the filling procedure
            if isequal(source_subspace,destin_subspace)
                % If the subspaces fully match, use the raw indices
                subspace_dim=size(source_subspace,1);
                L=L+sparse(source_subspace_index,destin_subspace_index,coeff(n)*ones(subspace_dim,1),spin_system.bas.nstates,spin_system.bas.nstates);
            else
                % Otherwise, use brute-force state-by-state matching
                [does_it_go_anywhere,where_it_goes_if_it_does]=ismember(source_subspace,destin_subspace,'rows');
                L=L+sparse(source_subspace_index(does_it_go_anywhere),destin_subspace_index(where_it_goes_if_it_does(does_it_go_anywhere)),...
                           coeff(n)*ones(nnz(does_it_go_anywhere),1),spin_system.bas.nstates,spin_system.bas.nstates);
            end
            
        end
        
    end
    
    % Add up the superoperators corresponding to the individual opspec lines
    A=A+L;
    
end

% If operator caching is enabled, update the cache
if isfield(spin_system.sys,'identifier')
    save(filename,'A');
end

end

% My philosophy, in essence, is the concept of man as a heroic being, with
% his own happiness as the moral purpose of his life, with productive
% achievement as his noblest activity, and reason as his only absolute.
%
% Ayn Rand

