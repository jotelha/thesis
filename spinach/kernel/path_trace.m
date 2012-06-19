% Liouvillian path tracing. Treats the user-supplied Liouvillian as the
% adjacency matrix of a graph, computes the weakly connected subgraphs of
% that graph and returns a cell array of projectors into the corresponding
% independently evolving subspaces. The projectors are to be used as
% follows:
%
%     L=P'*L*P;     rho=P'*rho;
%
% ilya.kuprov@oerc.ox.ac.uk

function projectors=path_trace(spin_system,L)

% Click forward for output
spin_system=click(spin_system,'forward');

% Get the connectivity matrix
G=sparse((abs(L)+abs(L'))/2>spin_system.tols.path_trace);

% Make sure isolated states do not get lost
G=or(G,speye(size(G)));

% Get the weakly connected subgraphs
member_nodes=scomponents(G);

% Determine the number of subspaces
n_subs=max(member_nodes);

% Return the subspace projectors
projectors=cell(1,n_subs);
subspace_dims=zeros(1,n_subs);
for n=1:n_subs
    node_index=find(member_nodes==n);
    subspace_dims(n)=length(node_index);
    projectors{n}=sparse(node_index,1:length(node_index),ones(1,length(node_index)),size(L,1),length(node_index));
end

% Display the subspace statistics
for n=unique(subspace_dims);
    report(spin_system,['path_trace: found ' num2str(nnz(subspace_dims==n)) ' independent subspace(s) of dimension ' num2str(n)]);
end

end

% "My dear fellow, who will let you?"
% "That's not the point. The point is, who will stop me?"
%
% Ayn Rand, "The Fountainhead"

