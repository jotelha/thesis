% Graph partitioning module. Analyzes the system connectivity graph and
% creates a list of all connected subgraphs of up to the user-specified
% size by crawling the graph in all available directions. Arguments:
%
%       connectivity_matrix - a matrix with 1 for connected spins and zero
%                             otherwise.
%
%       max_subgraph_size   - maximum connected subgraph size in the
%                             resulting list.
%
% ilya.kuprov@oerc.ox.ac.uk

function subgraphs=dfpt(connectivity_matrix,max_subgraph_size)

% Start at each spin in the network
subgraphs=eye(size(connectivity_matrix));

% Depth-first path tracing through the coupling graph
for n=1:(max_subgraph_size-1)
    
    % Determine immediate neigbours reachable from current subgraphs 
    neighbor_matrix=logical(subgraphs*connectivity_matrix);
    neighbor_matrix(logical(subgraphs))=0;
    
    % Determine the dimensions of the new subgraph descriptor array
    ngraphs=nnz(neighbor_matrix)+nnz(sum(neighbor_matrix,2)==0);
    nspins=size(connectivity_matrix,2);
    
    % Preallocate the new subgraph descriptor array
    new_subgraphs=zeros(ngraphs,nspins);
    
    % Grow each subgraph by one node in each available direction
    list_position=1;
    for k=1:size(subgraphs,1)
        spins_to_add=find(neighbor_matrix(k,:));
        if isempty(spins_to_add)
            new_subgraphs(list_position,:)=subgraphs(k,:);
            list_position=list_position+1;
        else
            for m=spins_to_add
                new_subgraphs(list_position,:)=subgraphs(k,:);
                new_subgraphs(list_position,m)=1;
                list_position=list_position+1;
            end
        end
    end
    
    % Remove repetitions
    new_subgraphs=logical(new_subgraphs);
    subgraphs=unique(new_subgraphs,'rows');
    
end

end

% Psychologists call this “cognitive dissonance”—the ability to make a
% compelling, heartfelt case for one thing while doing another. Being able
% to pull off this sort of trick is an essential skill in many professions:
% "even if his message bore no relation to his actions, it expressed
% precisely and succinctly what he should have been doing".
%
% http://www.economist.com/opinion/displayStory.cfm?Story_ID=10311230

