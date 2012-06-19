% Trajectory metric function. Returns a vector representing step-by-step
% "similarity" of the two state space trajectories. A trajectory is to be
% supplied as nstates x nsteps matrix. Methods:
%
% 'RSP'    - running scalar product. Computes scalar products between the
%            corresponding vectors of the two trajectories.
%
% 'NRSP'   - normalized running scalar product. Computes scalar products
%            between corresponding vectors of the two trajectories.
%            Each vector in both trajectories is pre-normalized.
%
% 'DVN'    - difference vector norms. The two trajectories are subtracted
%            and deviation norms returned.
%
% 'NDVN'   - normalized difference vector norms. The two trajectories are
%            subtracted and deviation norms returned. Each vector in both
%            trajectories is pre-normalized.
%
% 'SG-'    - prefix that switches on the state grouping. L(+) and L(-)
%            spin orders of every spin (standalone or in direct products
%            with other operators) would be deemed equivalent. For example,
%            'SG-RSP' requests state-grouped running scalar product.
%
% 'BSG-'   - prefix that switches on broad state grouping. All spin orders
%            of a given spin (standalone or in direct products with other
%            operators) would be deemed equivalent. For example, 'BSG-RSP'
%            requests broadly state-grouped running scalar product.
%
% The state grouping consists in summing the absolute squares of the 
% coefficients to be grouped and taking the square root. Because the
% grouped vectors can have unpredictable norms, normalized metrics, e.g.
% 'SG-NRSP' are preferable.
%
% The four-parameter call 
%
%        metric(trajectory_1,trajectory_2,method,statelist)
%
% assumes that the state spaces of the two trajectories are the same. The
% five-parameter call
%
%        metric(trajectory_1,trajectory_2,method,statelist_1,statelist_2)
%
% will compute the similarity metric within the common subspace of the two
% trajectories.
%
% Ilya Kuprov              ilya.kuprov@oerc.ox.ac.uk
% Konstantin Pervushin     p@trosy.com

function similarity=metric(trajectory_1,trajectory_2,method,statelist_1,statelist_2)

%% Common subspace analysis
if nargin==5
    
    % Make sure the state lists are valid
    if (numel(unique(statelist_1,'rows'))~=numel(statelist_1))||(numel(unique(statelist_2,'rows'))~=numel(statelist_2))
        error('metric: invalid state lists - repetitions found.');
    end
    if (size(trajectory_1,1)~=size(statelist_1,1))||(size(trajectory_2,1)~=size(statelist_2,1))
        error('metric: dimension mismatch between a trajectory and its statelist.');
    end
    if size(statelist_1,2)~=size(statelist_2,2)
        error('metric: different number of spins in the two state spaces.');
    end
    
    % Find the common state space
    disp(['metric: trajectory 1 state space dimension ' num2str(size(statelist_1,1))]);
    disp(['metric: trajectory 2 state space dimension ' num2str(size(statelist_2,1))]);
    disp('metric: mapping the common subspace...');
    [common_state_list,extraction_mask_1,extraction_mask_2]=intersect(statelist_1,statelist_2,'rows');
    disp(['metric: common subspace dimension ' num2str(size(common_state_list,1))]);
    
    % Project the trajectories into the common state space
    trajectory_1=trajectory_1(extraction_mask_1,:);
    trajectory_2=trajectory_2(extraction_mask_2,:);
    
elseif nargin==4

    % Assume that the statelists match if the user only supplied one.
    common_state_list=statelist_1;

end

%% State grouping
if strcmp(method(1:3),'SG-')||strcmp(method(1:3),'BSG')
    
    % Tell the user we're started
    disp('metric: collapsing equivalent subspaces...');
    
    if strcmp(method(1:3),'SG-')
        % Rename all L(+) states into L(-) states
        common_state_list(common_state_list==3)=1;
        % Update the method variable
        method=method(4:end);
    elseif strcmp(method(1:3),'BSG')
        % Rename all non-identity states into Lz
        common_state_list(common_state_list~=0)=2;
        % Update the method variable
        method=method(5:end);
    end
    
    % Index all unique and repeated states on the list
    [grouped_state_list,index_forward,index_backward]=unique(common_state_list,'rows');
    
    % Preallocate state-grouped trajectories
    grouped_trajectory_1=zeros(length(index_forward),size(trajectory_1,2));
    grouped_trajectory_2=zeros(length(index_forward),size(trajectory_2,2));
    
    % Group the trajectory lines corresponding to the states that are
    % flaged as identical in the indices (root-sum-square)
    for n=1:length(index_forward)
        grouped_trajectory_1(n,:)=sqrt(sum(abs(trajectory_1(index_backward==n,:)).^2,1));
        grouped_trajectory_2(n,:)=sqrt(sum(abs(trajectory_2(index_backward==n,:)).^2,1));
    end
    
    % Update the variables
    trajectory_1=grouped_trajectory_1;
    trajectory_2=grouped_trajectory_2;
    
    % Tell the user we're done
    disp(['metric: ' num2str(size(common_state_list,1)) ' states collected into ' num2str(size(grouped_state_list,1)) ' groups.']);
    
end

%% Metrics calculation

% Make sure the trajectory matrices have the same dimensions
if any(size(trajectory_1)~=size(trajectory_2))
    error('metric: trajectory matrix dimensions do not match.');
else
    trajectory_length=size(trajectory_1,2);
end

% Normalize the trajectories if necessary
if strcmp(method,'NRSP')||strcmp(method,'NDVN')
    for n=1:trajectory_length
        trajectory_1(:,n)=trajectory_1(:,n)/norm(trajectory_1(:,n));
    end
    for n=1:trajectory_length
        trajectory_2(:,n)=trajectory_2(:,n)/norm(trajectory_2(:,n));
    end
end

% Preallocate the result array
similarity=zeros(1,trajectory_length);

% Run the RSP procedure
if strcmp(method,'RSP')||strcmp(method,'NRSP')
    for n=1:trajectory_length
        similarity(n)=trajectory_1(:,n)'*trajectory_2(:,n);
    end
end

% Run the DVN procedure
if strcmp(method,'DVN')||strcmp(method,'NDVN')
    similarity=sqrt(sum(abs(trajectory_1-trajectory_2).^2,1));
end

end

% "How would we look for a new law? In general we look for a new law by the
% following process. First, we guess it. Then we... don't laugh. That's the
% damned truth. Then we compute the consequences of the guess... to see if
% this is right, to see if this law we guessed is right, to see what it
% would imply. And then we compare those computation results to nature. Or
% we say to compare it to experiment, or to experience. Compare it directly
% with observations to see if it works. If it disagrees with experiment,
% it's wrong. In that simple statement is the key to science. It doesn't
% make a difference how beautiful your guess is. It doesn't make a
% difference how smart you are, who made the guess or what his name is...
% If it disagrees with experiment, it's wrong. That's all there is to it."
%
% Richard Feynman

