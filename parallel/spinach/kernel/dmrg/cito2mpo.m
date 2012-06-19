% Computes an MPO representation for a CI tensor operator using sequential 
% singular value decompositions.
%
% Requirements: Tensor Toolbox (Sandia Labs).
%
% ilya.kuprov@oerc.ox.ac.uk

function mpo=cito2mpo(spin_system,cito,drop_tol)

% Get the number of legs and dimensions
nlegs=ndims(cito); dims=size(cito);

% Preallocate the result array
mpo=cell(nlegs,1);

% Process the first leg
A=tenmat(cito,1,2:nlegs); [U,S,V]=svd(double(A),'econ');
k=nnz(diag(S)>drop_tol); U=U(:,1:k); S=S(1:k,1:k); V=V(:,1:k);
backward_dim=1; physical_dim=dims(1); forward_dim=size(U,2);
mpo{1}=tensor(U,[backward_dim physical_dim forward_dim]);
cito=tensor(tenmat(S*V',1,2:nlegs,[forward_dim A.tsize(2:end)]));

% Process the middle legs
for n=2:(nlegs-1)

    % Shift a pair of indices into rows
    A=tenmat(cito,1:2,3:(nlegs-n+2));
    
    % Run the SVD
    [U,S,V]=svd(double(A),'econ');
    
    % Drop insignificant singular values
    if drop_tol<1
        % If tol<1, drop by tolerance
        k=nnz(diag(S)>drop_tol); k=max([1 k]);
    else
        % If tol>1, keep tol singular values
        k=drop_tol; k=min([numel(diag(S)) k]);
    end
    U=U(:,1:k); S=S(1:k,1:k); V=V(:,1:k);
    
    % Get the dimension data
    backward_dim=forward_dim; 
    physical_dim=dims(n);
    forward_dim=size(U,2);
    
    % Assign the MPO tensor
    mpo{n}=tensor(tenmat(U,1:2,3,[backward_dim physical_dim forward_dim]));
    
    % Fold away the dangling summation index
    cito=tensor(tenmat(S*V',1,2:(nlegs-n+1),[forward_dim A.tsize(3:end)]));
    
end

% Process the last leg
U=S*V'; mpo{nlegs}=tensor(U,[size(U) 1]);

% Put physical indices last
for n=1:nlegs
    mpo{n}=permute(mpo{n},[1 3 2]);
end

% Collect the legs in pairs
for n=1:(nlegs/2)
    mpo{2*n-1}=ttt(mpo{2*n-1},mpo{2*n},2,1);
    mpo{2*n-1}=permute(mpo{2*n-1},[1 3 2 4]);
end
mpo=mpo(1:2:end);

% Print some statistics
report(spin_system,'cito2mpo: MPO dimensions');
for n=1:numel(mpo)
    report(spin_system,['cito2mpo: ' num2str(size(mpo{n}),'%-5d ')]);
end

end

% Pronouncement of experts to the effect that something cannot be done has
% always irritated me.
%
% Leo Szilard

