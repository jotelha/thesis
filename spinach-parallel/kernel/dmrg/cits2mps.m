% Computes an MPS representation for a CI state tensor using sequential 
% singular value decompositions.
%
% Requirements: Tensor Toolbox (Sandia Labs).
%
% ilya.kuprov@oerc.ox.ac.uk

function mps=cits2mps(spin_system,cits,drop_tol)

% Get the number of legs
nlegs=ndims(cits); dims=size(cits);

% Preallocate the result array
mps=cell(nlegs,1);

% Process the first leg
A=tenmat(cits,1,2:nlegs); [U,S,V]=svd(double(A),'econ');
k=nnz(diag(S)>drop_tol); U=U(:,1:k); S=S(1:k,1:k); V=V(:,1:k);
backward_dim=1; physical_dim=dims(1); forward_dim=size(U,2);
mps{1}=tensor(U,[backward_dim physical_dim forward_dim]);
cits=tensor(tenmat(S*V',1,2:nlegs,[forward_dim A.tsize(2:end)]));

% Process the middle legs
for n=2:(nlegs-1)

    % Shift a pair of indices into rows
    A=tenmat(cits,1:2,3:(nlegs-n+2));
    
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
    
    % Assign the MPS tensor
    mps{n}=tensor(tenmat(U,1:2,3,[backward_dim physical_dim forward_dim]));
    
    % Fold away the dangling summation index
    cits=tensor(tenmat(S*V',1,2:(nlegs-n+1),[forward_dim A.tsize(3:end)]));
    
end

% Process the last leg
U=S*V'; mps{nlegs}=tensor(U,[size(U) 1]);

% Put physical indices last
for n=1:nlegs
    mps{n}=permute(mps{n},[1 3 2]);
end

% Print some statistics
report(spin_system,'cits2mps: MPS dimensions');
for n=1:nlegs
    report(spin_system,['cits2mps: ' num2str(size(mps{n}),'%-5d ')]);
end

end

% Niels Bohr, on Dirac's equation:
%
% "Simply put an explanation of the theory on a poster, tack it up on a tree
% in the jungle, and any elephant (a beast noted for its wisdom) that passed
% by would immediately become so engrossed trying to figure it out that it
% could be packed up and delivered to the Copenhagen Zoo before it realized
% anything had happened."

