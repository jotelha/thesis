% Compresses an MPS using sequential singular value decompositions.
%
% Requirements: Tensor Toolbox (Sandia Labs).
%
% ilya.kuprov@oerc.ox.ac.uk

function mps=mps_compress(spin_system,mps,drop_tol)

% Loop over the elements of MPS
for n=1:(numel(mps)-1)
    
    % Pull out a pair of matrices
    m_left=mps{n}; m_right=mps{n+1};
    
    % Shift the indices
    m_left_mat=tenmat(m_left,[1 3],2);
    m_right_mat=tenmat(m_right,1,[2 3]);
    
    % Run SVD on both
    [U_left,S_left,V_left]=svd(double(m_left_mat),'econ');
    [U_right,S_right,V_right]=svd(double(m_right_mat),'econ');
    
    % Get the middle matrix
    M=S_left*V_left'*U_right*S_right;
    
    % SVD the middle matrix
    [U,S,V]=svd(M,'econ');
    
    % Drop insignificant singular values
    if drop_tol<1
        % If tol<1, drop by tolerance
        k=nnz(diag(S)>drop_tol); k=max([1 k]);
    else
        % If tol>1, keep tol singular values
        k=drop_tol; k=min([numel(diag(S)) k]);
    end
    U=U(:,1:k); S=S(1:k,1:k); V=V(:,1:k);
    
    % Push to the sides
    U_left=U_left*U; V_right=(S*V'*V_right')';
    
    % Reassemble the pair of matrices
    mps{n}=tensor(tenmat(U_left,[1 3],2,[size(m_left,1) k size(m_left,3)]));
    mps{n+1}=tensor(tenmat(V_right',1,[2 3],[k size(m_right,2) size(m_right,3)]));
    
end

% Print some statistics
report(spin_system,'mps_compress: compressed MPS dimensions');
for n=1:numel(mps)
    report(spin_system,['mps_compress: ' num2str(size(mps{n}),'%-5d ')]);
end

end

% The most exciting phrase to hear in science, the one that heralds
% new discoveries, is not "Eureka!" but "That's funny...".
%
% Isaac Asimov

