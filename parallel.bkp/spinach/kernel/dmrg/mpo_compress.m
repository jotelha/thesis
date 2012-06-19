% Compresses an MPO using sequential singular value decompositions.
%
% Requirements: Tensor Toolbox (Sandia Labs).
%
% ilya.kuprov@oerc.ox.ac.uk

function mpo=mpo_compress(spin_system,mpo,drop_tol)

% Loop over the elements of MPO
for n=1:(numel(mpo)-1)
    
    % Pull out a pair of matrices
    m_left=mpo{n}; m_right=mpo{n+1};
    
    % Shift the indices
    m_left_mat=tenmat(m_left,[1 3 4],2);
    m_right_mat=tenmat(m_right,1,[2 3 4]);
    
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
    mpo{n}=tensor(tenmat(U_left,[1 3 4],2,[size(m_left,1) k size(m_left,3) size(m_left,4)]));
    mpo{n+1}=tensor(tenmat(V_right',1,[2 3 4],[k size(m_right,2) size(m_right,3) size(m_right,4)]));
    
end

% Print some statistics
report(spin_system,'compress_mpo: compressed MPO dimensions');
for n=1:numel(mpo)
    report(spin_system,['compress_mpo: ' num2str(size(mpo{n}),'%-5d ')]);
end

end

% She was twelve years old when she told Eddie Willers that she would run
% the railroad when they grew up. She was fifteen when it occurred to her
% for the first time that women did not run railroads and that people might
% object. To hell with that, she thought - and never worried about it again.
%
% Ayn Rand, "Atlas Shrugged"

