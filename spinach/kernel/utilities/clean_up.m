% Sparse matrix clean-up utility. Drops the non-zero elements with
% magnitude below a user-specified tolerance. 
%
% ilya.kuprov@oerc.ox.ac.uk

function A=clean_up(spin_system,A,nonzero_tol)

% Only consider sparse matrices
if issparse(A)
    
    % Respect the disable switch
    if ~any(strcmp(spin_system.sys.disable,'clean-up'))
        
        % Clean up the matrix
        A=A.*(abs(A)>nonzero_tol);
        
        % If the density exceeds the threshold, make the matrix dense
        if nnz(A)/numel(A)>spin_system.tols.dense_matrix
            A=full(A);
        end
        
    end
    
end

end

% "The most dangerous man to any government is the man who is able to think
%  things out for himself, without regard to the prevailing superstitions
%  and taboos." - H.L. Mencken

