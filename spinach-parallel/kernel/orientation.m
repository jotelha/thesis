% Anisotropic part of the Hamiltonian commutation superoperator
% for a specific spin system orientation. Arguments:
%
%     Q             -  rotational basis as returned by h_superop
%                      function.
%
%     euler_angles  -  a 1x3 vector or a vertical stack of 1x3 
%                      vectors specifying Euler angles (radians)
%                      relative to the input orientation.
%
% Output:
%
%     H    -    anisotropic part of the Hamiltonian commutation
%               superoperator (or a cell array thereof) for the 
%               specified Euler angles.
%
% This function may also be used for the corresponding Hilbert
% space operators, because the H -> [H, ] map is linear.
%
% ilya.kuprov@oerc.ox.ac.uk

function H=orientation(Q,euler_angles)

% Determine problem dimension
n_orientations=size(euler_angles,1); 
matrix_dim=size(Q{1,1},1);

% Preallocate the result array
H=cell(n_orientations,1);
for n=1:n_orientations
    H{n}=spalloc(matrix_dim,matrix_dim,10*matrix_dim);
end

% Compute Wigner matrices
D=cell(n_orientations,1);
for n=1:n_orientations
    D{n}=euler2wigner(euler_angles(n,1),...
                      euler_angles(n,2),...
                      euler_angles(n,3));
end

% Compute Hamiltonian commutation superoperators
for k=1:5
    for m=1:5
        for n=1:n_orientations
            H{n}=H{n}+Q{k,m}*D{n}(k,m);
        end
    end
end

% If there is only one operator to return, remove the cell
if numel(H)==1, H=H{1}; end

end

% 1f y0u c4n r34d 7h15, y0u r34||y n33d 70 637 |41d

