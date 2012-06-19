% Converts a CI tensor operator into a Spinach superoperator. Currently
% requires the Spinach calculation to have a complete basis set.
%
% Requirements: Tensor Toolbox (Sandia Labs).
%
% ilya.kuprov@oerc.ox.ac.uk

function superop=cito2spinach(spin_system,cito)

% Set the index order
row_index=zeros(1,ndims(cito)/2);
col_index=zeros(1,ndims(cito)/2);
for n=1:spin_system.comp.nspins
    row_index(n)=ndims(cito)-2*n+1;
    col_index(n)=ndims(cito)-2*n+2;
end

% Unroll the CI tensor operator
superop=double(tenmat(cito,row_index,col_index));

end

% The speed of light sucks.
%
% John Carmack

