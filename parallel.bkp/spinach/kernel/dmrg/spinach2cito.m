% Converts a Spinach suproperator into a CI tensor operator. Currently 
% requires the Spinach calculation to have a complete basis set. Odd
% indices of the resulting tensor point in the bra direction and even
% indices point in the ket direction.
%
% Requirements: Tensor Toolbox (Sandia Labs).
%
% ilya.kuprov@oerc.ox.ac.uk

function cito=spinach2cito(spin_system,superop)

% Roll up the superoperator
cito=tensor(superop,[spin_system.comp.mults.^2 spin_system.comp.mults.^2]);

% Get tensor indices into zipper order
index_order=zeros(ndims(cito),1);
for n=1:spin_system.comp.nspins
    index_order(2*n-1)=spin_system.comp.nspins-n+1;
    index_order(2*n)=2*spin_system.comp.nspins-n+1;
end
cito=permute(cito,index_order);

end

% I can, therefore I am.
%
% Simone Weil

