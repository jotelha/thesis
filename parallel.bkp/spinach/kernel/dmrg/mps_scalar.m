% Computes a scalar product of two matrix product states.
%
% Requirements: Tensor Toolbox (Sandia Labs).
%
% ilya.kuprov@oerc.ox.ac.uk

function scalar_product=mps_scalar(mps_a,mps_b)

% Bomb out if there is a dimension mismatch
if numel(mps_a)~=numel(mps_b)
    error('mpsprod: MPS dimensions do not match.');
end

% Compute the zipper
mpp=cell(size(mps_a));
for n=1:numel(mps_a)
    mpp{n}=ttt(tenfun(@conj,mps_a{n}),mps_b{n},3);
end

% Fold up the zipper
scalar_product=permute(mpp{1},[1 3 2 4]);
for n=2:numel(mpp)
    scalar_product=ttt(scalar_product,mpp{n},[3 4],[1 3]);
end
scalar_product=double(squeeze(scalar_product));

end

% Part of the inhumanity of the computer is that, once it is competently
% programmed and working smoothly, it is completely honest.
% 
% Isaac Asimov

