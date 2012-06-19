% Computes the result of an action by an MPO on an MPS.
%
% Requirements: Tensor Toolbox (Sandia Labs).
%
% ilya.kuprov@oerc.ox.ac.uk

function mpp=mpo_x_mps(spin_system,mpo,mps)

% Bomb out if there is a dimension mismatch
if numel(mpo)~=numel(mps)
    error('mpo_x_mps: MPS and MPO dimensions do not match.');
end

% Compute the zipper
mpp=cell(size(mpo));
for n=1:numel(mpo)
    mpp{n}=ttt(mpo{n},mps{n},4,3); a=tenmat(mpp{n},[1 4 2 5],3);
    mpp{n}=tensor(tenmat(a.data,[1 2],3,[size(mpp{n},1)*size(mpp{n},4) size(mpp{n},2)*size(mpp{n},5) size(mpp{n},3)]));
end

% Print some statistics
report(spin_system,'mpo_x_mps: MPO-MPS product dimensions');
for n=1:numel(mpp)
    report(spin_system,['mpo_x_mps: ' num2str(size(mpp{n}),'%-5d ')]);
end

end

% She knew the general doctrine on sex, held by people in one form or
% another, the doctrine that sex was an ugly weakness of man's lower
% nature, to be condoned regretfully. She experienced an emotion of
% chastity that made her shrink not from the desires of her body, but 
% from any contact with the minds who held this doctrine.
%
% Ayn Rand, "Atlas Shrugged"

