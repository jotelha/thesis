% Contracts the bond dimensions of all tensors in a matrix product
% state and returns the corresponding CI tensor state.
%
% Requirements: Tensor Toolbox (Sandia Labs).
%
% ilya.kuprov@oerc.ox.ac.uk

function cits=mps2cits(mps)

% Set the starting point
cits=mps{1};

% Loop over the spins
for n=2:numel(mps)
    cits=ttt(cits,mps{n},n,1);
end

% Squeeze out singleton dimensions
cits=squeeze(cits);

end

% Mid-2011, MPLS senior staff meeting, Oxford.
%
% Tim Softley: "...and we must, of course, mention our students' role in
% this success - after all, it is they that do most of the work in our
% labs..."
%
% Incredulous voice from the audience: "Speak for yourself!"

