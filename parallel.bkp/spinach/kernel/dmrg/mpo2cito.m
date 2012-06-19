% Contracts the indirect dimensions of all tensors in a matrix product
% operator and returns the corresponding CI tensor operator.
%
% Requirements: Tensor Toolbox (Sandia Labs).
%
% ilya.kuprov@oerc.ox.ac.uk

function cito=mpo2cito(mpo)

% Set the starting point
cito=mpo{1};

% Loop over the spins
for n=2:numel(mpo)
    cito=ttt(cito,mpo{n},2*(n-1),1);
end

% Squeeze out singleton dimensions
cito=squeeze(cito);

end

% We saw that we'd been given a law to live by, a moral law, they called
% it, which punished those who observed it - for observing it. The more you
% tried to live up to it, the more you suffered; the more you cheated it,
% the bigger reward you got.
%
% Ayn Rand, "Atlas Shrugged"

