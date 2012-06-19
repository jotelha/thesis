% Converts a Spinach state vector into a CI tensor state. Currently 
% requires the Spinach calculation to have a complete basis set. The
% resulting tensor state corresponds to a ket and should be contracted
% with the even indices of any CI tensor operator.
%
% Requirements: Tensor Toolbox (Sandia Labs).
%
% ilya.kuprov@oerc.ox.ac.uk

function cits=spinach2cits(spin_system,state_vector)

% Roll up the state vector
cits=tensor(state_vector,spin_system.comp.mults.^2);

% Reverse the index order
cits=permute(cits,spin_system.comp.nspins:-1:1);

end

% We do not hear the term "compassionate" applied to business executives or
% entrepreneurs, certainly not when they are engaged in their normal work.
% Yet in terms of results in the measurable form of jobs created, lives
% enriched, communities built, living standards raised, and poverty healed,
% a handful of capitalists has done infinitely more for mankind than all
% the self-serving politicians, academics, social workers, and religionists
% who march under the banner of "compassion".
%
% Nathaniel Branden

