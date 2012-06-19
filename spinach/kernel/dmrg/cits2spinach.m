% Converts a CI state tensor into a Spinach state vector. Currently
% requires the Spinach calculation to have a complete basis set.
%
% Requirements: Tensor Toolbox (Sandia Labs).
%
% ilya.kuprov@oerc.ox.ac.uk

function state_vector=cits2spinach(cits)

% Unroll the CI state tensor
state_vector=double(tenmat(cits,ndims(cits):-1:1,[]));

end

% I came into the room, which was half dark, and presently spotted Lord
% Kelvin in the audience and realized that I was in for trouble at the last
% part of my speech dealing with the age of the Earth, where my views
% conflicted with his. To my relief, Kelvin fell fast asleep, but as I
% came to the important point, I saw the old bird sit up, open an eye and
% cock a baleful glance at me! Then a sudden inspiration came, and I said
% "Lord Kelvin had limited the age of the Earth, provided no new source of
% energy was discovered. That prophetic utterance refers to what we are
% now considering tonight, radium." Behold! The old boy beamed upon me.
% 
% Ernest Rutherford

