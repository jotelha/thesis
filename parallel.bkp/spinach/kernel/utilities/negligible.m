% An aux function determining whether a given object deserves attention
% given the tolerance specified. Used in the internal decision making 
% performed by Spinach kernel functions.
%
% ilya.kuprov@oerc.ox.ac.uk

function answer=negligible(something,type,tolerance)

answer=~significant(something,type,tolerance);

end

% It's too bad that stupidity isn't painful.
%
% Anton Szandor LaVey