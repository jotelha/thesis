% Hadamard matrix product. Useful as a replacement for trace(A'*B) because
% trace(A'*B)=hadm(conj(A),B) and the latter only needs O(n^2) multiplica-
% tions as opposed to O(n^3) for trace(A'*B).
%
% ilya.kuprov@oerc.ox.ac.uk

function h=hadm(a,b)

    h=sum(sum(a.*b));

end

% An infinite number of mathematicians walk into a bar. The first one
% orders a pint of beer, the second one half a pint, the third one a
% quarter... "Bastards..." says the bartender and pours two pints.

