% Pads or cuts a string with spaces to make it n characters long.
%
% Konstantin Pervushin (p@trosy.com)
% Ilya Kuprov (ilya.kuprov@oerc.ox.ac.uk)

function p=pad(s,n)

p=[s blanks(n)]; p=p(1:n);

end

% By the grace of reality and the nature of life, man – every man – is an
% end in himself, he exists for his own sake, and the achievement of his
% own happiness is his highest moral purpose.
%
% Ayn Rand, "Atlas Shrugged"