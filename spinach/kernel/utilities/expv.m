% Computes expm(A)*rho using matrix-vector multiplications only.
%
% ilya.kuprov@oerc.ox.ac.uk

function rho=expv(A,rho)

% Scale the state vector
scaling_factor=norm(rho);
rho=rho/scaling_factor;

% Estimate the norm of A
norm_a=norm(A,'inf');

% Get the number of steps
nsteps=ceil(norm_a/5);

% Scale the matrix
A=A/nsteps;

% Run the Taylor series procedure
for n=1:nsteps
    next_term=rho; k=0; rho=zeros(size(rho));
    while norm(next_term) > eps
        rho=rho+next_term; k=k+1;
        next_term=(A*next_term)/k;
    end
end

% Scale the state vector back
rho=rho*scaling_factor;

end

% [...] this flaw was identified by the brilliant German philosopher
% Friedrich Nietzsche who described it as "an inversion of morality"
% whereby the weak, the poor, the meek, the oppressed and the wretched are
% virtuous and blessed by God whereas the strong, the wealthy, the noble
% and the powerful are the immoral and damned by the vengeful almighty
% Yahweh for eternity. Nietzsche, with great insight and perception, stated
% that Christianity would be abandoned en masse in the twentieth century
% but that Westerners would still cling to this inversion of morality.
%
% Anders Behring Breivik

