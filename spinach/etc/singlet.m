% Two-spin singlet state.
%
% ilya.kuprov@oerc.ox.ac.uk

function rho=singlet(spin_system,spin_a,spin_b)

unit=state(spin_system,{'E','E'},{spin_a,spin_b});
LzSz=state(spin_system,{'Lz','Lz'},{spin_a,spin_b});
LmSp=state(spin_system,{'L-','L+'},{spin_a,spin_b});
LpSm=state(spin_system,{'L+','L-'},{spin_a,spin_b});

rho=unit/4-(LzSz+0.5*(LpSm+LmSp)); 

end

% Physics is, hopefully, simple. Physicists are not.
% 
% Edward Teller 