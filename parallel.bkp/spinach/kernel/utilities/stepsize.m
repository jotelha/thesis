% Optimal step for time propagation.
%
% ilya.kuprov@oerc.ox.ac.uk

function [timestep,nsteps]=stepsize(L,interval)

% Catch incorrect calls
if (nargin==1)&&(nargout==2)
    error('stepsize: time interval is required to compute the number of steps.');
end

% Adapt to the output style
switch nargout
    
    case 1
        
        % Get the optimal time step
        timestep=1/norm(L,'inf');
    
    case 2
               
        % Get the optimal integer number of steps (accommodating negative and
        % imaginary intervals)
        nsteps=ceil(abs(interval)*norm(L,'inf'));
        
        % Get the time step
        timestep=interval/nsteps;
        
end

end

% Whenever anyone accuses some person of being 'unfeeling' he means that
% that person is just. He means that that person has no causeless emotions
% and will not grant him a feeling which he does not deserve... justice is
% the opposite of charity.
%
% Ayn Rand, "Atlas Shrugged"

