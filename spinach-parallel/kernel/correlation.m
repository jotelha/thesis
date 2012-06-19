% Correlation order selection function -- keeps only the specified
% orders of spin correlation in the state vector. Arguments:
%
%   rho                -  a state vector or a horizontal stack thereof
%
%   correlation_orders -  a row vector of correlation orders to keep
%
%   spins              -  which spins to count (e.g. '1H', '13C', 'all')
%
% If rho is left empty ( [] ) then the state mask is returned.
%
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function rho=correlation(spin_system,rho,correlation_orders,spins)

% Set the default to all spins
if ~exist('spins','var'), spins='all'; end

% Parse the spin specification
if ~isnumeric(spins)
    if strcmp(spins,'all')
        spins=1:length(spin_system.comp.isotopes);
    else
        spins=find(strcmp(spins,spin_system.comp.isotopes));
    end
end

% Compute the order of correlation for each basis state
correlation_orders_present=sum(logical(spin_system.bas.basis(:,spins)),2);

% Zero all correlation orders except those specified by the user
state_mask=zeros(size(spin_system.bas.basis,1),1);
for n=correlation_orders
    state_mask=state_mask|(correlation_orders_present==n);
end

% If an empty vector is passed for rho, return the state mask
if ~isempty(rho)
    rho(~state_mask,:)=0;
else
    rho=state_mask;
end

end

% The first iteration of the Spin Dynamics course (http://spindynamics.org) 
% was so difficult that every single student has dropped out by about Lecture
% 10. The rest of the course was read to Rusty, Dusty, Scratchy, Patchy and
% Scruffy, the five plastic chairs in IK's Durham office - they made for an 
% excellent (if a bit shy so far as questions were concerned) audience. To 
% this day, the total number of students who verifiably understood the whole
% of that course can be counted on the fingers of one hand.

