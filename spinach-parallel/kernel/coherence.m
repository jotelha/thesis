% Coherence selection function -- keeps only the specified coherence orders
% in the state vector. Arguments:
%
%     rho               -  a state vector or a horizontal stack thereof
%
%     coherence_orders  -  a row vector of coherence orders to keep
%
%     spins             -  which spins to count (e.g. '1H', '13C', 'all')
%
% If rho is left empty ( [] ) then the state mask is returned.
%
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function rho=coherence(spin_system,rho,coherence_orders,spins)

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

% Compute the projection quantum numbers for spins in each basis state
[~,M]=lin2lm(spin_system.bas.basis);

% Determine the total coherence order of each basis state
coherence_orders_present=sum(M(:,spins),2);

% Zero all coherence orders except those specified by the user
state_mask=zeros(size(spin_system.bas.basis,1),1);
for n=coherence_orders
    state_mask=state_mask|(coherence_orders_present==n);
end

% If an empty vector is passed for rho, return the state mask
if ~isempty(rho)
    rho(~state_mask,:)=0;
else
    rho=state_mask;
end

end

% "There was a time when men were afraid that somebody would reveal some
% secret of theirs that was unknown to their fellows. Nowadays, they're
% afraid that somebody will name what everybody knows. Have you practical
% people ever thought that that's all it would take to blast your whole,
% big, complex structure, with all your laws and guns—just somebody naming
% the exact nature of what it is you're doing?"
%
% Ayn Rand, "Atlas Shrugged"

