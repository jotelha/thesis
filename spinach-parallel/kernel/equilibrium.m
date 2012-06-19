% Thermal equilibrium state. Arguments:
%
%     contributions - can be set to:
%                          
%                     'isotropic' - uses the isotropic part of the
%                                   Hamiltoniian to compute the thermal
%                                   equilibrium state (default).
%
%                     'full'      - uses the complete Hamiltonian at input
%                                   orientation to compute the thermal 
%                                   equilibrium state.
% 
% If the temperature is set to zero during the call to create.m, returns the 
% high-temperature approximation to the thermal equilibrium state.
%
% If the temperature is specified, returns the accurate equilibrium state at 
% that temperature.
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function rho=equilibrium(spin_system,contributions)

% Set the defaults
if nargin==1
    contributions='isotropic';
end

% Click forward for output
spin_system=click(spin_system,'forward');

% Inform the user
report(spin_system,'equilibrium: computing the thermal equilibrium state...');

% Decide on the Hamiltonian to use
switch contributions
    case 'isotropic'
        report(spin_system,'equilibrium: using the isotropic part of the Hamiltonian...');
        H_right=h_superop(secularity(spin_system,'keep_all'),'right');
    case 'full'
        report(spin_system,'equilibrium: using the full Hamiltonian at the input orientation...');
        [H_right,Q]=h_superop(secularity(spin_system,'keep_all'),'right');
        H_right=H_right+orientation(Q,[0 0 0]);
    otherwise
        error('equilibrium: unknown contribution specification.');
end

% Get the unit state
unit=sparse(1,1,1,spin_system.bas.nstates,1,1);

% Decide between the high-temperature approximation and the accurate calculation
if spin_system.rlx.temperature==0
    
    % Warn the user
    report(spin_system,'equilibrium: WARNING - high temperature approximation.');
    
    % Use the high-temperature approximation
    rho=unit+H_right(:,1)/max(abs(H_right(:,1)));

else
    
    % Get the temperature factor
    beta=spin_system.tols.hbar/(spin_system.tols.kbol*spin_system.rlx.temperature);
    
    % Get the thermal equilibrium state
    rho=expv(beta*H_right,unit);
     
    % Normalize the thermal equilibrium state
    rho=rho/rho(1);

end
    
end

% Compassion is a wonderful thing. It's what one feels when one looks at a
% squashed caterpillar. An elevating experience. One can let oneself go and
% spread –- you know, like taking a girdle off. You don't have to hold your
% stomach, your heart or your spirit up –- when you feel compassion. All you
% have to do is look down. It's much easier. When you look up, you get a
% pain in the neck. Compassion is the greatest virtue. It justifies
% suffering. There's got to be  suffering in the world, else how would we
% be virtuous and feel compassion? ...Oh, it has an antithesis –- but such a
% hard, demanding one... Admiration, Mrs. Jones, admiration. But that takes
% more than a girdle. So I say that anyone for whom we can't feel sorry is
% a vicious person.
%
% Ayn Rand, "The Fountainhead"

