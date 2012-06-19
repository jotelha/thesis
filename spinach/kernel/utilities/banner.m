% Prints the banners.
%
% ilya.kuprov@oerc.ox.ac.uk

function banner(spin_system,identifier)

switch identifier
    case 'version_banner'
        report(spin_system,' ');
        report(spin_system,'===========================================');
        report(spin_system,'=                                         =');
        report(spin_system,'=            SPINACH 1.1.1054             =');
        report(spin_system,'=                                         =');
        report(spin_system,'=       Ilya Kuprov, Hannah Hogben,       =');
        report(spin_system,'=    Luke Edwards, Matthew Krzystyniak    =');
        report(spin_system,'=       Peter Hore, Gareth Charnock       =');
        report(spin_system,'=                                         =');
        report(spin_system,'=        Oxford e-Research Centre         =');
        report(spin_system,'=          University of Oxford           =');
        report(spin_system,'=                                         =');
        report(spin_system,'=         GNU Public License v2.5         =');
        report(spin_system,'=                                         =');
        report(spin_system,'===========================================');
        report(spin_system,' ');
    case 'spin_system_banner'
        report(spin_system,' ');
        report(spin_system,'===========================================');
        report(spin_system,'=                                         =');
        report(spin_system,'=               SPIN SYSTEM               =');
        report(spin_system,'=                                         =');
        report(spin_system,'===========================================');
        report(spin_system,' ');
    case 'basis_banner'
        report(spin_system,' ');
        report(spin_system,'===========================================');
        report(spin_system,'=                                         =');
        report(spin_system,'=                BASIS SET                =');
        report(spin_system,'=                                         =');
        report(spin_system,'===========================================');
        report(spin_system,' ');
    case 'sequence_banner'
        report(spin_system,' ');
        report(spin_system,'===========================================');
        report(spin_system,'=                                         =');
        report(spin_system,'=              PULSE SEQUENCE             =');
        report(spin_system,'=                                         =');
        report(spin_system,'===========================================');
        report(spin_system,' ');
    otherwise
        error('banner: unknown banner.');
end

end

% "The free man will ask neither what his country can do for him, nor what
%  he can do for his country." - Milton Friedman

