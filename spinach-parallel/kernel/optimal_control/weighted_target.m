% Applies weights to the target state amplitudes.
%
% ilya.kuprov@oerc.ox.ac.uk

function target=weighted_target(spin_system,state_type,weight_array)

% Preallocate the array
target=zeros(spin_system.bas.nstates,1);

% Assemble the target
switch state_type
    case 'Lx'
        for n=1:spin_system.comp.nspins
            L_plus=state(spin_system,'L+',n); L_minus=state(spin_system,'L-',n);
            target=target+weight_array(n)*(L_plus+L_minus)/2;
        end
    case 'Ly'
        for n=1:spin_system.comp.nspins
            L_plus=state(spin_system,'L+',n); L_minus=state(spin_system,'L-',n);
            target=target+weight_array(n)*(L_plus-L_minus)/2i;
        end
    case 'Lz'
        for n=1:spin_system.comp.nspins
            target=target+weight_array(n)*state(spin_system,'Lz',n);
        end
end
        
end

% "Everyone is entitled to his own opinion, but not his own facts."
%
% Daniel Patrick Moynihan

