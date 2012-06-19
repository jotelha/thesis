% Printing offset for Spinach output.
%
% ilya.kuprov@oerc.ox.ac.uk

function spin_system=click(spin_system,direction)

switch direction
    
    case 'forward'
        
        % Add one white space character in front of the output
        if isfield(spin_system,'click')
            spin_system.click=[spin_system.click ' '];
        else
            spin_system.click=' ';
        end
        
    case 'backward'
        
        % Remove one white space character from the front of the output
        if isfield(spin_system,'click')
            if any(spin_system.click)
                spin_system.click=spin_system.click(2:end);
            else
                spin_system.click='';
            end
        else
            spin_system.click='';
        end

end

% She wanted to injure him by her contempt - but he could not be injured,
% unless he respected her judgement. [...] An issue of quilt, he thought,
% had to rest on his own acceptance of the code of justice that pronounced
% him guilty. He did not accept it; he never had.
%
% Ayn Rand, "Atlas Shrugged"