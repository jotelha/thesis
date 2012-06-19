% An aux function determining whether a given object deserves attention
% given the tolerance specified. Used in the internal decision making 
% performed by Spinach kernel functions.
%
% ilya.kuprov@oerc.ox.ac.uk

function answer=significant(something,type,tolerance)

switch type
    
    case {'tensor','vector','scalar','operator'}
        
        if isempty(something)
            % Empty array is not significant
            answer=false;
        elseif nnz(something)==0
            % All-zero array is not significant
            answer=false;
        elseif max(abs(nonzeros(something)))<=tolerance
            % Anything containing only small elements is not significant
            answer=false;
        else
            % Everything else is significant
            answer=true;
        end
        
    otherwise
        error('significant: unknown object type.');
        
end

end

% There is a beast in man that needs to be exercised, not exorcised.
%
% Anton Szandor LaVey

