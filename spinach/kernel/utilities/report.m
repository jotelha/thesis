% Writes a log message to the console or an ACSII file.
%
% ilya.kuprov@oerc.ox.ac.uk

function report(spin_system,report_string)

switch spin_system.sys.output
    
    case 'console'

        % Dump the report string to the console
        if isfield(spin_system,'click')
            disp([spin_system.click report_string]);
        else
            disp(report_string);
        end
        
    case 'file'
        
        % Dump the report string to the file
        file_id=fopen(spin_system.sys.logfile,'a');
        if isfield(spin_system,'click')
            fprintf(file_id,'%s\n',[spin_system.click report_string]);
        else
            fprintf(file_id,'%s\n',report_string);
        end
        fclose(file_id);
        
    case 'hush'
        
        % Ignore the report string
        
    otherwise
        
        % Bomb out if anything unexpected is received
        error('report: output device not specified');
        
end

end

% "All parts should go together without forcing. You must remember that
% the parts you are reassembling were disassembled by you. Therefore, if
% you can’t get them together again, there must be a reason. By all means,
% do not use a hammer."
%
% IBM Manual, 1925 

