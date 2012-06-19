% Collecting the data of split jlh_outsourced_pulse_acquire() node package
% processes after their successfull completion and writing .dat spectrum
% output file
%       
% johannes.hoermann@googlemail.com
function jlh_sum_results(input_file_name)
    
    %loading output file names of node package processes into workspace
    load(input_file_name,'spin_system','parameters','output_file_name','logfile_name','spectrum_file_name');
    %creating own banner -- modified banner.m needed!!!
    banner(spin_system,'jlh_sum_results_banner');

    %preallocating fid array
    fid=zeros(parameters.npoints,1);

    for i=1:parameters.nodes
        %checking, whether each process wrote ouput file correctly
        if exist(output_file_name{i}, 'file')
            tmp = load( output_file_name{i}, 'fid');
            fid = fid + tmp.fid;
        else
           %reporting errror otherwise
           report(spin_system, [sprintf('jlh_sum_results: Error: jlh_outsourced_pulse_acquire instance %d did not write output file "', i) output_file_name{i} '"!']);
        end
        %print node package process lofiles on stdout
        logfile_content = fileread(logfile_name{i});
        fprintf(logfile_content);
    end
    
    % finalizing everything 
    
    % Apodization
    fid=apodization(fid,'crisp-1d');

    % Perform Fourier transform
    spectrum=fftshift(fft(fid,parameters.zerofill));


    %Compute the derivative if necessary
    if isfield(parameters,'derivative')
        spectrum=fft(ifft(spectrum).*fftdiff(parameters.derivative,length(spectrum),1)');
    end

    ax=axis_1d(spin_system,parameters);
    data = cat(2, transpose(ax), real(spectrum));
    save(spectrum_file_name,'data', '-ASCII');
end