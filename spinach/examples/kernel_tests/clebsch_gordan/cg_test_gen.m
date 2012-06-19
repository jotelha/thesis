% Generates a large number of physically valid index combinations
% to be used for testing of Clebsch-Gordan coefficients.
%
% Index order in the resulting file: L,M,L1,M1,L2,M2.
%
% ilya.kuprov@oerc.ox.ac.uk

function cg_test_gen()

% User-selectable parameters
n_index_lines=10000;
max_multiplicity=1000;

% Open the file for writing
file_id=fopen('cg_indices.txt','a');

% Loop over the index lines
for n=1:n_index_lines
    
    % Generate random L1 and L2 ranks
    L1=(ceil(max_multiplicity*rand(1))-1)/2;
    L2=(ceil(max_multiplicity*rand(1))-1)/2;
    
    % Generate random physically valid projections
    M1_range=-L1:L1; M2_range=-L2:L2;
    M1=M1_range(ceil(numel(M1_range)*rand(1)));
    M2=M2_range(ceil(numel(M2_range)*rand(1)));
    M=M1+M2;
    
    % Generate a random physically valid L rank
    L_range=max([abs(L1-L2) abs(M)]):(L1+L2);
    L=L_range(ceil(numel(L_range)*rand(1)));
    
    % Write the file
    fprintf(file_id,'%6.1f   %6.1f   %6.1f   %6.1f   %6.1f   %6.1f  \n',[L M L1 M1 L2 M2]);

end

% Close the file
fclose(file_id);

end

