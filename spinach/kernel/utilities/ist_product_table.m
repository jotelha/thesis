% Structure coefficient tables for the envelopes of su(mult) algebras. Disk
% cache is used and updated automatically.
%
% hannah.hogben@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function [product_table_left,product_table_right]=ist_product_table(spin_system,mult)

% Generate cache record name
table_file=[spin_system.sys.root_dir spin_system.sys.slash 'kernel' ...
                                     spin_system.sys.slash 'cache'  ...
                                     spin_system.sys.slash 'ist_product_table_' num2str(mult) '.mat'];

% Check the cache
if exist(table_file,'file')
    
    % Lift data from the cache if the file is already available
    load(table_file);
    
else
    
    % Get the irreducible spherical tensors
    T=irr_sph_ten(mult);
    
    % Preallocate the arrays
    product_table_left=zeros(mult^2,mult^2,mult^2);
    product_table_right=zeros(mult^2,mult^2,mult^2);
    
    % Get the structure coefficients
    for m=1:mult^2
        for k=1:mult^2
            normalization=sqrt(trace(T{k}*T{k}')*trace(T{m}*T{m}'));
            for n=1:mult^2
                product_table_left(n,m,k)=trace(T{n}*T{m}*T{k}')/normalization;
                product_table_right(n,m,k)=trace(T{m}*T{n}*T{k}')/normalization;
            end
        end
    end

    % Save the table to the cache
    save(table_file,'product_table_left','product_table_right');

end

end

% According to Oxford Chemistry folklore, Peter Atkins has once asked the
% following question at an interview for a Lecturer post: "What is it that
% you have done that a technician would not do?". The candidate produced a
% reasonable answer to the effect that a technician would not have the
% required skills. A better answer was suggested by a member of the
% Interview Board some time later: "And what is it that you, Atkins, have
% done that a journalist would not do?"

