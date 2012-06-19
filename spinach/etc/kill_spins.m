% Function to eliminate unwanted spins from sys and inter. An example of
% its use is in removing protons that exchange rapidly with the solvent and
% so don't affect the wanted NMR spectrum when input has been parsed from
% a log file. Only accepts the .matrix and .scalar fields of the inter
% struct.
%
% luke.edwards@chem.ox.ac.uk

function [sys,inter]=kill_spins(sys,inter,spins)

%% sys
% isotopes
sys.isotopes(spins)=[];

% warning 
if isfield(sys,'sym_group')
    error('kill_spins: call this function before defining symmetries')
end

%% inter
% zeeman
if isfield(inter.zeeman,'scalar')
    inter.zeeman.scalar(spins)=[];
elseif isfield(inter.zeeman,'matrix')
    inter.zeeman.matrix(spins)=[];
else
    disp('kill_spins: no Zeeman terms found. Only .matrix and .scalar accepted.')
end

% coupling
if isfield(inter.coupling,'scalar')
    inter.coupling.scalar(spins,:)=[];
    inter.coupling.scalar(:,spins)=[];
elseif isfield(inter.coupling,'matrix')
    inter.coupling.matrix(spins,:)=[];
    inter.coupling.matrix(:,spins)=[];
else
    disp('kill_spins: no coupling terms found. Only .matrix and .scalar accepted.')
end

% coordinates
if isfield(inter,'coordinates')
    inter.coordinates(spins)=[];
end

end