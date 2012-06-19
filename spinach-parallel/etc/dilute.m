% A function that splits a simulation involving dilute nuclei into
% a number of simulations each involving one of the nuclei. This avoids
% unphysical magnetization transfer between the different systems and
% reduces the matrix dimensions. This function is very rough and it is
% easy to find input that breaks it, but it does well enough for now.
%
% luke.edwards@chem.ox.ac.uk

function [sys,inter]=dilute(sys_0,inter_0,dilute_spins)

% Convert coupling specifications to one specific type
if isfield(inter_0.zeeman,'scalar')
    inter_0.zeeman.matrix=cell(size(inter_0.zeeman.scalar));
    for k=1:length(inter_0.zeeman.scalar)
        if ~isempty(inter.zeeman_0.scalar{k})
            inter.zeeman_0.matrix{k}=eye(3)*inter_0.zeeman.scalar{k};
        else
            inter.zeeman_0.matrix{k}=[];
        end
    end
    inter_0.zeeman=rmfield(inter_0.zeeman,'scalar');
end
if isfield(inter_0.zeeman,'eigs')
    error('dilute: not currently implemented');
end
if isfield(inter_0.coupling,'scalar')
    inter_0.coupling.matrix=cell(size(inter_0.coupling.scalar));
    for k=1:numel(inter_0.coupling.scalar)
        if ~isempty(inter_0.coupling.scalar{k})
            inter_0.coupling.matrix{k}=eye(3)*inter_0.coupling.scalar{k};
        else
            inter_0.coupling.matrix{k}=[];
        end
    end
    inter_0.coupling=rmfield(inter_0.coupling,'scalar');
end
if isfield(inter_0.coupling,'eigs')
    error('dilute: not currently implemented');
end
if ~isfield(inter_0.zeeman,'matrix')&&~isfield(inter_0.coupling,'matrix')
    error('dilute: couplings not specified in appropriate form.')
end

% Extract positions of dilute spins.
if ~isnumeric(dilute_spins)
    dilute_spins=find(strcmp(sys_0.isotopes,dilute_spins));
end

other_spins=[];
for k=1:length(sys_0.isotopes)
    if ~any(k==dilute_spins)
        other_spins=[other_spins,k];
    end
end
other_isotopes=sys_0.isotopes(other_spins);
other_zeeman=inter_0.zeeman.matrix(other_spins);

% Symmetry
if isfield(sys_0,'sym_spins')
    for k=1:length(sys_0.sym_spins)
        for n=1:length(sys_0.sym_spins{k})
            sys_0.sym_spins{k}(n)=find(sys_0.sym_spins{k}(n)==other_spins);
        end
    end
end

% Extract couplings between non-dilute spins.
other_coupling=cell(length(other_spins)+1);
for n=1:length(other_spins)
    for m=1:length(other_spins)
        other_coupling{n,m}=inter_0.coupling.matrix{other_spins(n),other_spins(m)};
    end
end

if isfield(inter_0,'coordinates')
    coordinates=cell(1,length(other_spins)+1);
    coordinates(1:end-1)=inter_0.coordinates(other_spins);
    coord_dilute=inter_0.coordinates(dilute_spins);
end
    

% Add in couplings to dilute spins.
n=1;
for k=dilute_spins
    sys(n)=sys_0;
    sys(n).isotopes=[other_isotopes,sys_0.isotopes(k)];
    
    inter(n)=inter_0;
    
    inter(n).zeeman.matrix=[other_zeeman;inter_0.zeeman.matrix(k)];
    
    inter(n).coupling.matrix=other_coupling;
    for m=1:length(other_spins)
        inter(n).coupling.matrix{length(other_spins)+1,m}=...
            inter_0.coupling.matrix{k,other_spins(m)};
        inter(n).coupling.matrix{m,length(other_spins)+1}=...
            inter_0.coupling.matrix{other_spins(m),k};
    end
    
    if exist('coordinates','var')
        inter(n).coordinates=coordinates;
        inter(n).coordinates{end}=coord_dilute{n};
    end
    n=n+1;
end

end