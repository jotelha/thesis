% Forms Spinach data structures from Gaussian parsing output.
% 
% Returns the data structures appropriate for NMR or ESR/SC simulations
% depending on the spins selected for inclusion. Syntax:
%
%       [sys,inter]=g03_to_spinach(props,nuclei,references,options)
%
% where props is the output of g03_parse, nuclei is a cell array of the
% following form:
%
%       {{'H','1H'},{'N','15N'}...}
%
% giving the list of elements and isotopes that should be included. All
% other spins will be ignored.
%
% References is a vector of reference chemical shifts in ppm. It is usually
% necessary to run Gaussian on, say, TMS with the same method, and the
% resulting chemical shift goes here. Computed absolute isotropic shielding
% values for TMS in vacuum:
% 
% GIAO                  13C       1H 
% B3LYP/6-31G*        189.6621  32.1833 
% B3LYP/6-311+G(2d,p) 182.4485  31.8201 
% HF/6-31G*           199.9711  32.5957 
% HF/6-311+G(2d,p)    192.5828  32.0710
%
% CSGT                  13C       1H 
% B3LYP/6-31G*        188.5603  29.1952 
% B3LYP/6-311+G(2d,p) 182.1386  31.7788 
% HF/6-31G*           196.8670  29.5517 
% HF/6-311+G(2d,p)    192.5701  31.5989 
%
% If the isotope list described above involves an electron, for example
%
%       {{'E','E'},{'H','1H'}...}
%
% then EPR mode is assumed -- chemical shielding and scalar couplings are
% ignored, but g-tensor and hyperfine couplings are included. A spin-1/2
% electron shell is assumed in this case.
%
% The following options are currently available:
%
%     options.min_j    -  scalar coupling threshold in Hz. J-couplings
%                         smaller than this value will be ignored in the
%                         NMR mode.
%
%     options.min_hfc  -  hyperfine coupling threshold in Hz. Hyperfine
%                         tensors with a Frobenius norm smaller than 
%                         this value will be ignored in the EPR mode.
%
%     options.purge    -  if set to 'on' in EPR mode, purges the spins
%                         with hyperfine coupling below options.min_hfc
%                         from the spin system.
%
%     options.no_xyz   -  if set to 1, causes the function to 
%                         ignore the coordinate information from the
%                         Gaussian log.
%
%     options.dilute   -  all couplings between spins specified in this
%                         cell array are erased.
%
% The following parameters (if found in the log) are currently returned:
%
%     sys.isotopes           Nspins x 1 cell array of strings
%
%     inter.coordinates      Nspins x 3 dense matrix, Angstrom. Not
%                            returned if there is an electron in the
%                            isotope list (in the EPR case it is not
%                            a good idea to use the molecular 
%                            coordinates for spins).
%
%     inter.zeeman.matrix    Nspins x 1 cell array of 3x3 matrices,
%                            ppm for nuclei, g-tensor for electrons.
%                            Zero interactions have zero matrices.
%
%     inter.coupling.matrix  Nspins x Nspins cell array of 3x3 matrices,
%                            all in Hz. Zero interactions have zero 
%                            matrices.
%
%     inter.coupling.scalar  Nspins x Nspins cell array of scalar
%                            couplings, all in Hz. Zero couplings are
%                            returned as zeros.
%
% ilya.kuprov@oerc.ox.ac.uk

function [sys,inter]=g03_to_spinach(props,nuclei,references,options)

% Index the nuclei to include
sys.isotopes={}; index=[]; ref_index=[];
for n=1:length(props.symbols)
    for k=1:length(nuclei)
        if strcmp(props.symbols{n},nuclei{k}{1})
            sys.isotopes=[sys.isotopes nuclei{k}{2}];
            index=[index n]; %#ok<AGROW>
            ref_index=[ref_index k]; %#ok<AGROW>
        end
    end
end
nspins=length(index);

% Decide whether to include coordinates
if nargin<4
    include_xyz=true();
else
    if ~isfield(options,'no_xyz')
        include_xyz=true();
    elseif ~options.no_xyz
        include_xyz=true();
    else
        include_xyz=false();
    end
end

% Process the coordinates
if include_xyz
    inter.coordinates=props.std_geom(index,:);
    inter.coordinates=num2cell(inter.coordinates,2);
end

% Decide which magnetic parameters to return
switch any(strcmp([nuclei{:}],'E'))
    
    % EPR parameterization
    case 1
        
        % Add the electron as the last spin in the isotope list
        sys.isotopes=[sys.isotopes, 'E']; nspins=nspins+1;
        
        % All Zeeman tensors are zero except for the g-tensor for the electron
        inter.zeeman.matrix=mat2cell(zeros(3*nspins,3),3*ones(nspins,1),3);
        inter.zeeman.matrix{nspins}=props.g_tensor.matrix;
        
        % All couplings are zero except for the hyperfine couplings to the electron
        inter.coupling.matrix=mat2cell(zeros(3*nspins,3*nspins),3*ones(nspins,1),3*ones(nspins,1));
        for n=1:(nspins-1)
            inter.coupling.matrix{n,end}=1e6*gauss2mhz(props.hfc.full.matrix{index(n)}/2);
            inter.coupling.matrix{end,n}=1e6*gauss2mhz(props.hfc.full.matrix{index(n)}/2);
        end
        
        % Kill the small hyperfine couplings
        if (nargin==4)&&isfield(options,'min_hfc')
            for n=1:nspins
                for k=1:nspins
                    if norm(inter.coupling.matrix{n,k},'fro')<options.min_hfc
                        inter.coupling.matrix{n,k}=[];
                    end
                end
            end
        end
        
        % Kill the nuclei that are not coupled to the electron
        if (nargin==4)&&isfield(options,'purge')&&strcmp(options.purge,'on')
            killing_pattern=cellfun(@isempty,inter.coupling.matrix(:,nspins)); killing_pattern(end)=0;
            inter.coupling.matrix(:,killing_pattern)=[];
            inter.coupling.matrix(killing_pattern,:)=[];
            inter.zeeman.matrix(killing_pattern)=[];
            sys.isotopes(killing_pattern)=[];
        end
        
    % NMR parameterization    
    case 0
        
        % Reference to bare nuclei if the user did not supply any reference values.
        if (nargin<3)||(~any(references))
            references=zeros(size(nuclei));
        end
        
        % Absorb Zeeman tensors
        inter.zeeman.matrix=mat2cell(zeros(3*nspins,3),3*ones(nspins,1),3);
        for n=1:nspins
            inter.zeeman.matrix{n}=-props.cst{index(n)}+eye(3)*references(ref_index(n));
        end
        
        % Absorb scalar couplings
        inter.coupling.scalar=zeros(nspins);
        if isfield(props,'j_couplings')
            inter.coupling.scalar=props.j_couplings(index,index)/2;
        end
        
        % Kill the small scalar couplings
        if (nargin==4)&&isfield(options,'min_j')
            inter.coupling.scalar=inter.coupling.scalar.*(abs(inter.coupling.scalar)>options.min_j);
        end
        
        % Return a cell array
        inter.coupling.scalar=num2cell(inter.coupling.scalar);
        
end

% Process the dilute spins
if (nargin==4)&&isfield(options,'dilute')
    for n=1:length(sys.isotopes)
        for k=1:length(sys.isotopes)
            if isfield(inter.coupling,'matrix')
                % Destroy coupling if both spins are dilute
                if any(strcmp(sys.isotopes{n},options.dilute))&&any(strcmp(sys.isotopes{k},options.dilute))
                    inter.coupling.matrix{n,k}=[];
                end
            end
            if isfield(inter.coupling,'scalar')
                % Destroy coupling if both spins are dilute
                if any(strcmp(sys.isotopes{n},options.dilute))&&any(strcmp(sys.isotopes{k},options.dilute))
                    inter.coupling.scalar{n,k}=[];
                end
            end
        end
    end
end

end

% ACHTUNG! ALLES LOOKENSPEEPERS
% Das Computermachine ist nicht fur gefingerpoken und mittengrabben.
% Ist easy schnappen der Springenwerk, blowenfusen und poppencorken
% mit Spitzensparken. Ist nicht fur gewerken bei das Dumpkopfen.
% Das rubbernecken Sichtseeren keepen Hands in die Pockets muss,
% relaxen und watchen die Blinkenlichten.

