% Prints various summaries.
%
% ilya.kuprov@oerc.ox.ac.uk

function summary(spin_system,topic,header)

switch topic
    
    case 'zeeman'
        report(spin_system,header);
        report(spin_system,'========================================================================================================================');
        report(spin_system,'   Spin  2S+1             Eigenvalues ( A | B | C )                    Eigenvectors ( A | B | C )           Isotropic   ');
        report(spin_system,'------------------------------------------------------------------------------------------------------------------------');
        first_line=true();
        for n=1:spin_system.comp.nspins
            if significant(spin_system.inter.zeeman.matrix{n},'tensor',0)
                if first_line
                    first_line=false();
                else
                    report(spin_system,' ');
                end
                [eigvals,dcm,iso]=tensor_analysis(spin_system,spin_system.inter.zeeman.matrix{n});
                report(spin_system,[pad(num2str(n),4) pad(spin_system.comp.isotopes{n},6) pad(num2str(spin_system.comp.mults(n)),6)...
                                    num2str(eigvals','%+0.5e   ') '   ' num2str(dcm(1,:),'%+0.3e   ') '   ' num2str(iso,'%+0.5e   ')]);
                report(spin_system,['                                                                ' num2str(dcm(2,:),'%+0.3e   ')]);
                report(spin_system,['                                                                ' num2str(dcm(3,:),'%+0.3e   ')]);
            end
        end
        report(spin_system,'========================================================================================================================');
    
    case 'coordinates'
        report(spin_system,header);
        report(spin_system,'======================================');
        report(spin_system,'N    Spin     X         Y         Z   ');
        report(spin_system,'--------------------------------------');
        for n=1:spin_system.comp.nspins
            report(spin_system,[strjust([num2str(n) blanks(3-length(num2str(n)))],'left') ' '...
                                strjust([spin_system.comp.isotopes{n} blanks(5-length(spin_system.comp.isotopes{n}))],'center') '  '...
                                num2str(spin_system.inter.coordinates{n},'%+5.3f   ') '  ' spin_system.comp.labels{n}]);
        end
        report(spin_system,'======================================');
        
    case 'distances'
        report(spin_system,header);
        report(spin_system,repmat('-',1,6*spin_system.comp.nspins+4));
        report(spin_system,['N   ' blanks(length(num2str(spin_system.comp.nspins))) num2str(1:spin_system.comp.nspins,['%d ' blanks(5-length(num2str(spin_system.comp.nspins)))])]);
        report(spin_system,repmat('-',1,6*spin_system.comp.nspins+4));
        report(spin_system,[num2str((1:spin_system.comp.nspins)','%d') blanks(spin_system.comp.nspins)'...
                            blanks(spin_system.comp.nspins)' num2str(cellfun(@norm,spin_system.inter.distmatrix),'%5.2f ')]);
        report(spin_system,repmat('-',1,6*spin_system.comp.nspins+4));
        
    case 'dipolar'
        report(spin_system,header);
        report(spin_system,'  (angstroms, degrees, angular frequencies, RDC scaling factors)  ');
        report(spin_system,'==================================================================');
        report(spin_system,'L    S      R       Theta      Phi        D/2           RDCSF     ');
        report(spin_system,'------------------------------------------------------------------');
        [from, to]=find(spin_system.inter.proxmatrix);
        for n=1:length(from)
            [phi,theta,dist]=cart2sph(spin_system.inter.distvectors{from(n),to(n)}(1),...
                                      spin_system.inter.distvectors{from(n),to(n)}(2),...
                                      spin_system.inter.distvectors{from(n),to(n)}(3));
            theta=pi/2-theta;
            A=0.5*spin_system.comp.gamma(from(n))*spin_system.comp.gamma(to(n))*1.054571628e-34*1e-7/(dist*1e-10)^3;
            ort=spin_system.inter.distvectors{from(n),to(n)}/norm(spin_system.inter.distvectors{from(n),to(n)});
            rdc_scaler=ort*spin_system.inter.order_matrix*ort';
            report(spin_system,[strjust([num2str(from(n)) blanks(3-length(num2str(from(n))))],'left') '  '...
                                strjust([num2str(to(n)) blanks(3-length(num2str(to(n))))],'left') '  '...
                                num2str(dist,'%0.3f  ') '   '...
                                blanks(7-length(num2str(180*theta/pi,'%0.2f'))) num2str(180*theta/pi,'%0.2f')  '   '...
                                blanks(7-length(num2str(180*phi/pi,'%0.2f'))) num2str(180*phi/pi,'%0.2f') '   '...
                                num2str(A,'%+0.4e') '   ' num2str(rdc_scaler,'%+0.4e')]);
        end
        report(spin_system,'==================================================================');
        report(spin_system,['create: proximity cut-off ' num2str(spin_system.tols.prox_cutoff,'%5.3f') ' Angstrom.']);
        
    case 'couplings'
        report(spin_system,header);
        report(spin_system,'========================================================================================================');
        report(spin_system,'L   S          Eigenvalues ( A | B | C )                 Eigenvectors ( A | B | C )          Isotropic  ');
        report(spin_system,'--------------------------------------------------------------------------------------------------------');
        first_line=true();
        for n=1:spin_system.comp.nspins
            for k=1:spin_system.comp.nspins
                if significant(spin_system.inter.coupling.matrix{n,k},'tensor',spin_system.tols.inter_cutoff)
                    if first_line
                        first_line=false();
                    else
                        report(spin_system,' ');
                    end
                    [eigvals,dcm,iso]=tensor_analysis(spin_system,spin_system.inter.coupling.matrix{n,k});
                    report(spin_system,[pad(num2str(n),4) pad(num2str(k),4)...
                                        num2str(eigvals','%+0.3e   ') '   ' num2str(dcm(1,:),'%+0.3e   ') '   ' num2str(iso,'%+0.3e   ')]);
                    report(spin_system,['                                                  ' num2str(dcm(2,:),'%+0.3e   ')]);
                    report(spin_system,['                                                  ' num2str(dcm(3,:),'%+0.3e   ')]);
                end
            end
        end
        report(spin_system,'========================================================================================================');
        
    case 'chemistry'
        report(spin_system,header);
        report(spin_system,'============================');
        report(spin_system,' N(from)  N(to)   Rate(Hz)  ');
        report(spin_system,'----------------------------');
        for n=1:length(spin_system.chem.rate)
            report(spin_system,[' ' strjust([num2str(spin_system.chem.from(n)) blanks(3-length(num2str(spin_system.chem.from(n))))],'left') '      '...
                                    strjust([num2str(spin_system.chem.to(n)) blanks(3-length(num2str(spin_system.chem.to(n))))],'left') '   '...
                                    num2str(spin_system.chem.rate(n),'%+0.3e')]);
        end
        report(spin_system,'============================');

    case 'rlx_rates'
        report(spin_system,header);
        report(spin_system,'========================================');
        report(spin_system,'N    Spin        R1             R2   ');
        report(spin_system,'----------------------------------------');
        for n=1:spin_system.comp.nspins
            report(spin_system,[strjust([num2str(n) blanks(3-length(num2str(n)))],'left') ' '...
                                strjust([spin_system.comp.isotopes{n} blanks(5-length(spin_system.comp.isotopes{n}))],'center') '  '...
                                num2str(spin_system.rlx.r1(n),'%+0.5e   ') '  '...
                                num2str(spin_system.rlx.r2(n),'%+0.5e   ') '  '...
                                spin_system.comp.labels{n}]);
        end
        report(spin_system,'========================================');
        
    case 'symmetry'
        report(spin_system,header);
        report(spin_system,'=====================');
        report(spin_system,' Group    Spins      ');
        report(spin_system,'---------------------');
        for n=1:length(spin_system.comp.sym.spins)
            report(spin_system,['  ' spin_system.comp.sym.group{n} '     ' num2str(spin_system.comp.sym.spins{n})]);
        end
        report(spin_system,'=====================');
        
    case 'basis'
        if spin_system.bas.nstates > spin_system.tols.basis_hush
            report(spin_system,['basis: over ' num2str(spin_system.tols.basis_hush) ' states in the basis - printing suppressed.']);
        else
            report(spin_system,'basis: final basis set summary (L,M quantum numbers in irreducible spherical tensor products).')
            report(spin_system,['N       ' blanks(length(num2str(spin_system.comp.nspins))) num2str(1:spin_system.comp.nspins,['%d ' blanks(7-length(num2str(spin_system.comp.nspins)))])]);
            for n=1:spin_system.bas.nstates
                current_line=blanks(7+8*spin_system.comp.nspins); spin_number=num2str(n);
                current_line(1:length(spin_number))=spin_number;
                for k=1:spin_system.comp.nspins
                    [L,M]=lin2lm(spin_system.bas.basis(n,k));
                    current_line(7+8*(k-1)+1)='(';
                    current_line(7+8*(k-1)+2)=num2str(L);
                    current_line(7+8*(k-1)+3)=',';
                    proj=num2str(M);
                    switch length(proj)
                        case 1
                            current_line(7+8*(k-1)+4)=proj;
                            current_line(7+8*(k-1)+5)=')';
                        case 2
                            current_line(7+8*(k-1)+4)=proj(1);
                            current_line(7+8*(k-1)+5)=proj(2);
                            current_line(7+8*(k-1)+6)=')';
                    end
                end
                report(spin_system,current_line);
            end
            report(spin_system,' ');
        end
        report(spin_system,['basis: state space dimension ' num2str(spin_system.bas.nstates) ' (' num2str(100*spin_system.bas.nstates/(prod(spin_system.comp.mults)^2)) '% of the full state space).']);
        
    otherwise
        error('summary: unknown topic.');
        
end

end

% 1: Blessed are the strong, for they shall possess the earth -- cursed are
% the weak, for they shall inherit the yoke.
%
% 2: Blessed are the powerful, for they shall be reverenced among men --
% cursed are the feeble, for they shall be blotted out.
%
% 3: Blessed are the bold, for they shall be masters of the world -- cursed
% are the humble, for they shall be trodden under hoofs.
%
% 4: Blessed are the victorious, for victory is the basis of right --
% cursed are the vanquished, for they shall be vassals forever.
%
% 5: Blessed are the iron-handed, for the unfit shall flee before them --
% cursed are the poor in spirit, for they shall be spat upon.
%
% 6: Blessed are the death-defiant, for their days shall be long in the
% lands -- cursed are the gazers toward a richer life beyond the grave, for
% they shall perish amidst plenty.
%
% 7: Blessed are the destroyers of false hope, for they are true Messiahs --
% cursed are the God-adorers, for they shall be shorn sheep.
%
% 8: Blessed are the valiant, for they shall obtain great treasure --
% cursed are the believers in good and evil, for they are frightened by
% shadows.
%
% 9: Blessed are those who believe in what is best for them, for never
% shall their minds be terrorized -- cursed are the "lambs of God", for
% they shall be bled whiter than snow.
%
% 10: Blessed is the man who has a sprinkling of enemies, for they shall
% make him a hero -- cursed is he who doeth good unto others who sneer upon
% him in return, for he shall be despised.
%
% 11: Blessed are the mighty-minded, for they shall ride the whirlwinds --
% cursed are they who teach lies for truth and truth for lies, for they are
% an abomination.
%
% 12: Thrice cursed are the weak whose insecurity makes them vile, for they
% shall serve and suffer.
%
% 13: The angel of self-deceit is camped in the souls of the "righteous" --
% the eternal flame of power through joy dwelleth within the flesh of a
% Satanist.
%
% Anton Szandor LaVey, "Satanic Bible"


 
