% Database of multiplicities and magnetogyric ratios.
%
% Uses data from the NMR Periodic Table (Rider University NMR Facility)
% http://arrhenius.rider.edu/nmr/NMR_tutor/periodic_table/nmr_pt_frameset.html
%
% matthew.krzystyniak@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function [gamma,multiplicity]=spin(name)

if regexp(name,'^E\d')
    
    % Treat high-spin electrons as a special case
    multiplicity=str2double(name(2:end));
    gamma=-1.760859770e11;
    
else
    
    % All other spins come from the database
    switch char(name)
        case 'E'
            multiplicity=2;
            gamma=-1.760859770e11;
        case '1H'
            multiplicity=2;
            gamma=26.7522205e7;
        case '2H'
            multiplicity=3;
            gamma=4.10662914e7;
        case '3H'
            multiplicity=2;
            gamma=28.5350e7;
        case '3He'
            multiplicity=2;
            gamma=-20.3801e7;
        case '6Li'
            multiplicity=3;
            gamma=3.9371e7;
        case '7Li'
            multiplicity=4;
            gamma=10.3976e7;
        case '9Be'
            multiplicity=4;
            gamma=-3.7606e7;
        case '11B'
            multiplicity=4;
            gamma=8.5847e7;
        case '13C'
            multiplicity=2;
            gamma=6.728286e7;
        case '14N'
            multiplicity=3;
            gamma=1.9337798e7;
        case '15N'
            multiplicity=2;
            gamma=-2.7126188e7;
        case '17O'
            multiplicity=6;
            gamma=-3.62808e7;
        case '19F'
            multiplicity=2;
            gamma=25.18148e7;
        case '21Ne'
            multiplicity=4;
            gamma=-2.1130e7;
        case '23Na'
            multiplicity=4;
            gamma=7.0704e7;
        case '25Mg'
            multiplicity=6;
            gamma=-1.6389e7;
        case '27Al'
            multiplicity=6;
            gamma=6.971e7;
        case '29Si'
            multiplicity=2;
            gamma=-5.314e7;
        case '31P'
            multiplicity=2;
            gamma=10.8394e7;
        case '33S'
            multiplicity=4;
            gamma=2.0557e7;
        case '35Cl'
            multiplicity=4;
            gamma=2.6242e7;
        case '37Cl'
            multiplicity=4;
            gamma=2.1844e7;
        case '39Ar'
            multiplicity=8;
            gamma=-1.78e7;
        case '39K'
            multiplicity=4;
            gamma=1.2499e7;
        case '40K'
            multiplicity=9;
            gamma=-1.75544e7;
        case '41K'
            multiplicity=4;
            gamma=0.6861e7;
        case '41Ca'
            multiplicity=8;
            gamma=-1.78e7;
        case '43Ca'
            multiplicity=8;
            gamma=-1.8028e7;
        case '45Sc'
            multiplicity=8;
            gamma=6.5088e7;
        case '47Ti'
            multiplicity=6;
            gamma=-1.5106e7;
        case '49Ti'
            multiplicity=8;
            gamma=-1.5110e7;
        case '50V'
            multiplicity=13;
            gamma=2.6721e7;
        case '51V'
            multiplicity=8;
            gamma=7.0492e7;
        case '53Cr'
            multiplicity=4;
            gamma=-1.5077e7;
        case '55Mn'
            multiplicity=6;
            gamma=6.6453e7;
        case '57Fe'
            multiplicity=2;
            gamma=0.8687e7;
        case '59Co'
            multiplicity=8;
            gamma=6.3015e7;
        case '61Ni'
            multiplicity=4;
            gamma=-2.3943e7;
        case '63Cu'
            multiplicity=4;
            gamma=7.1088e7;
        case '65Cu'
            multiplicity=4;
            gamma=7.6104e7;
        case '67Zn'
            multiplicity=6;
            gamma=1.6778e7;
        case '69Ga'
            multiplicity=4;
            gamma=6.4389e7;
        case '71Ga'
            multiplicity=4;
            gamma=8.1812e7;
        case '73Ge'
            multiplicity=10;
            gamma=-0.93660e7;
        case '75As'
            multiplicity=4;
            gamma=4.5961e7;
        case '77Se'
            multiplicity=2;
            gamma=5.1214e7;
        case '79Br'
            multiplicity=4;
            gamma=6.7256e7;
        case '83Br'
            multiplicity=10;
            gamma=-1.0331e7;
        case '85Rb'
            multiplicity=6;
            gamma=2.5923e7;
        case '87Rb'
            multiplicity=4;
            gamma=8.7851e7;
        case '87Sr'
            multiplicity=10;
            gamma=-1.1635e7;
        case '89Y'
            multiplicity=2;
            gamma=-1.3163e7;
        case '91Zr'
            multiplicity=6;
            gamma=-2.4975e7;
        case '93Nb'
            multiplicity=10;
            gamma=6.5674e7;
        case '95Mo'
            multiplicity=6;
            gamma=-1.7514e7;
        case '97Mo'
            multiplicity=6;
            gamma=-1.7884e7;
        case '99Tc'
            multiplicity=10;
            gamma=6.0503e7;
        case '99Ru'
            multiplicity=6;
            gamma=-1.2286e7;
        case '101Ru'
            multiplicity=6;
            gamma=-1.3773e7;
        case '103Rh'
            multiplicity=2;
            gamma=-0.8468e7;
        case '105Pd'
            multiplicity=6;
            gamma=-1.2305e7;
        case '107Ag'
            multiplicity=2;
            gamma=-1.0878e7;
        case '109Ag'
            multiplicity=2;
            gamma=-1.2519e7;
        case '111Cd'
            multiplicity=2;
            gamma=-5.7046e7;
        case '113Cd'
            multiplicity=2;
            gamma=-5.9609e7;
        case '113In'
            multiplicity=10;
            gamma=5.8845e7;
        case '115In'
            multiplicity=10;
            gamma=5.8971e7;
        case '115Sn'
            multiplicity=2;
            gamma=-8.8014e7;
        case '117Sn'
            multiplicity=2;
            gamma= -9.589e7;
        case '119Sn'
            multiplicity=2;
            gamma= -10.0318e7;
        case '121Sb'
            multiplicity=6;
            gamma= 6.4442e7;
        case '123Sb'
            multiplicity=8;
            gamma= 3.4904e7;
        case '123Te'
            multiplicity=2;
            gamma= -7.0576e7;
        case '125Te'
            multiplicity=2;
            gamma= -8.5087e7;
        case '127I'
            multiplicity=6;
            gamma= 5.3896e7;
        case '129Xe'
            multiplicity=2;
            gamma= -7.4521e7;
        case '131Xe'
            multiplicity=4;
            gamma= 2.2091e7;
        case '133Cs'
            multiplicity=8;
            gamma= 3.5339e7;
        case '135Ba'
            multiplicity=4;
            gamma= 2.657e7;
        case '137Ba'
            multiplicity=4;
            gamma= 2.973e7;
        case '138La'
            multiplicity=11;
            gamma= 3.5575e7;
        case '139La'
            multiplicity=8;
            gamma= 3.8085e7;
        case '139Ce'
            multiplicity=6;
            gamma= 2.906e7;
        case '141Pr'
            multiplicity=6;
            gamma= 7.765e7;
        case '143Nd'
            multiplicity=8;
            gamma= -1.474e7;
        case '145Nd'
            multiplicity=8;
            gamma= -0.913e7;
        case '147Pm'
            multiplicity=8;
            gamma= 3.613e7;
        case '147Sm'
            multiplicity=8;
            gamma= -1.1124e7;
        case '149Sm'
            multiplicity=8;
            gamma= -0.9175e7;
        case '151Eu'
            multiplicity=6;
            gamma= 6.5477e7;
        case '153Eu'
            multiplicity=6;
            gamma= 2.9371e7;
        case '155Gd'
            multiplicity=4;
            gamma= -0.8273e7;
        case '157Gd'
            multiplicity=4;
            gamma= -1.0792e7;
        case '159Tb'
            multiplicity=4;
            gamma= 6.4306e7;
        case '161Dy'
            multiplicity=6;
            gamma= -0.9206e7;
        case '163Dy'
            multiplicity=6;
            gamma= 1.275e7;
        case '156Ho'
            multiplicity=8;
            gamma=5.487e7;
        case '167Er'
            multiplicity=8;
            gamma= -0.7752e7;
        case '169Tm'
            multiplicity=2;
            gamma= -2.1376e7;
        case '171Yb'
            multiplicity=2;
            gamma= 4.7348e7;
        case '173Yb'
            multiplicity=6;
            gamma= -1.3025e7;
        case '175Lu'
            multiplicity=8;
            gamma= 3.0589e7;
        case '177Hf'
            multiplicity=8;
            gamma= 1.086e7;
        case '179Hf'
            multiplicity=10;
            gamma= -0.6821e7;
        case '181Ta'
            multiplicity=8;
            gamma= 3.2445e7;
        case '183W'
            multiplicity=2;
            gamma= 1.1283e7;
        case '185Re'
            multiplicity=6;
            gamma= 6.1057;
        case '187Re'
            multiplicity=6;
            gamma= 6.1682e7;
        case '187Os'
            multiplicity=2;
            gamma= 0.6193e7;
        case '189Os'
            multiplicity=4;
            gamma= 2.1072e7;
        case '191Ir'
            multiplicity=4;
            gamma= 0.4665e7;
        case '193Ir'
            multiplicity=4;
            gamma= 0.508e7;
        case '195Pt'
            multiplicity=2;
            gamma= 5.8383e7;
        case '197Au'
            multiplicity=4;
            gamma= 0.4692e7;
        case '199Hg'
            multiplicity=2;
            gamma= 4.8458e7;
        case '201Hg'
            multiplicity=4;
            gamma= -1.7888e7;
        case '203Tl'
            multiplicity=2;
            gamma= 15.5394e7;
        case '205Tl'
            multiplicity=2;
            gamma= 15.6922e7;
        case '207Pb'
            multiplicity=2;
            gamma= 5.6264e7;
        case '209Bi'
            multiplicity=10;
            gamma= 4.3752e7;
        case '209Po'
            multiplicity=2;
            gamma= 7.4e7;
        case '227Ac'
            multiplicity=4;
            gamma= 3.5e7;
        case '229Th'
            multiplicity=6;
            gamma= 0.4e7;
        case '231Pa'
            multiplicity=4;
            gamma= 0.4e7;
        case '235U'
            multiplicity=8;
            gamma=-0.4926e7;
        case '237Np'
            multiplicity=6;
            gamma=3.1e7;
        case '239Pu'
            multiplicity=2;
            gamma=0.972e7;
        case '243Am'
            multiplicity=6;
            gamma=1.54e7;
        case '247Cm'
            multiplicity=10;
            gamma=0.2e7;
        otherwise
            disp(['spin: ' name ' - no data available in the current NMR literature. See the spin.m function.']);
            error('spin: unknown spin.');
    end
    
end

% I mean that there is no way to disarm any man except through guilt.
% Through that which he himself has accepted as guilt. If a man has ever
% stolen a dime, you can impose on him the punishment intended for a bank
% robber and he will take it. He'll bear any form of misery, he'll feel
% that he deserves no better. If there’s not enough guilt in the world, we
% must create it. If we teach a man that it's evil to look at spring
% flowers and he believes us and then does it – we'll be able to do
% whatever we please with him. He won't defend himself. He won't feel he's
% worth it. He won't fight. But save us from the man who lives up to his
% own standards. Save us from the man of clean conscience. He's the man
% who'll beat us.
%
% Ayn Rand, "Atlas Shrugged"

