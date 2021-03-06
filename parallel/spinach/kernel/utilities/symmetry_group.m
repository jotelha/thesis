% Symmetry group database.
%
% ilya.kuprov@oerc.ox.ac.uk

function group=symmetry_group(group_name)

switch group_name
    case {'S2','Ci','C2','Cs'}
        group.name='S2 permutation group, aka Ci, aka C2, aka Cs (Abelian)';
        group.order=2;
        group.nclasses=2;
        group.class_sizes=[1 1];
        group.class{1}=[1 2];
        group.class{2}=[2 1];
        group.n_irreps=2;
        group.irrep_dims=[1 1];
        group.class_characters=[1  1;
                                1 -1];
    case {'S3','C3v','D3'}
        group.name='S3 permutation group, aka C3v, aka D3 (non-Abelian)';
        group.order=6;
        group.nclasses=3;
        group.class_sizes=[1 2 3];
        group.class{1}=[1 2 3];
        group.class{2}=[2 3 1; 3 1 2];
        group.class{3}=[1 3 2; 3 2 1; 2 1 3];
        group.n_irreps=3;
        group.irrep_dims=[1 1 2];
        group.class_characters=[1  1  1;
                                1  1 -1;
                                2 -1  0];
    case {'S4','Td'}
        group.name='S4 permutation group, aka Td (non-Abelian)';
        group.order=24;
        group.nclasses=5;
        group.class_sizes=[1 6 8 6 3];
        group.class{1}=[1 2 3 4];
        group.class{2}=[2 3 4 1; 2 4 1 3; 3 1 4 2; 3 4 2 1; 4 1 2 3; 4 3 1 2];
        group.class{3}=[1 3 4 2; 1 4 2 3; 4 2 1 3; 3 2 4 1; 2 4 3 1; 4 1 3 2; 3 1 2 4; 2 3 1 4];
        group.class{4}=[1 2 4 3; 1 3 2 4; 1 4 3 2; 2 1 3 4; 4 2 3 1; 3 2 1 4];
        group.class{5}=[2 1 4 3; 3 4 1 2; 4 3 2 1];
        group.n_irreps=5;
        group.irrep_dims=[1 1 2 3 3];
        group.class_characters=[1  1  1  1  1;...
                                1 -1  1 -1  1;...
                                2  0 -1  0  2;...
                                3  1  0 -1 -1;...
                                3 -1  0  1 -1];                     
    case 'S5'
        group.name='S5 permutation group (non-Abelian)';
        group.order=120;
        group.nclasses=7;
        group.class_sizes=[1 10 20 15 30 20 24];
        group.class{1}=[1 2 3 4 5];
        group.class{2}=[1 2 3 5 4; 1 2 4 3 5; 1 2 5 4 3; 1 3 2 4 5; 1 4 3 2 5;
                        1 5 3 4 2; 2 1 3 4 5; 3 2 1 4 5; 4 2 3 1 5; 5 2 3 4 1];
        group.class{3}=[1 2 4 5 3; 1 2 5 3 4; 1 3 4 2 5; 1 3 5 4 2; 1 4 2 3 5;
                        1 4 3 5 2; 1 5 2 4 3; 1 5 3 2 4; 2 3 1 4 5; 2 4 3 1 5;
                        2 5 3 4 1; 3 1 2 4 5; 3 2 4 1 5; 3 2 5 4 1; 4 1 3 2 5;
                        4 2 1 3 5; 4 2 3 5 1; 5 1 3 4 2; 5 2 1 4 3; 5 2 3 1 4];
        group.class{4}=[1 3 2 5 4; 1 4 5 2 3; 1 5 4 3 2; 2 1 3 5 4; 2 1 4 3 5;
                        2 1 5 4 3; 3 2 1 5 4; 3 4 1 2 5; 3 5 1 4 2; 4 2 5 1 3;
                        4 3 2 1 5; 4 5 3 1 2; 5 2 4 3 1; 5 3 2 4 1; 5 4 3 2 1];
        group.class{5}=[1 3 4 5 2; 1 3 5 2 4; 1 4 2 5 3; 1 4 5 3 2; 1 5 2 3 4;
                        1 5 4 2 3; 2 3 4 1 5; 2 3 5 4 1; 2 4 1 3 5; 2 4 3 5 1;
                        2 5 1 4 3; 2 5 3 1 4; 3 1 4 2 5; 3 1 5 4 2; 3 2 4 5 1;
                        3 2 5 1 4; 3 4 2 1 5; 3 5 2 4 1; 4 1 2 3 5; 4 1 3 5 2;
                        4 2 1 5 3; 4 2 5 3 1; 4 3 1 2 5; 4 5 3 2 1; 5 1 2 4 3; 
                        5 1 3 2 4; 5 2 1 3 4; 5 2 4 1 3; 5 3 1 4 2; 5 4 3 1 2];
        group.class{6}=[2 1 4 5 3; 2 1 5 3 4; 2 3 1 5 4; 2 4 5 1 3; 2 5 4 3 1;
                        3 1 2 5 4; 3 4 1 5 2; 3 4 5 2 1; 3 5 1 2 4; 3 5 4 1 2;
                        4 1 5 2 3; 4 3 2 5 1; 4 3 5 1 2; 4 5 1 3 2; 4 5 2 1 3; 
                        5 1 4 3 2; 5 3 2 1 4; 5 3 4 2 1; 5 4 1 2 3; 5 4 2 3 1];
        group.class{7}=[2 3 4 5 1; 2 3 5 1 4; 2 4 1 5 3; 2 4 5 3 1; 2 5 1 3 4; 
                        2 5 4 1 3; 3 1 4 5 2; 3 1 5 2 4; 3 4 2 5 1; 3 4 5 1 2;
                        3 5 2 1 4; 3 5 4 2 1; 4 1 2 5 3; 4 1 5 3 2; 4 3 1 5 2; 
                        4 3 5 2 1; 4 5 1 2 3; 4 5 2 3 1; 5 1 2 3 4; 5 1 4 2 3;
                        5 3 1 2 4; 5 3 4 1 2; 5 4 1 3 2; 5 4 2 1 3];
        group.n_irreps=7;
        group.irrep_dims=[1 1 4 4 6 5 5];
        group.class_characters=[1  1  1  1  1  1  1;
                                1 -1  1  1 -1 -1  1;
                                4  2  1  0  0 -1 -1;
                                4 -2  1  0  0  1 -1;
                                6  0  0 -2  0  0  1;
                                5  1 -1  1 -1  1  0;
                                5 -1 -1  1  1 -1  0];
                            
    case 'S6'
        group.name='S6 permutation group (non-Abelian)';
        group.order=720;
        group.nclasses=11;
        group.class_sizes=[1 15 40 45 90 120 144 15 40 90 120];
        group.class{1}=[1 2 3 4 5 6];
        group.class{2}=[1 2 3 4 6 5; 1 2 3 5 4 6; 1 2 3 6 5 4; 1 2 4 3 5 6; 1 2 5 4 3 6; 1 2 6 4 5 3; 1 3 2 4 5 6; 1 4 3 2 5 6; 1 5 3 4 2 6; 1 6 3 4 5 2;
                        2 1 3 4 5 6; 3 2 1 4 5 6; 4 2 3 1 5 6; 5 2 3 4 1 6; 6 2 3 4 5 1];
        group.class{3}=[1 2 3 5 6 4; 1 2 3 6 4 5; 1 2 4 5 3 6; 1 2 4 6 5 3; 1 2 5 3 4 6; 1 2 5 4 6 3; 1 2 6 3 5 4; 1 2 6 4 3 5; 1 3 4 2 5 6; 1 3 5 4 2 6;
                        1 3 6 4 5 2; 1 4 2 3 5 6; 1 4 3 5 2 6; 1 4 3 6 5 2; 1 5 2 4 3 6; 1 5 3 2 4 6; 1 5 3 4 6 2; 1 6 2 4 5 3; 1 6 3 2 5 4; 1 6 3 4 2 5;
                        2 3 1 4 5 6; 2 4 3 1 5 6; 2 5 3 4 1 6; 2 6 3 4 5 1; 3 1 2 4 5 6; 3 2 4 1 5 6; 3 2 5 4 1 6; 3 2 6 4 5 1; 4 1 3 2 5 6; 4 2 1 3 5 6;
                        4 2 3 5 1 6; 4 2 3 6 5 1; 5 1 3 4 2 6; 5 2 1 4 3 6; 5 2 3 1 4 6; 5 2 3 4 6 1; 6 1 3 4 5 2; 6 2 1 4 5 3; 6 2 3 1 5 4; 6 2 3 4 1 5];
        group.class{4}=[1 2 4 3 6 5; 1 2 5 6 3 4; 1 2 6 5 4 3; 1 3 2 4 6 5; 1 3 2 5 4 6; 1 3 2 6 5 4; 1 4 3 2 6 5; 1 4 5 2 3 6; 1 4 6 2 5 3; 1 5 3 6 2 4;
                        1 5 4 3 2 6; 1 5 6 4 2 3; 1 6 3 5 4 2; 1 6 4 3 5 2; 1 6 5 4 3 2; 2 1 3 4 6 5; 2 1 3 5 4 6; 2 1 3 6 5 4; 2 1 4 3 5 6; 2 1 5 4 3 6;
                        2 1 6 4 5 3; 3 2 1 4 6 5; 3 2 1 5 4 6; 3 2 1 6 5 4; 3 4 1 2 5 6; 3 5 1 4 2 6; 3 6 1 4 5 2; 4 2 3 1 6 5; 4 2 5 1 3 6; 4 2 6 1 5 3;
                        4 3 2 1 5 6; 4 5 3 1 2 6; 4 6 3 1 5 2; 5 2 3 6 1 4; 5 2 4 3 1 6; 5 2 6 4 1 3; 5 3 2 4 1 6; 5 4 3 2 1 6; 5 6 3 4 1 2; 6 2 3 5 4 1;
                        6 2 4 3 5 1; 6 2 5 4 3 1; 6 3 2 4 5 1; 6 4 3 2 5 1; 6 5 3 4 2 1];
        group.class{5}=[1 2 4 5 6 3; 1 2 4 6 3 5; 1 2 5 3 6 4; 1 2 5 6 4 3; 1 2 6 3 4 5; 1 2 6 5 3 4; 1 3 4 5 2 6; 1 3 4 6 5 2; 1 3 5 2 4 6; 1 3 5 4 6 2; 
                        1 3 6 2 5 4; 1 3 6 4 2 5; 1 4 2 5 3 6; 1 4 2 6 5 3; 1 4 3 5 6 2; 1 4 3 6 2 5; 1 4 5 3 2 6; 1 4 6 3 5 2; 1 5 2 3 4 6; 1 5 2 4 6 3; 
                        1 5 3 2 6 4; 1 5 3 6 4 2; 1 5 4 2 3 6; 1 5 6 4 3 2; 1 6 2 3 5 4; 1 6 2 4 3 5; 1 6 3 2 4 5; 1 6 3 5 2 4; 1 6 4 2 5 3; 1 6 5 4 2 3;
                        2 3 4 1 5 6; 2 3 5 4 1 6; 2 3 6 4 5 1; 2 4 1 3 5 6; 2 4 3 5 1 6; 2 4 3 6 5 1; 2 5 1 4 3 6; 2 5 3 1 4 6; 2 5 3 4 6 1; 2 6 1 4 5 3;
                        2 6 3 1 5 4; 2 6 3 4 1 5; 3 1 4 2 5 6; 3 1 5 4 2 6; 3 1 6 4 5 2; 3 2 4 5 1 6; 3 2 4 6 5 1; 3 2 5 1 4 6; 3 2 5 4 6 1; 3 2 6 1 5 4; 
                        3 2 6 4 1 5; 3 4 2 1 5 6; 3 5 2 4 1 6; 3 6 2 4 5 1; 4 1 2 3 5 6; 4 1 3 5 2 6; 4 1 3 6 5 2; 4 2 1 5 3 6; 4 2 1 6 5 3; 4 2 3 5 6 1; 
                        4 2 3 6 1 5; 4 2 5 3 1 6; 4 2 6 3 5 1; 4 3 1 2 5 6; 4 5 3 2 1 6; 4 6 3 2 5 1; 5 1 2 4 3 6; 5 1 3 2 4 6; 5 1 3 4 6 2; 5 2 1 3 4 6; 
                        5 2 1 4 6 3; 5 2 3 1 6 4; 5 2 3 6 4 1; 5 2 4 1 3 6; 5 2 6 4 3 1; 5 3 1 4 2 6; 5 4 3 1 2 6; 5 6 3 4 2 1; 6 1 2 4 5 3; 6 1 3 2 5 4; 
                        6 1 3 4 2 5; 6 2 1 3 5 4; 6 2 1 4 3 5; 6 2 3 1 4 5; 6 2 3 5 1 4; 6 2 4 1 5 3; 6 2 5 4 1 3; 6 3 1 4 5 2; 6 4 3 1 5 2; 6 5 3 4 1 2];
        group.class{6}=[1 3 2 5 6 4; 1 3 2 6 4 5; 1 3 4 2 6 5; 1 3 5 6 2 4; 1 3 6 5 4 2; 1 4 2 3 6 5; 1 4 5 2 6 3; 1 4 5 6 3 2; 1 4 6 2 3 5; 1 4 6 5 2 3;
                        1 5 2 6 3 4; 1 5 4 3 6 2; 1 5 4 6 2 3; 1 5 6 2 4 3; 1 5 6 3 2 4; 1 6 2 5 4 3; 1 6 4 3 2 5; 1 6 4 5 3 2; 1 6 5 2 3 4; 1 6 5 3 4 2;
                        2 1 3 5 6 4; 2 1 3 6 4 5; 2 1 4 5 3 6; 2 1 4 6 5 3; 2 1 5 3 4 6; 2 1 5 4 6 3; 2 1 6 3 5 4; 2 1 6 4 3 5; 2 3 1 4 6 5; 2 3 1 5 4 6;
                        2 3 1 6 5 4; 2 4 3 1 6 5; 2 4 5 1 3 6; 2 4 6 1 5 3; 2 5 3 6 1 4; 2 5 4 3 1 6; 2 5 6 4 1 3; 2 6 3 5 4 1; 2 6 4 3 5 1; 2 6 5 4 3 1;
                        3 1 2 4 6 5; 3 1 2 5 4 6; 3 1 2 6 5 4; 3 2 1 5 6 4; 3 2 1 6 4 5; 3 2 4 1 6 5; 3 2 5 6 1 4; 3 2 6 5 4 1; 3 4 1 5 2 6; 3 4 1 6 5 2;
                        3 4 5 2 1 6; 3 4 6 2 5 1; 3 5 1 2 4 6; 3 5 1 4 6 2; 3 5 4 1 2 6; 3 5 6 4 2 1; 3 6 1 2 5 4; 3 6 1 4 2 5; 3 6 4 1 5 2; 3 6 5 4 1 2;
                        4 1 3 2 6 5; 4 1 5 2 3 6; 4 1 6 2 5 3; 4 2 1 3 6 5; 4 2 5 1 6 3; 4 2 5 6 3 1; 4 2 6 1 3 5; 4 2 6 5 1 3; 4 3 2 5 1 6; 4 3 2 6 5 1;
                        4 3 5 1 2 6; 4 3 6 1 5 2; 4 5 1 3 2 6; 4 5 2 1 3 6; 4 5 3 1 6 2; 4 5 3 6 2 1; 4 6 1 3 5 2; 4 6 2 1 5 3; 4 6 3 1 2 5; 4 6 3 5 1 2;
                        5 1 3 6 2 4; 5 1 4 3 2 6; 5 1 6 4 2 3; 5 2 1 6 3 4; 5 2 4 3 6 1; 5 2 4 6 1 3; 5 2 6 1 4 3; 5 2 6 3 1 4; 5 3 2 1 4 6; 5 3 2 4 6 1; 
                        5 3 4 2 1 6; 5 3 6 4 1 2; 5 4 1 2 3 6; 5 4 2 3 1 6; 5 4 3 2 6 1; 5 4 3 6 1 2; 5 6 1 4 3 2; 5 6 2 4 1 3; 5 6 3 1 4 2; 5 6 3 2 1 4; 
                        6 1 3 5 4 2; 6 1 4 3 5 2; 6 1 5 4 3 2; 6 2 1 5 4 3; 6 2 4 3 1 5; 6 2 4 5 3 1; 6 2 5 1 3 4; 6 2 5 3 4 1; 6 3 2 1 5 4; 6 3 2 4 1 5; 
                        6 3 4 2 5 1; 6 3 5 4 2 1; 6 4 1 2 5 3; 6 4 2 3 5 1; 6 4 3 2 1 5; 6 4 3 5 2 1; 6 5 1 4 2 3; 6 5 2 4 3 1; 6 5 3 1 2 4; 6 5 3 2 4 1];
        group.class{7}=[1 3 4 5 6 2; 1 3 4 6 2 5; 1 3 5 2 6 4; 1 3 5 6 4 2; 1 3 6 2 4 5; 1 3 6 5 2 4; 1 4 2 5 6 3; 1 4 2 6 3 5; 1 4 5 3 6 2; 1 4 5 6 2 3;
                        1 4 6 3 2 5; 1 4 6 5 3 2; 1 5 2 3 6 4; 1 5 2 6 4 3; 1 5 4 2 6 3; 1 5 4 6 3 2; 1 5 6 2 3 4; 1 5 6 3 4 2; 1 6 2 3 4 5; 1 6 2 5 3 4; 
                        1 6 4 2 3 5; 1 6 4 5 2 3; 1 6 5 2 4 3; 1 6 5 3 2 4; 2 3 4 5 1 6; 2 3 4 6 5 1; 2 3 5 1 4 6; 2 3 5 4 6 1; 2 3 6 1 5 4; 2 3 6 4 1 5; 
                        2 4 1 5 3 6; 2 4 1 6 5 3; 2 4 3 5 6 1; 2 4 3 6 1 5; 2 4 5 3 1 6; 2 4 6 3 5 1; 2 5 1 3 4 6; 2 5 1 4 6 3; 2 5 3 1 6 4; 2 5 3 6 4 1; 
                        2 5 4 1 3 6; 2 5 6 4 3 1; 2 6 1 3 5 4; 2 6 1 4 3 5; 2 6 3 1 4 5; 2 6 3 5 1 4; 2 6 4 1 5 3; 2 6 5 4 1 3; 3 1 4 5 2 6; 3 1 4 6 5 2; 
                        3 1 5 2 4 6; 3 1 5 4 6 2; 3 1 6 2 5 4; 3 1 6 4 2 5; 3 2 4 5 6 1; 3 2 4 6 1 5; 3 2 5 1 6 4; 3 2 5 6 4 1; 3 2 6 1 4 5; 3 2 6 5 1 4; 
                        3 4 2 5 1 6; 3 4 2 6 5 1; 3 4 5 1 2 6; 3 4 6 1 5 2; 3 5 2 1 4 6; 3 5 2 4 6 1; 3 5 4 2 1 6; 3 5 6 4 1 2; 3 6 2 1 5 4; 3 6 2 4 1 5;
                        3 6 4 2 5 1; 3 6 5 4 2 1; 4 1 2 5 3 6; 4 1 2 6 5 3; 4 1 3 5 6 2; 4 1 3 6 2 5; 4 1 5 3 2 6; 4 1 6 3 5 2; 4 2 1 5 6 3; 4 2 1 6 3 5;
                        4 2 5 3 6 1; 4 2 5 6 1 3; 4 2 6 3 1 5; 4 2 6 5 3 1; 4 3 1 5 2 6; 4 3 1 6 5 2; 4 3 5 2 1 6; 4 3 6 2 5 1; 4 5 1 2 3 6; 4 5 2 3 1 6;
                        4 5 3 2 6 1; 4 5 3 6 1 2; 4 6 1 2 5 3; 4 6 2 3 5 1; 4 6 3 2 1 5; 4 6 3 5 2 1; 5 1 2 3 4 6; 5 1 2 4 6 3; 5 1 3 2 6 4; 5 1 3 6 4 2;
                        5 1 4 2 3 6; 5 1 6 4 3 2; 5 2 1 3 6 4; 5 2 1 6 4 3; 5 2 4 1 6 3; 5 2 4 6 3 1; 5 2 6 1 3 4; 5 2 6 3 4 1; 5 3 1 2 4 6; 5 3 1 4 6 2;
                        5 3 4 1 2 6; 5 3 6 4 2 1; 5 4 1 3 2 6; 5 4 2 1 3 6; 5 4 3 1 6 2; 5 4 3 6 2 1; 5 6 1 4 2 3; 5 6 2 4 3 1; 5 6 3 1 2 4; 5 6 3 2 4 1; 
                        6 1 2 3 5 4; 6 1 2 4 3 5; 6 1 3 2 4 5; 6 1 3 5 2 4; 6 1 4 2 5 3; 6 1 5 4 2 3; 6 2 1 3 4 5; 6 2 1 5 3 4; 6 2 4 1 3 5; 6 2 4 5 1 3; 
                        6 2 5 1 4 3; 6 2 5 3 1 4; 6 3 1 2 5 4; 6 3 1 4 2 5; 6 3 4 1 5 2; 6 3 5 4 1 2; 6 4 1 3 5 2; 6 4 2 1 5 3; 6 4 3 1 2 5; 6 4 3 5 1 2;
                        6 5 1 4 3 2; 6 5 2 4 1 3; 6 5 3 1 4 2; 6 5 3 2 1 4];
        group.class{8}=[2 1 4 3 6 5; 2 1 5 6 3 4; 2 1 6 5 4 3; 3 4 1 2 6 5; 3 5 1 6 2 4; 3 6 1 5 4 2; 4 3 2 1 6 5; 4 5 6 1 2 3; 4 6 5 1 3 2; 5 3 2 6 1 4;
                        5 4 6 2 1 3; 5 6 4 3 1 2; 6 3 2 5 4 1; 6 4 5 2 3 1; 6 5 4 3 2 1];
        group.class{9}=[2 3 1 5 6 4; 2 3 1 6 4 5; 2 4 5 1 6 3; 2 4 6 1 3 5; 2 5 4 6 1 3; 2 5 6 3 1 4; 2 6 4 5 3 1; 2 6 5 3 4 1; 3 1 2 5 6 4; 3 1 2 6 4 5;
                        3 4 5 6 1 2; 3 4 6 5 2 1; 3 5 4 1 6 2; 3 5 6 2 4 1; 3 6 4 1 2 5; 3 6 5 2 1 4; 4 1 5 2 6 3; 4 1 6 2 3 5; 4 3 5 6 2 1; 4 3 6 5 1 2;
                        4 5 1 3 6 2; 4 5 2 6 3 1; 4 6 1 3 2 5; 4 6 2 5 1 3; 5 1 4 6 2 3; 5 1 6 3 2 4; 5 3 4 2 6 1; 5 3 6 1 4 2; 5 4 1 6 3 2; 5 4 2 3 6 1;
                        5 6 1 2 3 4; 5 6 2 1 4 3; 6 1 4 5 3 2; 6 1 5 3 4 2; 6 3 4 2 1 5; 6 3 5 1 2 4; 6 4 1 5 2 3; 6 4 2 3 1 5; 6 5 1 2 4 3; 6 5 2 1 3 4];
        group.class{10}=[2 1 4 5 6 3; 2 1 4 6 3 5; 2 1 5 3 6 4; 2 1 5 6 4 3; 2 1 6 3 4 5; 2 1 6 5 3 4; 2 3 4 1 6 5; 2 3 5 6 1 4; 2 3 6 5 4 1; 2 4 1 3 6 5;
                         2 4 5 6 3 1; 2 4 6 5 1 3; 2 5 1 6 3 4; 2 5 4 3 6 1; 2 5 6 1 4 3; 2 6 1 5 4 3; 2 6 4 3 1 5; 2 6 5 1 3 4; 3 1 4 2 6 5; 3 1 5 6 2 4;
                         3 1 6 5 4 2; 3 4 1 5 6 2; 3 4 1 6 2 5; 3 4 2 1 6 5; 3 4 5 2 6 1; 3 4 6 2 1 5; 3 5 1 2 6 4; 3 5 1 6 4 2; 3 5 2 6 1 4; 3 5 4 6 2 1;
                         3 5 6 1 2 4; 3 6 1 2 4 5; 3 6 1 5 2 4; 3 6 2 5 4 1; 3 6 4 5 1 2; 3 6 5 1 4 2; 4 1 2 3 6 5; 4 1 5 6 3 2; 4 1 6 5 2 3; 4 3 1 2 6 5;
                         4 3 2 5 6 1; 4 3 2 6 1 5; 4 3 5 1 6 2; 4 3 6 1 2 5; 4 5 1 6 2 3; 4 5 2 1 6 3; 4 5 6 1 3 2; 4 5 6 2 1 3; 4 5 6 3 2 1; 4 6 1 5 3 2;
                         4 6 2 1 3 5; 4 6 5 1 2 3; 4 6 5 2 3 1; 4 6 5 3 1 2; 5 1 2 6 3 4; 5 1 4 3 6 2; 5 1 6 2 4 3; 5 3 1 6 2 4; 5 3 2 1 6 4; 5 3 2 6 4 1;
                         5 3 4 6 1 2; 5 3 6 2 1 4; 5 4 1 2 6 3; 5 4 2 6 1 3; 5 4 6 1 2 3; 5 4 6 2 3 1; 5 4 6 3 1 2; 5 6 1 3 4 2; 5 6 2 3 1 4; 5 6 4 1 3 2;
                         5 6 4 2 1 3; 5 6 4 3 2 1; 6 1 2 5 4 3; 6 1 4 3 2 5; 6 1 5 2 3 4; 6 3 1 5 4 2; 6 3 2 1 4 5; 6 3 2 5 1 4; 6 3 4 5 2 1; 6 3 5 2 4 1;
                         6 4 1 2 3 5; 6 4 2 5 3 1; 6 4 5 1 3 2; 6 4 5 2 1 3; 6 4 5 3 2 1; 6 5 1 3 2 4; 6 5 2 3 4 1; 6 5 4 1 2 3; 6 5 4 2 3 1; 6 5 4 3 1 2];
        group.class{11}=[2 3 4 5 6 1; 2 3 4 6 1 5; 2 3 5 1 6 4; 2 3 5 6 4 1; 2 3 6 1 4 5; 2 3 6 5 1 4; 2 4 1 5 6 3; 2 4 1 6 3 5; 2 4 5 3 6 1; 2 4 5 6 1 3;
                         2 4 6 3 1 5; 2 4 6 5 3 1; 2 5 1 3 6 4; 2 5 1 6 4 3; 2 5 4 1 6 3; 2 5 4 6 3 1; 2 5 6 1 3 4; 2 5 6 3 4 1; 2 6 1 3 4 5; 2 6 1 5 3 4;
                         2 6 4 1 3 5; 2 6 4 5 1 3; 2 6 5 1 4 3; 2 6 5 3 1 4; 3 1 4 5 6 2; 3 1 4 6 2 5; 3 1 5 2 6 4; 3 1 5 6 4 2; 3 1 6 2 4 5; 3 1 6 5 2 4;
                         3 4 2 5 6 1; 3 4 2 6 1 5; 3 4 5 1 6 2; 3 4 5 6 2 1; 3 4 6 1 2 5; 3 4 6 5 1 2; 3 5 2 1 6 4; 3 5 2 6 4 1; 3 5 4 2 6 1; 3 5 4 6 1 2;
                         3 5 6 1 4 2; 3 5 6 2 1 4; 3 6 2 1 4 5; 3 6 2 5 1 4; 3 6 4 2 1 5; 3 6 4 5 2 1; 3 6 5 1 2 4; 3 6 5 2 4 1; 4 1 2 5 6 3; 4 1 2 6 3 5;
                         4 1 5 3 6 2; 4 1 5 6 2 3; 4 1 6 3 2 5; 4 1 6 5 3 2; 4 3 1 5 6 2; 4 3 1 6 2 5; 4 3 5 2 6 1; 4 3 5 6 1 2; 4 3 6 2 1 5; 4 3 6 5 2 1;
                         4 5 1 2 6 3; 4 5 1 6 3 2; 4 5 2 3 6 1; 4 5 2 6 1 3; 4 5 6 2 3 1; 4 5 6 3 1 2; 4 6 1 2 3 5; 4 6 1 5 2 3; 4 6 2 3 1 5; 4 6 2 5 3 1;
                         4 6 5 2 1 3; 4 6 5 3 2 1; 5 1 2 3 6 4; 5 1 2 6 4 3; 5 1 4 2 6 3; 5 1 4 6 3 2; 5 1 6 2 3 4; 5 1 6 3 4 2; 5 3 1 2 6 4; 5 3 1 6 4 2;
                         5 3 4 1 6 2; 5 3 4 6 2 1; 5 3 6 1 2 4; 5 3 6 2 4 1; 5 4 1 3 6 2; 5 4 1 6 2 3; 5 4 2 1 6 3; 5 4 2 6 3 1; 5 4 6 1 3 2; 5 4 6 3 2 1;
                         5 6 1 2 4 3; 5 6 1 3 2 4; 5 6 2 1 3 4; 5 6 2 3 4 1; 5 6 4 1 2 3; 5 6 4 2 3 1; 6 1 2 3 4 5; 6 1 2 5 3 4; 6 1 4 2 3 5; 6 1 4 5 2 3;
                         6 1 5 2 4 3; 6 1 5 3 2 4; 6 3 1 2 4 5; 6 3 1 5 2 4; 6 3 4 1 2 5; 6 3 4 5 1 2; 6 3 5 1 4 2; 6 3 5 2 1 4; 6 4 1 3 2 5; 6 4 1 5 3 2;
                         6 4 2 1 3 5; 6 4 2 5 1 3; 6 4 5 1 2 3; 6 4 5 3 1 2; 6 5 1 2 3 4; 6 5 1 3 4 2; 6 5 2 1 4 3; 6 5 2 3 1 4; 6 5 4 1 3 2; 6 5 4 2 1 3];
        group.n_irreps=11;
        group.irrep_dims=[1 1 5 5 10 10 9 9 5 5 16];
        group.class_characters=[1   1  1  1  1  1  1  1  1  1  1;
                                1  -1  1  1 -1 -1  1 -1  1  1 -1;
                                5   3  2  1  1  0  0 -1 -1 -1 -1;
                                5  -3  2  1 -1  0  0  1 -1 -1  1;
                                10  2  1 -2  0 -1  0 -2  1  0  1;
                                10 -2  1 -2  0  1  0  2  1  0 -1;
                                9   3  0  1 -1  0 -1  3  0  1  0;
                                9  -3  0  1  1  0 -1 -3  0  1  0;
                                5   1 -1  1 -1  1  0 -3  2 -1  0;
                                5  -1 -1  1  1 -1  0  3  2 -1  0;
                                16  0 -2  0  0  0  1  0 -2  0  0];
                                            
      case 'D2h'
        group.name='D2h point group (Abelian)';
        group.order=8;
        group.nclasses=8;
        group.class_sizes=[1 1 1 1 1 1 1 1];
        group.class{1}=[1 2 3 4]; % E
        group.class{2}=[4 3 2 1]; % C2z
        group.class{3}=[2 1 4 3]; % C2y
        group.class{4}=[3 4 1 2]; % C2x
        group.class{5}=[4 3 2 1]; % i
        group.class{6}=[1 2 3 4]; % sigma_xy
        group.class{7}=[3 4 1 2]; % sigma_xz
        group.class{8}=[2 1 4 3]; % sigma_yz
        group.n_irreps=8;
        group.irrep_dims=[1 1 1 1 1 1 1 1];
        group.class_characters=[1  1  1  1  1  1  1  1;
                                1  1 -1 -1  1  1 -1 -1;
                                1 -1  1 -1  1 -1  1 -1;
                                1 -1 -1  1  1 -1 -1  1;
                                1  1  1  1 -1 -1 -1 -1;
                                1  1 -1 -1 -1 -1  1  1;
                                1 -1  1 -1 -1  1 -1  1;
                                1 -1 -1  1 -1  1  1 -1];
    otherwise
        error(['group: symmetry group ' group_name ' not implemented.']);
end

% Get the list of group elements
group.elements=vertcat(group.class{:});

% Transform the class-wise character list into an element-wise list
group.characters=zeros(length(group.class_sizes),sum(group.class_sizes));
for n=1:length(group.class_sizes)
    for k=1:group.class_sizes(n)
        group.characters(:,sum(group.class_sizes(1:(n-1)))+k)=group.class_characters(:,n);
    end
end

end

% It's so wonderful to see a great, new, crucial idea which is not mine!
%
% Ayn Rand, "Atlas Shrugged"

