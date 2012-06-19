% Propagating a certain range lower_boundary to upper_boundray of 
% orientations, using predefined spin system with prepared isotropic 
% Liouvillian from .mat-file input_file_name, writing resulting fid to
% .mat-file output_file_name
%       
% johannes.hoermann@googlemail.com

function fid=jlh_outsourced_pulse_acquire(lower_boundary, upper_boundary, input_file_name, output_file_name)

    %Loading variable dump from common input file
    tmp = load(input_file_name);
    %variables have to be defined explicitely, otherwise conflicting with
    %use of parfor:
    spin_system = tmp.spin_system;
    Iso = tmp.Iso;
    Q = tmp.Q;
    Ly = tmp.Ly;
    rho = tmp.rho;
    phi = tmp.phi;
    theta = tmp.theta;
    weight = tmp.weight;
    coil = tmp.coil;
    parameters = tmp.parameters;
    
    %preallocating fid array
    fid=zeros(parameters.npoints,1);
   
    report(spin_system, sprintf('jlh_outsourced_pulse_acquire: Starting to calculate orientation %d to %d.',lower_boundary, upper_boundary));

    
    parforTic = tic;
    parfor n=lower_boundary:upper_boundary %jlh only certain range
        %jlh - benchmarking
        singleOrientationTic = tic;
        msg1 = sprintf('jlh_outsourced_pulse_acquire: Orientation %d', n);
        %report(spin_system, msg);

        % Get the current orientation
        % jlh - the only line differing from for to for
        L=Iso+orientation(Q,[phi(n) theta(n) 0]);

        %jlh - report the size of the Louvillian
        sizeL = size(L);
        s1 = sizeL(1);
        s2 = sizeL(2);
        stot = s1*s2;
        mbytes = double(stot*8)/(1024*1024);
        msg2 =  [ sprintf('jlh_outsourced_pulse_acquire: dim(L) = %d x %d = %d ~ %.2f MB', s1,s2,stot,mbytes) ' ; class(L) = ' class(L) ];

        % Apply the pulse
        rho_pulsed=step(spin_system,Ly,rho,pi/2);

        % Run the simulation
        fid=fid+weight(n)*evolution(spin_system,L,coil,rho_pulsed,1/parameters.sweep,parameters.npoints-1,'observable');

        singleOrientationTime = toc(singleOrientationTic);
        msg3 = sprintf('jlh_outsourced_pulse_acquire: orientation computation time: %f s', singleOrientationTime);
        report(spin_system, msg1);
        report(spin_system, msg2);
        report(spin_system, msg3);
    end
    parforTime = toc(parforTic);
    timeSpentMsg = sprintf('jlh_outsourced_pulse_acquire: parfor computation time: %f s', parforTime);
    report(spin_system, timeSpentMsg);
    %dump results to output_file_name
    save(output_file_name,'fid');
end