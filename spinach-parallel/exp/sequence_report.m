% Experiment summaries.
%
% matthew.krzystyniak@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function sequence_report(spin_system,parameters)

% Click forward for output
spin_system=click(spin_system,'forward');

if isfield(parameters,'nuclei_f1')
    report(spin_system,['sequence_report: F1 channel nuclei: ' parameters.nuclei_f1]);
end

if isfield(parameters,'nuclei_f2')
    report(spin_system,['sequence_report: F2 channel nuclei: ' parameters.nuclei_f2]);
end

if isfield(parameters,'nuclei')
    report(spin_system,['sequence_report: nuclei: ' parameters.nuclei]);
end

if isfield(parameters,'sweep_f1')
    report(spin_system,['sequence_report: F1 channel sweep width: ' num2str(parameters.sweep_f1)  ' Hz']);
end

if isfield(parameters,'sweep_f2')
    report(spin_system,['sequence_report: F2 channel sweep width: ' num2str(parameters.sweep_f2)  ' Hz']);
end

if isfield(parameters,'sweep')
    report(spin_system,['sequence_report: sweep width: ' num2str(parameters.sweep)  ' Hz']);
end

if isfield(parameters,'npoints_f1')
    report(spin_system,['sequence_report: number of points in F1 channel: ' num2str(parameters.npoints_f1)]);
end

if isfield(parameters,'npoints_f2')
    report(spin_system,['sequence_report: number of points in F2 channel: ' num2str(parameters.npoints_f2)]);
end

if isfield(parameters,'npoints')
    report(spin_system,['sequence_report: number of points: ' num2str(parameters.npoints)]);
end

if isfield(parameters,'zerofill_f1')
    report(spin_system,['sequence_report: F1 channel zerofilled to: ' num2str(parameters.zerofill_f1)]);
end

if isfield(parameters,'zerofill_f2')
    report(spin_system,['sequence_report: F2 channel zerofilled to: ' num2str(parameters.zerofill_f2)]);
end

if isfield(parameters,'zerofill')
    report(spin_system,['sequence_report: zerofilled to: ' num2str(parameters.zerofill)]);
end

if isfield(parameters,'offset_f1')
    report(spin_system,['sequence_report: F1 offset: ' num2str(parameters.offset_f1) ' Hz']);
end

if isfield(parameters,'offset_f2')
    report(spin_system,['sequence_report: F2 offset: ' num2str(parameters.offset_f2) ' Hz']);
end

if isfield(parameters,'offset')
    report(spin_system,['sequence_report: offset: ' num2str(parameters.offset) ' Hz']);
end

if isfield(parameters,'J')
    report(spin_system,['sequence_report: J-coupling assumed in magnetization transfer stages: ' num2str(parameters.J) ' Hz']);
end

if isfield(parameters,'J_near')
    report(spin_system,['sequence_report: near-range J-coupling assumed in magnetization transfer stages: ' num2str(parameters.J_near) ' Hz']);
end

if isfield(parameters,'J_far')
    report(spin_system,['sequence_report: far-range J-coupling assumed in magnetization transfer stages: ' num2str(parameters.J_far) ' Hz']);
end

if isfield(parameters,'tmix')
    report(spin_system,['sequence_report: mixing time: ' num2str(parameters.tmix) ' seconds.']);
end

if isfield(parameters,'duration')
    report(spin_system,['sequence_report: trace duration: ' num2str(parameters.duration) ' seconds.']);
end

if isfield(parameters,'electrons')
    report(spin_system,['sequence_report: electrons: ' parameters.electrons]);
end

end

% Virtue is rewarded in this world, remember. Natural law makes no false
% judgments. Its decisions are true and just even when dreadful. The victor
% gets the gold and the land every time. He also gets the fairest maidens,
% the glory tributes. And -- why should it be otherwise? Why should the
% delights of life go to failures and cowards? Why should the spoils of
% battle belong to the unwarlike? That would be insanity, utterly unnatural
% and immoral.
%
% Ragnar Redbeard, "Might Is Right"


