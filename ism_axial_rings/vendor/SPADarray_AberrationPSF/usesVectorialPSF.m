function tf = usesVectorialPSF(sim)
%USESVECTORIALPSF True when the simulation requests vectorial diffraction.

    tf = false;
    if isfield(sim, 'includesVectorialPolarization') && ...
            ~isempty(sim.includesVectorialPolarization)
        tf = logical(sim.includesVectorialPolarization);
    end
    if isfield(sim, 'diffractionModel') && ~isempty(sim.diffractionModel)
        model = lower(char(sim.diffractionModel));
        tf = tf || contains(model, 'vector');
        tf = tf && ~contains(model, 'scalar');
    end
end
