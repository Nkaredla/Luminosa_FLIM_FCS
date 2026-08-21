function intensity = psfBesselAirInterface(sim, coeffs, wavelength, stageZ, emitterHeight)
%PSFBESSELAIRINTERFACE Diffraction-model dispatch for air-interface PSFs.

    if usesVectorialPSF(sim)
        intensity = vectorialPSFBesselAirInterface( ...
            sim, coeffs, wavelength, stageZ, emitterHeight);
    else
        intensity = scalarPSFBesselAirInterface( ...
            sim, coeffs, wavelength, stageZ, emitterHeight);
    end
end
