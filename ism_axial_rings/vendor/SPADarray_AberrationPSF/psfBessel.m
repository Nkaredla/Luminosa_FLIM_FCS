function psf = psfBessel(sim, coeffs, wavelength)
%PSFBESSEL Diffraction-model dispatch for homogeneous PSFs.

    if usesVectorialPSF(sim)
        psf = vectorialPSFBessel(sim, coeffs, wavelength);
    else
        psf = scalarPSFBessel(sim, coeffs, wavelength);
    end
end
