function activeOn = updateBlinking(activeOn, bleached, opts)
    import membrane_tracking.focused_ism.*

    totalRate = opts.blinkOffRateS + opts.blinkOnRateS;
    if totalRate == 0
        activeOn(bleached) = false;
        return;
    end
    transitionMass = 1 - exp(-totalRate * opts.dtS);
    pOff = opts.blinkOffRateS / totalRate * transitionMass;
    pOn = opts.blinkOnRateS / totalRate * transitionMass;

    wasOn = activeOn & ~bleached;
    wasOff = ~activeOn & ~bleached;

    activeOn(wasOn) = rand(sum(wasOn), 1) >= pOff;
    activeOn(wasOff) = rand(sum(wasOff), 1) < pOn;
    activeOn(bleached) = false;
end
