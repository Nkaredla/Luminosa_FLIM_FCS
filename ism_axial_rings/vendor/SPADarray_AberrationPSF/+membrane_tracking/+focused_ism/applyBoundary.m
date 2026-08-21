function pos = applyBoundary(pos, boxSizeUm, boundaryCondition)
    import membrane_tracking.focused_ism.*

    for dim = 1:2
        L = boxSizeUm(dim);
        lo = -L/2;
        period = 2 * L;

        switch boundaryCondition
            case 'reflecting'
                wrapped = mod(pos(:,dim) - lo, period);
                overHalf = wrapped > L;
                wrapped(overHalf) = period - wrapped(overHalf);
                pos(:,dim) = lo + wrapped;

            case 'periodic'
                pos(:,dim) = lo + mod(pos(:,dim) - lo, L);
        end
    end
end
