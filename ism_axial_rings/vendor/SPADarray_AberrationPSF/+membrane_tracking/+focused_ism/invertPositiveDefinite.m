function [matrixInverse, ok] = invertPositiveDefinite(matrix)
    import membrane_tracking.focused_ism.*

    matrix = 0.5 * (matrix + matrix.');
    [factor, flag] = chol(matrix, 'lower');
    ok = flag == 0 && all(isfinite(factor(:))) && ...
        all(diag(factor) > sqrt(eps) * max(diag(factor)));
    if ok
        matrixInverse = factor.' \ (factor \ eye(size(matrix)));
        matrixInverse = 0.5 * (matrixInverse + matrixInverse.');
        ok = all(isfinite(matrixInverse(:))) && ...
            all(diag(matrixInverse) > 0);
    else
        matrixInverse = nan(size(matrix));
    end
end
