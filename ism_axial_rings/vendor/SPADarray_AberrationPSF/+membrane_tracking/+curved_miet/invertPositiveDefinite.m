function [matrixInverse, ok] = invertPositiveDefinite(matrix)
    import membrane_tracking.curved_miet.*

    matrix = double(matrix);
    matrix = 0.5 * (matrix + matrix.');
    matrixInverse = nan(size(matrix));
    ok = all(isfinite(matrix(:))) && size(matrix,1) == size(matrix,2);
    if ~ok
        return;
    end
    [upper, flag] = chol(matrix);
    if flag ~= 0 || rcond(matrix) < 1e-12
        ok = false;
        return;
    end
    matrixInverse = upper \ (upper.' \ eye(size(matrix)));
    matrixInverse = 0.5 * (matrixInverse + matrixInverse.');
    ok = all(isfinite(matrixInverse(:)));
end
