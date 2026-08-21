function covariance = localizationCovarianceFromTable(tableData, row)
    import membrane_tracking.focused_ism.*

    covariance = [tableData.crbVarXUm2(row), ...
        tableData.crbCovXYUm2(row); ...
        tableData.crbCovXYUm2(row), tableData.crbVarYUm2(row)];
    covariance = 0.5 * (covariance + covariance.');
end
