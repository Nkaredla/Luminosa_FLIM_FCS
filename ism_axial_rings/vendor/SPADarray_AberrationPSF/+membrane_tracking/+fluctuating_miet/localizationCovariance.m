function covariance = localizationCovariance(tableData, row)
    import membrane_tracking.fluctuating_miet.*

    covariance = [tableData.crbVarXUm2(row), tableData.crbCovXYUm2(row); ...
        tableData.crbCovXYUm2(row), tableData.crbVarYUm2(row)];
    covariance = 0.5*(covariance + covariance.');
end
