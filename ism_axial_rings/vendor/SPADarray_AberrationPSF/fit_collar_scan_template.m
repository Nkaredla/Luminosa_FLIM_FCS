%% fit_collar_scan_template.m
% Template workflow for calibrating the equivalent Zernike trajectory of a
% water-immersion objective correction collar from detector-resolved 5x5 bead scans.
%
% This file assumes that each bead stack at each collar setting has already
% been reduced to a fitted coefficient vector a_i by your single-stack
% inverse solver. Replace the synthetic_example() block with real fitted
% coefficients from the raw-25-image bead fit.
%
% Coefficient order:
%   [defocus, astig_x, astig_y, coma_x, coma_y, spherical]

clear; clc; close all;

[collar_mm, coeffs] = synthetic_example();
model = fit_quadratic_collar_model(collar_mm, coeffs, 1e-6);

disp('Reference collar c0 ='); disp(model.c0);
disp('a0 ='); disp(model.a0);
disp('g1 ='); disp(model.g1);
disp('g2 ='); disp(model.g2);

plot_calibration(collar_mm, coeffs, model);

%% -------- local functions --------

function model = fit_quadratic_collar_model(collar_mm, coeffs, ridge)
    if nargin < 3, ridge = 1e-6; end
    collar_mm = collar_mm(:);
    c0 = median(collar_mm);
    u = collar_mm - c0;
    X = [ones(size(u)), u, u.^2];
    beta = (X.'*X + ridge*eye(3)) \ (X.'*coeffs);
    model.c0 = c0;
    model.a0 = beta(1,:);
    model.g1 = beta(2,:);
    model.g2 = beta(3,:);
end

function pred = predict_model(model, collar_mm)
    u = collar_mm(:) - model.c0;
    pred = model.a0 + u.*model.g1 + (u.^2).*model.g2;
end

function plot_calibration(collar_mm, coeffs, model)
    names = {'defocus','astig_x','astig_y','coma_x','coma_y','spherical'};
    dense = linspace(min(collar_mm), max(collar_mm), 200).';
    pred = predict_model(model, dense);

    figure('Color','w');
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    for j = 1:size(coeffs,2)
        nexttile;
        plot(collar_mm, coeffs(:,j), 'o', 'LineWidth', 1.2); hold on;
        plot(dense, pred(:,j), '-', 'LineWidth', 1.4);
        grid on;
        xlabel('Collar setting (mm)');
        ylabel('Coeff. (waves RMS)');
        title(names{j});
    end
end

function [collar_mm, coeffs] = synthetic_example()
    collar_mm = [0.13; 0.15; 0.17; 0.19; 0.21];
    u = collar_mm - 0.17;

    coeffs = [...
        -0.12*u + 0.50*u.^2, ...   % defocus
        0.00*u + 0.06*u.^2, ...    % astig_x
        0.03*u - 0.03*u.^2, ...    % astig_y
        0.02*u, ...                % coma_x
        -0.01*u, ...               % coma_y
        0.35*u + 1.10*u.^2];       % spherical

    rng(4);
    coeffs = coeffs + 0.01*randn(size(coeffs));
end
