function checkDetectorShiftSign(sim)
    stack = normalizedStack(sim, struct());

    xc = zeros(size(sim.detXY,1),1);
    yc = zeros(size(sim.detXY,1),1);

    [X, Y] = meshgrid(sim.x, sim.y);

    for k = 1:size(stack,3)
        img = stack(:,:,k);
        s = sum(img(:));
        xc(k) = sum(X(:).*img(:)) / s;
        yc(k) = sum(Y(:).*img(:)) / s;
    end

    figure;
    quiver(sim.detXY(:,1), sim.detXY(:,2), xc, yc, 0);
    axis equal;
    grid on;
    xlabel('detector x');
    ylabel('detector y');
    title('Channel-image centroids vs detector offsets');
end