function [feld, phi, rr, zz] = FocusImage3D(rho,z,vol,flag,tsh,maxangle,clfflag,varargin)
%--------------------------------------------------------------------------
% FocusImage3D   Render or reconstruct a 3-D azimuthally-symmetric field.
%
% TWO USAGE MODES
% ───────────────
% (A) Azimuthal-harmonic input  (ndims(rho) <= 2)
%     rho  : [Nr×Nz] radial coordinate grid (or pass rho(:,1), z(1,:))
%     z    : [Nr×Nz] axial coordinate grid
%     vol  : [Nr×Nz×(2M+1)] azimuthal Fourier coefficients
%              vol(:,:,1)        = m=0 (DC)
%              vol(:,:,1+j)      = cos(j·φ) amplitude
%              vol(:,:,(end+1)/2+j) = sin(j·φ) amplitude
%     flag : (unused / ignored in this mode)
%
% (B) Pre-computed 3-D Cartesian / cylindrical grid  (ndims(rho) > 2)
%     rho  : [Ny×Nx×Nz] radial distance grid  r = sqrt(x²+y²)
%     z    : [Ny×Nx×Nz] axial coordinate grid
%     vol  : [Ny×Nx×Nz] field/intensity values on the same grid
%     flag : [Ny×Nx×Nz] azimuthal angle grid  φ = atan2(y,x)
%
%     For Cartesian data build the grids with:
%         [X,Y,Z] = meshgrid(xvec, yvec, zvec);
%         rho = sqrt(X.^2+Y.^2);  flag = atan2(Y,X);  z = Z;
%
% DISPLAY (nargout == 0)
%     Renders nested isosurfaces at levels  tsh * max(vol).
%     Uses 'material metal' + Gouraud lighting for a polished look.
%
% OPTIONAL INPUTS
%     tsh      : isosurface levels as fraction of peak
%                [default: 1./exp(1:3) ≈ [0.368, 0.135, 0.050]]
%     maxangle : azimuthal extent in radians  [default: 2π]
%     clfflag  : pass any non-empty value to suppress adding new lights
%                (useful when calling inside a subplot loop)
%     varargin : extra Name/Value pairs forwarded to  set(gca, ...)
%
% OUTPUTS (when called with output arguments)
%     feld : reconstructed 3-D field array  [Ny×Nx×Nz(×nCh)]
%     phi  : azimuthal angle grid
%     rr   : radial coordinate grid
%     zz   : axial coordinate grid
%--------------------------------------------------------------------------

if nargin<5 || isempty(tsh)
    tsh = 1./exp(1:3);
end
if nargin<6 || isempty(maxangle)
    maxangle = 2*pi;
end

if ndims(rho)>2
    % ── Mode B: 3-D grid already provided ────────────────────────────────
    rr   = rho;
    zz   = z;
    feld = vol;
    phi  = flag;
else
    % ── Mode A: azimuthal-harmonic reconstruction ────────────────────────
    maxphi = 2^7;
    [zz,rr,phi] = meshgrid(z(1,:), rho(:,1), (0:maxphi)/maxphi*maxangle);

    if nargin<4 || isempty(flag)
        % Reconstruct field from azimuthal harmonics
        if nargout>0
            feld = zeros(size(rr,1), size(rr,2), size(rr,3), size(vol,4));
        else
            feld = 0*rr;
        end
        for k = 1:size(rr,3)
            psi = squeeze(phi(1,1,k));
            tmp = vol(:,:,1,:);
            for j = 1:(size(vol,3)-1)/2
                tmp = tmp + vol(:,:,j+1,:)*cos(j*psi) ...
                          + vol(:,:,(end+1)/2+j,:)*sin(j*psi);
            end
            if nargout==0
                if isreal(vol)
                    feld(:,:,k)   = sum(tmp,4);
                else
                    feld(:,:,k)   = sum(abs(tmp).^2, 4);
                end
            else
                feld(:,:,k,:) = tmp;
            end
        end
    else
        feld = sum(vol,4);
    end
end

% ── Rendering ─────────────────────────────────────────────────────────────
if nargout==0
    tsh = sort(tsh);                             % ascending → outer first
    if length(tsh)>1
        al = fliplr(1./(1:length(tsh)));         % outer = most transparent
        co = flipud([ones(length(tsh),1), ...
                     (0:length(tsh)-1)'/length(tsh), ...
                     zeros(length(tsh),1)]);      % orange→red palette
    else
        al = 0.5;
        co = [1 0.8 0.5];
    end

    if nargin<7 || isempty(clfflag)
        % clf;                                   % (commented out: caller manages clf)
    end

    % Render isosurface layers
    %   rr.*cos(phi) = x,  rr.*sin(phi) = y,  zz = z
    for j = 1:length(tsh)
        patch(isosurface(rr.*cos(phi), rr.*sin(phi), zz, feld, ...
              tsh(j)*max(feld(:))), ...
              'FaceColor', co(j,:), 'EdgeColor', 'none', 'FaceAlpha', al(j));
    end

    axis image
    view(3);
    material metal
    if nargin<7 || isempty(clfflag)
        light('position', [0  1  0])
        light('position', [1 -1  0])
        light('position',[-1 -1  0])
    end
    lighting gouraud
    axis vis3d
    box on
    if nargin>7
        set(gca, varargin{:});
    end
    cameratoolbar

    xlabel('\itx\rm (\mum)');
    ylabel('\ity\rm (\mum)');
    zlabel('\itz\rm (\mum)')

    clear feld phi
end

% mpcolor(squeeze(rr(:,50,:).*cos(phi(:,50,:))),squeeze(rr(:,50,:).*sin(phi(:,50,:))),squeeze(feld(:,50,:)))
