function [ afftrans, ll, h ] = spm12_maff ( mri, sampling, fwhm, tpm, afftrans, regtype )
% Affine registration to MNI space using mutual information
% FORMAT [M,ll,h] = spm12_maff8(P,samp,fwhm,tpm,M0,regtyp)
% P       - filename or structure handle of image
% samp    - distance between sample points (mm).  Small values are
%           better, but things run more slowly.
% fwhm    - smoothness estimate for computing a fudge factor.  Estimate
%           is a full width at half maximum of a Gaussian (in mm). 
% tpm     - data structure encoding a tissue probability map, generated
%           via spm12_load_priors8.m.
% M0      - starting estimates for the affine transform (or [] to use
%           default values).
% regtype - regularisation type
%           'mni'     - registration of European brains with MNI space
%           'eastern' - registration of East Asian brains with MNI space
%           'rigid'   - rigid(ish)-body registration
%           'subj'    - inter-subject registration
%           'none'    - no regularisation
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm12_maff8.m 6421 2015-04-23 16:54:25Z john $


% Generates a 4 voxels FWHM gaussian kernel.
% krn       = exp ( - ( -256: 256 ) .^ 2 / 8 * 1.3 )';
krn       = exp ( - ( -256: 256 ) .^ 2 / 8 )';
krn       = krn / sum ( krn );


% Gets the MRI volume and metadata.
volume    = mri.anatomy;
mridim    = mri.dim ( 1: 3 );
mritrans  = mri.transform;
vxsize    = sqrt ( sum ( mritrans ( 1: 3, 1: 3 ) .^2 ) );


% Rounds the sampling factor to the nearest voxel.
sampling  = max ( [ 1 1 1 ], round ( sampling * [ 1 1 1 ] ./ vxsize ) );

% Calculates the fudge factor for the given Gaussian filter.
dummy     = ( fwhm + mean ( vxsize ) ) / sqrt ( 8 * log (2) );
ffactor   = sqrt ( prod ( 4 * pi * ( dummy ./ vxsize ./ sampling ) .^ 2 + 1 ) );

% Initializes the matrix.
[ a, b ]  = spm12_affine_priors ( regtype );
mu        = cat ( 1, zeros ( 6, 1 ), a );
Alpha0    = blkdiag ( eye (6) * 1e-5, b ) * ffactor;


% Extracts the sampling points.
xs ( :, 1, 1 ) = 1: sampling (1): mridim (1);
ys ( 1, :, 1 ) = 1: sampling (2): mridim (2);
zs ( 1, 1, : ) = 1: sampling (3): mridim (3);

% Samples the MRI volume by the desired factor.
volume    = volume ( xs, ys, zs, : );

% Normalizes between 0 and 255 ignoring higher and lower 0.05% of voxels.
values    = sort ( volume ( isfinite ( volume (:) ) & ( volume (:) ~= 0 ) & ( volume (:) ~= 3024 ) ) );
minval    = values ( round ( numel ( values ) * 0.0005 ) );
maxval    = values ( round ( numel ( values ) * 0.9995 ) );

volume    = ( volume - minval ) / ( maxval - minval ) * 255;
volume    = uint8 ( volume );

% Masks out not valid values (NaN, Inf, 0, -3924 and -1500).
msk       = isfinite ( volume ) & volume ~= 0 & volume ~= -3924 & volume ~= -1500;
msk       = any ( msk, 4 );

% Initializes the a-priori probabilities structure.
approb    = cell ( 1, length ( zs ) );


if ~isempty ( afftrans )
    sol  = M2P ( afftrans );
else
    sol  = mu;
end

ll   = -Inf;

h1       = ones(256,numel(tpm.dat));
nsearch  = 12;
for iter=1:200

    stepsize = 1;
    for search = 1: nsearch,
        if iter>1 
            sol1    = sol - stepsize*dsol;
        else
            sol1    = sol;
        end
        penalty = 0.5*(sol1-mu)'*Alpha0*(sol1-mu);
        
        % Combines the current transformation with the native one.
        afftrans       = P2M(sol1);
        T       = tpm.M \ afftrans * mritrans;
        
        % Calculates the real-world coordinates of each voxel.
        xsr    = bsxfun ( @plus, bsxfun ( @plus, T (1,1) * xs, T (1,2) * ys ), T (1,3) * zs ) + T (1,4);
        ysr    = bsxfun ( @plus, bsxfun ( @plus, T (2,1) * xs, T (2,2) * ys ), T (2,3) * zs ) + T (2,4);
        zsr    = bsxfun ( @plus, bsxfun ( @plus, T (3,1) * xs, T (3,2) * ys ), T (3,3) * zs ) + T (3,4);
        
        % Masks out the voxels whose z dimension is bellow 1mm.
        zmsk   = msk & zsr > 1;
        
        % Applies the mask to the volume.
        mvol   = double ( volume ) + 1;
        mvol ( ~zmsk ) = NaN;
        
        % Goes through each slice.
        for i = 1: length ( zs )
            
            % Gets the current slice.
            xr    = xsr  ( :, :, i );
            yr    = ysr  ( :, :, i );
            zr    = zsr  ( :, :, i );
            
            % Uses only the masked valules.
            xr    = xr ( zmsk ( :, :, i ) );
            yr    = yr ( zmsk ( :, :, i ) );
            zr    = zr ( zmsk ( :, :, i ) );
            
            % Ignores blank slices.
            if ~numel ( xr ), continue; end
            
            % Finds the a-priori probability of each not-zero voxel.
            approb {i} = spm12_sample_priors ( tpm, xr, yr, zr );
        end
        

        ll1 = 0;
        for subit=1:32
            h  = zeros(256,numel(tpm.dat))+eps;
            if ~rem(subit,4),
                ll0 = ll1;
                ll1 = 0;
            end
            
            for i = 1: length ( zs )
                
                % Gets the voxels inclded in the mask.
                gm = mvol ( :, :, i );
                gm = gm ( ~isnan ( gm ) );
                
                % Ignores blank slices.
                if ~numel ( gm ), continue, end
                
                q  = h1 ( gm, : ) .* cat ( 2, approb {i} {:} );
                sq = sum ( q, 2 ) + eps;
                
                if ~rem(subit,4),
                    ll1 = ll1 + sum ( log ( sq ) );
                end
                for k = 1: size ( h, 2 ),
                    h(:,k) = h(:,k) + accumarray(gm,q(:,k)./sq,[256 1]);
                end
            end
            
            % Convolves H and the Gaussian kernel and normalizes.
            h1  = conv2 ( h + eps, krn, 'same' ) / sum ( h (:) );
            h1  = h1 ./ ( sum ( h1, 2 ) * sum ( h1, 1 ) );
            
            % If the variation is negligible exits.
            if ~rem ( subit, 4 ) && ( ll1 - ll0 ) / sum ( h (:) ) < 1e-5
                break
            end
        end
        
        

        ssh   = sum(h(:));
        ll1   = (sum(sum(h.*log(h1))) - penalty)/ssh/log(2);

        if iter==1, break; end               % No need for search
        if abs(ll1-ll)<1e-4, return; end     % Seems to have converged
        if (ll1<ll)
            stepsize = stepsize*0.5;        % Worse solution. Try again
            if search==nsearch, return; end; % Give up trying
        else
            break;                           % Better solution.  Carry on to next GN step.
        end
    end
    ll    = ll1;
    sol   = sol1;

    % Computes 1st and approximate second derivatives for computing Gauss-Newton update
    Alpha = zeros ( 12 );    % 2nd derivatives with respect to an affine transform
    Beta  = zeros ( 12, 1 ); % 1st derivatives with respect to an affine transform
    
    
    
    % Goes through each slice.
    for i = 1: length ( zs )
        
        % Gets the current slice.
        xr    = xsr  ( :, :, i );
        yr    = ysr  ( :, :, i );
        zr    = zsr  ( :, :, i );
        vol   = mvol ( :, :, i );
        
        % Uses only the masked valules.
        xr    = xr   ( zmsk ( :, :, i ) );
        yr    = yr   ( zmsk ( :, :, i ) );
        zr    = zr   ( zmsk ( :, :, i ) );
        vol   = vol  ( zmsk ( :, :, i ) );
        
        % Ignores blank slices.
        if ~numel ( vol ), continue; end
        
        % Gets the a priori probabilities and its derivatives.
        [ b, db1, db2, db3 ] = spm12_sample_priors ( tpm, xr, yr, zr );
        
        mi   = sum ( h1 ( vol, : ) .* cat ( 2, b   {:} ), 2 );
        dmi1 = sum ( h1 ( vol, : ) .* cat ( 2, db1 {:} ), 2 ) ./ mi;
        dmi2 = sum ( h1 ( vol, : ) .* cat ( 2, db2 {:} ), 2 ) ./ mi;
        dmi3 = sum ( h1 ( vol, : ) .* cat ( 2, db3 {:} ), 2 ) ./ mi;
        
        % Gets the coordinates of each masked voxel.
        [ x, y ] = find ( zmsk ( :, :, i ) );
        x  = xs (x);
        y  = ys (y);
        z  = zs (i);
        
        % Convert from derivatives w.r.t. displacements at each voxel to
        % derivatives w.r.t. affine basis functions (at each voxel).
        A = cat ( 2, ...
            dmi1 .* x (:), dmi2 .* x (:), dmi3 .* x (:), ...
            dmi1 .* y (:), dmi2 .* y (:), dmi3 .* y (:), ...
            dmi1 .* z (:), dmi2 .* z (:), dmi3 .* z (:), ...
            dmi1,          dmi2,          dmi3 );
        
        % Updates the derivative matrices.
        Alpha = Alpha + A'*A;
        Beta  = Beta  - sum(A,1)';
    end
    
    % Convert from derivatives w.r.t. affine matrix, to derivatives w.r.t. parameters
    R     = derivs(tpm.M,sol,mritrans);
    Alpha = R'*Alpha*R;
    Beta  = R'*Beta;

    % Gauss-Newton increment direction
    dsol  = ((Alpha+Alpha0)\(Beta+Alpha0*(sol-mu)));
end

afftrans = P2M ( sol );


%==========================================================================
% function P = M2P(M)
%==========================================================================
function P = M2P(M)
% Polar decomposition parameterisation of affine transform,
% based on matrix logs
J  = M(1:3,1:3);
V  = sqrtm(J*J');
R  = V\J;

lV = logm(V);
lR = -logm(R);
if sum(sum(imag(lR).^2))>1e-6, error('Rotations by pi are still a problem.'); end;
P       = zeros(12,1);
P(1:3)  = M(1:3,4);
P(4:6)  = lR([2 3 6]);
P(7:12) = lV([1 2 3 5 6 9]);
P       = real(P);


%==========================================================================
% function M = P2M(P)
%==========================================================================
function M = P2M(P)
% Polar decomposition parameterisation of affine transform,
% based on matrix logs

% Translations
D      = P(1:3);
D      = D(:);

% Rotation part
ind    = [2 3 6];
T      = zeros(3);
T(ind) = -P(4:6);
R      = expm(T-T');

% Symmetric part (zooms and shears)
ind    = [1 2 3 5 6 9];
T      = zeros(3);
T(ind) = P(7:12);
V      = expm(T+T'-diag(diag(T)));

M      = [V*R D ; 0 0 0 1];


%==========================================================================
% function R = derivs(MF,P,MG)
%==========================================================================
function R = derivs(MF,P,MG)
% Numerically compute derivatives of Affine transformation matrix w.r.t.
% changes in the parameters.
R  = zeros(12,12);
M0 = MF\P2M(P)*MG;
M0 = M0(1:3,:);
for i=1:12
    dp     = 0.0000001;
    P1     = P;
    P1(i)  = P1(i) + dp;
    M1     = MF\P2M(P1)*MG;
    M1     = M1(1:3,:);
    R(:,i) = (M1(:)-M0(:))/dp;
end
