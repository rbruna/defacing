function [ tx, ty, tz ] = spm12_warpaffine ( mridim, affine, dswarp, sam2ori )

% Defines the voxel coordinates.
x = permute ( 1: mridim (1), [ 2 1 3 ] );
y = permute ( 1: mridim (2), [ 1 2 3 ] );
z = permute ( 1: mridim (3), [ 1 3 2 ] );

% Extends the coordinates to a 3D grid.
[ x, y, z ] = ndgrid ( x, y, z );


% Applies the warping, if provided.
if nargin > 2
    
    % Sanitizes the warping matrix.
    dswarp  = double ( dswarp );
    
    % Transforms the voxels coordinates to sampled space.
    if nargin > 3
        xs = x * sam2ori ( 1, 1 ) + sam2ori ( 1, 4 );
        ys = y * sam2ori ( 2, 2 ) + sam2ori ( 2, 4 );
        zs = z * sam2ori ( 3, 3 ) + sam2ori ( 3, 4 );
    else
        xs = x;
        ys = y;
        zs = z;
    end
    
    % Upsamples the warping matrix using cubic splines.
    BSwarpX = spm12_bsplinc ( dswarp ( :, :, :, 1 ), 3 );
    warpX   = spm12_bsplins ( BSwarpX, xs, ys, zs, 3 );
    BSwarpY = spm12_bsplinc ( dswarp ( :, :, :, 2 ), 3 );
    warpY   = spm12_bsplins ( BSwarpY, xs, ys, zs, 3 );
    BSwarpZ = spm12_bsplinc ( dswarp ( :, :, :, 3 ), 3 );
    warpZ   = spm12_bsplins ( BSwarpZ, xs, ys, zs, 3 );
    
    % Applies the nonlinear warp.
    x = x + warpX;
    y = y + warpY;
    z = z + warpZ;
end


% Applies the affine transformation matrix.
tx = affine ( 1, 1 ) * x + affine ( 1, 2 ) * y + affine ( 1, 3 ) * z + affine ( 1,4 );
ty = affine ( 2, 1 ) * x + affine ( 2, 2 ) * y + affine ( 2, 3 ) * z + affine ( 2,4 );
tz = affine ( 3, 1 ) * x + affine ( 3, 2 ) * y + affine ( 3, 3 ) * z + affine ( 3,4 );

% If only one output argument, concatenates the three 3D matrices.
if nargout == 1
    tx = cat ( 4, tx, ty, tz );
end
