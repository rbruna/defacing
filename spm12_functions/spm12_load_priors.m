function tpm = spm12_load_priors ( tpmname, tpmvols )

% Based on SPM12 functions:
% * spm12_load_priors8 by John Ashburner.


% Load the tissue probability maps for segmentation
% FORMAT tpm = spm12_load_priors8(V)
% V   - structures of image volume information (or filenames)
% tpm - a structure for tissue probabilities
%
% This function is intended to be used in conjunction with spm12_sample_priors.
% V = spm12_vol(P);
% T = spm12_load_priors(V);
% B = spm12_sample_priors(T,X,Y,Z);
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm12_load_priors8.m 5962 2014-04-17 12:47:43Z spm $


% Loads the TPM template.
tpmmri    = myft_read_mri ( tpmname );


% Selects all the volumes, if none provided.
if nargin < 2
    tpmvols   = 1: size ( tpmmri.anatomy, 4 );
end

% Gets the selected volumes.
vols      = tpmmri.anatomy ( :, :, :, tpmvols );
vols      = double ( vols );


% Stores the data and metadata.
tpm       = [];
tpm.name  = tpmname;
tpm.mri   = tpmmri;
tpm.vols  = tpmvols;
tpm.dat   = [];
tpm.M     = tpmmri.transform;

% Extracts the backgrounds (average of the first and last slice).
tpm.bg1   = squeeze ( mean ( mean ( vols ( :, :, 1,   : ) ) ) );
tpm.bg2   = squeeze ( mean ( mean ( vols ( :, :, end, : ) ) ) )';

% Defines the degree of the B-spline and the zero value.
tpm.deg   = 2;
tpm.tiny  = 1e-4;


% Intializes the B-spline cell array.
tpm.dat   = cell ( size ( tpmvols ) );

% Goes through each volume.
for vindex = 1: numel ( tpmvols )
    
    % Sets the degree.
    deg       = tpm.deg - 1;
    
    % Gets the logarithm of the volume (avoiding zeroes).
    logvol    = log ( vols ( :, :, :, vindex ) + tpm.tiny );
    
    % Gets the B-spline coefficents of the current volume.
%     tpm.dat { vindex } = spm_bsplinc ( logvol, [ deg deg deg  0 0 0 ] );
%     tpm.dat { vindex } = spm12_bsplinc ( logvol, [ deg deg deg ] );
    tpm.dat { vindex } = spm12_bsplinc ( logvol, deg );
end
