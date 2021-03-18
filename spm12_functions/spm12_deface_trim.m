function mri = spm12_deface_trim ( mri )

% The MRI can be a 4D matrix with several channels.
% All 3D matrices must share dimensions and transformations.

% fprintf ( 1, '  Loading the TPMs as a-priori structures.\n' );

% Loads the tissue probability maps.
tpm          = spm12_load_priors ( 'NYH_TPM.nii.gz' );
tpm_face     = spm12_load_priors ( 'NYH_TPM_face.nii.gz' );


% fprintf ( 1, '  Linearly transforming the MRI to the TPM space.\n' );

% Makes sure that the MRI is in double precission.
mri.anatomy  = double ( mri.anatomy );

% Performs a two-step affine registration from native to TPM space.
nat2tpm      = eye (4);
nat2tpm      = spm12_maff ( mri_mni, 3, 16, tpm, nat2tpm, 'mni' );
nat2tpm      = spm12_maff ( mri_mni, 3,  0, tpm, nat2tpm, 'mni' );

% Gets the linear transformation MRI voxels to TPM voxels.
vox2tpmvox   = tpm.M \ nat2tpm * mri_mni.transform;


% fprintf ( 1, '  Non-linearly transforming the MRI to the TPM space.\n' );
    
% Finds a non-linear transformation to/from native from/to TPM space.
cfg          = [];
cfg.fwhm     = 0;
cfg.samp     = 3;
cfg.biasreg  = 1e-3;
cfg.biasfwhm = 60;
cfg.lkp      = [ 1 2 3 3 4 4 4 5 5 5 5 6 6 ];
cfg.reg      = [ 0 1e-3 .5 .05 .2 ];

preprocinfo  = spm12_preproc ( mri_mni, tpm, nat2tpm, cfg );

% Gets the non-linear warping and the upsampling transformation.
dswarp       = preprocinfo.Twarp;
sam2ori      = inv ( preprocinfo.MT );


% fprintf ( 1, '  Transforming the MNI face mask to subject space.\n' );

% Gets the non-linear correspondence of each voxel.
[ x, y, z ]  = spm12_warpaffine ( mri.dim, vox2tpmvox, dswarp, sam2ori );

% Calculates the a-priori probability for each voxel.
tpmpriors    = spm12_sample_priors ( tpm_face, x, y, z );

% Generates the face mask.
facemask     = tpmpriors {1} > 0.1;


% fprintf ( 1, '  Generating an anonymized (defaced) version of the data.\n' );

% Trims the face.
mri.anatomy ( facemask ) = 0;
