clc
clear
close all

% Defines the files locations.
config.path.mri  = '../../data/anatomy/T1/';
config.path.def  = '../../data/anatomy/T1-def/';
config.path.trim = '../../data/anatomy/T1-deftrim/';
config.path.patt = '*.nii.gz';

% Compress (or not) the output file.
config.gzip      = true;


% Adds the 'functions' folder to the path.
addpath ( sprintf ( '%s/spm12_functions/', pwd ) );
addpath ( sprintf ( '%s/functions/', pwd ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path

% Adds the FT toolboxes that will be required.
ft_hastoolbox ( 'freesurfer', 1, 1 );


% Loads the tissue probability maps.
tpm           = spm12_load_priors ( 'NYH_TPM.nii.gz' );
tpm_face      = spm12_load_priors ( 'NYH_TPM_face.nii.gz' );


% Gets the list of input MRIs.
files         = dir ( sprintf ( '%s%s', config.path.mri, config.path.patt ) );

% Goes through all the files.
for findex = 1: numel ( files )
    
    % Gets the base file name.
    mrifile       = files ( findex ).name;
    basename      = regexprep ( files ( findex ).name, '\.nii(.gz)?$', '' );
    
    fprintf ( 1, 'Working with file %s.\n', mrifile );
    
    % Reads the MRI file.
    mri            = myft_read_mri ( sprintf ( '%s%s', config.path.mri, mrifile ) );
    
    % Tries to load the AC-PC landmarks.
    if exist ( sprintf ( '%s%s.mat', config.path.mri, basename ), 'file' )
        fprintf ( 1, '  Transforming the MRI to AC-PC coordinates.\n' );
        
        % Loads the landmark definition.
        landmark      = load ( sprintf ( '%s%s', config.path.mri, basename ) );
        
        % Transforms the MRI to AC-PC coordinates.
        cfg           = [];
        cfg.method    = 'fiducial';
        cfg.coordsys  = 'acpc';
        cfg.fiducial  = landmark;
        
%         mri_acpc       = ft_volumerealign ( cfg, mri );
        [ ~, mri_acpc ]  = evalc ( 'ft_volumerealign ( cfg, mri );' );
        
        % Gets the transformation matrix from native to AC-PC.
        nat2acpc      = mri_acpc.transform / mri.transform;
        
        % Uses this transformation as starting point.
        nat2tpm       = nat2acpc;
        
    % Otherwise tries to load the transformation to Talairach space.
    elseif exist ( sprintf ( '%s%s_nat2tal.xfm', config.path.mri, basename ), 'file' )
        
        % Loads the transformation to Talairach space.
        nat2tal       = myfs_read_xfm ( sprintf ( '%s%s_nat2tal.xfm', config.path.mri, basename ) );
        
        % Uses this transformation as starting point.
        nat2tpm       = nat2tal;
        
    % As a last resort tries to perform a rigid body transformation.
    % Courtesy of Rik Henson.
    else
        
        fprintf ( 1, '  Doing 6-param rigid-body transform to TPM space.\n' );
        
        % Adds SPM12 lite to the path.
        ft_hastoolbox ( 'spm12', 1, 1 );
        
        % Gets the MRI file.
        mri_file      = sprintf ( '%s%s', config.path.mri, mrifile );
        
        % Uncompress the MRI file (not required?).
        if strcmp ( mri_file ( ( end - 2 ): end ), '.gz' )
            mri_file      = gunzip ( mri_file );
            mri_file      = mri_file {1};
            rmmri         = true;
        else
            rmmri         = false;
        end
        
        % Loads the MRI file and the MNI template.
        natmri        = spm_vol ( mri_file );
        mnimri        = spm_vol ( which ( 'single_subj_T1.nii' ) );
        
        % Co-registrates the MRI file to the MNI template.
        flags         = struct ();
        %flags.cost_fun = 'nmi';
        flags.cost_fun = 'ncc';
        
        dummy         = spm_coreg ( natmri, mnimri, flags );
        
        % Gets the native-to-MNI transformation in matrix form.
        nat2mni       = spm_matrix ( dummy (:)' );
        
        % Uses this transformation as starting point.
        nat2tpm       = nat2mni;
        
        % Removes the uncompressed MRI, if required.
        if rmmri, delete ( mri_file ), end
    end
    
    
    % Reslices the volume into a temporary AC-PC-like shape.
    mri_acpc       = mri;
    mri_acpc.transform = nat2tpm * mri.transform;
    
%     mri_acpc     = ft_volumereslice ( [], mri_acpc );
    [ ~, mri_acpc ] = evalc ( 'ft_volumereslice ( [], mri_acpc );' );
    
    % Performs a two-step affine registration from AC-PC to TPM space.
    acpc2tpm       = eye (4);
    acpc2tpm       = spm12_maff ( mri_acpc, 3, 16, tpm, acpc2tpm, 'mni' );
    acpc2tpm       = spm12_maff ( mri_acpc, 3,  0, tpm, acpc2tpm, 'mni' );
    
    % Combines both transformations.
    nat2tpm        = nat2tpm * acpc2tpm;

    
%     % Performs a two-step affine registration from native to TPM space.
%     nat2tpm        = spm12_maff ( mri, 3, 16, tpm, nat2tpm, 'mni' );
%     nat2tpm        = spm12_maff ( mri, 3,  0, tpm, nat2tpm, 'mni' );
    
    % Gets the linear transformation MRI voxels to TPM voxels.
    vox2tpmvox     = tpm.M \ nat2tpm * mri.transform;
    
    
    fprintf ( 1, '  Calculating a non-linear transformation to TPM space.\n' );
    
    % Finds a non-linear transformation to/from native from/to TPM space.
    cfg            = [];
    cfg.fwhm       = 0;
    cfg.samp       = 3;
    cfg.biasreg    = 1e-3;
    cfg.biasfwhm   = 60;
    cfg.lkp        = [ 1 2 3 3 4 4 4 5 5 5 5 6 6 ];
    cfg.reg        = [ 0 1e-3 .5 .05 .2 ];
    
    preprocinfo    = spm12_preproc ( mri, tpm, nat2tpm, cfg );
    
    % Gets the non-linear warping and the upsampling transformation.
    dswarp         = preprocinfo.Twarp;
    sam2ori        = inv ( preprocinfo.MT );
    
    
    fprintf ( 1, '  Non-linearly mapping the face tissue of the TPM to the MRI.\n' );
    
    % Gets the non-linear correspondence of each voxel.
    [ wx, wy, wz ] = spm12_warpaffine ( mri.dim, vox2tpmvox, dswarp, sam2ori );
    
    % Calculates the a-priori probability for each voxel.
    tpmpriors      = spm12_sample_priors ( tpm_face, wx, wy, wz );
    
    % Generates the face mask.
    facemask       = tpmpriors {1} > 0.1;
    
    
    fprintf ( 1, '  Blurring the face of the image.\n' );
    
    % Generates an anonymazed anatomy by blurring the original data.
    mriblur        = imgaussfilt3 ( mri.anatomy, 30 );
    
    % Makes a copy of the original MRI data and replaces the face.
    mri_def        = mri;
    mri_def.anatomy ( facemask ) = mriblur ( facemask );
    
    
    fprintf ( 1, '  Removing the face of the image.\n' );
    
    % Makes a copy of the original MRI data and removes the face.
    mri_trim       = mri;
    mri_trim.anatomy ( facemask ) = 0;
    
    
    fprintf ( 1, '  Saving the de-faced MRI data.\n' )
    
    % Saves the de-faced data in NIfTI format.
    ft_write_mri ( sprintf ( '%s%s.nii', config.path.def, basename ), mri_def.anatomy, 'transform', mri.transform, 'dataformat', 'nifti' )
    ft_write_mri ( sprintf ( '%s%s.nii', config.path.trim, basename ), mri_trim.anatomy, 'transform', mri.transform, 'dataformat', 'nifti' )
    
    % Gzips the data, if requested.
    if config.gzip
        
        % Gzips the data.
        gzip ( sprintf ( '%s%s.nii', config.path.def, basename ) )
        gzip ( sprintf ( '%s%s.nii', config.path.trim, basename ) )
        
        % Removes the uncompressed NIfTI files.
        delete ( sprintf ( '%s%s.nii', config.path.def, basename ) )
        delete ( sprintf ( '%s%s.nii', config.path.trim, basename ) )
    end
end
