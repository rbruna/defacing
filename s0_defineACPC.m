clc
clear
close all

% Defines the files locations.
config.path.mri  = '../../data/anatomy/T1/';
config.path.patt = '*.nii.gz';


% Adds the 'functions' folders to the path.
addpath ( sprintf ( '%s/functions/', pwd ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path

% Adds the FT toolboxes that will be required.
ft_hastoolbox ( 'freesurfer', 1, 1 );


% Gets the file list.
files         = dir ( sprintf ( '%s%s', config.path.mri, config.path.patt ) );

% Goes through all the files.
for findex = 1: numel ( files )
    
    % Gets the base file name.
    mrifile       = files ( findex ).name;
    basename      = regexprep ( mrifile, '\.nii(.gz)?$', '' );
    
    fprintf ( 1, 'Loading NIfTI file %s.\n', mrifile );
    
    % Gets the MRI file.
    mri           = myft_read_mri ( sprintf ( '%s%s', config.path.mri, mrifile ) );
    mri.coordsys  = 'ras';
    
    
    % Asks for the AC-PC landmarks.
    cfg           = [];
    cfg.coordsys  = 'acpc';
    
    dummy         = ft_volumerealign ( cfg, mri );
    landmark      = dummy.cfg.fiducial;
    drawnow
    
    
    % Saves the landmarks.
    save ( '-v6', sprintf ( '%s%s.mat', config.path.mri, basename ), '-struct', 'landmark' )
end
