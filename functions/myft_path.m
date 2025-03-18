function myft_path

% Search for a version of FieldTrip in the path.
if isempty ( which ( 'ft_defaults' ) )

    % Looks for FieldTrip in the current folder.
    ft_path = dir ( sprintf ( '%s/fieldtrip-*', pwd ) );
    ft_path = ft_path ( [ ft_path.isdir ] );

    % If not found, looks for FieldTrip in the parent folder.
    if ~numel ( ft_path )
        ft_path = dir ( sprintf ( '%s/fieldtrip-*', fileparts ( pwd ) ) );
        ft_path = ft_path ( [ ft_path.isdir ] );
    end
    
    % If not found, rises an error.
    if ~numel ( ft_path )
        error ( ...
            'preproc:NoFT', [ ...
            'FieldTrip not found neither in your path, the current folder, nor its parent folder.\n' ...
            'This script cannot continue.' ] )
    end

    % Otherwise adds FieldTrip to the path.
    addpath ( sprintf ( '%s/%s/', ft_path ( end ).folder, ft_path ( end ).name ) )
end

% Initializes FieldTrip.
ft_defaults

% Disables FieldTrip feedback.
global ft_default;
ft_default.showcallinfo = 'no';
ft_default.checkconfig  = 'silent';
