function xfmmat = myfs_read_xfm ( xfmfile )
% Reads FreeSurfer transformation files (*.xfm).
% Only works with linear transformations.

% Reads the content of the file as a string.
fid = fopen ( xfmfile, 'rt' );
xfmraw = fread ( fid, [ 1 Inf ], '*char' );
fclose ( fid );


% Gets the lines corresponding to the transformation matrix.
xfmmat = regexp ( xfmraw, 'Linear_Transform = (.*);', 'tokens' );

if isempty ( xfmmat )
    error ( 'No linear transformation matrix found in the file.' )
end


% Parses the transformation matrix.
xfmmat = xfmmat {1} {1};
xfmmat = str2num ( xfmmat ); %#ok<ST2NM>

% Checks the matrix size.
if size ( xfmmat, 1 ) ~= 3 || size ( xfmmat, 2 ) ~= 4
    error ( 'Incorrect transformation matrix size.' )
end


% Expands the transformation matrix to 4 x 4.
xfmmat ( 4, 4 ) = 1;
