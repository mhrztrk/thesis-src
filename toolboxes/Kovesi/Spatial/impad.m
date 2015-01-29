% IMPAD - adds zeros to the boundary of an image
%
% Usage:  paddedim = impad(im, b)
%
% Arguments:     im - Image to be padded (greyscale or colour)
%                 b - Width of padding boundary to be added
%
% Returns: paddedim - Padded image of size rows+2*b x cols+2*b
%
% See also: IMTRIM, IMSETBORDER

% Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% http://www.csse.uwa.edu.au/~pk/research/matlabfns/
% June  2010

function pim = impad(im, b)

    assert(b >= 1, 'Padding size must be >= 1')

    b = round(b);     % Ensure integer
    [rows, cols, channels] = size(im);
    pim = zeros(rows+2*b, cols+2*b, channels);
    pim(1+b:rows+b, 1+b:cols+b, 1:channels) = im;