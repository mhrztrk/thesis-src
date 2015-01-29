% STRENDSWITH - tests if a string ends with a specified substring
%
% Usage: b = strendswith(str, substr)
%
%  Arguments:
%        str   - string to be tested
%     substr   - ending of string that we are hoping to find
%
%   Returns: true/false. Note that case of strings is ignored
% 
% See also: STRSTARTSWITH

% Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% http://www.csse.uwa.edu.au/~pk/research/matlabfns/
% June  2010

function b = strendswith(str, substr)
    
    % Compute index of character in str that should match with the the first
    % character of str
    s = length(str) - length(substr) + 1;
    
    % True if s > 0 and all appropriate characters match (ignoring case)
    b =  s > 0 && strcmpi(str(s:end), substr);
    
