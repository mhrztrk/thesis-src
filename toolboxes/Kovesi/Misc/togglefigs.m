% TOGGLEFIGS - toggles figures to aid comparison
%
% Usage:  togglefigs(figs)
%
% Argument:  figs - figure numbers, separated by commas, to toggle
%
% Example if you have 3 figures you wish to compare manually drag them until
% they are perfectly overlaid, then use
% >> togglefigs(1, 2, 3)

% PK March 2010

function togglefigs(varargin)
    
    figs = getfigs(varargin(:));
    
    fprintf('Hit any key to toggle images, ''X'' to exit\n'); 
    while 1
        for n = 1:length(figs)
            figure(figs(n));
            
            pause;
            ch = get(gcf,'CurrentCharacter');
            if lower(ch)=='x'
                return
            end
        end
    end
    
%------------------------------------------
function figs = getfigs(arg)
    
    figs = zeros(size(arg));
    for n = 1:length(arg)
        figs(n) = arg{n};
    end
    