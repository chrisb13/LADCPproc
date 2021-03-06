function h = pcolorn(x,y,c)
%   map to imagesc
%   M. Visbeck 2004

if nargin < 1
    error('Too few input arguments.');
elseif nargin > 4
    error('Too many input arguments.')
end

if nargin == 1
    hh = imagesc(flipud(x))
elseif nargin == 3
    hh = imagesc(x,-y,c);
else
    error('Must have one or three input arguments.')
end
