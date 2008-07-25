function om_save_full(data,filename,format)

% OM_SAVE_FULL   Save full Matrix
%
%   Save full Matrix
%
%   SYNTAX
%       OM_SAVE_FULL(DATA,FILENAME,FORMAT)
%
%       FORMAT : can be 'ascii' or 'binary' (default)
%

% $Id$
% $LastChangedBy$
% $LastChangedDate$
% $Revision$

me = 'OM_SAVE_FULL';

if nargin == 0
    eval(['help ',lower(me)])
    return
end

if nargin < 3
    format = 'binary';
end

switch format
case 'binary'
    disp(['Saving file ',filename])
    file = fopen(filename,'w');
    dims = size(data);
    fwrite(file,dims,'uint32','ieee-le');
    fwrite(file,data(:),'double','ieee-le');
    fclose(file);
case 'ascii'
    data = double(data);
    save(filename,'data','-ASCII','-double','-v6')
otherwise
    error([me,' : Unknown file format'])
end

