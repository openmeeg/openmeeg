function save_bin(data,filename)

% SAVE_BIN   Save binary matrix
%
%   Save binary matrix
%
%   SYNTAX
%       SAVE_BIN(DATA,FILENAME)
%

%
%   Created by Alexandre Gramfort on 2007-05-22.
%

file = fopen(filename,'w');
dims = size(data);
fwrite(file,dims,'uint32');
fwrite(file,data(:),'double');
fclose(file);
