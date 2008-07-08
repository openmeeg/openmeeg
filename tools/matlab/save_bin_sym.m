function save_bin_sym(data,filename)

% SAVE_BIN_SYM   Save binary symmetric Matrix
%
%   Save binary symmetric Matrix
%   
%   SYNTAX
%       SAVE_BIN_SYM(DATA,FILENAME)
%   

%
%   Created by Alexandre Gramfort on 2007-09-26.
%

dims = size(data);
assert(dims(1) == dims(2),'Matrix non square')
assert(isempty(find(data ~= data')),'Matrix non symmetric')

file = fopen(filename,'w');
dim = dims(1);
fwrite(file,dim,'uint32');
fwrite(file,data(tril(ones(dim,dim)) > 0),'double');
fclose(file);

