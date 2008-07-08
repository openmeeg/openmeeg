function save_txt_sym(data,filename)

% SAVE_TXT_SYM   Save ascii symmetric Matrix
%
%   Save ascii symmetric Matrix
%   
%   SYNTAX
%       SAVE_TXT_SYM(DATA,FILENAME)
%   

%
%   Created by Alexandre Gramfort on 2007-09-26.
%

dims = size(data);
assert(dims(1) == dims(2),'Matrix non square')
assert(isempty(find(data ~= data')),'Matrix non symmetric')

dim=dims(1);
for i=1:dim
    if i == 1
        dlmwrite(filename, data(i,i:end), 'delimiter', '\t','precision',18);
    else
        dlmwrite(filename, data(i,i:end), 'delimiter', '\t','-append','precision',18);
    end
end
