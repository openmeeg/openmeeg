function save_txt(data,filename)

% SAVE_TXT   Save ascii matrix
%
%   Save ascii matrix
%
%   SYNTAX
%       SAVE_TXT(DATA,FILENAME)
%

%
%   Created by Alexandre Gramfort on 2007-05-22.
%

data = double(data);
save(filename,'data','-ASCII','-double','-v6')
