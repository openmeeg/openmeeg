% =============================
% = Testing standard matrices =
% =============================

data = randn(4,9);

save_bin(data,'test.bin');
save_txt(data,'test.txt');

data_txt = load_txt('test.txt');
data_bin = load_bin('test.bin');

norm(data_txt - data_bin)
norm(data_txt - data)
norm(data_bin - data)

delete 'test.txt'
delete 'test.bin'

% =============================
% = Testing symmetric matrices =
% =============================

data = randn(5,5);
data = randn(2,2);
data = (data+data')/2;

save_bin_sym(data,'test.bin');
save_txt_sym(data,'test.txt');

data_txt = load_txt_sym('test.txt');
data_bin = load_bin_sym('test.bin');

norm(data_txt - data_bin)
norm(data_txt - data)
norm(data_bin - data)

delete 'test.txt'
delete 'test.bin'

