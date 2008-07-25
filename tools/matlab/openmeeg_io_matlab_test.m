% =============================
% = Testing standard matrices =
% =============================

data = randn(4,9);

om_save_full(data,'test.txt','ascii');
om_save_full(data,'test.bin','binary');

data_txt = om_load_full('test.txt','ascii');
data_bin = om_load_full('test.bin','binary');

norm(data_txt - data_bin)
norm(data_txt - data)
norm(data_bin - data)

delete 'test.txt'
delete 'test.bin'

% =============================
% = Testing symmetric matrices =
% =============================

randn('seed',0);
data = randn(5,5);
data = (data+data')/2;

om_save_sym(data,'test.txt','ascii');
om_save_sym(data,'test.bin','binary');

data_txt = om_load_sym('test.txt','ascii');
data_bin = om_load_sym('test.bin','binary');

norm(data_txt - data_bin)
norm(data_txt - data)
norm(data_bin - data)

delete 'test.txt'
delete 'test.bin'

% =============================
% = Testing sparse matrices =
% =============================

data = sprand(5,5,0.5);

om_save_sparse(data,'test.txt','ascii');
om_save_sparse(data,'test.bin','binary');

data_txt = om_load_sparse('test.txt','ascii');
data_bin = om_load_sparse('test.bin','binary');

norm(full(data_txt - data_bin))
norm(full(data_txt - data))
norm(full(data_bin - data))

delete 'test.txt'
delete 'test.bin'

