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
