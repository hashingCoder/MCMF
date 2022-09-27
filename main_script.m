clear all

i = 2;

datasets={ 'blogcatalog','PPI','wiki','flickr'};
data_dir = 'data/network';

dataset = datasets{i};

load(sprintf('%s/netmf_%s.mat', data_dir, dataset))

dim = 128;

n = length(net);

epsilon = 1e-6;
max_iter = 100;
gamma = 0.08;

random_B = randi([0,1],n,dim);
random_B = 2*random_B - ones(n,dim);
init_b = random_B;

tic;
[b] = MCMF(net,init_b,gamma,epsilon,max_iter);
t = toc;

bit_blance = b'*ones(n,1);