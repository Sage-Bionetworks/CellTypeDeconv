%Clustering on COIL20
clear;
load('COIL20.mat');	
nClass = length(unique(gnd));

%Normalize each data vector to have L2-norm equal to 1  
fea = NormalizeFea(fea);

%Clustering in the original space
rand('twister',5489);
label = litekmeans(fea,nClass,'Replicates',20);
MIhat = MutualInfo(gnd,label);
disp(['Clustering in the original space. MIhat: ',num2str(MIhat)]);
%Clustering in the original space. MIhat: 0.7386

%NMF learning
options = [];
options.maxIter = 100;
options.alpha = 0;
%when alpha = 0, GNMF boils down to the ordinary NMF.
rand('twister',5489);
[U,V] = GNMF(fea',nClass,[],options); %'

%Clustering in the NMF subspace
rand('twister',5489);
label = litekmeans(V,nClass,'Replicates',20);
MIhat = MutualInfo(gnd,label);
disp(['Clustering in the NMF subspace. MIhat: ',num2str(MIhat)]);
%Clustering in the NMF subspace. MIhat:  0.74361

%GNMF learning
options = [];
options.WeightMode = 'Binary';  
W = constructW(fea,options);
options.maxIter = 100;
options.alpha = 100;
rand('twister',5489);
[U,V] = GNMF(fea',nClass,W,options); %'

%Clustering in the GNMF subspace
rand('twister',5489);
label = litekmeans(V,nClass,'Replicates',20);
MIhat = MutualInfo(gnd,label);
disp(['Clustering in the GNMF subspace. MIhat: ',num2str(MIhat)]);
%Clustering in the GNMF subspace. MIhat: 0.87599