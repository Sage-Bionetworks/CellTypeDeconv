cd e:\sagedocs\CellTypeDeconv\NMF_manifold\
load('SynMouseESprog_1000.mat')


A = zeros(100);
for i = 1:99
A(i,i+1) = 1;
A(i+1,i) = 1;
end
A = sparse(A);


options = [];
options.WeightMode = 'Binary';
options.maxIter = 100;
options.maxIter = 100;
options.alpha = .1;


[U2,V2] = nnmf(Dat,2);
V2 = V2*inv(diag(sum(V2)));

[U,V] = GNMF(full(Dat),2,A,options,U2,V2'); %'
V = V';
V = V*inv(diag(sum(V)));
plot(1:100,V(1,1:100))
hold on
plot(1:100,Xact(1,1:100),'r')
plot(Xact(1,1:100),V(1,1:100))


corr([1:100]',V2(1,1:100)','Type','Spearman')
corr([1:100]',V(1,1:100)','Type','Spearman')

plot(Xact(1,1:100),Xact(1,1:100),'g',Xact(1,1:100),V(1,1:100),'b',Xact(1,1:100),V2(1,1:100))
