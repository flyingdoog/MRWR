function [ave, oave] = RWRM(beta,mute)

tic
load DBLP_4AREA\AC_Nets_v4.mat
S_Net = A_Net;
D_Net = C_Net;
D_Class = C_Label;
A_Net = CA_Net';

indexes = find(D_Class);


AT = A_Net';
AT = AT ./(sum(AT)+1e-5);


D_Net = D_Net./sum(D_Net);
A = A_Net./(sum(A_Net)+1e-5);
% A(isnan(A)) = 0;
% AT(isnan(AT)) = 0;
S_Net = S_Net ./ sum(S_Net);

D_Net = D_Net*beta;
AT = AT*beta;
A = A*(1-beta);
S_Net = S_Net*(1-beta);

P = [D_Net,AT; A, S_Net];
P = P ./sum(P);
P(isnan(P)) = 0;
deg = sum(P~=0, 2);
isum =0;
count=0;
osum =0;

psum=0;
rsum=0;

opsum=0;
orsum=0;


tsize = size(D_Net,1);
vsize = size(P,1);

for index = indexes'

    if count>100
        break
    end
    
    T = find(D_Class == D_Class(index));
    len = length(T)/length(indexes)*tsize;
    
    [bestset, pindex,pscores, bestcond,bestcut,bestvol] = pprgrow(sparse(P),index,'maxexpand', 10000);
    pr = zeros(vsize,1);
    pr(pindex) = pscores;

    [V, I] = sort(pr, 'descend');
    I(I>tsize)=[];
    S = I(1:int16(len));
    S = intersect(S,indexes);

    
    T = find(D_Class == D_Class(index));

    f1 = 2*length(intersect(S, T))/(length(S)+length(T));
    
    
    isum = isum+f1; 
    
    
    if mute
        count = count+1;
        ave = isum/count;
    else
        count = count+1
        ave = isum/count
    end
end 
toc
end 