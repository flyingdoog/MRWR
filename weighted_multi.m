function [ave] = wRWRM(alpha,mlevel,mute)
tic

% load Baidu_SDdata\DS_Nets.mat;
% if mlevel==1
%     D_Class = D_Lv1_Class;
% else
%     D_Class = D_Lv2_Class;
% end

% load Gene-Disease\DG_Nets_v2.mat
% if mlevel==1
%     D_Net = D_Net;
%     D_Class = D_Label;
%     S_Net = G_Net;
%     A_Net = A_Net';
% else
%     S_Net = D_Net;
%     D_Net = G_Net;
%     D_Class = G_Label;
% end

load DBLP_4AREA\AC_Nets_v4.mat
if mlevel==1
    D_Net = A_Net;
    D_Class = A_Label';
    S_Net = C_Net;
    A_Net = CA_Net;
else
    S_Net = A_Net;
    D_Net = C_Net;
    D_Class = C_Label;
    A_Net = CA_Net';
end
% 
% load synthetic\1m.mat
% if mlevel==1
%     D_Net = A_Net;
%     D_Class = A_Label;
%     S_Net = B_Net;
%     A_Net = AB;
% else
%     S_Net = A_Net;
%     D_Net = B_Net;
%     D_Class = B_Label;
%     A_Net = AB';
% end



indexes = find(D_Class);


AT = A_Net';
AT = AT ./(sum(AT)+1e-5);


D_Net = D_Net./sum(D_Net);
A = A_Net./(sum(A_Net)+1e-5);
% A(isnan(A)) = 0;
% AT(isnan(AT)) = 0;
S_Net = S_Net ./ sum(S_Net);

D_Net = D_Net*alpha;
AT = AT*alpha;
A = A*(1-alpha);
S_Net = S_Net*(1-alpha);

P = [D_Net,AT; A, S_Net];
P = P ./sum(P);
P(isnan(P)) = 0;
deg = sum(P~=0, 2);
isum =0;
count=0;
osum =0;

tsize = size(D_Net,1);
vsize = size(P,1);

for index = indexes'
%length(find(D_Class))
%unique(D_Class)

    if count>1000
        break
    end
    
    T = find(D_Class == D_Class(index));
    len = length(T)/length(indexes)*tsize;
    
    [bestset, pindex,pscores, bestcond,bestcut,bestvol] = weighted_pprgrow(sparse(P),index,'maxexpand', 10000);
    pr = zeros(vsize,1);
    pr(pindex) = pscores;
    pr = pr./deg;
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