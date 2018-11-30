clear
load Baidu_SDdata\DS_Nets.mat;
indexes = find(D_Lv2_Class);
P = D_Net*diag(1./sum(D_Net));
deg = sum(P~=0, 2)
isum =0;
count=0;
osum =0;

for index = indexes'
%length(find(D_Lv2_Class))
%unique(D_Lv2_Class)
    count = count+1
    
    T = find(D_Lv2_Class == D_Lv2_Class(index));
    len = length(T)/length(indexes)*9721; 
    
    [bestset,obestset, pindex,pscores,opindex, opscores, bestcond,bestcut,bestvol] = pprgrow(sparse(P),index,'maxexpand', 1000);
    pr = zeros(9721,1);
    pr(pindex) = pscores;
    pr = pr./deg;
    [V, I] = sort(pr, 'descend');
    S = I(1:int16(len));
    S = intersect(S,indexes);


    opr = zeros(9721,1);
    opr(opindex) = opscores;
    opr = opr./deg;
    [oV, oI] = sort(opr, 'descend');
    oS = oI(1:int16(len));
    oS = intersect(oS,indexes);
    
    T = find(D_Lv2_Class == D_Lv2_Class(index));
    f1 = 2*length(intersect(S, T))/(length(S)+length(T));
    of1 = 2*length(intersect(oS, T))/(length(oS)+length(T));
    
    
    isum = isum+f1;
    ave = isum/count
    
    osum = osum+of1;
    oave = osum/count
end 