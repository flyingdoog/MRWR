clear
parfor i = [1:10] 
    if i~=6
       [ave6,oave6] = multi(i/10.0,2,1);
        ave(i)=ave6
        rave(i)=oave6
        disp([i,ave6,oave6])
        disp('completed')
    end
end
disp('CA')
% matlabpool('close')