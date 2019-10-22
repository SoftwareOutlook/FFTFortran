function [res1,res2,res3,res4,res5]=dataperm(res1,res2,res3,res4,res5)

[n,m] = size(res1);
res = zeros(n,11);
%resi = res;

% Permute rows to correct order
for i=1:5
    switch i
        case 1
            r = res1;
        case 2
            r = res2;
        case 3
            r = res3;
        case 4
            r = res4;
        case 5
            r = res5;
    end
    [I,J] = sort(r(:,4));
    r = r(J,:);
    [I,J] = sort(r(:,3));
    r = r(J,:);
    [I,J] = sort(r(:,2));
    r = r(J,:);
    [I,J] = sort(r(:,1));
    r = r(J,:);
    switch i
        case 1
            res1=r;
        case 2
            res2=r;
        case 3
            res3=r;
        case 4
            res4=r;
        case 5
            res5=r;
    end    
end
            
