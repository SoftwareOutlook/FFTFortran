
ps=[1,2,4,8,12,16,24];
thrs=[1,2,4,8,12,16,24];

Tbl = -1*ones(7*4,7);

for i=1:7
    for j=1:7
        p=ps(i)
        thr=thrs(j)
        
        if (p*thr <25)
            
            I = find(R(:,1)==lib & R(:,2)==p & R(:,3)==thr & R(:,4)==n)
            
            if (I >0)
                %I = find(R(:,1)==lib & R(:,4)==n);
                
                Tbl(4*(i-1)+1  ,j) = (R(I,end))';
                if (i==1 & j==1)
                    Dval = R(I,end);
                end
                
                Tbl(4*(i-1)+2  ,j) = Tbl(4*(i-1)+1,j)/Dval;
                
                I = find(Ri(:,1)==lib & Ri(:,2)==p & Ri(:,3)==thr & Ri(:,4)==n);
                %I = find(Ri(:,1)==lib & Ri(:,4)==n);
                
                Tbl(4*(i-1)+3  ,j) = (Ri(I,end))';
                
                if (i==1 & j==1)
                    Ival = Ri(I,end);
                end
                Tbl(4*(i-1)+4  ,j) = Tbl(4*(i-1)+3  ,j)/Ival;
            end
        end
    end
end

Tbl