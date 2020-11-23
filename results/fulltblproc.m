
ps=[1,2,4,8,12,16,24];
thrs=[1,2,4,8,12,16,24];

Tbl = -1*ones(7*4,7);
Rt=R;
Rit=Ri;

for i=1:7
    for j=1:7
        p=ps(i);
        thr=thrs(j);
        
        if (p*thr <25)
            
            I = find(R(:,1)==lib & R(:,2)==p & R(:,3)==thr & R(:,4)==n);
            It = find(Rt(:,1)==lib & Rt(:,2)==p & Rt(:,3)==thr & Rt(:,4)==n);
            %It = find(Rt(:,1)==lib & Rt(:,2)==thr & Rt(:,3)==n)
            
            if (I >0)
                %I = find(R(:,1)==lib & R(:,4)==n);
                
                Tbl(4*(i-1)+1  ,j) = (R(I,end))';
                if (i==1 & j==1)
                    Dval = Rt(It,end);
                end
                
                Tbl(4*(i-1)+2  ,j) = Tbl(4*(i-1)+1,j)/Dval;
                
                I = find(Ri(:,1)==lib & Ri(:,2)==p & Ri(:,3)==thr & Ri(:,4)==n);
                It = find(Rit(:,1)==lib & Rit(:,2)==p & Rit(:,3)==thr & Rit(:,4)==n);
                %It = find(Rit(:,1)==lib & Rit(:,2)==thr & Rit(:,3)==n);
                %I = find(Ri(:,1)==lib & Ri(:,4)==n);
                
                Tbl(4*(i-1)+3  ,j) = (Ri(I,end))';
                
                if (i==1 & j==1)
                    Ival = Rit(It,end);
                end
                Tbl(4*(i-1)+4  ,j) = Tbl(4*(i-1)+3  ,j)/Ival;
            end
        end
    end
end

strtbl = cell(28,21);

for i=1:28
    for j=1:9
        strtbl{i,2*j} = " & ";
    end
    for j=1:7
        strtbl{i,2*j+5} = "    -    ";
    end
    strtbl{i,20} = " \\ ";
    strtbl{i,21} = " ";
    strtbl{i,1} = "    ";
    strtbl{i,3} = "    ";
end

for i=1:7
    strtbl{4*(i-1)+1,1} = ps(i);
    strtbl{4*(i-1)+1,3} = rc;
    strtbl{4*(i-1)+1,5} = "ftime";
    strtbl{4*(i-1)+2,5} = "fratio";
    strtbl{4*(i-1)+3,5} = "itime";
    strtbl{4*(i-1)+4,5} = "iratio";    
    strtbl{4*(i-1)+4,21} = "\hline";  
end

for i=1:14
    for j=1:7
        if (Tbl(2*(i-1)+1,j) ~= -1)
            strtbl{2*(i-1)+1,2*j+5} = compose("%3.2e",Tbl(2*(i-1)+1,j));
        end
        if (Tbl(2*(i-1)+2,j) ~= -1)
            strtbl{2*(i-1)+2,2*j+5} = compose("%4.2f",Tbl(2*(i-1)+2,j));
        end
        
    end
end

fulltbl = cell(28,1);
for i=1:28
   str = strtbl{i,1};
   for j=2:21
       str = str + strtbl{i,j};
   end
   fulltbl{i,1} = str;
end

fulltblstr = string(fulltbl)