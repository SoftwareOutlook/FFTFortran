

I = find(R(:,1)==lib & R(:,2)==1 & R(:,4)==n);
%I = find(R(:,1)==lib & R(:,4)==n);

Dval = (R(I,end))'
Drto = Dval/Dval(1)

I = find(Ri(:,1)==lib & Ri(:,2)==1 & Ri(:,4)==n);
%I = find(Ri(:,1)==lib & Ri(:,4)==n);

Ival = (Ri(I,end))'
Irto = Ival/Ival(1)