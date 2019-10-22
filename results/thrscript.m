for i=[1,4,16]
    I = find(R(:,1)==lib & R(:,2)==1 & R(:,3)==i );
    switch i
        case 1
           ll='o'; cc=[0 0.4470 0.7410];
        case 4
            ll='x'; cc=[0.8500 0.3250 0.0980];
        case 8
            ll='s'; cc=[0.9290 0.6940 0.1250];
        case 16
            ll='d'; cc=[0.4940 0.1840 0.5560];
        case 24
            ll='p'; cc=[0.4660 0.6740 0.1880];
        
    end
    ll = strcat(lls,ll);
    loglog(R(I,4),R(I,11),ll,'linewidth',2,'MarkerSize',8.0,'color',cc)
    hold on
end
