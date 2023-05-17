function yout = selfpooldata(xaug,LibraryType)
switch LibraryType
    case 4
        if size(xaug,1)==4   
        xaug = xaug';
        end
        u = xaug(:,1);   
        v = xaug(:,2);
        r = xaug(:,3);
        d = xaug(:,4);
        %Silvestre model
        yout = [u.*abs(u),u.*v,u.*r,v.*v,v.*r,v.*abs(v),r.*r,r.*abs(r),u.*u.*d,r.*r.*r,u.*u.*d.*d,1./u];
        
end

end

