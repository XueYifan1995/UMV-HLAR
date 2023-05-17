function yout = selfpooldatalist(ahat,LibraryType)

switch LibraryType
         
    case 4
        newout(1) = {''};  
        newout{1,2} = ['udot'];
        newout{1,3} = ['vdot'];
        newout{1,4} = ['rdot'];
        for k=1:size(ahat,1)   
            for j=1:size(ahat,2)
                newout{k+1,1+j} = ahat(k,j);
            end
        end
        newout{2,1} = ['u|u|'];  
        newout{3,1} = ['uv'];
        newout{4,1} = ['ur'];
        newout{5,1} = ['vv'];
        newout{6,1} = ['vr'];
        newout{7,1} = ['v|v|'];
        newout{8,1} = ['rr'];
        newout{9,1} = ['r|r|'];
        newout{10,1} = ['uud'];
        newout{11,1} = ['rrr'];
        newout{12,1} = ['uudd'];
        newout{13,1} = ['1/u'];
        
        yout = newout 
        
end
