function [xref,yref,mref] = reflection(x3,y3,m,lineimg)
mref = -m;
if(x3==1)
    if((1-y3)/mref+1>0 && (1-y3)/mref+1<size(lineimg,1))
        yref = 1;
        xref = (1-y3)/mref+x3;
    elseif((mref*(size(lineimg,1) - x3)+y3)<size(lineimg,2) && (mref*(size(lineimg,1) - x3)+y3)>0)
        xref = size(lineimg,1);
        yref = mref*(xref - x3)+y3;
    elseif(y3+mref*(1-x3)>0 && y3+mref*(1-x3)<size(lineimg,2))
        yref = size(lineimg,2);
        xref = x3+((size(lineimg,2)-y3)/mref);
    end
elseif(y3 == 1)
    if(x3+(size(lineimg,2)-y3)/mref>0 && x3+(size(lineimg,2)-y3)/mref<size(lineimg,1))
        yref = size(lineimg,2);
        xref = x3+((size(lineimg,2)-y3)/mref);
    elseif((mref*(size(lineimg,1) - x3)+y3)<size(lineimg,2) && (mref*(size(lineimg,1) - x3)+y3)>0)
        xref = size(lineimg,1);
        yref = mref*(xref - x3)+y3;
    elseif(y3+mref*(1-x3)>0 && y3+mref*(1-x3)<size(lineimg,2))
        xref = 1;
        yref = mref*(1-x3)+y3;
    end
 elseif(x3==size(lineimg,1))
    if((1-y3)/mref+x3>0 && (1-y3)/mref+x3<size(lineimg,1))
        yref = 1;
        xref = (1-y3)/mref+x3;
    elseif(x3+((size(lineimg,2)-y3)/mref)<size(lineimg,1) && x3+((size(lineimg,2)-y3)/mref)>0)
        yref = size(lineimg,2);
        xref = x3+((size(lineimg,2)-y3)/mref);
    elseif(y3+mref*(1-x3)>0 && y3+mref*(1-x3)<size(lineimg,2))
        xref = 1;
        yref = mref*(1-x3)+y3;
    end
elseif(y3 == size(lineimg,2))
    if((1-y3)/mref+x3>0 && (1-y3)/mref+x3<size(lineimg,1))
        yref = 1;
        xref = (1-y3)/mref+x3;
    elseif((mref*(size(lineimg,1) - x3)+y3)<size(lineimg,2) && (mref*(size(lineimg,1) - x3)+y3)>0)
        xref = size(lineimg,1);
        yref = mref*(xref - x3)+y3;
    elseif(y3+mref*(1-x3)>0 && y3+mref*(1-x3)<size(lineimg,2))
        xref = 1;
        yref = mref*(1-x3)+y3;
    end
 end
