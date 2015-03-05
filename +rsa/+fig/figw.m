function h=figw(figI)

    if ~exist('figI','var')||isempty(figI), figI=gcf; end
    if figI
        h=figure(figI); 
    else
        h=figure; 
    end
    set(h,'Color','w');

end%function
