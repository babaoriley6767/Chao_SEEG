function cellarray=stringunpad(string,ns,padlen)

cellarray=cell(ns,1);

for i=1:ns
    j=padlen*(i-1)+1;
    cellarray{i}=cellstr(string(j:j+padlen-1)');
end

end
