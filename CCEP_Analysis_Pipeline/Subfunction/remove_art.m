function y = remove_art(a,x)     % a: data x: stim time
    t1 = a(x-15:x-5);t2 = a(x+5:x+15);
    r1 = flip(t1);r2 = flip(t2);
    x1 = 1;x2 = 0;tt=1/10;
      for j=1:11 
          t1(j) = r1(j)*x1;
          t2(j) = r2(j)*x2;
          x1 = x1-tt;          
          x2 = x2+tt;
      end
    y=a;
    for j =1:11
        y(x+j-1-5) = t1(j)+t2(j);
end

