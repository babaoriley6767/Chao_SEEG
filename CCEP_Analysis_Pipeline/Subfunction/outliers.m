function y =outliers(B)

  Q1=prctile(B,25);
  Q3=prctile(B,75);
  IQR = Q3-Q1;
  B(B>1.5*IQR+Q3|B<Q1-1.5*IQR)=[];
  y = mean(B);
end