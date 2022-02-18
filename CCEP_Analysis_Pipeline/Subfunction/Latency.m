function loc = Latency( a )
  b=-1*a;
  [locs pks]=peakseek(b,100);
  loc=locs;
end


