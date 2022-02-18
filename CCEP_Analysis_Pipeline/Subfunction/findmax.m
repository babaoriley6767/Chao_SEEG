function [amp latency] =findmax(x,st,ed)
  amp = -10000;
  latency = 0;
  for i = st:ed
      if abs(x(i))>amp
          amp = abs(x(i));
          latency = i;
      end
  end
  