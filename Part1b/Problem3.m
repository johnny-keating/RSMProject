%Problem 3
m=3;
n=6;
A=[2	2	-4	1	1	0
  -1	3	7	1	1	1
   5	-2	1	0	-1	1 ];
c= [ 4 -5 -3 -4  1  1]';
b = [10  15 19 ]';
[z, x, pie, indices, exitflag] = fullsimplex(A, b, c, m, n)
