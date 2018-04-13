/* this module implements the SpringRank routine of Larrimore, based on his
   github matlab code.  Note that as implemented it assumes a connected
   component.
  
   Author: Moody
   Source: https://github.com/cdebacco/SpringRank
   Date: 4.10.2018

   Input: An adjacency matrix for a connected component.
   output: a 1-dimensional score of position in the hierarchy, relative to 
     node n which has value 0.

*/

start springrank(a);
  n=nrow(a);
  dout=a[,+];
  din=a[+,];
  din=din`;

  *print dout din;
  dNout=dout[n];
  dNin=din[n];

  dout=diag(dout[1:N-1]);
  din=diag(din[1:n-1]);

*out-degree and in-degree for vertex N;
dNout = dout[N];
dNin = din[N];
at=a`;

B=Dout+Din - a[1:n-1,1:n-1] - a[1:n-1,1:n-1]` - repeat(a[n,1:n-1],n-1,1) 
    -repeat(at[n,1:n-1] ,n-1,1);

lb = vecdiag(Dout)-vecdiag(Din)+dNout-dNin;
t = solve(B, lb);
*print b lb t;

s=t//0;
return(s);

finish;
