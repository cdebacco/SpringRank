# SpringRank

SpringRank is a method for ranking the nodes in a directed network or from a set of pairwise comparisons. 

In this repo, we provide three implementations of the SpringRank model described in:

* [1] C. De Bacco, D. B. Larremore and C. Moore, *A physical model for efficient ranking in networks*, Science Advances, Vol 4, **7**, eaar8260, 2018. [Science Advances PDF](http://advances.sciencemag.org/content/4/7/eaar8260). [arXiv PDF](https://arxiv.org/abs/1709.09002) and [here](http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf).  

If you use this code please cite [1].


## What's included:
- `python` : Python code and a test script.
- `matlab` : MATLAB code and a test script.
- `r` : R code and a test script.
- `iml_sas:` SAS/IML code.
- `data` : Contains sample adjacency files to test the code and a sample result.

## Other implementations:
- `Haskell` : by Jon Zingale [github.com/jonzingale/Haskell/](https://github.com/jonzingale/Haskell/tree/master/SpringRank).
- `Python` : Fast Sparse alternative code by Larremore Lab on PyPI [github.com/LarremoreLab/SpringRank](https://github.com/LarremoreLab/SpringRank)

## Python Notes:
Need to make a directory called `data` at the same level of the `python` folder. 
To make one, just type from the command line, inside that folder: 
* `mkdir data`

#### Input format.
The directed adjacency matrix should be formatted as an edge list with 3 columns:

`node1 node2 3 `

The first and second columns are the source and target nodes of that edge, respectively; the third is the edge weigth (must be integer). In this example the edge node1 --> node2 exists with weight 3.

#### Output.
One file will be generated inside the `data` folder containg the SpringRank scores ordered from highest to lowest. The output file will be inside `data` folder with names:
- `networkname_SpringRank_a0.0_l0_1.0_l1_1.0.dat`  for the case `alpha=0` and `l_0=l_1=1.0`. 

The first column is the node id, the second column is the SpringRank score.




Copyright (c) 2017 [Caterina De Bacco](http://cdebacco.com) and [Daniel B Larremore](https://larremorelab.github.io)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
