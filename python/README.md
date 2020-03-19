- `SpringRank/` folder containing the main function to calculate SpringRank scores. To use it in your code, add it by writing on top:
    - `from SpringRank import SpringRank`  
    see  `test_calculate_SpringRank.py` for an example.

- `SpringRank_tools.py` containes the function to generate a network according to the SpringRank generative model.

- `tools.py` containes auxiliary functions needed to process the input and output.

- `test_calculate_SpringRank.py` is a sample script for testing the code. It gives an example of how to use the code: it calculates the SpringRank scores for the sample network contained in `../data`.  

- `test_spring_rank_dense_vs_sparse`  is a sample script for testing the efficency using sparse vs dense matrices. Use always sparse matrices when you can! 
