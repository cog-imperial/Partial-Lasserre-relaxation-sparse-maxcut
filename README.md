# Partial-Lasserre-relaxation-sparse-maxcut
Implements the partial and augmented relaxations described in "Partial Lasserre relaxation for sparse Max-Cut", by Campos et al. (2020).
When using this code for any publications please reference this software package as:
```
@article{Campos2020Partial,
  author = {Juan S. Campos, Ruth Misener, and Panos Parpas},
  title = {{Partial Lasserre relaxation for sparse Max-Cut}},
  note = {Available on http://www.optimization-online.org},
  year = {2020}
  month = {10} 
}
```

## Description
Given a sparse graph, the partial and augmented relaxations strengthen the first order sparse Lasserre relaxation of the Max-Cut problem described in "Sums of squares and semidefinite program relaxations for polynomial optimization problems with structured sparsity" by Waki et al. (2006), by including some of the second order positive semidefinite constraints. The tightness of the partial relaxation is determined by limiting the size of the second order positive semidefinite constraints that can be included. For more details see: http://www.optimization-online.org/DB_FILE/2020/10/8063.pdf

## Requirements
The code has been tested in a Linux machine with Ubuntu 18.04 and using:
- MATLAB 2018a
- MOSEK 8.1 Fusion API for C++
- Boost 1.58

## Files and folders
1. Makefile
2. gen_cliques.m: MATLAB file that calculates the maximal cliques of the graph.\
This code was taken from SparsePOP version 3.01 by Hayato Waki, Sunyoung Kim, Masakazu Kojima, Masakazu Muramatsu, Hiroshi Sugimoto and Makoto Yamashita. Specifically, it corresponds to lines 97 to 135 of the file genClique.m contained in SparsePOP301/subPrograms/Mfiles. \
Availability: https://sparsepop.sourceforge.io/. \
<ins>**Refer to the License in the file SparsePOP301/readme.txt for usage of this function.**</ins> 
3. select_cliques_and_subsets.m: MATLAB file that calculates the final set of maximal cliques and subsets that will be included in the partial or augmented relaxations, as well as the order of the positive semidefinite constraint (i.e., if they must be included as 1st or 2nd order positive semidefinite constraints).
4. maxcut_partial.cpp: implements and solves the partial or augmented relaxations using MOSEK.
5. graphs.tar.gz: compressed folder containing all the graphs used for the numerical experiments in Campos et al. (2020). The format of each .txt file is: first row has the number of nodes of the graph and number of edges, row k has the 3 elements (i,j,val) corresponding to the nonzero weight val for the edge between node i and j. Three types of graphs:
	* Randoms graphs: created using the SparsePOP version 3.01 function randomUnconst.m. The name of the graph has the form W_size_u_i_w.txt, where size corresponds to the number of nodes of the graph (300 or 500), u is an upper bound for the size of the maximal cliques of the randomly created chordal graphs used to construct the graph (4,6,8,10), w is the type of weight used (1 for graphs with nonzero weights in {-1, 1}, 2 for weights in {1, 2, ..., 10} and 3 for weights in {-10, -9, ..., -1, 1, 2, ..., 10}), and i for the total number of random graphs created of class size, u, w (10 instances in total from 0 to 9). See Section 4.1 in Campos et al. (2020) for more details.
	* Physics graphs: toroidal grid graphs (see Frauke Liers (2004)). These graphs have names starting with t2g or t3g (e.g., t2g10_5555.txt) and were taken from the Max-Cut instances from applications in statistical physics of the Biq Mac library (http://biqmac.aau.at/biqmaclib.html accessed 20/03/2020).
	* Sparse random graphs from CS-TSSOS: random sparse graphs taken from https://wangjie212.github.io/jiewang/code.html (accessed 25/06/2020). These graphs have names g20.txt, g40.txt, ..., g180.txt, g200.txt (see Wang et al. (2020)).
6. graphs: folder where the graphs files must be saved to be used by select_cliques_and_subsets.m and the executable maxcut_partial.
	

## Installation
The installation described here is for a Linux machine along with MOSEK 8.1 Fusion API for C++ and Boost 1.58. For a different operating system and/or MOSEK/Boost versions, it might be necessary to make further changes to the Makefile.
1. Download the folder Partial_relaxation.
1. Replace /path/to/Mosek on the second line of the Makefile by the path to the location of MOSEK on your machine.
2. From the terminal compile the file maxcut_partial.cpp using **make maxcut_partial** 

## Instructions
* To use the partial or augmented relaxations on a sparse graph W, 3 parameters must be selected (see Section 3.1 in Campos et al. (2020)):	
1. r: Max size of the maximal cliques to be included as 2nd order SDP constraints. 
2. p: For the augmented relaxation, this paramenter corresponds to the number of subsets to be included as 2nd order SDP constraints for the maximal cliques with more than r elements. Set this parameter to 0 to use only the partial relaxation.
3. H: Type of heuristic to select the p subsets for the augmented relaxation. Set this parameter to H1 to use only the partial relaxation (the code will ignore this parameter when p=0, i.e., when the partial relaxation is used).

* After selecting these parameters the calculation of the partial relaxation is done in 2 steps:
1. Calculate the maximal cliques of the graph, and the p subsets for each maximal clique larger than r. This is done by using the MATLAB function select_cliques_and_subsets.m (this function uses gen_cliques.m to calculate the maximal cliques following the code in SparsePOP version 3.01 by Hayato Waki, Sunyoung Kim, Masakazu Kojima, Masakazu Muramatsu, Hiroshi Sugimoto and Makoto Yamashita).\
This function will create 2 .txt files: 
	* clique_W_r_p_H.txt. File with the maximal cliques (and subsets if p>0): each row contains a maximal clique (or subset) with 1 in position i if the ith node is part of the maximal clique (or subset) or zero otherwise.
	* clique_aux_W_r_p_H.txt. File with a column vector with 1 in position i if the set in position i in clique_W_r_p_H.txt needs to be included as 1st order SDP constraint, or 2 if needs to be included as 2nd order.
	
2. After the two .txt files are created, the executable maxcut_partial is used to calculate the partial relaxation: ./maxcut_partial r p H W

### Example
To calculate the augmented relaxation for the graph W_300_4_0_1.txt in the folder Partial_relaxation/graphs, with r=5, p=1, H=H5:
* Use MATLAB to run the function: select_cliques_and_subsets('W_300_4_0_1', 5, 1, 'H5'); 
* From the terminal use: ./maxcut_partial 5 1 H5 W_300_4_0_1

To calculate the partial relaxation for the same graph and r then set r=5, p=0, H=H1:
* Use MATLAB to run the function: select_cliques_and_subsets('W_300_4_0_1', 5, 0, 'H1'); 
* From the terminal use: ./maxcut_partial 5 0 H1 W_300_4_0_1

	
**IMPORTANT. It is assumed that the graphs are contained in the folder Partial_relaxation/graphs, with extension .txt and have the format described in numeral 5 of the section Files and folders.**

## Author
Juan S. Campos

## License
All the codes contained in this package, with the exception of gen_cliques.m, are released under the BSD 3-Clause License. Please refer to the LICENSE file for details.\
**For the usage of gen_cliques.m please refer to the License in the file SparsePOP301/readme.txt, available in https://sparsepop.sourceforge.io/.**

## Acknowledgements
This work was funded by an Engineering & Physical Sciences Research Council Research Fellowship [Grant Number EP/P016871/1].

## References
- Juan S. Campos, Ruth Misener, and Panos Parpas. Partial Lasserre relaxation for sparse Max-Cut". http://www.optimization-online.org/DB_FILE/2020/10/8063.pdf, 2020.
- Frauke Liers. Contributions to Determining Exact Ground States Of Ising Spin Glasses And To Their Physics. PhD thesis, Verlag nicht ermittelbar, 2004.
- Frauke Liers, Michael Junger, Gerhard Reinelt, and Giovanni Rinaldi. Computing exact ground states of hard Ising spin glass problems by branch-and-cut. New optimization algorithms in physics, 50(47-68):6, 2004.
- Hayato Waki, Sunyoung Kim, Masakazu Kojima, and Masakazu Muramatsu. Sums of squares and semidefinite program relaxations for polynomial optimization problems with
structured sparsity. SIAM Journal on Optimization, 17(1):218–242, 2006.
- Jie Wang, Victor Magron, Jean B Lasserre, and Ngoc Hoang Anh Mai. CS-TSSOS: Correlative and term sparsity for large-scale polynomial optimization. arXiv preprint
arXiv:2005.02828, 2020.
