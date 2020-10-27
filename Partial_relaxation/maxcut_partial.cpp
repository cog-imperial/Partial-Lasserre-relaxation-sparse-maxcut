/*
BSD 3-Clause License
Copyright (c) 2020, Juan S. Campos
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
This code calculates the partial and augmented relaxations as described in
"Partial Lasserre relaxation for sparse Max-Cut" by Campos, J.S., Misener R., and Parpas, P. 2020.
*/


#include <iostream>
#include "fusion.h"
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <typeinfo>
#include <chrono>
#include <math.h> 
#include <algorithm>
#include <vector>
#include <iomanip>
#include <stdlib.h> 
#include <time.h>    
#include <boost/functional/hash.hpp>
#include <unordered_map>
#define BILLION 1e9


using namespace std;
using namespace mosek::fusion;
using namespace monty;

struct timespec my_start;


void read_matrix(const char *, Matrix::t&); 
void read_matrix_sparse_format(const char *, Matrix::t&); 

template <typename Container>
struct container_hash
{
	size_t operator()(Container const& c) const {
		return boost::hash_range(c.begin(), c.end());
	}
};

int main(int argc, char ** argv)
{
 	 // Order of arguments: max_size_cliques (r), no_subsets to add (p), type of heuristic (H1 to H5), name  of the instance 
	char aux_string[300];	
	// Read weight matrix W
	// First row must contain the number of nodes and edges
	Matrix::t W;
	sprintf(aux_string,"graphs/%s.txt",argv[4]);
	read_matrix_sparse_format(aux_string, W);	
    
   
	int n = W->numColumns();		
	int max_size_clique = atoi(argv[1]);// 0,2,3 ...
	int no_subsets = atoi(argv[2]);// 0,1,2,3 ...
	char heuristic [300]; 
	sprintf(heuristic,"%s",argv[3]);// H1, H2, ..., H5
	
	int total_cliques;	
	double time_creating_mosek = 0;
	double time_solving_SDP = 0;
	double time_creating_cliques = 0;
	double time_generating_LB = 0;
	double time_total = 0;
	
	auto start = std::chrono::high_resolution_clock::now();	
	auto start_creating = std::chrono::high_resolution_clock::now();	
	
    // Read maximal cliques and subset generated from matlab
    Matrix::t cliques;
    Matrix::t first_second_order;  
    
   	sprintf(aux_string,"clique_%s_%d_%d_%s.txt",argv[4], max_size_clique, no_subsets,heuristic);
	read_matrix(aux_string, cliques);
    sprintf(aux_string,"clique_aux_%s_%d_%d_%s.txt",argv[4], max_size_clique, no_subsets,heuristic);
	read_matrix(aux_string, first_second_order);
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout<<"time reading W, cliques, and first_second_order: "<<elapsed.count()<<endl;
	
	start = std::chrono::high_resolution_clock::now();
	total_cliques = (cliques->numRows());
	
	vector <int> size_cliques(total_cliques);
	vector <vector<int>> idx_cliques(total_cliques);
	vector <int> size_blk(total_cliques);	
	int total_var_all_blocks = 0;
	for (int i = 0; i < total_cliques; i++)
    {
		int aux_count = 0;
		
	    for(int j=0; j < n; j++)
	    {
			if(cliques->get(i,j)>0)
			{
				aux_count = aux_count + 1;
				idx_cliques[i].push_back(j);
			}
		}
		size_cliques[i] = aux_count;
		if((first_second_order->get(i,0) == 1) | (aux_count == 1))
		{
			size_blk[i] = aux_count + 1;			
		}
		else
		{
			size_blk[i] = (aux_count + 2)*(aux_count + 1)/2 - aux_count;			
		}	
		total_var_all_blocks = total_var_all_blocks + (size_blk[i]*size_blk[i]);
    }
    finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start;
	cout<<"time creating size_cliques and idx_cliques: "<<elapsed.count()<<endl;	
	
    // Generate unordered_map containing the monomials needed for the cliques 
    start = std::chrono::high_resolution_clock::now();
    vector <int> alpha (4,0);
    unordered_map<vector<int>, int, container_hash<vector<int>>> map; 
     
    for(int i=0; i< n; i++)
    {
		alpha[0] = i;
		map.insert({alpha,i});
	} 
	int value_map = n;
	
    for (int i = 0; i < total_cliques; i++)
    {				
		int degree = 2;
		if(first_second_order->get(i,0) == 2)
		{
			degree = 4;
		}
		
		for(int j=2; j<=degree; j++)
		{
			if(j <= size_cliques[i])
			{								
				vector <int> perm(size_cliques[i], 0);
                fill_n(perm.begin(), j, 1);
                int count_aux = 0;
                do
                {
                    for (int k = 0; k < size_cliques[i]; k++)
                    {
                        if (perm[k] == 1)
                        {
                            alpha[count_aux] = idx_cliques[i][k];
                            count_aux = count_aux + 1;
                        }
                    }
                    count_aux = 0;
                    if(map.find(alpha) == map.end())
                    {
                        map.insert({alpha,value_map});
                        value_map = value_map + 1;
                    }
                    fill(alpha.begin(), alpha.end(),0);
                } while (prev_permutation(perm.begin(), perm.end()));
			}			
		}
	}
	
	int total_y = map.size();
	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start;
	cout<<"time creating unordered_map: "<<elapsed.count()<<endl;
          
    // Generating SDP problem: max b*y, st. Aty - C = S, S>=0
    
    // Creating mosek model
    Model::t M = new Model("sdo1"); 
    // Creating SDP blocks and lifted monomials variable and calculating norm_C
    double norm_C = 0;
    Variable::t S[total_cliques];    
    for(int i=0; i<total_cliques; i++)
    {
		sprintf(aux_string,"S_%d",i);			
	    S[i] = M->variable(aux_string, Domain::inPSDCone(size_blk[i]));	 
	    norm_C = max(norm_C, sqrt(double(size_blk[i])));  
	}
	
	norm_C = (max(1.0,norm_C));
	
	// Creating b
	auto bsubi = new_array_ptr<int, 1>(n*(n-1)*0.5 + n);
	auto bsubj = new_array_ptr<int, 1>(n*(n-1)*0.5 + n);
	auto bcof  = new_array_ptr<double, 1>(n*(n-1)*0.5 + n);
	double obj_constant = 0;
	
	fill(alpha.begin(), alpha.end(),0);
	double norm_b = 0;
	vector<double> pert(n);
	srand(time(NULL));
	double perturbation = 0;//////////////////////////////
	for(int i = 0; i < n; i++)
	{
		pert[i] = perturbation*(double(rand())/RAND_MAX);
		norm_b = norm_b + pert[i];
		
		for(int j = i+1; j < n; j++)
		{
			norm_b = 0.25*4*W->get(i,j)*W->get(i,j) + norm_b;
		}
	}	
	norm_b = max(1.0,sqrt(norm_b));	
	int aux_count = 0;	
	for(int i = 0; i < n; i++)
	{
		alpha[0] = i;
		alpha[1] = 0;
		(*bsubi)[aux_count] = map[alpha];
		(*bsubj)[aux_count] = 0;
		(*bcof)[aux_count] = pert[i]/norm_b;
		aux_count = aux_count + 1;		
		for(int j = i+1; j < n; j++)
		{
			alpha[0] = i;
			alpha[1] = j;
			(*bsubi)[aux_count] = map[alpha];
			(*bsubj)[aux_count] = 0;
			(*bcof)[aux_count] = -0.25*2*(W->get(i,j))/norm_b;
			obj_constant = obj_constant + W->get(i,j);
			aux_count = aux_count + 1;			
		}
	}	
	obj_constant = 2*obj_constant*0.25;	
	Matrix::t b = Matrix::sparse(total_y, 1, bsubi, bsubj, bcof); 
	
	// Generating A 
    start = std::chrono::high_resolution_clock::now();
    int total_constraints = total_var_all_blocks;
    
	Expression::t Ax_b;
	int count_const;
	vector<double> norm_A(total_cliques);	
	for(int i=0; i < total_cliques; i++)
    {   
		int size_vec_blk_i = size_blk[i]*size_blk[i];
		int size_blk_i = size_blk[i];
		int size_clique_i = size_cliques[i];
		int total_constraints = size_blk_i*(size_blk_i - 1);
		count_const = 0;
		auto msubi = new_array_ptr<int, 1>(total_constraints);
	    auto msubj = new_array_ptr<int, 1>(total_constraints);
	    auto mcof  = new_array_ptr<double, 1>(total_constraints);
	    norm_A[i] = max(1.0, sqrt(double(total_constraints)));
	    
		fill(alpha.begin(), alpha.end(),0);
		double const_coef = 1/norm_A[i];
		for(int j = 1; j <= size_clique_i; j++)
		{			
			alpha[0] = idx_cliques[i][j-1];			
			(*msubi)[count_const] = j;
			(*msubj)[count_const] = map[alpha];				
			(*mcof)[count_const] = const_coef;
			(*msubi)[count_const+1] = j*(size_blk_i);
			(*msubj)[count_const+1] = map[alpha];
			(*mcof)[count_const+1] = const_coef;
			count_const = count_const + 2;
		} 	
		
		for(int j = 0; j < size_clique_i; j++)
		{			
			for(int k = j + 1; k < size_cliques[i]; k++)
			{
				alpha[0] = idx_cliques[i][j];
				alpha[1] = idx_cliques[i][k];
				(*msubi)[count_const] = (j+1)*(size_blk[i]) + k + 1;
				(*msubj)[count_const] = map[alpha];				
				(*mcof)[count_const] = const_coef;
				(*msubi)[count_const+1] = (k+1)*(size_blk[i]) + j + 1;
				(*msubj)[count_const+1] = map[alpha];
				(*mcof)[count_const+1] = const_coef;
				count_const = count_const + 2;
			}
		} 	
				
		if((first_second_order->get(i,0) == 2) & (size_clique_i > 1))
		{
			fill(alpha.begin(), alpha.end(),0);
			int count_aux = 1;
		    for(int k = 0; k < size_clique_i; k++)
		    {
				for(int l = k+1; l < size_clique_i; l++)
				{			        
			        alpha[0] = idx_cliques[i][k];
				    alpha[1] = idx_cliques[i][l];
				    (*msubi)[count_const] = (size_clique_i) + count_aux;
				    (*msubj)[count_const] = map[alpha];				
				    (*mcof)[count_const] = const_coef;
				    (*msubi)[count_const+1] = ((size_clique_i+count_aux)*size_blk_i);
				    (*msubj)[count_const+1] = map[alpha];
				    (*mcof)[count_const+1] = const_coef;
				    count_const = count_const + 2;
				    count_aux = count_aux + 1;
			    }
			}
		  
		   		
		    for(int j = 0; j < size_clique_i; j++)
		    {			
				count_aux = 1;
			    for(int k = 0; k < size_clique_i-1; k++)
			    {
					for(int l = k+1; l < size_clique_i; l++)
					{
				        fill(alpha.begin(), alpha.end(),0);
				        if(j == k)
				        {
							alpha[0] = idx_cliques[i][l];
						}
						else if(j == l)
						{
							alpha[0] = idx_cliques[i][k];
						}
						else if(j<k)
						{
							alpha[0] = idx_cliques[i][j];
							alpha[1] = idx_cliques[i][k];
							alpha[2] = idx_cliques[i][l];
						}  
						else if (j<l)
						{
							alpha[0] = idx_cliques[i][k];
							alpha[1] = idx_cliques[i][j];
							alpha[2] = idx_cliques[i][l];
						} 
						else 
						{
							alpha[0] = idx_cliques[i][k];
							alpha[1] = idx_cliques[i][l];
							alpha[2] = idx_cliques[i][j];
						} 
				        
				        (*msubi)[count_const] = (j+1)*size_blk_i + size_clique_i + count_aux;				       
				        (*msubj)[count_const] = map[alpha];		
				        (*mcof)[count_const] = const_coef;				       
				        (*msubi)[count_const+1] = (size_clique_i + count_aux)*size_blk_i + j + 1;
				        (*msubj)[count_const+1] = map[alpha];
				        (*mcof)[count_const+1] = const_coef;
				        count_const = count_const + 2;
				        count_aux = count_aux + 1;				        
			        }
				}
		    } 	
		    
		   
		    count_aux = 0;
		    int count_aux2;
		    for(int j = 0; j < size_clique_i-1; j++)
		    {			
				for(int k = j+1; k < size_clique_i; k++)
			    {
					count_aux = count_aux + 1;
					count_aux2 = 0;
					for(int l = 0; l < size_clique_i-1; l++)
					{
						for(int m = l+1; m < size_clique_i; m++)
					    {
							count_aux2 = count_aux2 + 1;			
							if((j==l) & (k==m))
							{
								continue;
							}
							fill(alpha.begin(), alpha.end(),0);
				            if(j == l)
				            {
								if(k<m)
								{
								    alpha[0] = idx_cliques[i][k];
								    alpha[1] = idx_cliques[i][m];
								}
								else
								{
									alpha[0] = idx_cliques[i][m];
								    alpha[1] = idx_cliques[i][k];
								}							    
						    }
						    else if(k == m)
						    {
							    if(j<l)
								{
								    alpha[0] = idx_cliques[i][j];
								    alpha[1] = idx_cliques[i][l];
								}
								else
								{
									alpha[0] = idx_cliques[i][l];
								    alpha[1] = idx_cliques[i][j];
								}
						    }
						    else if(j == m)
						    {							
							    alpha[0] = idx_cliques[i][l];
							    alpha[1] = idx_cliques[i][k];
								
						    }
						    else if(k == l)
						    {							
							    alpha[0] = idx_cliques[i][j];
							    alpha[1] = idx_cliques[i][m];
								
						    }
						    else if(k<l)
						    {
							    alpha[0] = idx_cliques[i][j];
							    alpha[1] = idx_cliques[i][k];
							    alpha[2] = idx_cliques[i][l];
							    alpha[3] = idx_cliques[i][m];							    
						    }  
						    else if(m<j)
						    {
							    alpha[0] = idx_cliques[i][l];
							    alpha[1] = idx_cliques[i][m];
							    alpha[2] = idx_cliques[i][j];
							    alpha[3] = idx_cliques[i][k];							    
						    }
						    else if((j<l) & (m<k))
						    {
							    alpha[0] = idx_cliques[i][j];
							    alpha[1] = idx_cliques[i][l];
							    alpha[2] = idx_cliques[i][m];
							    alpha[3] = idx_cliques[i][k];							    
						    } 
						    else if((j<l) & (k<m))
						    {
							    alpha[0] = idx_cliques[i][j];
							    alpha[1] = idx_cliques[i][l];
							    alpha[2] = idx_cliques[i][k];
							    alpha[3] = idx_cliques[i][m];							    
						    } 
						    else if((l<j) & (m<k))
						    {
							    alpha[0] = idx_cliques[i][l];
							    alpha[1] = idx_cliques[i][j];
							    alpha[2] = idx_cliques[i][m];
							    alpha[3] = idx_cliques[i][k];							    
						    }  
						    else if((l<j) & (k<m))
						    {
							    alpha[0] = idx_cliques[i][l];
							    alpha[1] = idx_cliques[i][j];
							    alpha[2] = idx_cliques[i][k];
							    alpha[3] = idx_cliques[i][m];							    
						    } 
				             
				            (*msubi)[count_const] = (size_clique_i + count_aux)*size_blk_i + size_clique_i + count_aux2;				           
				            (*msubj)[count_const] = map[alpha];	
				            (*mcof)[count_const] = const_coef;		                   
				            count_const = count_const + 1;				            
						}
			        }
				}
		    }
	    }	    
	    
	   
	    auto msubid = new_array_ptr<int, 1>(size_blk_i);
	    auto msubjd = new_array_ptr<int, 1>(size_blk_i);
	    auto mcofd  = new_array_ptr<double, 1>(size_blk_i);
	    for (int j = 0; j < size_blk[i]; j++)
	    {
			(*msubid)[j] = (j*size_blk_i)+j;
			(*msubjd)[j] = 0;
			(*mcofd)[j] = 1.0/(norm_C*norm_A[i]);
		}
	 
	    Matrix::t A = Matrix::sparse(total_y,size_vec_blk_i, msubj, msubi, mcof);  
	    
	    if (i==0)
	    {
			Ax_b = Expr::mul(A, Expr::reshape(Expr::mul(S[i],-1),size_vec_blk_i));
		}
		else
		{
	        Ax_b = Expr::add(Ax_b, Expr::mul(A, Expr::reshape(Expr::mul(S[i],-1),size_vec_blk_i)));
	    }
	   	
	}
	Ax_b = Expr::sub(Ax_b,b);
	Constraint::t con = M->constraint(Ax_b, Domain::equalsTo(0.0));
	cout<<"size map: "<<map.size()<<endl;	
	cout<<"Total_constraints: "<<b->numRows()<<endl;
	
	
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout<<"Time creating constraints: "<<elapsed.count()<<endl;
    
    
    // Creating objective function
    
    Expression::t obj = Expr::constTerm(0);
    for (int i=0;i<total_cliques;i++)
    {
		int size_blk_i = size_blk[i];
	    auto msubid = new_array_ptr<int, 1>(size_blk_i);
	    auto msubjd = new_array_ptr<int, 1>(size_blk_i);
	    auto mcofd  = new_array_ptr<double, 1>(size_blk_i);	    
	    for (int j = 0; j < size_blk_i; j++)
	    {
			(*msubid)[j] = j;
			(*msubjd)[j] = j;
			(*mcofd)[j] = 1.0/(norm_C*norm_A[i]);
		}
		
		Matrix::t D_ones = Matrix::sparse(size_blk_i,size_blk_i, msubid, msubjd, mcofd);		
		obj = Expr::add(obj, Expr::dot(S[i],D_ones));	
		
	}
	
    auto finish_creating = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_creating = finish_creating - start_creating;
    //Solve SDP
  	M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );			
	M->objective(ObjectiveSense::Minimize, obj);		
	M->setSolverParam("optimizerMaxTime", 10800);
	
	auto start_solving_SDP = std::chrono::high_resolution_clock::now();	
	
	M->solve();	 			
	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start_solving_SDP;	
	time_solving_SDP = elapsed.count();
	
		
	std::cout << std::fixed;
    std::cout << std::setprecision(4);
	cout<<"time creating problem: "<<elapsed_creating.count()<<endl;
	cout<<"time solving SDP: "<<elapsed.count()<<endl;
	cout<<"Obj_value_partial_relaxation: "<<M->primalObjValue()*norm_b*norm_C + obj_constant<<endl;
	return 0;
}


void read_matrix(const char *filename, Matrix::t& result)  
{
	int row, col;    
	int i = 0,j = 0;	
    ifstream fin(filename); 
	if (!fin.is_open())
	{
		cout << "Error opening file: "<<filename;
		exit(1);
	}

	string line, aux;
	// First line must have in the first two positions the number of rows and columns
	getline(fin, line);
	stringstream stream(line);
	getline(stream,aux,',');
	row = stoi(aux);		
	getline(stream,aux,',');
	col = stoi(aux);
	
	auto msubi = new_array_ptr<int, 1>(row*col);
	auto msubj = new_array_ptr<int, 1>(row*col);
	auto mcof  = new_array_ptr<double, 1>(row*col);
	int aux2 = 0;
	while (getline(fin, line))
   	{	
		stringstream stream(line);
		while(stream.good())
		{
			getline(stream,aux,',');
			(*msubi)[aux2] = i;
			(*msubj)[aux2] = j;
			(*mcof)[aux2] = stod(aux);									
			j++;
			aux2 = aux2 + 1;
			if (j == col)
			{
				i = i+1;
				j = 0;
			}
		}
	}
	result = Matrix::sparse(row, col, msubi, msubj, mcof);
}

void read_matrix_sparse_format(const char *filename, Matrix::t& result)  
{
	// reads symmetric  matrix in sparse format given as 
	// n nnz
	// i_1 j_1 val_1
	// i_2 j_2 val_2
	// ...
	// i1_nnz j1_nnz val_nnz
	// where n is the dimension of the matrix, and nnz the non-zero elements.
	// It gives only the upper triangular values.
	
	int row, col;    	
	string val, i, j;	
    ifstream fin(filename); 
	if (!fin.is_open())
	{
		cout << "Error opening file: "<<filename;
		exit(1);
	}

	string line, aux;
	// First line must have the dimension of the square matrix and the number of non-zeros
	getline(fin, line);
	stringstream stream(line);
	getline(stream,aux,' ');
	row = stoi(aux);
	col = row;
	getline(stream,aux,' ');
	int nnz = stoi(aux);
	auto msubi = new_array_ptr<int, 1>(2*nnz);
	auto msubj = new_array_ptr<int, 1>(2*nnz);
	auto mcof  = new_array_ptr<double, 1>(2*nnz);
	int aux2 = 0;
	while(fin >> i >> j >> val)
	{
		(*msubi)[aux2] = stod(i)-1;
		(*msubj)[aux2] = stod(j)-1;
		(*msubi)[aux2+1] = stod(j)-1;	
		(*msubj)[aux2+1] = stod(i)-1;
		(*mcof) [aux2] = stod(val);
		(*mcof) [aux2+1] = stod(val);
		aux2 = aux2+2;			
	}
	
	result = Matrix::sparse(row, col, msubi, msubj, mcof);
}
