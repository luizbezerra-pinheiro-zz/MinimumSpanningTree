#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <iterator>
#include <mpi.h>

#include "graph.hpp"
#include "edge.hpp"
#include "unionFind.hpp"


#define INFINITY __DBL_MAX__

  //variables
  double* allEdges;
  double* mst;
  void start_communication_measure();
  void end_communication_measure();
  int nProcessors;          /* Number of processors */
  int rank;                 /* MPI rank */
  int root = 0;             /* Rank of root process (which holds the final result) */
  int EProcess;             /* Number of vertices that each process will get */
  int extra_edges;
  double total_distance;    /* Total distance of MST (relevant only in root process) */
  std::list<edge> l;

  FILE * f;                 /* Input file containing graph. */

   double *minEdges;
   double *minEdgesRecv;
   int edges_added_to_mst;   /* Number of edges currently added to MST */

  double startTime;          /* Start time of computation */
  double totalTime;          /* End time of computation */

  /*************************************************************************/
  //Initialise the MPI
  void init_MPI(int argc, char** argv) {

	/* Parsing graph size */
	//f = fopen(argv[1], "rb");
	//fread(&nVertices, sizeof(nVertices), 1, f);

	//starting_vertex = atoi(argv[2]);

	/* Initialize MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* Create 'edge' type (u, v, distance) */
	//MPI_Type_contiguous(3, MPI_INT, &edge);
	//MPI_Type_commit(&edge);

	/* Create reduce operation */
	//MPI_Op_create((MPI_User_function *)reduceMin, 1, &reduce_op);

	//total_comm_time = 0;
  }

    void transform_list(mst){
    edge e;
    for(int i =0 ; i < edges_added_to_mst; i += 3){
        e = new edge(mst[i], mst[i+1], mst[i+3]);
        l.push_back(e);
    }
}

double* transform_graph(graph *g, int E){
    std::list<edge> l = g->getAllEdges();
    double* transformedGraph;

    transformedGraph= (double*)malloc(3*E*sizeof(double));
    int count = 0;
    for (it = l->begin(); it != l->end(); it++){        
        transformedGraph[count*3] = (double)(it->get_source());    
        transformedGraph[count*3+1] = (double)(it->get_target());
        transformedGraph[count*3+2] = it->weight();
        count++;
    }
    return transformedGraph;
}

  void initData(graph* g) {

	//check_conditions();

	//Variables initialisation

	/* Create data structures */
	EProcess = nVertices / nProcessors;
	extra_edges = nVertices - EProcess * nProcessors;

	int i;
	/* Dynamic init of edges*/
	allEdges = transform_graph(g);

	/* Dynamic init of minEdges*/
	if(rank == nProcessors-1){
        minEdges = (double*)malloc((EProcess + extra_edges) * 3*sizeof(double));
        for(int i = 0; i < EProcess + extra_edges; i++)
            minEdges[3*i + 2] = INFINITY;
	}
	
	else{
        minEdges = (double*)malloc(EProcess * 3);
        for(int i = 0; i < EProcess; i++)
            minEdges[3*i + 2] = INFINITY;
	}
	

	if (rank == root) {
		mst = (double *) malloc(3*(nVertices-1) * sizeof(double));
		edges_added_to_mst = 0;
	}

	//inputData();
} 

void GarbageCollector() {

	/* Free memory used by dynamic data structures */
    free(allEdges);
    free(minEdges);
    free(minEdgesRecv);
	/* Shut down MPI */
	MPI_Finalize();
}

  /*************************************************************************/
  
  int main(int argc, char*argv[]){  
      
    //Data check and inicialisation
    const char *infile = argv[1];
    const char *weightFile = argv[2];

    int size;
    std::ifstream inFile(infile);
    std::ifstream weightFile(weightfile);

    fscanf(inFile, "%d", &size);
    
    graph* g = new  graph (size, inFile, weightFile);
    
    int nVertices = g->get_size();   /* Number of vertices in the graph */
    unionFind * parent = new unionFind(nVertices);
    int nEdges;
    //Initialise
    //initialise MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);
	MPI_Status st;
    MPI_Barrier(MPI_COMM_WORLD);
 
    startTime = MPI_Wtime();
    
    //Initialise Data structures
    //getting nb edges
    std::list<edge>::iterator it;
    nEdges =0;
    for (it = l->begin(); it != l->end(); it++)
        nEdges++;
	allEdges = transform_graph(g, nEdges);

	//Variables initialisation
	//divide the edges among the processors
    EProcess = nEdges / nProcessors;
	extra_edges = nEdges - EProcess * nProcessors;

	/* Dynamic init of minEdges*/
	minEdges = (double*)malloc(nVertices * 3*sizeof(double));
        
	if (rank == root) {
		mst = (double *) malloc(3*(nVertices-1) * sizeof(double));
		minEdgesRecv = (double*)malloc((3 * nVertices) * nProcessors*sizeof(double));
        edges_added_to_mst = 0;
	}

    int start = rank * EProcess*3;
    int end = (rank+1)* EProcess*3;
    if(rank == nProcessors-1){
        end += extra_edges*3*;
    }
    
    while(edges_added_to_mst < nVertices-1){
		
        //get min edge for a processor
		for(int i = 0; i < nVertices; i++)
            minEdges[3*i + 2] = INFINITY;

        for(int i = start; i < end;i+=3){
            //is it a valid edge?
            if(parent->find((int)allEdges[i]) != parent->find((int)allEdges[i+1])){
             //is it a better edge?   
                if(allEdges[i+2] < minEdges[(parent->find(i))+2]){ 
                    //refresh
                    minEdges[(parent->find(i))] = parent->find((int)allEdges[i]);
                    minEdges[(parent->find(i))+1] = parent->find((int)allEdges[i+1]);
                    minEdges[(parent->find(i))+2] = allEdges[i+2];
                }
                if(allEdges[i+2] < minEdges[(parent->find(i+1))+2]){ 
                    //refresh
                    minEdges[(parent->find(i+1))] = parent->find((int)allEdges[i+1]);
                    minEdges[(parent->find(i+1))+1] = parent->find((int)allEdges[i]);
                    minEdges[(parent->find(i+1))+2] = allEdges[i+2];
                }
            }
        }
        //here we share ours min edges with root
        MPI_Gather(minEdges, nVertices * 3, MPI_DOUBLE,minEdgesRecv , 3 * V * nProcessors, MPI_DOUBLE, root, MPI_COMM_WORLD);
        
        //We do a manually MPI_reduce, it means, we take the minedge among the min edge s received 
		if(rank == root){
            for(int i = 0; i < nVertices; i++){
                for(int t = 0; t < nProcessors; t++){
                    if(minEdgesRecv[3*i+2 + t*(3*V)] < minEdges[3*i+2]){
                        //update
                        minEdges[3*i] = parent->find((int)minEdgesRecv[3*i + t*(3*V)]);
                        minEdges[3*i + 1] = parent->find((int)minEdgesRecv[3*i + 1 + t*(3*V)]);
                        minEdges[3*i + 2] = minEdgesRecv[3*i + 2 + t*(3*V)];
                    }
                }
            }
            //now minEdges contains the min edge going out of ifro rank == root
            //we broadcast this edge value so that everyone can update your tree
		    MPI_Bcast(minEdges, 3 * V , MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        }
        
		// add new edges to MST
		for (int i = 1; i <= nVertices; i++)
		{
			if (minEdges[i * 3 + 2] != DBL_MAX)
			{
				if (parent->find(minEdges[3*i]) != parent->find(minEdges[3*i+1])
				{
                    //add to the union
                    parent->merge(minEdges[3*i], minEdges[3*i+1]);
                    //add edge to the mst
                    mst[edges_added_to_mst] = minEdges[3*i];
                    mst[edges_added_to_mst+1] = minEdges[3*i+1];
                    mst[edges_added_to_mst+2] = minEdges[3*i+2];
                    //update nb of edges added
                    edges_added_to_mst++;
				}
			}
		}
	}
    //Results
    printf("excution time = %lf\n",MPI_Wtime() - startTime);
    
    transform_list(mst);
	
    GarbageCollector();
    
    return 0;
  }