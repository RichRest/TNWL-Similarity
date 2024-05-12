// main.cpp
// A main reading in temporal graphs and performing the TNWL labeling, kernel and similarity measure
//
// Richard Restetzki


#include "input_graph.h"
#include "tg.h"


// Random seed used within the plots of the thesis (used in the functionality output_rd_similarities() )
extern const uint32_t seed = 123123;

// Set up the container for the Voronoi Tessellations (used for xyz data)
double x_min=0,x_max=18.801999;
double y_min=0,y_max=18.801999;
double z_min=0,z_max=18.801999;
extern const bool x_periodic=true, y_periodic=true, z_periodic=true;
// Specify the number of blocks the container is split up in (relevant for runtime optimization on xyz data)
extern const int n_x=6,n_y=6,n_z=6;
// Restrict a simulation to the first n timesteps (instance size control for xyz data)
int max_sim_time=1000;
// Restrict the out degree of vertices in a temporal delaunay graph (instance size and randomness control for xyz data)
extern const int max_sim_out_deg=6;


using namespace std;
using namespace tglib;

int main(int argc, char** argv) {

    // Specify the input data for the graph and the respective format
    string data_1 = "MIM_OAc_n";
    string input_format_1 = "xyz";
    string data_2 = "MIM_OAc_i";
    string input_format_2 = "xyz";
    // In Color_Mapping the compression function g from (TNWL3) is stored
    Color_Mapping Color_Mapping_1;
    // Create the temporal graphs
    IncidentLists graph_1;
    create_TG(Color_Mapping_1, graph_1, data_1, input_format_1);
    IncidentLists graph_2;
    create_TG(Color_Mapping_1, graph_2, data_2, input_format_2);
    // Perform TNWL and a label analysis on the graphs
    TNWL(Color_Mapping_1, {&graph_1, &graph_2}, 2, 1, false, 1000);
    label_analysis(Color_Mapping_1, graph_1, data_1, graph_2, data_2);
    graph_1.clear();
    graph_2.clear();
    Color_Mapping_1.clear();


    // Output random similarities
    max_sim_time = 1000;
    create_TG(Color_Mapping_1, graph_1, data_1, input_format_1);
    output_rd_similarities(Color_Mapping_1, graph_1, 51, 1000);
    graph_1.clear();
    Color_Mapping_1.clear();


    // Output interval similarities
    max_sim_time = numeric_limits<int>::max();
    x_max = 15.64;
    y_max = 15.64;
    z_max = 15.64;
    string data_3 = "H2O";
    string input_format_3 = "xyz";
    create_TG(Color_Mapping_1, graph_1, data_3, input_format_3);
    TNWL(Color_Mapping_1, {&graph_1}, 10, 1, false, 1000);
    output_similarities(Color_Mapping_1, graph_1);
    graph_1.clear();
    Color_Mapping_1.clear();


    // Output gram matrix of a TUdataset
    string data_4 = "infectious_ct1";
    string input_format_4 = "TU";
    int nr_graphs = 200;
    // The temporal graphs are represented through the event label on the edges, "graph_1" is the graph sum of all the temporal graphs
    create_TG(Color_Mapping_1, graph_1, data_4, input_format_4);
    TNWL(Color_Mapping_1, {&graph_1}, 10, 2, true, 50);
    print_gram_of_events(Color_Mapping_1, graph_1, nr_graphs);
}