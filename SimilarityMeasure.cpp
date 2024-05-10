// SimilarityMeasure.cpp
// A main reading in chemical simulations and testing the WL-neighborhood similarity measure on it
//
// Richard Restetzki

#include "input_graph.h"
#include "tg.h"

extern const uint32_t seed = 123123;

using namespace std;

int main(int argc, char** argv) {

    // Read the data points from a given file
    if(argc > 1){
        const char* edge_file = "TG_edges.txt";
        tglib::NodeId num_nodes = create_edges(argv[1], argv[2], edge_file);
        //tglib::NodeId num_nodes = create_edges(argv[1], edge_file);
        if(num_nodes != 0){
            tglib::IncidentLists temp_graph(num_nodes);
            tglib::Color_Mapping Color_Mapping;
            //create_TG(Color_Mapping, temp_graph, argv[1], edge_file);
            create_TG(Color_Mapping, temp_graph, argv[1], edge_file, false);
            cout << "Edges in the graph: " << temp_graph.getNumberOfEdges() << "\n";
            temp_graph.initialize_Event(argv[3]);
            temp_graph.colorEdges_VertBased(Color_Mapping, argv[4]);

            //num_nodes = create_edges(argv[2], edge_file);
            //tglib::IncidentLists graph2(num_nodes);
            //create_TG(Color_Mapping, graph2, argv[2], edge_file);
            //cout << "Edges in the graph: " << graph2.getNumberOfEdges() << "\n";

            vector<tglib::IncidentLists *> graphs;
            graphs.emplace_back(&temp_graph);
            //graphs.emplace_back(&graph2);
            cout << "number of graphs: " << graphs.size() << "\n";
            TNWL(Color_Mapping, graphs, 5, 1, true, 500);
            //tglib::IncidentLists graph2 = temp_graph;
            print_gram_of_events(Color_Mapping, temp_graph, 200);

            //print_color_densities(argv[1], Color_Mapping, temp_graph);
            //Color_Mapping.Bucket_collisions();
            //output_similarities(Color_Mapping, temp_graph);
            //cout << "The similarity of both graphs is: " << TNWL_similarity(temp_graph, graph2, Color_Mapping) << "\n";
            //label_analysis(Color_Mapping, temp_graph, graph2);
            //output_rd_similarities(Color_Mapping, temp_graph, 51, 1000);
        }
    }
}