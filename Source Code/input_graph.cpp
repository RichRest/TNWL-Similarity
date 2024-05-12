// chem_input.cpp
// Source file defining the functionalities to process chemical simulations
//
// Richard Restetzki

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include "voro++/voro++.hh"
#include "tg.h"
#include "input_graph.h"

using namespace std;
using namespace voro;

int sgn(double val){
    return (0<val) - (val<0);
}

double eucl_dist(vector<double>& point_1, vector<double>& point_2){
    double squared_dist = 0;
    if(point_1.size()!=point_2.size()){
        cout << "ERROR euclidean distance of two points with differing dimensions can't be computed!\n";
    }
    else{
        for(int i=0; i<point_1.size(); i++){
            squared_dist += pow(point_1[i]-point_2[i],2);
        }
    }
    return sqrt(squared_dist);
}

// Performs a Voronoi Tesselation on the spatial positions read from the "read_file" and writes the edges corresponding to the Delaunay graph in (v,w,t) format to the "output_file"
int Tesselate_slice(ifstream& read_file, int partcount, int timepoint, FILE* output_file, size_t degree_limit, double radial_limit) {
    // Create a container with the geometry given above
    container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,x_periodic,y_periodic,z_periodic,4);
    istringstream linestream;
    string line, word;
    int partindex, argcount, wordcount, linecount = 0;
    vector<vector<double>> Nodes_pos(partcount, vector<double>(3));
    vector<Voro_neighbor> neighbors_dist;

    // Loops over all lines of the current timeslice
    while( (read_file.peek() != EOF) && linecount < partcount){
        getline(read_file, line);
        argcount = 0;
        wordcount = 0;
        linestream.clear();
        linestream.str(line);

        // Loops over all the words within the current line (stored in line)
        while(getline(linestream, word, ' ')){
            // Exclude empty words occurring between duplicate spaces
            if(!word.empty()){
                if(wordcount==0){
                    // Set the index of the particle from this line
                    partindex = linecount;
                    //partindex = stoi(word);
                }
                else if( (wordcount==1) || (argcount>0&&argcount<3) ){
                    // Store the coordinates for the particle
                    Nodes_pos[linecount][argcount] = stod(word);
                    // If the xyz file holding the simulation is created such that molecules stick together over a periodic boundary, then it might happen that
                    // the spatial position of an atom is outside the boundaries. In this case, those atoms need to be folded back into the box which is done
                    // with the switch on the current spatial dimension that is read from "linestream" if the value is outside the border
                    switch (argcount) {
                        case 0:
                            if (Nodes_pos[linecount][argcount] < x_min) Nodes_pos[linecount][argcount] += x_max - x_min;
                            else if (Nodes_pos[linecount][argcount] > x_max)
                                Nodes_pos[linecount][argcount] -= x_max - x_min;
                            break;
                        case 1:
                            if (Nodes_pos[linecount][argcount] < y_min) Nodes_pos[linecount][argcount] += y_max - y_min;
                            else if (Nodes_pos[linecount][argcount] > y_max)
                                Nodes_pos[linecount][argcount] -= y_max - y_min;
                            break;
                        case 2:
                            if (Nodes_pos[linecount][argcount] < z_min) Nodes_pos[linecount][argcount] += z_max - z_min;
                            else if (Nodes_pos[linecount][argcount] > z_max)
                                Nodes_pos[linecount][argcount] -= z_max - z_min;
                            break;
                        default:
                            break;
                    }
                    argcount++;
                }
                wordcount++;
            }
        }
        con.put(partindex, Nodes_pos[linecount][0], Nodes_pos[linecount][1], Nodes_pos[linecount][2]);
        linecount++;
    }
    if (linecount != partcount){
        cout << "Invalid data format!\n";
        return -1;
    }
    else{
        // Output connecting lines between neighboring particles
        int curr_particle, nr_edges=0;
        double dim_x = x_max-x_min, dim_y = y_max-y_min, dim_z = z_max-z_min;
        double min_dim = min({dim_x,dim_y,dim_z});
        double x,y,z,r;
        vector<double> point;
        vector<int> neighbors;
        c_loop_all elt_in(con);
        voronoicell_neighbor c;

        // Loop over all particles (curr_particle) and its neighbors (neighbor_id)
        // This will create every edge twice, once in each direction which is intended
        if(elt_in.start()) do if (con.compute_cell(c, elt_in))  {
                    neighbors_dist.clear();
                    elt_in.pos(curr_particle,x,y,z,r);
                    // ATTENTION: !Testing showed that for some particles the neighbors vector might hold duplicates!
                    // This can happen because of the periodic boundary which is computed using
                    // copies of vertices outside the actual boundary.
                    // For large enough simulations it is very unlikely or even impossible that this
                    // happens, but the code takes care of it anyway when reading the data from the file.
                    c.neighbors(neighbors);
                    for (auto neighbor_id:neighbors){
                        point = {x,y,z};
                        double dist_to_neighbor = eucl_dist(point, Nodes_pos[neighbor_id]);
                        // The Euclidean distance might not be correct if the boundaries are periodic this is checked
                        if(dist_to_neighbor>min_dim/2){
                            // Only in this case the distance over a periodic border could be less than the direct distance
                            double periodic_neigh_dist;
                            for(int i=1; i<8; i++){
                                // This switch evaluates all due to periodicity possible foldings and checks if they result in a shorter distance.
                                switch (i) {
                                    case 1:
                                        if(z_periodic){
                                            point = {x,y,z + sgn(Nodes_pos[neighbor_id][2]-z)*dim_z};
                                            periodic_neigh_dist = eucl_dist(point, Nodes_pos[neighbor_id]);
                                            if(dist_to_neighbor>periodic_neigh_dist){
                                                dist_to_neighbor = periodic_neigh_dist;
                                            }
                                        }
                                        break;
                                    case 2:
                                        if(y_periodic) {
                                            point = {x,y + sgn(Nodes_pos[neighbor_id][1]-y)*dim_y,z};
                                            periodic_neigh_dist = eucl_dist(point, Nodes_pos[neighbor_id]);
                                            if(dist_to_neighbor>periodic_neigh_dist){
                                                dist_to_neighbor = periodic_neigh_dist;
                                            }
                                        }
                                        break;
                                    case 3:
                                        if(y_periodic && z_periodic) {
                                            point = {x,y + sgn(Nodes_pos[neighbor_id][1]-y)*dim_y,z + sgn(Nodes_pos[neighbor_id][2]-z)*dim_z};
                                            periodic_neigh_dist = eucl_dist(point, Nodes_pos[neighbor_id]);
                                            if(dist_to_neighbor>periodic_neigh_dist){
                                                dist_to_neighbor = periodic_neigh_dist;
                                            }
                                        }
                                        break;
                                    case 4:
                                        if(x_periodic) {
                                            point = {x + sgn(Nodes_pos[neighbor_id][0]-x)*dim_x,y,z};
                                            periodic_neigh_dist = eucl_dist(point, Nodes_pos[neighbor_id]);
                                            if(dist_to_neighbor>periodic_neigh_dist){
                                                dist_to_neighbor = periodic_neigh_dist;
                                            }
                                        }
                                        break;
                                    case 5:
                                        if(x_periodic && z_periodic) {
                                            point = {x + sgn(Nodes_pos[neighbor_id][0]-x)*dim_x,y,z + sgn(Nodes_pos[neighbor_id][2]-z)*dim_z};
                                            periodic_neigh_dist = eucl_dist(point, Nodes_pos[neighbor_id]);
                                            if(dist_to_neighbor>periodic_neigh_dist){
                                                dist_to_neighbor = periodic_neigh_dist;
                                            }
                                        }
                                        break;
                                    case 6:
                                        if(x_periodic && y_periodic) {
                                            point = {x + sgn(Nodes_pos[neighbor_id][0]-x)*dim_x,y + sgn(Nodes_pos[neighbor_id][1]-y)*dim_y,z};
                                            periodic_neigh_dist = eucl_dist(point, Nodes_pos[neighbor_id]);
                                            if(dist_to_neighbor>periodic_neigh_dist){
                                                dist_to_neighbor = periodic_neigh_dist;
                                            }
                                        }
                                        break;
                                    case 7:
                                        if(x_periodic && y_periodic && z_periodic) {
                                            point = {x + sgn(Nodes_pos[neighbor_id][0]-x)*dim_x,y + sgn(Nodes_pos[neighbor_id][1]-y)*dim_y,z + sgn(Nodes_pos[neighbor_id][2]-z)*dim_z};
                                            periodic_neigh_dist = eucl_dist(point, Nodes_pos[neighbor_id]);
                                            if(dist_to_neighbor>periodic_neigh_dist){
                                                dist_to_neighbor = periodic_neigh_dist;
                                            }
                                        }
                                        break;
                                    default:
                                        break;
                                }
                            }
                        }
                        if(dist_to_neighbor<=radial_limit){
                            neighbors_dist.push_back({neighbor_id, dist_to_neighbor});
                        }
                    }
                    if(neighbors_dist.size()>degree_limit){
                        // Sort ensures that only the closest "degree_limit" neighbors get connected
                        sort(neighbors_dist.begin(),neighbors_dist.end());
                    }
                    // Write the first "degree_limit" edges into the edge file in (v,w,t) format
                    for (int i=0; i < min(neighbors_dist.size(), degree_limit); i++){
                        fprintf(output_file,"%d %d %d\n",curr_particle, neighbors_dist[i].partindex, timepoint);
                        nr_edges++;
                    }
                } while (elt_in.inc());
        return nr_edges;
    }
}

// Takes the edges in (v,w,t) format and adds them to the temporal graph
void initialize_edges(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, const string& edge_file, bool time_sorted) {
    ifstream myfile;
    tglib::Time time, min_time=std::numeric_limits<tglib::Time>::max(), max_time=std::numeric_limits<tglib::Time>::min();
    myfile.open(edge_file.c_str());
    if(myfile.is_open()) {
        string line, word;
        int wordcount;
        tglib::NodeId inputedge[2];
        istringstream linestream;
        tglib::EdgeId label = -1;
        while(myfile.peek() != EOF) {
            getline(myfile, line);
            linestream.clear();
            linestream.str(line);
            wordcount = 0;
            while(getline(linestream, word, ' ')){
                // Loops over all the words within the current line (stored in line)
                if(wordcount<2) inputedge[wordcount] = stoi(word);
                else if(wordcount == 2){
                    time = stol(word);
                }
                else if(wordcount == 3){
                    // Only if 4 columns are present, the edges are labeled
                    if(Color_Mapping.Initial_Color.find(word) == Color_Mapping.Initial_Color.end()) Color_Mapping.Initial_Color[word] = Color_Mapping.Initial_Color.size();
                    label = Color_Mapping.Initial_Color[word];
                }
                wordcount++;
            }
            if(time < min_time) min_time = time;
            if(time > max_time) max_time = time;
            tglib::TemporalEdge Edge {inputedge[0],inputedge[1],time};
            if(label != -1) Edge.colors.push_back(label);
            if(time_sorted){
                temp_graph.addEdge(Edge);
            }
            else temp_graph.insertEdge(Edge);
        }

        if(min_time != 0){
            temp_graph.shift_time_labels(-min_time);
            max_time -= min_time;
            min_time = 0;
        }
        temp_graph.set_ti(min_time, max_time);
        std::cout << "The created graph has edges ranging from time " << temp_graph.getTimeInterval().first << " to " << temp_graph.getTimeInterval().second << "\n";
    }
}

void label_Vertices(const string& input_file, tglib::IncidentLists& temp_graph) {
    ifstream myfile;
    myfile.open(input_file.c_str());

    if (myfile.is_open()) {
        string line, word;
        size_t part_count = temp_graph.getNumberOfNodes();
        istringstream linestream;

        // Ignore the first two lines
        getline(myfile, line);
        getline(myfile, line);

        // Read the first column of each line which is the atom type and set it as label
        for (int i = 0; i < part_count; i++) {
            getline(myfile, line);
            linestream.clear();
            linestream.str(line);
            getline(linestream, word, ' ');
            temp_graph.setNodeLabel(i, word.c_str());
        }
    }
    else cout << "Unable to open file" << "\n";
}

// For every time step in the simulation, reads the spatial atom positions from the xyz format, performs a Voronoi Tesselation on it, and writes the edges
// corresponding to the Delaunay graph in (v,w,t) format to the edge_file
void Voro_graph_edges(const string& input_file, const string& edge_file, int& partcount, tglib::Time limit, size_t degree_limit, double radial_limit) {
    // Writes all edges of the resulting graphs to the edge_file
    // If degree_limit is an active constraint, the resulting graph might not be undirected!
    ifstream myfile;
    myfile.open(input_file.c_str());
    if(myfile.is_open()) {
        string line;
        tglib::Time timepoint = 0;

        // The first line holds the number of particles/atoms within the simulation
        getline(myfile, line);
        partcount = stoi(line);

        // Ignore the second line
        getline(myfile, line);

        // Read the first timeslice and write the edges into a file
        FILE *output= fopen(edge_file.c_str(), "w");
        Tesselate_slice(myfile, partcount, timepoint, output, degree_limit, radial_limit);
        timepoint++;

        // Read the rest
        while(myfile.peek() != EOF && timepoint < limit){
            // Ignore the first two lines
            getline(myfile,line);
            getline(myfile,line);

            // Read the next timeslice
            Tesselate_slice(myfile, partcount, timepoint, output, degree_limit, radial_limit);

            timepoint++;
        }
        fclose(output);
    }
    else {
        cout << "Unable to open file\n";
    }
}

// Reads the file containing the adjacency matrix and the file containing the temporal information and writes the temporal edges in (v,w,t) format to the edge file
void TU_edges(const string& input_file, const string& time_file, const string& edge_file, int& partcount){
    ifstream input, time_info;
    input.open(input_file.c_str());
    time_info.open(time_file.c_str());
    if(input.is_open() && time_info.is_open()) {
        tglib::NodeId max_vertex=-1;
        string input_line, time_info_line, word;
        FILE *output= fopen(edge_file.c_str(), "w");
        istringstream linestream;
        int v,w,timepoint;
        while(input.peek() != EOF && time_info.peek() != EOF){
            getline(input, input_line);
            linestream.str(input_line);
            getline(linestream, word, ',');
            v = stoi(word)-1;
            getline(linestream, word, ' ');
            getline(linestream, word, ' ');
            w = stoi(word)-1;
            linestream.clear();

            getline(time_info, time_info_line);
            timepoint = stoi(time_info_line);

            if(v>max_vertex) max_vertex = v;
            if(w>max_vertex) max_vertex = w;
            fprintf(output,"%d %d %d\n", v, w, timepoint);
        }
        fclose(output);
        partcount = max_vertex+1;
    }
}

// Takes a file as input with the columns v|w|t|...|label with the label column at position "label_pos"
// Maps the vertex indices to the numbers 0,...,n-1 and writes the edges corresponding to the new vertices in (v,w,t) format to the edge file
void v_w_t_edges_relabel(const string& input_file, const string& edge_file, int& partcount, int label_pos){
    ifstream myfile;
    myfile.open(input_file.c_str());
    if(myfile.is_open()) {
        string line, word;
        istringstream linestream;
        int wordcount;
        unordered_map<string, size_t > Vert_relabel;
        string in_Node_id;
        FILE *output= fopen(edge_file.c_str(), "w");
        while(myfile.peek() != EOF){
            wordcount = 0;
            getline(myfile, line);
            replace(line.begin(), line.end(), ',', ' ');
            linestream.clear();
            linestream.str(line);
            // Loop over all the words within the current line (stored in line)
            while(getline(linestream, word, ' ')){
                // Exclude empty words occurring between duplicate spaces
                if(!word.empty()){
                    if(wordcount<2){
                        // The Node corresponding to "word" is mapped if it is encountered for the first time
                        if(Vert_relabel.find(word) == Vert_relabel.end()) Vert_relabel[word] = Vert_relabel.size();
                        fprintf(output, "%lu ", Vert_relabel[word]);
                    }
                    else if(wordcount == 2){
                        fprintf(output, "%ld", stol(word));
                    }
                    else if(wordcount+1 == label_pos){
                        fprintf(output, " %s\n", word.c_str());
                    }
                    wordcount++;
                }
            }
        }
        partcount = Vert_relabel.size();
    }
    else {
        cout << "Unable to open file\n";
    }
}

void create_TG(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, const string& input_file, const string& filetype){
    string edge_file = "TG_edges.txt";
    tglib::NodeId num_nodes = 0;

    if(filetype == "xyz"){
        string sim_input = "datasets/" + input_file + "/" + input_file + ".xyz";
        Voro_graph_edges(sim_input, edge_file, num_nodes, max_sim_time, max_sim_out_deg);
        if(num_nodes > 0) {
            temp_graph.addNodes(num_nodes);
            initialize_edges(Color_Mapping, temp_graph, edge_file);
            cout << "# Edges in the graph resulting from " << input_file << ": " << temp_graph.getNumberOfEdges() << "\n";
            // Set labels originating from vertices and compute the initial edge colors from them
            label_Vertices(sim_input, temp_graph);
            temp_graph.colorEdges_VertBased(Color_Mapping);
        }
        else cout << "Temporal graph could not be created!\n";
    }
    else if(filetype == "TU"){
        TU_edges("datasets/" + input_file+"/"+input_file+"_A.txt", "datasets/" + input_file+"/"+input_file+"_edge_attributes.txt", edge_file, num_nodes);
        if(num_nodes > 0) {
            temp_graph.addNodes(num_nodes);
            initialize_edges(Color_Mapping, temp_graph, edge_file, false);
            temp_graph.initialize_Event( ("datasets/" + input_file+"/"+input_file+"_graph_indicator.txt").c_str() );
            temp_graph.colorEdges_VertBased(Color_Mapping, "datasets/" + input_file+"/"+input_file+"_node_labels.txt");
        }
        else cout << "Temporal graph could not be created!\n";
    }
    else if(filetype =="vwt5"){
        v_w_t_edges_relabel(input_file, edge_file, num_nodes, 5);
        if(num_nodes > 0) {
            initialize_edges(Color_Mapping, temp_graph, edge_file);
            cout << "# Edges in the graph resulting from " << input_file << ": " << temp_graph.getNumberOfEdges() << "\n";
        }
    }
    else cout << "Filetype is not known. ";
}