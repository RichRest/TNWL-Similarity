#ifndef CHEM_INPUT_H
#define CHEM_INPUT_H

// chem_input.h
// Header file declaring functionalities to process chemical simulations
//
// Richard Restetzki

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "voro++/voro++.hh"
#include "tg.h"

extern double x_min,x_max,y_min,y_max,z_min,z_max;
extern const bool x_periodic,y_periodic,z_periodic;
extern const int n_x,n_y,n_z,max_sim_out_deg;
extern int max_sim_time;

using namespace std;
using namespace voro;

struct Voro_neighbor {
    // WL-Iteration depth corresponding to the color
    int partindex;
    // Unique color_id within for the respective depth
    double distance;

    inline bool operator<(const Voro_neighbor &other) const {
        if(this->distance<other.distance){
            return true;
        }
        return false;
    }
};

// Basic functionalities
int sgn(double val);
double eucl_dist(vector<double>& point_1, vector<double>& point_2);
int Tesselate_slice(ifstream& read_file, int partcount, int timepoint, FILE* output_file, size_t degree_limit = std::numeric_limits<size_t>::max(), double radial_limit = std::numeric_limits<double>::max());
void initialize_edges(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, const string& edge_file, bool time_sorted = true);
void label_Vertices(const string& input_file, tglib::IncidentLists& temp_graph);

// Input conversions into edge file
void Voro_graph_edges(const string& input_file, const string& edge_file, int& partcount, tglib::Time time_limit = tglib::inf, size_t degree_limit = std::numeric_limits<size_t>::max(), double radial_limit = std::numeric_limits<double>::max());
void TU_edges(const string& input_file, const string& time_file, const string& edge_file, int& partcount);
void v_w_t_edges_relabel(const string& input_file, const string& edge_file, int& partcount, int label_pos=0);

// Temporal graph creation function
void create_TG(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, const string& input_file, const string& filetype);

#endif // CHEM_INPUT_H