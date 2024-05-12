// tg.cpp
// Source file defining the functions of the graph class "IncidentLists" as well as related functions
//
// Richard Restetzki

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <functional>
#include <filesystem>
#include <queue>
#include <numeric>
#include <cmath>
#include <string>
#include <iomanip>
#include "tg.h"

using namespace std;

namespace tglib{

    void Edge_NH::sort() {
        if(Side_2 < Side_1){
            std::vector<color> temp = Side_1;
            Side_1 = Side_2;
            Side_2 = temp;
        }
    }

    void Color_Mapping::WL_color_reset() {
        Color_Mappings.clear();
    }

    void Color_Mapping::clear() {
        depth = 0;
        temp_range = 0;
        Initial_Color.clear();
        Color_Mappings.clear();
    }

    // analyzes the collisions within the hashmaps
    void Color_Mapping::Bucket_collisions() {
        {
            int counter = 0, empty_buckets;
            size_t collisions = 0;
            for(std::unordered_map<Edge_NH, EdgeId>& map: Color_Mappings){
                empty_buckets = 0;
                for(size_t bucket=0; bucket < map.bucket_count(); bucket++){
                    if(map.bucket_size(bucket)>1){
                        collisions += map.bucket_size(bucket)-1;
                        //std::cout << "The bucket " << bucket << " has size " << map.bucket_size(bucket) << "\n";
                        for(auto iter = map.begin(bucket); iter != map.end(bucket); iter++){
                            //std::cout << "Element is: " << iter->second << "\n";
                        }
                        //std::cout << "\n";
                    }
                    else{
                        empty_buckets++;
                    }
                }
                std::cout << "The mapping from WL-depth " << counter << " to WL-depth " << counter+1 << " has " << collisions << " collisions.\n";
                std::cout << "And " << map.bucket_count() << " many buckets of which " << empty_buckets << " are empty\n";
                counter++;
                collisions = 0;
            }
        }
    }

    // Constructs a graph with n nodes and no edges
    IncidentLists::IncidentLists(NodeId n) {
    addNodes(n);
    ti = {inf, 0};
    }

    // Returns n
    [[nodiscard]] size_t IncidentLists::getNumberOfNodes() const{
        return nodes.size();
    }

    // Returns m
    [[nodiscard]] EdgeId IncidentLists::getNumberOfEdges() const{
        return num_edges;
    }

    // Returns the timespan of the graph
    [[nodiscard]] TimeInterval IncidentLists::getTimeInterval() const {
        return ti;
    }

    // Returns a reference to the node with id nid
    TGNode& IncidentLists::getNode(size_t nid) {
        return nodes[nid];
    }

    // Returns a reference to the node vector
    std::vector<TGNode> &IncidentLists::getNodes() {
        return nodes;
    }

    // Returns an iterator pointing to the in-version of the given edge
    vector<TemporalEdge>::iterator IncidentLists::findEdge_in(const TGEdge& e){
        // This inserts all edges after the first - take care that the edges all have the same tail and same time!!
        auto pointer = std::lower_bound(nodes[e.v].inEdges.begin(), nodes[e.v].inEdges.end(), e, pred);
        while( (pointer != nodes[e.v].inEdges.end()) && (pointer->t <= e.t) ){
            if(pointer->u == e.u && pointer->t == e.t){
                return pointer;
            }
            pointer++;
        }
        return nodes[e.v].inEdges.end();
    }

    // Returns an iterator pointing to the out-version of the given edge
    vector<TemporalEdge>::iterator IncidentLists::findEdge_out(const TGEdge& e){
        // This inserts all edges after the first - take care that the edges all have the same tail and same time!!
        auto pointer = std::lower_bound(nodes[e.u].outEdges.begin(), nodes[e.u].outEdges.end(), e, pred);
        while( (pointer != nodes[e.u].outEdges.end()) && (pointer->t <= e.t) ){
            if(pointer->v == e.v && pointer->t == e.t){
                return pointer;
            }
            pointer++;
        }
        return nodes[e.u].outEdges.end();
    }

    // Stores pointers to all out-versions of edges with the specified color in "edges"
    void IncidentLists::get_color_outEdges(vector<TGEdge*>& edges, color color){
        for(TGNode& node: nodes){
            for(TemporalEdge& edge: node.outEdges) {
                if(edge.colors.size()-1==color.depth && edge.colors.back()==color.color_id){
                    edges.push_back(&edge);
                }
            }
        }
    }

    void IncidentLists::reset_Edges() {
        for(TGNode& node: nodes){
            for(TemporalEdge& edge: node.outEdges){
                edge.traversed = false;
            }
            for(TemporalEdge& edge: node.inEdges){
                edge.traversed = false;
            }
        }
    }

    // Clears all colors added by WL
    void IncidentLists::WL_color_reset() {
        for(TGNode& node: nodes){
            for(TemporalEdge& edge: node.outEdges){
                edge.colors.erase(edge.colors.begin()+1, edge.colors.end());
            }
            for(TemporalEdge& edge: node.inEdges){
                edge.colors.erase(edge.colors.begin()+1, edge.colors.end());
            }
        }
    }

    void IncidentLists::clear() {
        nodes.clear();
        num_edges = 0;
        ti = {numeric_limits<Time>::max(), numeric_limits<Time>::min()};
    }

    // Adds num new nodes
    void IncidentLists::addNodes(NodeId num) {
        NodeId size = nodes.size();
        for (NodeId i = 0; i < num; ++i) {
            TGNode node;
            node.id = i + size;
            nodes.push_back(node);
        }
    }

    // Sets the node label of a vertex
    void IncidentLists::setNodeLabel(NodeId nid, const char* label) {
        nodes[nid].label = label;
    }

    // Calculate the neighborhood at a vertex
    void IncidentLists::Vertex_neighborhood(NodeId node, std::vector<color>& NH, int iteration, Time start_time, Time end_time, bool out_going){
        TemporalEdge refedge{0,0,start_time};
        NH.clear();
        vector<tglib::TGEdge>::iterator pointer, end;
        if(out_going){
            pointer = std::lower_bound(nodes[node].outEdges.begin(), nodes[node].outEdges.end(), refedge, pred);
            end = nodes[node].outEdges.end();
        }
        else{
            pointer = std::lower_bound(nodes[node].inEdges.begin(), nodes[node].inEdges.end(), refedge, pred);
            end = nodes[node].inEdges.end();
        }
        while( (pointer != end) && (pointer->t <= end_time) ){
            // Here, the colors vector of the edges could have gained a new entry already
            // This is the case if "colors.size()==iteration+1"
            // -> Take entry [iteration-1] if that applies, otherwise use colors.back()
            if( pointer->colors.size()==iteration+1 ){
                NH.push_back( {static_cast<Depth>(iteration),pointer->colors[iteration-1]} );
            }
            else {
                NH.push_back({static_cast<Depth>(pointer->colors.size()), pointer->colors.back()});
            }
            pointer++;
        }
        // Sort the neighborhood
        sort(NH.begin(), NH.end());
    }

    // Adds a new edge
    void IncidentLists::addEdge(const TGEdge& e) {
        nodes[e.u].outEdges.push_back(e);
        nodes[e.v].inEdges.push_back(e);
        num_edges++;
    }

    void IncidentLists::insertEdge(const TGEdge& e){
        auto pointer = std::lower_bound(nodes[e.u].outEdges.begin(), nodes[e.u].outEdges.end(), e, pred);
        nodes[e.u].outEdges.insert(pointer, e);
        pointer = std::lower_bound(nodes[e.v].inEdges.begin(), nodes[e.v].inEdges.end(), e, pred);
        nodes[e.v].inEdges.insert(pointer, e);
    }

    // Sets the edge to be traversed
    int IncidentLists::traversed_edge(const TemporalEdge& edge){
        int count = 0;
        auto pointer = findEdge_out(edge);
        if(pointer != nodes[edge.u].outEdges.end()){
            pointer->traversed = true;
            count++;
        }
        pointer = findEdge_in(edge);
        if(pointer != nodes[edge.v].inEdges.end()){
            pointer->traversed = true;
            count++;
        }
        return count;
    }

    // A function inserting the neighborhoods of the edge w.r.t. the given parameters
    void IncidentLists::Edge_neighborhood(const TemporalEdge& edge, Edge_NH& NH, int iteration, Time temp_range, bool directed){
        // The first edge in the neighborhood of u having the same label as "edge" will be ignored as "edge" itself is not considered
        // This is technically not necessary but lowers the size of the neighborhood and thus might reduce collisions within the map
        // Only applicable if the edge can be found, e.g. not in the strict_forward case
        // Store the neighborhood at u
        Vertex_neighborhood(edge.u, NH.Side_1, iteration, edge.t+1, edge.t + temp_range, true);
        // Store the neighborhood at v
        // If undirected, the neighborhood at v will also be searched like u and then NH is sorted. If directed, the incoming edges of the tail of "edge" will be considered.
        if(directed){
            Vertex_neighborhood(edge.v, NH.Side_2, iteration, edge.t-temp_range, edge.t-1, false);
        }
        else{
            Vertex_neighborhood(edge.v, NH.Side_2, iteration, edge.t+1, edge.t + temp_range, true);
            NH.sort();
        }
        NH.self_color.depth = iteration;
        NH.self_color.color_id = edge.colors.back();
    }

    // Sets the timespan of a graph
    void IncidentLists::set_ti(Time start, Time end) {
        ti.first = start;
        ti.second = end;
    }

    // Shift all edges within the graph by "shift"
    void IncidentLists::shift_time_labels(Time shift){
        for(TGNode& node: nodes){
            for(TemporalEdge& edge: node.outEdges){
                edge.t += shift;
            }
            for(TemporalEdge& edge: node.inEdges){
                edge.t += shift;
            }
        }
    }

    void IncidentLists::initialize_Event(const char* const& input_file) {
        int event, count=0;
        ifstream input;
        input.open(input_file);
        if(input.is_open()){
            string line;
            while(input.peek() != EOF){
                getline(input, line);
                event = stoi(line);
                for(TemporalEdge& edge: nodes[count].outEdges){
                    edge.event = event;
                }
                for(TemporalEdge& edge: nodes[count].inEdges){
                    edge.event = event;
                }
                count++;
            }
        }
    }

    // Initializes the edges with a color using the node labels of the incident vertices
    void IncidentLists::colorEdges_VertBased(Color_Mapping& Color_Mapping, bool direction_sensitive) {
        std::string edge_label;
        for(TGNode& node: nodes){
            for(TemporalEdge& edge: node.outEdges){
                // Color all the edges that do not have an initial color yet
                if(edge.colors.empty()){
                    if(node.label < nodes[edge.v].label || direction_sensitive){
                        edge_label = node.label + nodes[edge.v].label;
                    }
                    else{
                        edge_label = nodes[edge.v].label + node.label;
                    }

                    if(Color_Mapping.Initial_Color.find(edge_label) == Color_Mapping.Initial_Color.end()){
                        Color_Mapping.Initial_Color[edge_label] = Color_Mapping.Initial_Color.size();
                    }
                    edge.colors.push_back(Color_Mapping.Initial_Color.at(edge_label));
                    auto pointer = findEdge_in(edge);
                    pointer->colors.push_back(Color_Mapping.Initial_Color.at(edge_label));
                }
            }
        }
        std::cout << "# of initial colors: " << Color_Mapping.Initial_Color.size() << "\n";
    }

    void IncidentLists::colorEdges_VertBased(Color_Mapping& Color_Mapping, const string& input_file){
        ifstream input;
        input.open(input_file.c_str());
        if(input.is_open()){
            vector<int> activation_times(getNumberOfNodes(), numeric_limits<int>::max());
            int count=0;
            string line, word;
            istringstream linestream;
            while(input.peek() != EOF){
                getline(input, line);
                linestream.str(line);
                getline(linestream, word, ',');
                getline(linestream, word, ',');
                if(getline(linestream, word, ' ')){
                    getline(linestream, word, ',');
                    activation_times[count] = stoi(word);
                }
                linestream.clear();
                count++;
            }
            std::string edge_label;
            for(TGNode& node: nodes){
                for(TemporalEdge& edge: node.outEdges) {
                    // Color all the edges that do not have an initial color yet
                    if (edge.colors.empty()) {
                        edge_label = std::to_string(activation_times[edge.u]<=edge.t) + std::to_string(activation_times[edge.v]<=edge.t);
                        if (Color_Mapping.Initial_Color.find(edge_label) == Color_Mapping.Initial_Color.end()) {
                            Color_Mapping.Initial_Color[edge_label] = Color_Mapping.Initial_Color.size();
                        }
                        edge.colors.push_back(Color_Mapping.Initial_Color.at(edge_label));
                        auto pointer = findEdge_in(edge);
                        pointer->colors.push_back(Color_Mapping.Initial_Color.at(edge_label));
                    }
                }
            }
            std::cout << "# of initial colors: " << Color_Mapping.Initial_Color.size() << "\n";
        }
    }

    void IncidentLists::colorEdges_uniformly(Color_Mapping& Color_Mapping){
        for(TGNode& node: nodes) {
            for (TemporalEdge &edge: node.outEdges) {
                edge.colors.push_back(0);
            }
            for (TemporalEdge &edge: node.inEdges) {
                edge.colors.push_back(0);
            }
        }
    }

    void IncidentLists::randomize_edges(Color_Mapping& Color_Mapping, mt19937& generator, int count){
        NodeId vertex;
        int counter, added=0, max_tries=1000;
        TemporalEdge edge;
        uniform_int_distribution<mt19937::result_type> rd_insert(0,1);
        uniform_int_distribution<mt19937::result_type> rd_vertex(0,getNumberOfNodes()-1);
        uniform_int_distribution<mt19937::result_type> rd_time(ti.first,ti.second);
        for(int i=0; i< count; i++){
            counter=0;
            vertex = rd_vertex(generator);
            edge = {vertex, static_cast<NodeId>(rd_vertex(generator)), static_cast<Time>(rd_time(generator))};
            // If the edge is already present within the graph, try with another random edge up to 10000 times
            while(findEdge_out(edge) == nodes[edge.u].outEdges.end() && counter < max_tries){
                edge = {vertex, static_cast<NodeId>(rd_vertex(generator)), static_cast<Time>(rd_time(generator))};
                counter++;
            }
            if(counter < max_tries){
                added++;
                insertEdge(edge);
            }
            vertex = rd_vertex(generator);
            uniform_int_distribution<mt19937::result_type> rd_edge(0,nodes[vertex].outEdges.size()-1);
            int edge_nr = rd_edge(generator);
            edge = nodes[vertex].outEdges[edge_nr];
            nodes[vertex].outEdges.erase(nodes[vertex].outEdges.begin() + edge_nr);
            auto pointer = findEdge_in(edge);
            if(pointer == nodes[edge.v].inEdges.end()){
                cout << "ERROR in_edge not existent!\n";
            }
            else{
                nodes[edge.v].inEdges.erase(pointer);
            }
        }
        colorEdges_VertBased(Color_Mapping);
    }

    void IncidentLists::remove_small_cells(tglib::Color_Mapping& Color_Mapping, std::vector<size_t>& hist_vec, int min_count) {
        // Remove all the edge colors for which "hist_vec" has counts below "min_count"
        for (tglib::TGNode &node: nodes) {
            auto pointer = node.outEdges.begin();
            while (pointer != node.outEdges.end()) {
                if (pointer->colors.size() == Color_Mapping.depth+1 && hist_vec[pointer->colors.back()]<min_count) {
                    pointer->colors.pop_back();
                }
                pointer++;
            }
            pointer = node.inEdges.begin();
            while (pointer != node.inEdges.end()) {
                if (pointer->colors.size() == Color_Mapping.depth+1 && hist_vec[pointer->colors.back()]<min_count) {
                    pointer->colors.pop_back();
                }
                pointer++;
            }
        }
    }

    bool IncidentLists::BFS(TemporalEdge& starting_edge, Time temp_range, NodeId depth, EdgeId& edges_count, std::function<void(tglib::TGEdge&)> func, bool color_search){
        // This function explores the graph starting from the starting edge via breadth first search and return the number of edges encountered through "edges_count"
        // The search is performed until the given "depth" is reached and if "func_arg" != -1 the function "func" is performed on all visited vertices
        // Only edges varying by "temp_range" time will be considered in each step
        // If "color_search" is true, the search will also be stopped as soon as the color of the starting edge is found
        // Returns true if BFS was successful and reached the given depth/ reached the starting color - else false
        tglib::EdgeId current_depth_edges, current_depth_traversed, next_depth_edges=0;
        edges_count=0;
        std::queue<tglib::TemporalEdge*> q;
        int same_color_found = 0;
        color edge_color{static_cast<Depth>(starting_edge.colors.size()-1), starting_edge.colors.back()};

        // Set the edge to be traversed
        if(traversed_edge(starting_edge)<2){
            std::cout << "Starting Edge in BFS not found!\n";
            return false;
        }

        // Initialize the queue with neighborhood of the starting_edge
        //next_depth_edges += queue_add_neighbors(q, nodes[starting_edge.u], {starting_edge.t-temp_range, starting_edge.t+temp_range});
        next_depth_edges += queue_add_neighbors(q, nodes[starting_edge.v], {starting_edge.t-temp_range, starting_edge.t+temp_range});
        NodeId current_depth=1;
        while(current_depth<=depth && same_color_found == 0 && !q.empty()){
            current_depth_traversed=0;
            current_depth_edges=next_depth_edges;
            next_depth_edges=0;
            while(!q.empty() && current_depth_traversed<current_depth_edges){
                // First, it is checked if the edge has the same color and the number of edges searched is always increased
                if(color_search && q.front()->colors.size() == edge_color.depth + 1 && q.front()->colors.back() == edge_color.color_id) {
                    same_color_found++;
                }
                func(*q.front());
                edges_count++;
                // Only add edges starting in v to the queue as the edges added always start in u and thus u at (a different) time interval from the last depth was already considered.
                // This ensures that there always exists a temporal path from the root to the leaves in the BFS tree with maximum (possibly negative) transition times of "temp_range"
                TGNode& v=nodes[q.front()->v];
                next_depth_edges += queue_add_neighbors(q, v, {q.front()->t-temp_range, q.front()->t+temp_range});
                q.pop();
                current_depth_traversed++;
            }
            current_depth++;
        }

        reset_Edges();
        if (color_search){
            // Adjust the edges encountered to be all the edges seen until "curr_depth"-1 plus "curr_depth_edges"/"same_color_found"
            // This ensures that the "edges_count" is set to the expected value of edges to be encountered until one reaches the same color given edges at the same depth are chosen randomly
            if (same_color_found==0){
                cout << "After " << edges_count << " edges the same color was not found!\n";
                return false;
            }
            edges_count = (edges_count-current_depth_edges) + (current_depth_edges/same_color_found);
        }
        return true;
    }

    // Computes the WL-colors of all edges in a graph for one iteration
    void IncidentLists::Edge_WL(Color_Mapping& Color_Mapping, Depth iteration, Time temp_range, bool directed) {
        Edge_NH NH;
        for(TGNode& node: nodes){
            for(TemporalEdge& edge: node.outEdges){
                // Only in this case a new color has to be computed - if there are not "iteration" many colors, the refinement was too sharp for this color previously
                // And the neighborhood for this color does not need to be evaluated again (a proof can be found in chapter 3 of the thesis).
                if(edge.colors.size()==iteration){
                    // Compute the Neighborhood of "edge"
                    Edge_neighborhood(edge, NH, iteration, temp_range, directed);
                    // Check if the NH was already seen, insert and color it otherwise
                    if(Color_Mapping.Color_Mappings[iteration-1].find(NH) == Color_Mapping.Color_Mappings[iteration-1].end()){
                        Color_Mapping.Color_Mappings[iteration-1][NH] = Color_Mapping.Color_Mappings[iteration-1].size();
                    }
                    // Add the new or known color to the color vector of both occurrences of the edge
                    edge.colors.push_back(Color_Mapping.Color_Mappings[iteration-1].at(NH));
                    auto pointer = findEdge_in(edge);
                    pointer->colors.push_back(Color_Mapping.Color_Mappings[iteration-1].at(NH));
                }
            }
        }
    }
} // tglib

size_t sieve(tglib::Color_Mapping& Color_Mapping, vector<tglib::IncidentLists *>& graphs, tglib::Depth depth, int min_count){
    vector<size_t> hist_vec(Color_Mapping.Color_Mappings[depth-1].size(),0);
    size_t sieved_colors = 0;
    for (auto graph: graphs){
        color_hist(hist_vec, Color_Mapping, *graph, depth);
    }
    for (auto i: hist_vec){
        if(i<min_count) sieved_colors++;
    }
    for (auto graph: graphs){
        graph->remove_small_cells(Color_Mapping, hist_vec, min_count);
    }
    return sieved_colors;
}

void assign_Moment(tglib::IncidentLists& temp_graph, tglib::TGEdge& edge, int event, vector<tglib::TGEdge*>& edges){
    if(edge.event!=event){
        edge.event = event;
        edges.push_back(&edge);
    }
}

void no_Assignment(tglib::TGEdge& edge){}

void create_Moment(tglib::IncidentLists& temp_graph, tglib::TGEdge* edge, int event, tglib::Depth depth){
    vector<tglib::TGEdge*> edges;
    tglib::EdgeId edges_count=0;
    std::function<void(tglib::TGEdge&)> assign_curr_Moment = bind(&assign_Moment, temp_graph, placeholders::_1, event, edges);
    temp_graph.BFS(*edge, 0, depth, edges_count, assign_curr_Moment, false);
    // Maybe store the amount of equal edges found and return it
}

bool pred(const tglib::TemporalEdge &a, const tglib::TemporalEdge &b){
    if(a.t < b.t) return true;
    return false;
}

bool comp(const vector<size_t> &a, const vector<size_t> &b){
    // this ordering guarantees not ascending order
    if(a[1] > b[1]) return true;
    return false;
}

double cos_distance(vector<vector<size_t>> &vec1, vector<vector<size_t>> &vec2){
    double similarity=0, norming_sum=0, norming_factor;

    // Compute the inner product (vec1,vec2)
    // Summing over all vector entries yield the same result as appending all vectors together
    for(int i=0; i<vec1.size(); i++){
        similarity += std::inner_product(vec1[i].begin(),vec1[i].end(),vec2[i].begin(),0.0);
    }

    // Norm the product s.t. the result will be a number within [0,1]
    // If [vec1] and [vec2] were just vectors the code would look like this:
    // norming_factor = std::inner_product(vec1.begin(),vec1.end(),vec1.begin(),0.0) * std::inner_product(vec2.begin(),vec2.end(),vec2.begin(),0.0);
    // Summing over all vector entries yield the same result as appending all vectors together
    for(auto i: vec1){
        norming_sum += std::inner_product(i.begin(),i.end(),i.begin(),0.0);
    }
    norming_factor = norming_sum;
    norming_sum=0;
    for(auto i: vec2){
        norming_sum += std::inner_product(i.begin(),i.end(),i.begin(),0.0);
    }
    norming_factor *= norming_sum;
    norming_factor = std::sqrt(norming_factor);

    return similarity/norming_factor;
}

double average_color_dist(tglib::IncidentLists& temp_graph, tglib::color color, tglib::Time temp_range){
    size_t dist_sum=0;
    int curr_count=0, color_count=0, BFS_count=0;
    tglib::EdgeId dist;

    for(const tglib::TGNode& node: temp_graph.getNodes()) {
        for(const tglib::TemporalEdge& edge: node.outEdges){
            if(edge.colors.size()==color.depth+1 && edge.colors.back()==color.color_id){
                color_count++;
            }
        }
    }
    int sample_int = max(1, color_count/1000);
    std::cout << "color_size: " << color_count << "\n";
    for(tglib::TGNode& node: temp_graph.getNodes()) {
        for(tglib::TemporalEdge& edge: node.outEdges){
            if(edge.colors.size()==color.depth+1 && edge.colors.back()==color.color_id){
                if(curr_count%sample_int==0){
                    if(!temp_graph.BFS(edge, temp_range, temp_graph.getNumberOfNodes(), dist, no_Assignment)){
                        std::cout << "same color not found!\n";
                    }
                    dist_sum += dist;
                    BFS_count++;
                }
                curr_count++;
            }
        }
    }
    return (double)dist_sum/BFS_count;
}

tglib::EdgeId queue_add_neighbors(std::queue<tglib::TemporalEdge*>& q, tglib::TGNode& Node, tglib::TimeInterval interval){
    // Adds all edges within "interval" that are adjacent to "Node" that are not traversed to the queue and sets them to traversed
    tglib::TemporalEdge refedge{0,0,interval.first};
    tglib::EdgeId edge_count=0;
    auto pointer = std::lower_bound(Node.outEdges.begin(), Node.outEdges.end(), refedge, pred);
    while( (pointer != Node.outEdges.end()) && (pointer->t <= interval.second) ){
        if ( !pointer->traversed ){
            // If the current edge is not the starting_edge or its inverse, add it to q
            q.push(& *pointer);
            edge_count++;
            pointer->traversed = true;
        }
        pointer++;
    }
    return edge_count;
}

void frequent_colors(vector<vector<size_t>>& color_nr_freq, tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, tglib::TimeInterval TimeInterval){
    vector<size_t> hist_vec(Color_Mapping.Color_Mappings[Color_Mapping.depth-1].size(), 0), color(2,0);
    color_hist(hist_vec, Color_Mapping, temp_graph, Color_Mapping.depth, TimeInterval);
    vector<vector<size_t>> color_freq(hist_vec.size());
    color_nr_freq.assign(hist_vec.size(), vector<size_t>(2,0));
    for(size_t i=0; i<color_nr_freq.size(); i++){
        color_nr_freq[i][0]=i;
        color_nr_freq[i][1]=hist_vec[i];
    }
    sort(color_nr_freq.begin(), color_nr_freq.end(), comp);
    std::cout << "frequent_colors successful!\n";
}

tglib::Time frequent_interval(std::vector<size_t>& occurrence_count, size_t& maximum, tglib::Time timespan){
    // Find the maximum of a sliding sum over the corresponding "timespan"-many consecutive time_points
    size_t occurrences=0;
    tglib::Time start_time=0;
    for(tglib::Time i=0; i<timespan; i++){
        occurrences += occurrence_count[i];
    }
    maximum = occurrences;
    for(tglib::Time i=timespan; i<occurrence_count.size(); i++){
        occurrences -= occurrence_count[i-timespan];
        occurrences += occurrence_count[i];
        if(occurrences>maximum){
            maximum = occurrences;
            start_time=i-timespan+1;
        }
    }
    return start_time;
}

void find_frequent_intervals(vector<tglib::Time>& start_times, vector<size_t>& maxima, tglib::IncidentLists& temp_graph, tglib::color color, vector<tglib::Time> timespans){
    tglib::Time graph_duration = temp_graph.getTimeInterval().second-temp_graph.getTimeInterval().first+1;
    if (any_of(timespans.begin(), timespans.end(), [graph_duration](tglib::Time timespan) {return timespan>graph_duration;})){
        std::cout << "timespans is not matching the graph time interval!\n";
        return;
    }
    else if (start_times.size() != maxima.size() || start_times.size() != timespans.size()) {
        std::cout << "sizes of start_times, maxima, and timespans do not match!\n";
        return;
    }
    else {
        std::vector<size_t> occurrence_count(graph_duration, 0);
        // Store the occurrences of the color at every time_point in "occurrence_count"
        for (const tglib::TGNode &node: temp_graph.getNodes()) {
            auto pointer = node.outEdges.begin();
            while (pointer != node.outEdges.end()) {
                if (pointer->colors.size() == color.depth+1 && pointer->colors.back() == color.color_id) {
                    occurrence_count[pointer->t]++;
                }
                pointer++;
            }
        }
        for (int i = 0; i < start_times.size(); i++){
            start_times[i] = frequent_interval(occurrence_count, maxima[i], timespans[i]);
        }
    }
}

// Analyzes how often neighborhoods change (implemented for chemical simulations)
int check_change(tglib::IncidentLists& temp_graph, tglib::NodeId nid) {
    const tglib::TGNode& Node = temp_graph.getNode(nid);
    vector<tglib::TemporalEdge> neighbors = Node.outEdges;
    vector<tglib::NodeId> previous, current;
    tglib::NodeId curr_time = 0, equal_steps = 0;
    tglib::EdgeId iter = 0;
    bool equal_step = true;

    while (iter < neighbors.size()) {
        if(curr_time != neighbors[iter].t) {
            // Compare the neighbors from the previous timepoint to the one just passed
            sort(current.begin(),current.end());
            if (previous.size() != current.size()){
                equal_step = false;
            }
            else{
                for (int i=0; i<previous.size(); i++){
                    if(previous[i]!=current[i]){
                        equal_step = false;
                        break;
                    }
                }
            }
            if (equal_step){
                equal_steps++;
            }

            // The current vector gets swapped to the previous and that one gets flushed
            previous = current;
            current.clear();

            // Adjust the timepoint and start collecting the neighbors for the new timepoint in the current vector
            curr_time++;
            equal_step = true;
            current.push_back(neighbors[iter].v);
        }
        else {
            // Keep reading all the neighbors at this timepoint
            current.push_back(neighbors[iter].v);
        }
        iter++;
    }
    return equal_steps;
}

void color_hist(std::vector<std::vector<size_t>>& hist_vec, tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, tglib::TimeInterval interval){
    // Count the occurrences of every color simultaneously while looping over all outgoing edges
    bool conf_format = true;
    if(hist_vec.size() == Color_Mapping.depth+1){
        if(Color_Mapping.Initial_Color.size() != hist_vec[0].size()) conf_format = false;
        for(int i=0; i<Color_Mapping.depth; i++){
            if(Color_Mapping.Color_Mappings[i].size() != hist_vec[i+1].size()) conf_format = false;
        }
    }
    else conf_format = false;

    if(!conf_format){
        cout << "wrong input vector format for color_hist!\n";
    }
    else{
        tglib::TemporalEdge refedge{0,0,interval.first};
        for(const tglib::TGNode& node: temp_graph.getNodes()) {
            auto pointer = std::lower_bound(node.outEdges.begin(), node.outEdges.end(), refedge, pred);
            while ((pointer != node.outEdges.end()) && (pointer->t <= interval.second)) {
                hist_vec[pointer->colors.size()-1][pointer->colors.back()]++;
                pointer++;
            }
        }
    }
}

void color_hist(std::vector<size_t>& hist_vec, tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, tglib::Depth depth, tglib::TimeInterval interval){
    // Count the occurrences of every color simultaneously while looping over all outgoing edges
    if(Color_Mapping.Color_Mappings[depth-1].size() != hist_vec.size()){
        cout << "wrong input vector format for color_hist at specific depth!\n";
        cout << Color_Mapping.Color_Mappings[depth-1].size() << " != " << hist_vec.size();
    }
    else{
        tglib::TemporalEdge refedge{0,0,interval.first};
        for(const tglib::TGNode& node: temp_graph.getNodes()) {
            auto pointer = std::lower_bound(node.outEdges.begin(), node.outEdges.end(), refedge, pred);
            while ((pointer != node.outEdges.end()) && (pointer->t <= interval.second)) {
                if(pointer->colors.size() == depth+1) hist_vec[pointer->colors.back()]++;
                pointer++;
            }
        }
    }
}

void event_color_hist(std::vector<std::vector<size_t>>& hist_vec, tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, int event){
    // Count the occurrences of the colors in subgraph "event" simultaneously while looping over all outgoing edges
    bool conf_format = true;
    if(hist_vec.size() == Color_Mapping.depth+1){
        if(Color_Mapping.Initial_Color.size() != hist_vec[0].size()) conf_format = false;
        for(int i=0; i<Color_Mapping.depth; i++){
            if(Color_Mapping.Color_Mappings[i].size() != hist_vec[i+1].size()) conf_format = false;
        }
    }
    else conf_format = false;

    if(!conf_format){
        cout << "wrong input vector format for color_hist!\n";
    }
    else{
        for(const tglib::TGNode& node: temp_graph.getNodes()) {
            auto pointer = node.outEdges.begin();
            while (pointer != node.outEdges.end()) {
                if(pointer->event==event) hist_vec[pointer->colors.size()-1][pointer->colors.back()]++;
                pointer++;
            }
        }
    }
}

double interval_similarity(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, tglib::TimeInterval interval1, tglib::TimeInterval interval2){
    vector<vector<size_t>> vec1;
    vec1.emplace_back(Color_Mapping.Initial_Color.size(),0);
    for(int i=0; i<Color_Mapping.depth; i++){
        vec1.emplace_back(Color_Mapping.Color_Mappings[i].size(),0);
    }
    vector<vector<size_t>> vec2=vec1;
    color_hist(vec1, Color_Mapping, temp_graph, interval1);
    color_hist(vec2, Color_Mapping, temp_graph, interval2);
    return cos_distance(vec1,vec2);
}

// Prints the colors having the highest ratio between the densest interval and the normal density
void print_color_densities(const char* const& input_file, tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, size_t min_occurrence, bool cut_degenerated_end){
    // Initialization
    std::cout << "print_color_densities_started\n";
    vector<tglib::Time> intervals{3,10,30,100};
    vector<size_t> max_density_color(4,0);
    vector<size_t> maxima(4,0);
    vector<double> max_density_color_vals(4,0);
    vector<tglib::Time> max_density_color_starts(4,0);
    vector<tglib::Time> curr_color_starts(4);
    vector<vector<size_t>> color_nr_freq;
    tglib::TimeInterval TimeInterval = temp_graph.getTimeInterval();
    size_t iter=0;
    if(cut_degenerated_end){
        TimeInterval.second -= Color_Mapping.depth*Color_Mapping.temp_range;
    }
    tglib::Time graph_duration = TimeInterval.second-TimeInterval.first;

    // Computing the frequent colors
    vector<vector<size_t>> color_freq(Color_Mapping.Color_Mappings.back().size(), vector<size_t> (2,0));
    frequent_colors(color_nr_freq, Color_Mapping, temp_graph, TimeInterval);
    FILE *output_file= fopen("Labels_Histogram.txt", "w");
    while(iter<color_nr_freq.size() && color_nr_freq[iter][1] >= min_occurrence){
        fprintf(output_file,"Color %zu, #%zu with density relations", color_nr_freq[iter][0], color_nr_freq[iter][1]);
        // find the intervals for color "color_nr_freq[iter][0]" having the most occurrences
        find_frequent_intervals(curr_color_starts, maxima, temp_graph, {Color_Mapping.depth, static_cast<tglib::EdgeId>(color_nr_freq[iter][0])}, intervals);
        for(int j=0; j<intervals.size(); j++){
            double occurrence_ratio = ((double)graph_duration*maxima[j])/((double)intervals[j]*color_nr_freq[iter][1]);
            fprintf(output_file, " %f", occurrence_ratio);
            // update the currently densest color if necessary
            if(max_density_color_vals[j]<occurrence_ratio){
                max_density_color_vals[j] = occurrence_ratio;
                max_density_color[j] = color_nr_freq[iter][0];
                max_density_color_starts[j] = curr_color_starts[j];
            }
        }
        if(iter<0){
            fprintf(output_file, " %f", average_color_dist(temp_graph, {Color_Mapping.depth, static_cast<tglib::EdgeId>(color_nr_freq[iter][0])}, 3));
            fprintf(output_file, " %f", average_color_dist(temp_graph, {Color_Mapping.depth, static_cast<tglib::EdgeId>(color_nr_freq[iter][0])}, 10));
        }
        fprintf(output_file, "\n");
        std::cout << "Currently iter is at " << iter << "\n";
        iter++;
        if(color_nr_freq.size()==iter){
            std::cout << "A break because no color has lower occurrences than min_occurrence " << min_occurrence << "!\n";
            break;
        }
    }

    // Plotting the most dense intervals with their respective color
    for(int i=0; i<intervals.size(); i++){
        plot_labels(input_file, temp_graph, {Color_Mapping.depth, static_cast<tglib::EdgeId>(max_density_color[i])}, max_density_color_starts[i], intervals[i]);
    }
}

void plot_labels(const char* const& input_file, tglib::IncidentLists& temp_graph, tglib::color color, tglib::Time start_time, tglib::Time timespan){
    string foldername = "timesteps_"+to_string(start_time)+"-"+to_string(start_time+timespan-1)+"_edge_type_"+to_string(color.depth)+","+to_string(color.color_id);
    std::filesystem::create_directory(foldername);
    vector<tglib::TemporalEdge> Edges;
    tglib::TemporalEdge refedge{0,0,start_time};
    // Store all the edges that have to be plotted
    for(const tglib::TGNode& node: temp_graph.getNodes()) {
        auto pointer = std::lower_bound(node.outEdges.begin(), node.outEdges.end(), refedge, pred);
        while (pointer != node.outEdges.end() && pointer->t < start_time+timespan) {
            if(pointer->colors.size() == color.depth+1 && pointer->colors.back() == color.color_id){
                Edges.push_back(*pointer);
            }
            pointer++;
        }
    }
    std::sort(Edges.begin(), Edges.end(), pred);
    FILE *edge_file= fopen( ("./" + foldername + "/edges_plot").c_str(), "w");
    FILE *node_file= fopen( ("./" + foldername + "/nodes_plot").c_str(), "w");
    FILE *node_plain_file= fopen( ("./" + foldername + "/nodes_plain").c_str(), "w");
    ifstream myfile;
    istringstream linestream;
    string line, word;
    tglib::NodeId partcount = temp_graph.getNumberOfNodes();
    tglib::Time curr_time;
    size_t edge_iter=0;
    vector<vector<double>> Nodes_pos(partcount, vector<double>(3));
    vector<bool> Node_covered(partcount, false);
    myfile.open(input_file);
    if(myfile.is_open()) {
        // Skip the input file up to the time_point "start_time"
        for(curr_time=0; curr_time< start_time; curr_time++){
            for(int j=0; j<partcount+2; j++){
                getline(myfile, line);
            }
        }
        while(curr_time < start_time+timespan && myfile.peek() != EOF){
            getline(myfile, line);
            getline(myfile, line);
            int argcount, wordcount;
            for(tglib::NodeId linecount=0; linecount<partcount; linecount++){
                if(myfile.peek() != EOF){
                    getline(myfile, line);
                    argcount = 0;
                    wordcount = 0;
                    linestream.clear();
                    linestream.str(line);
                    //Loop over all the words within the current line (stored in line)
                    while(getline(linestream, word, ' ')){
                        //Exclude empty words occurring between duplicate spaces
                        if(!word.empty()){
                            if( (wordcount==1) || (argcount>0&&argcount<3) ){
                                // Store the coordinates for the particle
                                Nodes_pos[linecount][argcount] = stod(word);
                                argcount++;
                            }
                            wordcount++;
                        }
                    }
                }
                else return;
            }
            // write all relevant edges at timepoint "start_time+time_iter" to the plot file
            while(Edges[edge_iter].t==curr_time){
                tglib::NodeId u=Edges[edge_iter].u;
                tglib::NodeId v=Edges[edge_iter].v;
                fprintf(edge_file,"%g %g %g %d\n", Nodes_pos[u][0], Nodes_pos[u][1], Nodes_pos[u][2], Edges[edge_iter].t);
                fprintf(edge_file,"%g %g %g %d\n\n\n", Nodes_pos[v][0], Nodes_pos[v][1], Nodes_pos[v][2], Edges[edge_iter].t);
                if(!Node_covered[u]){
                    fprintf(node_file, "%d %g %g %g\n", u, Nodes_pos[u][0], Nodes_pos[u][1], Nodes_pos[u][2]);
                    fprintf(node_plain_file, "%d\n", u);
                    Node_covered[u] = true;
                }
                if(!Node_covered[v]){
                    fprintf(node_file, "%d %g %g %g\n", v, Nodes_pos[v][0], Nodes_pos[v][1], Nodes_pos[v][2]);
                    fprintf(node_plain_file, "%d\n", v);
                    Node_covered[v] = true;
                }
                edge_iter++;
                if(Edges.size()==edge_iter){
                    break;
                }
            }
            curr_time++;
        }
    }
}

// Plots all edges corresponding to a specific event (only applicable if xyz format is input_file)
void plot_event(string input_file, tglib::IncidentLists& temp_graph, int event){
    string foldername = "event_"+to_string(event);
    std::filesystem::create_directory(foldername);
    vector<tglib::TemporalEdge> Edges;
    // Store all the edges that have to be plotted
    for(const tglib::TGNode& node: temp_graph.getNodes()) {
        for(const tglib::TGEdge& edge: node.outEdges){
            if(edge.event == event) Edges.push_back(edge);
        }
    }
    std::sort(Edges.begin(), Edges.end(), pred);
    FILE *edge_file= fopen( ("./" + foldername + "/edges_plot").c_str(), "w");
    FILE *node_file= fopen( ("./" + foldername + "/nodes_plot").c_str(), "w");
    FILE *node_plain_file= fopen( ("./" + foldername + "/nodes_plain").c_str(), "w");
    ifstream myfile;
    istringstream linestream;
    string line, word;
    tglib::NodeId partcount = temp_graph.getNumberOfNodes();
    tglib::Time curr_time, start_time = Edges[0].t, timespan = 3;
    size_t edge_iter=0;
    vector<vector<double>> Nodes_pos(partcount, vector<double>(3));
    vector<bool> Node_covered(partcount, false);
    myfile.open( (input_file + "/" + input_file + ".xyz").c_str() );
    if(myfile.is_open()) {
        // Skip the input file up to the time_point "start_time"
        for(curr_time=0; curr_time< start_time; curr_time++){
            for(int j=0; j<partcount+2; j++){
                getline(myfile, line);
            }
        }
        while(curr_time < start_time+timespan && myfile.peek() != EOF){
            getline(myfile, line);
            getline(myfile, line);
            int argcount, wordcount;
            for(tglib::NodeId linecount=0; linecount<partcount; linecount++){
                if(myfile.peek() != EOF){
                    getline(myfile, line);
                    argcount = 0;
                    wordcount = 0;
                    linestream.clear();
                    linestream.str(line);
                    //Loop over all the words within the current line (stored in line)
                    while(getline(linestream, word, ' ')){
                        //Exclude empty words occurring between duplicate spaces
                        if(!word.empty()){
                            if( (wordcount==1) || (argcount>0&&argcount<3) ){
                                // Store the coordinates for the particle
                                Nodes_pos[linecount][argcount] = stod(word);
                                argcount++;
                            }
                            wordcount++;
                        }
                    }
                }
                else return;
            }
            // write all relevant edges at timepoint "start_time+time_iter" to the plot file
            vector<double> centers = Nodes_pos[Edges[0].u];
            while(Edges[edge_iter].t==curr_time){
                tglib::NodeId u=Edges[edge_iter].u;
                tglib::NodeId v=Edges[edge_iter].v;
                // If the edges are going through the periodic border, the plot will not represent the actual structure
                // Thus, the spatial positions are folded where necessary
                if(!Node_covered[u]){
                    if(Nodes_pos[u][0]-centers[0] > x_max/2) Nodes_pos[u][0] -= x_max;
                    else if(Nodes_pos[u][0]-centers[0] < -x_max/2) Nodes_pos[u][0] += x_max;
                    if(Nodes_pos[u][1]-centers[1] > y_max) Nodes_pos[u][1] -= y_max;
                    else if(Nodes_pos[u][1]-centers[1] < -y_max) Nodes_pos[u][1] += y_max;
                    if(Nodes_pos[u][2]-centers[2] > z_max) Nodes_pos[u][2] -= z_max;
                    else if(Nodes_pos[u][2]-centers[2] < -z_max) Nodes_pos[u][2] += z_max;
                    fprintf(node_file, "%d %g %g %g\n", u, Nodes_pos[u][0], Nodes_pos[u][1], Nodes_pos[u][2]);
                    fprintf(node_plain_file, "%d\n", u);
                    Node_covered[u] = true;
                }
                if(!Node_covered[v]){
                    if(Nodes_pos[v][0]-centers[0] > x_max/2) Nodes_pos[v][0] -= x_max;
                    else if(Nodes_pos[v][0]-centers[0] < -x_max/2) Nodes_pos[v][0] += x_max;
                    if(Nodes_pos[v][1]-centers[1] > y_max) Nodes_pos[v][1] -= y_max;
                    else if(Nodes_pos[v][1]-centers[1] < -y_max) Nodes_pos[v][1] += y_max;
                    if(Nodes_pos[v][2]-centers[2] > z_max) Nodes_pos[v][2] -= z_max;
                    else if(Nodes_pos[v][2]-centers[2] < -z_max) Nodes_pos[v][2] += z_max;
                    fprintf(node_file, "%d %g %g %g\n", v, Nodes_pos[v][0], Nodes_pos[v][1], Nodes_pos[v][2]);
                    fprintf(node_plain_file, "%d\n", v);
                    Node_covered[v] = true;
                }
                fprintf(edge_file,"%g %g %g %d\n", Nodes_pos[u][0], Nodes_pos[u][1], Nodes_pos[u][2], Edges[edge_iter].t);
                fprintf(edge_file,"%g %g %g %d\n\n\n", Nodes_pos[v][0], Nodes_pos[v][1], Nodes_pos[v][2], Edges[edge_iter].t);
                edge_iter++;
                if(Edges.size()==edge_iter){
                    break;
                }
            }
            curr_time++;
        }
    }
}

void output_rd_similarities(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, int n, int randomization_step){
    int rd_operations = 0;
    mt19937 mt(seed);
    tglib::IncidentLists graph_copy = temp_graph;
    FILE *output_file= fopen("Rd_similarities", "w");
    double similarity;
    for(int i=0; i<n; i++){
        TNWL(Color_Mapping, {&temp_graph, &graph_copy}, 5, 1, false, 1000);
        similarity = TNWL_similarity(temp_graph, graph_copy, Color_Mapping);
        cout << "After " << rd_operations << " edge randomizations the similarity is: " << similarity << "\n";
        fprintf(output_file,"%d %g\n", rd_operations, similarity);
        Color_Mapping.WL_color_reset();
        temp_graph.WL_color_reset();
        graph_copy.WL_color_reset();
        graph_copy.randomize_edges(Color_Mapping, mt, randomization_step);
        rd_operations += randomization_step;
    }
}

void output_similarities(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, int precision){
    std::cout << "Similarities:\n";
    for(tglib::Time i=0; i<10000; i+=1000){
        for(tglib::Time j=0; j<10000; j+=1000){
            tglib::TimeInterval interval1{i,i+999};
            tglib::TimeInterval interval2{j,j+999};
            std::cout << std::setprecision(precision) << interval_similarity(Color_Mapping, temp_graph, interval1, interval2) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    for(tglib::Time i=0; i<10000; i+=2000){
        for(tglib::Time j=0; j<10000; j+=1000){
            tglib::TimeInterval interval1{i,i+1999};
            tglib::TimeInterval interval2{j,j+999};
            std::cout << std::setprecision(precision) << interval_similarity(Color_Mapping, temp_graph, interval1, interval2) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    for(tglib::Time i=0; i<10000; i+=5000){
        for(tglib::Time j=0; j<10000; j+=1000){
            tglib::TimeInterval interval1{i,i+4999};
            tglib::TimeInterval interval2{j,j+999};
            std::cout << std::setprecision(precision) << interval_similarity(Color_Mapping, temp_graph, interval1, interval2) << " ";
        }
        std::cout << "\n";
    }
    return;
}

// Computes WL labels up to a given depth of all graphs in a list
void TNWL(tglib::Color_Mapping& Color_Mapping, vector<tglib::IncidentLists *> graphs, tglib::Depth depth, tglib::Time temp_range, bool directed, int min_count){
    bool stop=false;
    Color_Mapping.temp_range = temp_range;
    for (int iteration=1; iteration<=depth && !stop; iteration++){
        Color_Mapping.Color_Mappings.emplace_back();
        Color_Mapping.depth = iteration;
        int counter = 1;
        for (auto graph: graphs){
            graph->Edge_WL(Color_Mapping, iteration, temp_range, directed);
            counter++;
        }

        size_t sieved_count = sieve(Color_Mapping, graphs, iteration, min_count);
        std::cout << "In iteration " << iteration << " " << Color_Mapping.Color_Mappings.back().size() << " colors were identified of which " << sieved_count << " are sieved." << "\n";
        // Set "stop" to true as soon as no color exceeds min_count within the current iteration
        if(sieved_count==Color_Mapping.Color_Mappings.back().size()){
            stop = true;
            Color_Mapping.depth = iteration-1;
        }
    }
    cout << "Final depth of WL: " << Color_Mapping.depth << "\n";
}

// prints the matrix holding the pairwise similarities of the events in the given temporal graph
void print_gram_of_events(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, int nr_events){
    vector<vector<size_t>> sample_vec;
    sample_vec.emplace_back(Color_Mapping.Initial_Color.size(),0);
    for(int i=0; i<Color_Mapping.depth; i++){
        sample_vec.emplace_back(Color_Mapping.Color_Mappings[i].size(),0);
    }
    std::vector<std::vector<std::vector<size_t>>> hist_vecs(nr_events, sample_vec);
    for(int i=0; i<nr_events; i++){
        event_color_hist(hist_vecs[i], Color_Mapping, temp_graph, i+1);
    }
    FILE *output= fopen("Gram_matrix_of_events", "w");

    // Similarity of 0,0 is 1 and is entered separately
    fprintf(output, "1");
    // Similarity of 0,j ; j>0 for the first row
    for(int j=1; j<nr_events; j++){
        fprintf(output, " %f", cos_distance(hist_vecs[0], hist_vecs[j]));
    }
    // Similarity of i,j ; i>0 for all other rows
    for(int i=1; i<nr_events; i++){
        fprintf(output, "\n %f", cos_distance(hist_vecs[i], hist_vecs[0]));
        for(int j=1; j<nr_events; j++){
            fprintf(output, " %f", cos_distance(hist_vecs[i], hist_vecs[j]));
        }
    }
}

// Analyzes the labels at the max depth only
void label_analysis(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& graph_1, string input_file_1, tglib::IncidentLists& graph_2, string input_file_2){
    tglib::Depth depth = Color_Mapping.depth;
    size_t nr_colors = Color_Mapping.Color_Mappings[depth-1].size();
    std::vector<size_t> hist_vec_1(nr_colors,0), hist_vec_2(nr_colors,0);
    color_hist(hist_vec_1, Color_Mapping, graph_1, depth);
    color_hist(hist_vec_2, Color_Mapping, graph_2, depth);

    size_t id_1=0, id_2=0, id_common=0;
    double min_1 = numeric_limits<double>::max(), min_2 = numeric_limits<double>::max(), min_common = numeric_limits<double>::max();
    double ratio_1div2, ratio_2div1;
    cout << "Total colors analyzed: " << Color_Mapping.Color_Mappings[depth-1].size() << "\n";
    // Store the colors for which the ratios are lowest and the color for which it is closest to 1
    for(size_t i=0; i<nr_colors; i++){
        if(hist_vec_1[i]==0){
            if(hist_vec_2[i]==0){
                ratio_1div2 = numeric_limits<double>::max();
                ratio_2div1 = numeric_limits<double>::max();
            }
            else{
                ratio_1div2 = 0;
                ratio_2div1 = numeric_limits<double>::max();
            }
        }
        else if(hist_vec_2[i]==0){
            ratio_1div2 = numeric_limits<double>::max();
            ratio_2div1 = 0;
        }
        else{
            ratio_1div2 = 1.0*hist_vec_1[i]/hist_vec_2[i];
            ratio_2div1 = 1.0*hist_vec_2[i]/hist_vec_1[i];
        }
        if(ratio_1div2<min_1){
            min_1 = ratio_1div2;
            id_1 = i;
        }
        if(ratio_2div1<min_2){
            min_2 = ratio_2div1;
            id_2 = i;
        }
        if(max(ratio_1div2, ratio_2div1)<min_common){
            min_common = max(ratio_1div2, ratio_2div1);
            id_common = i;
        }
    }
    // Due to taking the minimum the ids are swapped for the unique labels
    cout << "The most unique label for graph_1 is " << id_2 << " with " << hist_vec_1[id_2] << ", " << hist_vec_2[id_2] << " occurrences in the first and second graph, respectively\n";
    cout << "The most unique label for graph_2 is " << id_1 << " with " << hist_vec_1[id_1] << ", " << hist_vec_2[id_1] << " occurrences in the first and second graph, respectively\n";
    cout << "The most common label for both graphs is " << id_common << " with " << hist_vec_1[id_common] << ", " << hist_vec_2[id_common] << " occurrences in the first and second graph, respectively\n";

    // Search for the first edge with the respective label, set the event label of all edges in the "depth"-neighborhood to 1,2,3, respectively and plot that event
    std::vector<tglib::TGEdge*> edges;
    graph_1.get_color_outEdges(edges, {depth, static_cast<tglib::EdgeId>(id_2)});
    create_Moment(graph_1, edges[0], 1, depth);
    plot_event(input_file_1, graph_1, 1);

    edges.clear();
    graph_2.get_color_outEdges(edges, {depth, static_cast<tglib::EdgeId>(id_1)});
    create_Moment(graph_2, edges[0], 2, depth);
    plot_event(input_file_2, graph_2, 2);

    edges.clear();
    graph_1.get_color_outEdges(edges, {depth, static_cast<tglib::EdgeId>(id_common)});
    create_Moment(graph_1, edges[0], 3, depth);
    plot_event(input_file_1, graph_1, 3);
}

double TNWL_similarity(tglib::IncidentLists& graph_1, tglib::IncidentLists& graph_2, tglib::Color_Mapping& Color_Mapping){
    // For this function to work, both graphs must have had their colors assigned through "Color_Mapping"
    vector<vector<size_t>> vec1;
    vec1.emplace_back(Color_Mapping.Initial_Color.size(),0);
    for(int i=0; i<Color_Mapping.depth; i++){
        vec1.emplace_back(Color_Mapping.Color_Mappings[i].size(),0);
    }
    vector<vector<size_t>> vec2=vec1;
    color_hist(vec1, Color_Mapping, graph_1);
    color_hist(vec2, Color_Mapping, graph_2);
    return cos_distance(vec1,vec2);
}