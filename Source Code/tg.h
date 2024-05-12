#ifndef SIMILARITYMEASURE_TG_H
#define SIMILARITYMEASURE_TG_H

// tg.h
// Header file declaring the graph class "IncidentLists" as well as related functions
//
// Richard Restetzki

#include <cstdint>
#include <limits>
#include <utility>
#include <string>
#include <functional>
#include <random>
#include <tuple>
#include <queue>
#include <iostream>
#include <algorithm>
#include <set>

extern const u_int32_t seed;
extern double x_min,x_max,y_min,y_max,z_min,z_max;

using namespace std;

namespace tglib {

    using NodeId = int32_t;
    using EdgeId = int64_t;
    using Time = int32_t;
    using TimeInterval = std::pair<Time, Time>;
    using Depth = unsigned short;
    inline auto const inf = std::numeric_limits<Time>::max();

    struct TemporalEdge {
        // Tail of the edge
        NodeId u;

        // Head of the edge
        NodeId v;

        // The timepoint of the edge
        Time t;

        // The event an edge is assigned to
        int event;

        // Storing the information if an edge was already traversed by some method
        // This must always be set to false at the end of any function!
        bool traversed = false;

        // A vector of colors (integers) of the edge with the 0-th entry holding the starting color,
        // the i-th entry holding the color of the edge after the i-th iteration of WL
        std::vector<EdgeId> colors;
    };

    // A struct representing the unique identifier of a color (depth,color_id)
    struct color {
        // WL-Iteration depth corresponding to the color
        Depth depth;
        // Unique color_id within for the respective depth
        EdgeId color_id;

        inline bool operator==(const color &other) const {
            if((this->depth==other.depth) && (this->color_id==other.color_id)){
                return true;
            }
            return false;
        }

        inline bool operator<(const color &other) const {
            if(this->depth<other.depth || (this->depth==other.depth && this->color_id<other.color_id)){
                return true;
            }
            return false;
        }
    };

    // A container for edge neighborhood information
    struct Edge_NH {
         // Vector of incident edge color on side 1 - this should always be the smaller vector (size, then values)
        std::vector<color> Side_1;
         // Vector of incident edge color on side 2
        std::vector<color> Side_2;
        // The color of the edge itself
        color self_color;
        // Sort Side_1 and Side_2 s.t. Side_1 is smaller
        void sort();

        inline bool operator==(const Edge_NH &other) const {
            if ((this->Side_1 == other.Side_1) && (this->Side_2 == other.Side_2) && this->self_color == other.self_color) {
                return true;
            }
            return false;
        }
    };

    template<typename E>
    struct TGNodeT {
        NodeId id{};
        std::string label{};
        // A vector storing all the outgoing edges
        std::vector<E> outEdges;
        std::vector<E> inEdges;
    };

    // Setting the used node and edge type
    using TGNode = TGNodeT<TemporalEdge>;
    using TGEdge = TemporalEdge;

    // A struct containing all the color mappings
    struct Color_Mapping {
        // A map converting the labels (string) of the two incident vertices of an edge to an edge color (integer).
        // The map stores all currently seen label combinations.
        std::unordered_map<std::string, EdgeId> Initial_Color;

        // For every depth a map converting the color neighborhood (vector of integer numbers) of the edge to a new edge color (integer).
        // The map stores all currently seen label combinations.
        std::vector< std::unordered_map<Edge_NH, EdgeId> > Color_Mappings;

        // The maximum depth of the WL-iterations
        Depth depth;

        // The amount of timesteps at the end of a simulation which are cut off
        Time temp_range;

        void WL_color_reset();
        void clear();
        void Bucket_collisions();
    };

    // Representation of the temporal graph as the incident lists of all vertices
    // The incident lists are in temporally ascending order
    class IncidentLists {

    public:
        // Constructors
        IncidentLists() = default;
        explicit IncidentLists(NodeId n);

        // Access functions
        [[nodiscard]] size_t getNumberOfNodes() const;
        [[nodiscard]] EdgeId getNumberOfEdges() const;
        [[nodiscard]] TimeInterval getTimeInterval() const;
        TGNode& getNode(size_t nid);
        std::vector<TGNode> &getNodes();
        std::vector<TemporalEdge>::iterator findEdge_out(const TGEdge& e);
        std::vector<TemporalEdge>::iterator findEdge_in(const TGEdge& e);
        void get_color_outEdges(std::vector<TGEdge*>& edges, color color);

        // Reset functions
        void reset_Edges();
        void WL_color_reset();
        void clear();

        // Basic vertex functions
        void addNodes(NodeId num);
        void setNodeLabel(NodeId nid, const char* label);
        void Vertex_neighborhood(NodeId node, std::vector<color>& NH, int iteration, Time start_time, Time end_time, bool outgoing);

        // Basic edge functions
        void addEdge(const TGEdge& e);
        void insertEdge(const TGEdge& e);
        int traversed_edge(const TemporalEdge& edge);
        void Edge_neighborhood(const TemporalEdge& edge, Edge_NH& NH, int iteration, Time temp_range, bool directed = true);

        // Basic functionalities
        void set_ti(Time start, Time end);
        void shift_time_labels(Time shift);
        void initialize_Event(const char* const& input_file);
        void colorEdges_VertBased(Color_Mapping& Color_Mapping, bool direction_sensitive = false);
        void colorEdges_VertBased(Color_Mapping& Color_Mapping, const string& input_file);
        void colorEdges_uniformly(Color_Mapping& Color_Mapping);
        void randomize_edges(Color_Mapping& Color_Mapping, std::mt19937& generator, int count);
        void remove_small_cells(Color_Mapping& Color_Mapping, std::vector<size_t>& hist_vec, int min_count);

        // Main functionalities
        bool BFS(TemporalEdge& starting_edge, Time temp_range, NodeId depth, EdgeId& edges_count, std::function<void(tglib::TGEdge&)> func, bool color_search = true);
        void Edge_WL(Color_Mapping& Color_Mapping, Depth depth, Time temp_range, bool directed = true);

    private:
        std::vector<TGNode> nodes;
        EdgeId num_edges{};
        TimeInterval ti;
    };
} // tglib

// A hash function for Edge_NH enabling the use of an unorganized map
namespace std {

    template<>
    struct hash<tglib::Edge_NH>{
        inline size_t operator()(const tglib::Edge_NH& NH) const {
            // Create a key that is unique for the Edge_NH
            string key;
            key += std::to_string(NH.self_color.depth) + "-" + std::to_string(NH.self_color.color_id) + ":";
            for(auto i: NH.Side_1){
                key += std::to_string(i.depth) + "-" + std::to_string(i.color_id) + ";";
            }
            // Replace the last delimiter with a "&" indicating the switch to the second neighborhood
            key.pop_back();
            key += "&";
            for(auto i: NH.Side_2){
                key += std::to_string(i.depth) + "-" + std::to_string(i.color_id) + ";";
            }
            return hash<std::string>()(key);
        }
    };
} // std

// Basic functionalities
size_t sieve(tglib::Color_Mapping& Color_Mapping, std::vector<tglib::IncidentLists *> &graphs, tglib::Depth depth, int min_count);
void assign_Moment(tglib::IncidentLists& temp_graph, tglib::TGEdge& edge, int event, std::vector<tglib::TGEdge*>& edges);
void no_Assignment(tglib::TGEdge& edge);
void create_Moment(tglib::IncidentLists& temp_graph, tglib::TGEdge* edge, int event, tglib::Depth depth);
bool pred(const tglib::TemporalEdge &a, const tglib::TemporalEdge &b);
bool comp(const vector<size_t> &a, const vector<size_t> &b);
double cos_distance(std::vector<std::vector<size_t>> &vec1, std::vector<std::vector<size_t>> &vec2);
double average_color_dist(tglib::IncidentLists& temp_graph, tglib::color color, tglib::Time temp_range);
tglib::EdgeId queue_add_neighbors(std::queue<tglib::TemporalEdge*>& queue, tglib::TGNode& Node, tglib::TimeInterval interval);
void frequent_colors(std::vector<std::vector<size_t>> &color_nr_freq, tglib::Color_Mapping &Color_Mapping, tglib::IncidentLists &temp_graph, tglib::TimeInterval pair);
tglib::Time frequent_interval(std::vector<size_t>& occurrence_count, size_t& maximum, tglib::Time timespan);
void find_frequent_intervals(std::vector<tglib::Time>& start_times, std::vector<size_t>& maxima, tglib::IncidentLists& temp_graph, tglib::color color, std::vector<tglib::Time> timespans);
int check_change(tglib::IncidentLists& temp_graph, tglib::NodeId nid);
void color_hist(std::vector<std::vector<size_t>>& hist_vec, tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, tglib::TimeInterval interval = tglib::TimeInterval (std::numeric_limits<tglib::Time>::min(),std::numeric_limits<tglib::Time>::max()) );
void color_hist(std::vector<size_t>& hist_vec, tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, tglib::Depth depth, tglib::TimeInterval interval = tglib::TimeInterval (std::numeric_limits<tglib::Time>::min(),std::numeric_limits<tglib::Time>::max()) );
void event_color_hist(std::vector<std::vector<size_t>>& hist_vec, tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, int event);
double interval_similarity(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, tglib::TimeInterval interval1, tglib::TimeInterval interval2);

// Output functions
void print_color_densities(const char* const& input_file, tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, size_t min_occurrence = 1000, bool cut_degenerated_colors = true);
void plot_labels(const char* const& input_file, tglib::IncidentLists& temp_graph, tglib::color color, tglib::Time start_time, tglib::Time timespan);
void plot_event(string input_file, tglib::IncidentLists& temp_graph, int event);
void output_rd_similarities(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, int n, int randomization_step);
void output_similarities(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, int precision = 4);

// Main functionalities
void TNWL(tglib::Color_Mapping& Color_Mapping, std::vector<tglib::IncidentLists *> graphs, tglib::Depth depth, tglib::Time temp_range, bool directed = true, int min_count = 0);
void print_gram_of_events(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& temp_graph, int nr_events);
void label_analysis(tglib::Color_Mapping& Color_Mapping, tglib::IncidentLists& graph_1, string input_file_1, tglib::IncidentLists& graph_2, string input_file_2);
double TNWL_similarity(tglib::IncidentLists& graph_1, tglib::IncidentLists& graph_2, tglib::Color_Mapping& Color_Mapping);

#endif //SIMILARITYMEASURE_TG_H