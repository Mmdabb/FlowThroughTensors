#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>
#include <sstream>
#include "Json.hpp"

using json = nlohmann::json;

// CSV Parser from the provided code
class CDTACSVParser {
public:
    char Delimiter;
    // for DataHub CSV files

    std::ifstream inFile;
    std::string mFileName;

    std::vector<std::string> LineFieldsValue;
    std::vector<int> LineIntegerVector;
    std::vector<std::string> Headers;
    std::map<std::string, int> FieldsIndices;

    CDTACSVParser()
        : Delimiter{ ',' }
    {
    }

    ~CDTACSVParser()
    {
        if (inFile.is_open())
            inFile.close();
    }





// inline member functions
std::vector<std::string> GetHeaderVector() { return Headers; }
void CloseCSVFile() { inFile.close(); }

int ParserIntSequence(std::string str, std::vector<int>& vect)
{

    std::stringstream ss(str);

    int i;

    while (ss >> i)
    {
        vect.push_back(i);

        if (ss.peek() == ';')
            ss.ignore();
    }

    return vect.size();
}


// definitions of CDTACSVParser member functions
void ConvertLineStringValueToIntegers()
{
    LineIntegerVector.clear();
    for (unsigned i = 0; i < LineFieldsValue.size(); ++i)
    {
        std::string si = LineFieldsValue[i];
        int value = atoi(si.c_str());

        if (value >= 1)
            LineIntegerVector.push_back(value);
    }
}

bool OpenCSVFile(std::string fileName, bool b_required)
{
    mFileName = fileName;
    inFile.open(fileName.c_str());

    if (inFile.is_open())
    {
            std::string s;
            std::getline(inFile, s);
            std::vector<std::string> FieldNames = ParseLine(s);
            Headers = FieldNames;

            for (size_t i = 0; i < FieldNames.size(); i++)
            {
                std::string tmp_str = FieldNames.at(i);
                size_t start = tmp_str.find_first_not_of(" ");

                std::string name;
                if (start == std::string::npos)
                {
                    name = "";
                }
                else
                {
                    name = tmp_str.substr(start);
                }
                FieldsIndices[name] = (int)i;
            }

        return true;
    }
    else
    {
        if (b_required)
        {
            std::cerr << "Error: Could not open file " << fileName << std::endl;
        }
        return false;
    }
}

bool ReadRecord()
{
    LineFieldsValue.clear();

    if (inFile.is_open())
    {
        std::string s;
        std::getline(inFile, s);
        if (s.length() > 0)
        {
            LineFieldsValue = ParseLine(s);
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}



std::vector<std::string> ParseLine(std::string line)
{
    std::vector<std::string> SeperatedStrings;
    std::string subStr;

    if (line.length() == 0)
        return SeperatedStrings;

    std::istringstream ss(line);

    if (line.find_first_of('"') == std::string::npos)
    {
        while (std::getline(ss, subStr, Delimiter))
        {
            SeperatedStrings.push_back(subStr);
        }

        if (line.at(line.length() - 1) == ',')
        {
            SeperatedStrings.push_back("");
        }
    }
    else
    {
        while (line.length() > 0)
        {
            size_t n1 = line.find_first_of(',');
            size_t n2 = line.find_first_of('"');

            if (n1 == std::string::npos &&
                n2 == std::string::npos)  // last field without double quotes
            {
                subStr = line;
                SeperatedStrings.push_back(subStr);
                break;
            }

            if (n1 == std::string::npos && n2 != std::string::npos)  // last field with double
                // quotes
            {
                size_t n3 = line.find_first_of('"', n2 + 1);  // second double quote

                // extract content from double quotes
                subStr = line.substr(n2 + 1, n3 - n2 - 1);
                SeperatedStrings.push_back(subStr);

                break;
            }

            if (n1 != std::string::npos && (n1 < n2 || n2 == std::string::npos))
            {
                subStr = line.substr(0, n1);
                SeperatedStrings.push_back(subStr);
                if (n1 < line.length() - 1)
                {
                    line = line.substr(n1 + 1);
                }
                else  // comma is the last char in the line string, push an empty string to the back
                      // of vector
                {
                    SeperatedStrings.push_back("");
                    break;
                }
            }

            if (n1 != std::string::npos && n2 != std::string::npos && n2 < n1)
            {
                size_t n3 = line.find_first_of('"', n2 + 1);  // second double quote
                subStr = line.substr(n2 + 1, n3 - n2 - 1);
                SeperatedStrings.push_back(subStr);
                size_t idx = line.find_first_of(',', n3 + 1);

                if (idx != std::string::npos)
                {
                    line = line.substr(idx + 1);
                }
                else
                {
                    break;
                }
            }
        }
    }
    return SeperatedStrings;
}

bool GetStringValueByFieldName(std::string field_name,
    std::string& value,
    bool required_field)
{
    if (FieldsIndices.find(field_name) == FieldsIndices.end())
    {
        if (required_field)
        {
            std::cerr << "[ERROR] Field " << field_name << " in file " << mFileName << " does not exist. Please check the file." << std::endl;
        }
        return false;
    }
    else
    {
        if (LineFieldsValue.size() == 0)
        {
            return false;
        }

        unsigned int index = FieldsIndices[field_name];
        if (index >= LineFieldsValue.size())
        {
            return false;
        }
        std::string str_value = LineFieldsValue[index];

        if (str_value.length() <= 0)
        {
            return false;
        }

        value = str_value;
        return true;
    }
}

template <class T>
bool GetValueByFieldName(std::string field_name,
    T& value,
    bool required_field)
{
    if (FieldsIndices.find(field_name) == FieldsIndices.end())
    {
        if (required_field)
        {
            std::cerr << "[ERROR] Field " << field_name << " in file " << mFileName.c_str() << " does not exist. Please check the file." << std::endl;
        }
        return false;
    }
    else
    {
        if (LineFieldsValue.size() == 0)
        {
            return false;
        }

        int size = (int)(LineFieldsValue.size());
        if (FieldsIndices[field_name] >= size)
        {
            return false;
        }

        std::string str_value = LineFieldsValue[FieldsIndices[field_name]];

        if (str_value.length() <= 0)
        {
            return false;
        }

        std::istringstream ss(str_value);

        T converted_value;
        ss >> converted_value;

        if (ss.fail())
        {
            return false;
        }

        if (required_field)
        {
  //          if (converted_value < 0)
  //              converted_value = 0;
        }

        value = converted_value;
        return true;
    }
}
};

// Forward declarations
class Node;
class Link;
class Route;
class ODPair;
class Network;
class TrafficAssignment;

// Node class
class Node {
public:
    int node_id;
    int zone_id;
    double x_coord;
    double y_coord;

    Node(int id, int zone, double x, double y) : node_id(id), zone_id(zone), x_coord(x), y_coord(y) {}
};

// Link class with BPR function
class Link {
public:
    int link_id;
    int from_node_id;
    int to_node_id;
    double free_flow_time;
    int lanes;
    double capacity;
    double vdf_alpha;
    double vdf_beta;
    double volume;
    double current_time;

    Link(int id, int from, int to, double fftt, int num_lanes, double cap, double alpha, double beta)
        : link_id(id), from_node_id(from), to_node_id(to), free_flow_time(fftt),
        lanes(num_lanes), capacity(cap), vdf_alpha(alpha), vdf_beta(beta),
        volume(0.0), current_time(fftt) {}

    // BPR function to calculate travel time based on volume
    double calculateTravelTime() {
        current_time = free_flow_time * (1.0 + vdf_alpha * pow(volume / (capacity * lanes), vdf_beta));
        return current_time;
    }
};

// Route class
class Route {
public:
    int mode;
    int origin_zone;
    int dest_zone;
    int route_id;
    std::vector<int> link_ids;
    double flow;
    double cost;

    Route(int m, int o, int d, int id, const std::vector<int>& links)
        : mode(m), origin_zone(o), dest_zone(d), route_id(id), link_ids(links), flow(0.0), cost(0.0) {}

    // Calculate route cost based on link travel times
    double calculateCost(const std::unordered_map<int, Link>& links) {
        cost = 0.0;
        for (int link_id : link_ids) {
            auto it = links.find(link_id);
            if (it != links.end()) {
                cost += it->second.current_time;
            }
        }
        return cost;
    }
};

// Origin-Destination pair class
class ODPair {
public:
    int mode;
    int origin;
    int destination;
    double demand;
    std::vector<int> route_indices;

    ODPair(int m, int o, int d, double dem = 0.0)
        : mode(m), origin(o), destination(d), demand(dem) {}
};

// Network class to hold all network elements
class Network {
public:
    std::unordered_map<int, Node> nodes;
    std::unordered_map<int, Link> links;
    std::vector<Route> routes;
    std::vector<ODPair> od_pairs;

    int num_modes;
    int num_nodes;
    int num_zones;
    int num_links;

    Network() : num_modes(0), num_nodes(0), num_zones(0), num_links(0) {}

    // Load network from JSON file
    bool loadFromJSON(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        json data;
        try {
            file >> data;
        }
        catch (const std::exception& e) {
            std::cerr << "Error parsing JSON: " << e.what() << std::endl;
            return false;
        }

        // Read settings
        if (data.contains("settings")) {
            num_modes = data["settings"].value("num_modes", 1);
        }

        // Read nodes
        if (data.contains("nodes")) {
            for (const auto& node_data : data["nodes"]) {
                int id = node_data["node_id"];
                int zone = node_data["zone_id"];
                double x = node_data["x_coord"];
                double y = node_data["y_coord"];

                nodes.emplace(id, Node(id, zone, x, y));
            }
            num_nodes = nodes.size();

            // Determine number of zones
            num_zones = 0;
            for (const auto& node_pair : nodes) {
                if (node_pair.second.zone_id > num_zones) {
                    num_zones = node_pair.second.zone_id;
                }
            }
        }

        // Read links
        if (data.contains("links")) {
            for (const auto& link_data : data["links"]) {
                int id = link_data["link_id"];
                int from = link_data["from_node_id"];
                int to = link_data["to_node_id"];
                double fftt = link_data["fftt"];
                int lanes = link_data["lanes"];
                double capacity = link_data["capacity"];
                double alpha = link_data["vdf_alpha"];
                double beta = link_data["vdf_beta"];

                links.emplace(id, Link(id, from, to, fftt, lanes, capacity, alpha, beta));
            }
            num_links = links.size();
        }

        // Read routes
        if (data.contains("routes")) {
            for (const auto& route_data : data["routes"]) {
                int mode = route_data["mode"];
                int origin = route_data["o_zone_id"];
                int dest = route_data["d_zone_id"];
                int route_id = route_data["route_id"];
                std::vector<int> link_ids;

                for (const auto& link_id : route_data["link_ids"]) {
                    link_ids.push_back(link_id);
                }

                routes.emplace_back(mode, origin, dest, route_id, link_ids);

                // Check if OD pair exists
                bool found = false;
                for (auto& od : od_pairs) {
                    if (od.mode == mode && od.origin == origin && od.destination == dest) {
                        od.route_indices.push_back(routes.size() - 1);
                        found = true;
                        break;
                    }
                }

                // If OD pair doesn't exist, create it
                if (!found) {
                    ODPair od(mode, origin, dest);
                    od.route_indices.push_back(routes.size() - 1);
                    od_pairs.push_back(od);
                }
            }
        }

        // Read OD demands if available
        if (data.contains("od_demands")) {
            for (const auto& demand_data : data["od_demands"]) {
                int mode = demand_data["mode"];
                int origin = demand_data["origin"];
                int destination = demand_data["destination"];
                double demand_value = demand_data["demand"];

                // Find the corresponding OD pair and set demand
                for (auto& od : od_pairs) {
                    if (od.mode == mode && od.origin == origin && od.destination == destination) {
                        od.demand = demand_value;
                        break;
                    }
                }
            }
        }

        return true;
    }

    // Load network from CSV files using the robust parser
    bool loadFromCSV(const std::string& node_file, const std::string& link_file, const std::string& route_file) {
        // Load nodes
        if (!loadNodesFromCSV(node_file)) {
            return false;
        }

        // Load links
        if (!loadLinksFromCSV(link_file)) {
            return false;
        }

        // Load routes
        if (!loadRoutesFromCSV(route_file)) {
            return false;
        }

        return true;
    }

private:
    bool loadNodesFromCSV(const std::string& filename) {
        CDTACSVParser parser;
        if (!parser.OpenCSVFile(filename, true)) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        while (parser.ReadRecord()) {
            int node_id;
            int zone_id;
            double x_coord;
            double y_coord;

            if (parser.GetValueByFieldName("node_id", node_id, true) &&
                parser.GetValueByFieldName("zone_id", zone_id, true) &&
                parser.GetValueByFieldName("x_coord", x_coord, true) &&
                parser.GetValueByFieldName("y_coord", y_coord, true)) {

                nodes.emplace(node_id, Node(node_id, zone_id, x_coord, y_coord));
            }
        }

        num_nodes = nodes.size();

        // Determine number of zones
        num_zones = 0;
        for (const auto& node_pair : nodes) {
            if (node_pair.second.zone_id > num_zones) {
                num_zones = node_pair.second.zone_id;
            }
        }

        parser.CloseCSVFile();
        return true;
    }

    bool loadLinksFromCSV(const std::string& filename) {
        CDTACSVParser parser;
        if (!parser.OpenCSVFile(filename, true)) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        while (parser.ReadRecord()) {
            int link_id;
            int from_node_id;
            int to_node_id;
            double fftt;
            double lanes;
            double capacity;
            double vdf_alpha;
            double vdf_beta;

            if (parser.GetValueByFieldName("link_id", link_id, true) &&
                parser.GetValueByFieldName("from_node_id", from_node_id, true) &&
                parser.GetValueByFieldName("to_node_id", to_node_id, true) &&
                parser.GetValueByFieldName("fftt", fftt, true) &&
                parser.GetValueByFieldName("lanes", lanes, true) &&
                parser.GetValueByFieldName("capacity", capacity, true) &&
                parser.GetValueByFieldName("VDF_alpha", vdf_alpha, true) &&
                parser.GetValueByFieldName("VDF_beta", vdf_beta, true)) {

                links.emplace(link_id, Link(link_id, from_node_id, to_node_id, fftt, lanes, capacity, vdf_alpha, vdf_beta));
            }
        }

        num_links = links.size();
        parser.CloseCSVFile();
        return true;
    }

    bool loadRoutesFromCSV(const std::string& filename) {
        CDTACSVParser parser;
        if (!parser.OpenCSVFile(filename, true)) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }


         num_modes = 1; // Default to 1 mode if not specified

        while (parser.ReadRecord()) {
            int mode;
            int o_zone_id;
            int d_zone_id;
            int route_id;

            if (parser.GetValueByFieldName("mode", mode, true) &&
                parser.GetValueByFieldName("o_zone_id", o_zone_id, true) &&
                parser.GetValueByFieldName("d_zone_id", d_zone_id, true) &&
                parser.GetValueByFieldName("route_id", route_id, true)) {

                // For link_ids, we need to handle them specially since they're a sequence
                std::vector<int> link_ids;

                // Try to see if there's a field called "link_ids" that contains the sequence
                std::string link_ids_str;
                if (parser.GetStringValueByFieldName("link_ids", link_ids_str, false)) {

                    parser.ParserIntSequence(link_ids_str, link_ids);

                }

                if (!link_ids.empty()) {
                    routes.emplace_back(mode, o_zone_id, d_zone_id, route_id, link_ids);

                    // Check if OD pair exists
                    bool found = false;
                    for (auto& od : od_pairs) {
                        if (od.mode == mode && od.origin == o_zone_id && od.destination == d_zone_id) {
                            od.route_indices.push_back(routes.size() - 1);
                            found = true;
                            break;
                        }
                    }

                    // If OD pair doesn't exist, create it
                    if (!found) {
                        ODPair od(mode, o_zone_id, d_zone_id);
                        od.route_indices.push_back(routes.size() - 1);
                        od_pairs.push_back(od);
                    }
                }
                else {
                    std::cerr << "Warning: No link sequence found for route " << route_id
                        << " (mode=" << mode << ", o=" << o_zone_id << ", d=" << d_zone_id << ")" << std::endl;
                }
            }
        }

        parser.CloseCSVFile();
        return true;
    }
};

// Traffic Assignment class
class TrafficAssignment {
private:
    Network& network;
    double logit_theta;
    int max_iterations;
    double convergence_threshold;

public:
    TrafficAssignment(Network& net, double theta = 0.1, int max_iter = 100, double conv_threshold = 0.001)
        : network(net), logit_theta(theta), max_iterations(max_iter), convergence_threshold(conv_threshold) {}

    // Main assignment method
    bool assign() {
        // Initialize
        resetAssignment();

        for (int iter = 0; iter < max_iterations; iter++) {
            // 1. Calculate link travel times based on current volumes
            updateLinkTravelTimes();

            // 2. Calculate route costs
            updateRouteCosts();

            // 3. Calculate route flows using logit model
            std::vector<double> old_flows = getRouteFlows();
            updateRouteFlows();

            // 4. Update link volumes based on route flows
            updateLinkVolumes();

            // 5. Check convergence
            double gap = calculateGap(old_flows);
            std::cout << "Iteration " << iter << ", Gap: " << gap << std::endl;

            if (gap < convergence_threshold) {
                std::cout << "Converged after " << iter + 1 << " iterations." << std::endl;
                return true;
            }
        }

        std::cout << "Maximum iterations reached without convergence." << std::endl;
        return false;
    }

    // Reset assignment to initial state
    void resetAssignment() {
        // Reset link volumes
        for (auto& link_pair : network.links) {
            link_pair.second.volume = 0.0;
            link_pair.second.current_time = link_pair.second.free_flow_time;
        }

        // Reset route flows and costs
        for (auto& route : network.routes) {
            route.flow = 0.0;
            route.cost = 0.0;
        }
    }

    // Update link travel times based on current volumes
    void updateLinkTravelTimes() {
        for (auto& link_pair : network.links) {
            link_pair.second.calculateTravelTime();
        }
    }

    // Update route costs based on link travel times
    void updateRouteCosts() {
        for (auto& route : network.routes) {
            route.calculateCost(network.links);
        }
    }

    // Update route flows using logit model
    void updateRouteFlows() {
        for (auto& od : network.od_pairs) {
            double demand = od.demand;

            // Calculate denominator for logit model
            double exp_sum = 0.0;
            for (int route_idx : od.route_indices) {
                double cost = network.routes[route_idx].cost;
                exp_sum += exp(-logit_theta * cost);
            }

            // Calculate flow for each route
            for (int route_idx : od.route_indices) {
                double cost = network.routes[route_idx].cost;

                if (exp_sum > 0) {
                    network.routes[route_idx].flow = demand * exp(-logit_theta * cost) / exp_sum;
                }
                else {
                    // Distribute evenly if all costs are too high
                    network.routes[route_idx].flow = demand / od.route_indices.size();
                }
            }
        }
    }

    // Update link volumes based on route flows
    void updateLinkVolumes() {
        // Reset link volumes
        for (auto& link_pair : network.links) {
            link_pair.second.volume = 0.0;
        }

        // Accumulate flows
        for (const auto& route : network.routes) {
            for (int link_id : route.link_ids) {
                auto it = network.links.find(link_id);
                if (it != network.links.end()) {
                    it->second.volume += route.flow;
                }
            }
        }
    }

    // Get current route flows
    std::vector<double> getRouteFlows() {
        std::vector<double> flows;
        for (const auto& route : network.routes) {
            flows.push_back(route.flow);
        }
        return flows;
    }

    // Calculate gap between iterations
    double calculateGap(const std::vector<double>& old_flows) {
        double gap = 0.0;
        double total_flow = 0.0;

        for (size_t i = 0; i < old_flows.size(); i++) {
            gap += std::abs(network.routes[i].flow - old_flows[i]);
            total_flow += network.routes[i].flow;
        }

        return total_flow > 0 ? gap / total_flow : gap;
    }

    // API: Change link capacity
    void changeLinkCapacity(int link_id, double new_capacity) {
        auto it = network.links.find(link_id);
        if (it != network.links.end()) {
            it->second.capacity = new_capacity;
        }
    }

    // API: Change OD demand
    void changeODDemand(int mode, int origin, int destination, double new_demand) {
        for (auto& od : network.od_pairs) {
            if (od.mode == mode && od.origin == origin && od.destination == destination) {
                od.demand = new_demand;
                break;
            }
        }
    }

    // API: Output link volumes
    std::vector<std::pair<int, double>> getLinkVolumes() {
        std::vector<std::pair<int, double>> result;
        for (const auto& link_pair : network.links) {
            result.emplace_back(link_pair.first, link_pair.second.volume);
        }
        return result;
    }

    // API: Output route costs
    std::vector<std::tuple<int, int, int, int, double>> getRouteCosts() {
        std::vector<std::tuple<int, int, int, int, double>> result;
        for (const auto& route : network.routes) {
            result.emplace_back(route.mode, route.origin_zone, route.dest_zone, route.route_id, route.cost);
        }
        return result;
    }

};

    // Example usage
    int main() {
        // Create network
        Network network;

        // Load network (choose one method)
        // 1. Load from JSON
       // bool loaded = network.loadFromJSON("network.json");

        // 2. Load from CSV files
         bool loaded = network.loadFromCSV("node.csv", "link.csv", "route_assignment.csv");

        if (!loaded) {
            std::cerr << "Failed to load network." << std::endl;
            return 1;
        }

        // Create traffic assignment
        TrafficAssignment assignment(network, 0.1);

        // Run assignment
        assignment.assign();

        // Output results
        auto link_volumes = assignment.getLinkVolumes();
        auto route_costs = assignment.getRouteCosts();

        std::cout << "Link Volumes:" << std::endl;
        for (const auto& pair : link_volumes) {
            int link_id = pair.first;
            double volume = pair.second;
            std::cout << "Link " << link_id << ": " << volume << std::endl;
        }

        std::cout << "Route Costs:" << std::endl;
        for (const auto& tuple : route_costs) {
            int mode = std::get<0>(tuple);
            int origin = std::get<1>(tuple);
            int dest = std::get<2>(tuple);
            int route_id = std::get<3>(tuple);
            double cost = std::get<4>(tuple);

            std::cout << "Mode " << mode << ", O " << origin << ", D " << dest
                << ", Route " << route_id << ": " << cost << std::endl;
        }

        // Export results
      //  exportNetworkToJSON(network, "network_results.json");

        return 0;
    }
