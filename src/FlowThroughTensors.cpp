#include "Json.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

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

    CDTACSVParser() : Delimiter{','} {}

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
                    else  // comma is the last char in the line string, push an empty string to the
                          // back of vector
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

    bool GetStringValueByFieldName(std::string field_name, std::string& value, bool required_field)
    {
        if (FieldsIndices.find(field_name) == FieldsIndices.end())
        {
            if (required_field)
            {
                std::cerr << "[ERROR] Field " << field_name << " in file " << mFileName
                          << " does not exist. Please check the file." << std::endl;
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
    bool GetValueByFieldName(std::string field_name, T& value, bool required_field)
    {
        if (FieldsIndices.find(field_name) == FieldsIndices.end())
        {
            if (required_field)
            {
                std::cerr << "[ERROR] Field " << field_name << " in file " << mFileName.c_str()
                          << " does not exist. Please check the file." << std::endl;
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
                double fftt = link_data["VDF_fftt"];
                int lanes = link_data["lanes"];
                double capacity = link_data["capacity"];
                double alpha = link_data["VDF_alpha"];
                double beta = link_data["VDF_beta"];

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
        std::cout << "Loading nodes from file: " << filename << std::endl;

        CDTACSVParser parser;
        if (!parser.OpenCSVFile(filename, true)) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        int count = 0;
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
                count++;

                // Print progress every 100 nodes
                if (count % 100 == 0) {
                    std::cout << "  Loaded " << count << " nodes so far..." << std::endl;
                }
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
        std::cout << "Successfully loaded " << num_nodes << " nodes with " << num_zones << " zones" << std::endl;
        return true;
    }

    bool loadLinksFromCSV(const std::string& filename) {
        std::cout << "Loading links from file: " << filename << std::endl;
        CDTACSVParser parser;
        if (!parser.OpenCSVFile(filename, true)) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        int count = 0;

        // Print the headers
        std::cout << "\nFirst few links with detailed attributes:" << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << "Link ID | From Node | To Node | Free Flow Time | Lanes | Capacity | Alpha | Beta" << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;

        while (parser.ReadRecord()) {
            int from_node_id;
            int to_node_id;


            // Set default values first
            int link_id = count + 1; // Auto-generate link ID if not found
            double fftt = 1.0;       // Default free flow travel time (1 minute)
            double lanes = 1.0;      // Default 1 lane
            double capacity = 1800.0; // Default capacity (1800 vehicles per hour per lane)
            double vdf_alpha = 0.15;  // Default BPR alpha parameter
            double vdf_beta = 4.0;    // Default BPR beta parameter

            // Try to read link_id (not required)
            parser.GetValueByFieldName("link_id", link_id, false);

            // Required fields - only from_node_id and to_node_id
            if (parser.GetValueByFieldName("from_node_id", from_node_id, true) &&
                parser.GetValueByFieldName("to_node_id", to_node_id, true)) {

                // Non-required fields with defaults - using GMNS standard names
                parser.GetValueByFieldName("VDF_fftt", fftt, false);
                parser.GetValueByFieldName("lanes", lanes, false);
                parser.GetValueByFieldName("capacity", capacity, false);
                parser.GetValueByFieldName("VDF_alpha", vdf_alpha, false);
                parser.GetValueByFieldName("VDF_beta", vdf_beta, false);

                // Check for reasonable values and apply constraints
                if (fftt <= 0) fftt = 1.0;        // Ensure positive travel time
                if (lanes <= 0) lanes = 1.0;      // Ensure at least one lane
                if (capacity <= 0) capacity = 1800.0; // Ensure positive capacity
                if (vdf_alpha < 0) vdf_alpha = 0.15;  // Ensure non-negative alpha
                if (vdf_beta < 0) vdf_beta = 4.0;     // Ensure non-negative beta

                // Create the link
                links.emplace(link_id, Link(link_id, from_node_id, to_node_id, fftt, lanes, capacity, vdf_alpha, vdf_beta));

                // Print detailed info for first 5 links
                if (count < 5) {
                    std::cout << std::setw(7) << link_id << " | ";
                    std::cout << std::setw(9) << from_node_id << " | ";
                    std::cout << std::setw(7) << to_node_id << " | ";
                    std::cout << std::setw(14) << std::fixed << std::setprecision(2) << fftt << " | ";
                    std::cout << std::setw(5) << std::fixed << std::setprecision(0) << lanes << " | ";
                    std::cout << std::setw(8) << std::fixed << std::setprecision(0) << capacity << " | ";
                    std::cout << std::setw(5) << std::fixed << std::setprecision(2) << vdf_alpha << " | ";
                    std::cout << std::setw(4) << std::fixed << std::setprecision(2) << vdf_beta << std::endl;
                }

                // Print progress every 100 links
                if (count > 0 && count % 100 == 0) {
                    std::cout << "Loaded " << count << " links so far..." << std::endl;
                }

                count++;
            }
            else {
                std::cout << "Warning: Skipping link due to missing required fields at line "
                    << (count + 1) << std::endl;
            }



        }

        std::cout << "-----------------------------------------------------------------------" << std::endl;

        num_links = links.size();
        parser.CloseCSVFile();

        // Calculate total capacity stats
        double total_capacity = 0.0;
        double min_capacity = std::numeric_limits<double>::max();
        double max_capacity = 0.0;

        for (const auto& link_pair : links) {
            double link_capacity = link_pair.second.capacity * link_pair.second.lanes;
            total_capacity += link_capacity;
            min_capacity = std::min(min_capacity, link_capacity);
            max_capacity = std::max(max_capacity, link_capacity);
        }

        double avg_capacity = num_links > 0 ? total_capacity / num_links : 0;

        std::cout << "Successfully loaded " << num_links << " links" << std::endl;
        std::cout << "Capacity statistics:" << std::endl;
        std::cout << "  - Total: " << total_capacity << std::endl;
        std::cout << "  - Average: " << avg_capacity << " per link" << std::endl;
        std::cout << "  - Min: " << min_capacity << std::endl;
        std::cout << "  - Max: " << max_capacity << std::endl;

        return true;
    }
    bool loadRoutesFromCSV(const std::string& filename) {
        std::cout << "Loading routes from file: " << filename << std::endl;

        CDTACSVParser parser;
        if (!parser.OpenCSVFile(filename, true)) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        num_modes = 1; // Default to 1 mode if not specified
        int count = 0;
        int skipped = 0;
        double total_initial_volume = 0.0;

        // Print header for first few routes
        std::cout << "\nFirst few routes with detailed attributes:" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
        std::cout << "Route ID | Mode | Origin | Destination | # Links | Initial Volume | First 3 Links" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------------" << std::endl;

        while (parser.ReadRecord()) {
            // Default values
            int mode = 1;
            int o_zone_id = 0;
            int d_zone_id = 0;
            int route_id = count + 1;
            double initial_volume = 0.0;

            // Try to read values - only require origin and destination
            parser.GetValueByFieldName("mode", mode, false);
            bool has_origin = parser.GetValueByFieldName("o_zone_id", o_zone_id, false);
            bool has_dest = parser.GetValueByFieldName("d_zone_id", d_zone_id, false);
            parser.GetValueByFieldName("route_id", route_id, false);

            // Try to read initial volume if available
            parser.GetValueByFieldName("volume", initial_volume, false);
            parser.GetValueByFieldName("flow", initial_volume, false); // Alternative field name

            // Check if we have the minimum required fields
            if (has_origin && has_dest) {
                // For link_ids, we need to handle them specially since they're a sequence
                std::vector<int> link_ids;
                std::string link_ids_str;

                if (parser.GetValueByFieldName("link_ids", link_ids_str, false)) {
                    parser.ParserIntSequence(link_ids_str, link_ids);
                }

                if (!link_ids.empty()) {
                    // Create the route with the initial volume as flow
                    Route route(mode, o_zone_id, d_zone_id, route_id, link_ids);
                    route.flow = initial_volume; // Set initial flow/volume
                    routes.push_back(route);

                    total_initial_volume += initial_volume;

                    // Print detailed info for first 5 routes
                    if (count < 5) {
                        std::cout << std::setw(8) << route_id << " | ";
                        std::cout << std::setw(4) << mode << " | ";
                        std::cout << std::setw(7) << o_zone_id << " | ";
                        std::cout << std::setw(11) << d_zone_id << " | ";
                        std::cout << std::setw(7) << link_ids.size() << " | ";
                        std::cout << std::setw(14) << std::fixed << std::setprecision(2) << initial_volume << " | ";

                        // Print first 3 links (or fewer if route has fewer)
                        for (size_t i = 0; i < std::min(size_t(3), link_ids.size()); i++) {
                            std::cout << link_ids[i];
                            if (i < std::min(size_t(2), link_ids.size() - 1)) {
                                std::cout << ", ";
                            }
                        }
                        if (link_ids.size() > 3) {
                            std::cout << ", ...";
                        }
                        std::cout << std::endl;
                    }

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

                    // Print progress every 100 routes
                    if (count > 0 && count % 100 == 0) {
                        std::cout << "Loaded " << count << " routes so far..." << std::endl;
                    }

                    count++;
                }
                else {
                    std::cerr << "Warning: No link sequence found for route " << route_id
                        << " (mode=" << mode << ", o=" << o_zone_id << ", d=" << d_zone_id << ")" << std::endl;
                    skipped++;
                }
            }
            else {
                std::cerr << "Warning: Skipping route due to missing origin/destination at line "
                    << (count + skipped + 1) << std::endl;
                skipped++;
            }
        }

        std::cout << "-----------------------------------------------------------------------------------------" << std::endl;

        parser.CloseCSVFile();

        // Print summary statistics
        std::cout << "Successfully loaded " << count << " routes" << std::endl;
        std::cout << "Skipped " << skipped << " invalid routes" << std::endl;
        std::cout << "Total initial route volume: " << total_initial_volume << std::endl;
        std::cout << "Route statistics:" << std::endl;

        // Calculate some route statistics
        if (count > 0) {
            int min_links = INT_MAX;
            int max_links = 0;
            double total_links = 0;
            std::map<std::pair<int, int>, int> od_route_count;
            double max_initial_volume = 0.0;
            double routes_with_volume = 0;

            for (const auto& route : routes) {
                int route_links = route.link_ids.size();
                min_links = std::min(min_links, route_links);
                max_links = std::max(max_links, route_links);
                total_links += route_links;

                // Count routes per OD pair
                od_route_count[std::make_pair(route.origin_zone, route.dest_zone)]++;

                // Track volume statistics
                if (route.flow > 0) {
                    routes_with_volume++;
                    max_initial_volume = std::max(max_initial_volume, route.flow);
                }
            }

            double avg_links = total_links / count;
            int max_routes_per_od = 0;
            int od_with_one_route = 0;

            for (const auto& od_count : od_route_count) {
                max_routes_per_od = std::max(max_routes_per_od, od_count.second);
                if (od_count.second == 1) {
                    od_with_one_route++;
                }
            }

            std::cout << "  - Number of OD pairs: " << od_pairs.size() << std::endl;
            std::cout << "  - Avg. links per route: " << avg_links << std::endl;
            std::cout << "  - Min links in a route: " << min_links << std::endl;
            std::cout << "  - Max links in a route: " << max_links << std::endl;
            std::cout << "  - Max routes for a single OD pair: " << max_routes_per_od << std::endl;
            std::cout << "  - OD pairs with only one route: " << od_with_one_route
                << " (" << (100.0 * od_with_one_route / od_pairs.size()) << "%)" << std::endl;
            std::cout << "  - Routes with initial volume: " << routes_with_volume
                << " (" << (100.0 * routes_with_volume / count) << "%)" << std::endl;
            std::cout << "  - Maximum initial route volume: " << max_initial_volume << std::endl;
        }

        // After loading all routes, calculate OD demands based on initial volumes
        std::cout << "\nCalculating OD demands from initial route volumes..." << std::endl;
        int updated_od_pairs = 0;
        double total_demand = 0.0;

        for (auto& od : od_pairs) {
            double total_od_volume = 0.0;

            // Sum up all route volumes for this OD pair
            for (int route_idx : od.route_indices) {
                total_od_volume += routes[route_idx].flow;
            }

            // Update OD demand based on sum of route volumes
            if (total_od_volume > 0) {
                od.demand = total_od_volume;
                total_demand += total_od_volume;
                updated_od_pairs++;
            }
        }

        std::cout << "Updated demands for " << updated_od_pairs << " OD pairs based on initial route volumes." << std::endl;
        std::cout << "Total system demand from initial volumes: " << total_demand << std::endl;

        // Print first 5 OD pairs with demand for verification
        if (updated_od_pairs > 0) {
            std::cout << "\nSample OD pairs with demand:" << std::endl;
            std::cout << "------------------------------------------------------" << std::endl;
            std::cout << "Mode | Origin | Destination | Demand | Routes" << std::endl;
            std::cout << "------------------------------------------------------" << std::endl;

            int shown = 0;
            for (const auto& od : od_pairs) {
                if (od.demand > 0 && shown < 5) {
                    std::cout << std::setw(4) << od.mode << " | "
                        << std::setw(6) << od.origin << " | "
                        << std::setw(11) << od.destination << " | "
                        << std::setw(6) << std::fixed << std::setprecision(1) << od.demand << " | "
                        << od.route_indices.size() << std::endl;
                    shown++;
                }
            }
            std::cout << "------------------------------------------------------" << std::endl;
        }

        return true;
    }
};


// Add this to your TrafficAssignment class
class TrafficAssignmentLog {
private:
    Network& network;
    double logit_theta;
    int max_iterations;
    double convergence_threshold;

    // Logging options
    bool enable_logging;
    std::string log_directory;
    int log_level; // 0: minimal, 1: basic, 2: detailed, 3: verbose

    double previous_gap = std::numeric_limits<double>::max();
    double last_gap_improvement = 0.0;

public:
    TrafficAssignmentLog(Network& net, double theta = 0.1, int max_iter = 100, double conv_threshold = 0.001,
        bool logging = false, const std::string& log_dir = "./logs/", int verbose_level = 1)
        : network(net), logit_theta(theta), max_iterations(max_iter), convergence_threshold(conv_threshold),
        enable_logging(logging), log_directory(log_dir), log_level(verbose_level) {

        // Create log directory if it doesn't exist and logging is enabled
        if (enable_logging) {
            createLogDirectory();
        }
    }

    // Main assignment method
    bool assign() {
        // Initialize
        resetAssignment();

        // Initial logs
        if (enable_logging) {
            logNetworkState("initial");
        }

        for (int iter = 0; iter < max_iterations; iter++) {

            // 1. Update link volumes based on  route flows (initial route flow  at iteration 0)
            updateLinkVolumes();

            if (enable_logging && log_level >= 1) {
                logLinkVolumes(iter);
            }

            // 2. Calculate link travel times based on current volumes
            updateLinkTravelTimes();

            if (enable_logging && log_level >= 2) {
                logLinkTravelTimes(iter);
            }

            // 3. Calculate route costs
            updateRouteCosts();


            if (enable_logging && log_level >= 2) {
                logRouteCosts(iter);
            }

            // 4. Calculate route flows using logit model
            std::vector<double> old_flows = getRouteFlows();
            updateRouteFlows(iter);

            if (enable_logging && log_level >= 1) {
                logRouteFlows(iter);
                logODFlows(iter);
            }



            // 5. Check convergence
            double gap = calculateEquilibriumGap();
            std::cout << "Iteration " << iter << ", Gap: " << std::fixed << std::setprecision(2) << (gap * 100.0) << "%" << std::endl;

            if (enable_logging) {
                logIterationSummary(iter, gap);
            }

            if (gap < convergence_threshold) {
                std::cout << "Converged after " << iter + 1 << " iterations." << std::endl;

                // Final logs
                if (enable_logging) {
                    logNetworkState("final");
                    logSummaryReport(iter + 1);
                }

                return true;
            }
        }

        std::cout << "Maximum iterations reached without convergence." << std::endl;

        // Final logs even without convergence
        if (enable_logging) {
            logNetworkState("final_no_convergence");
            logSummaryReport(max_iterations);
        }

        return false;
    }

    // Enable or disable logging
    void setLogging(bool enable, const std::string& log_dir = "", int verbose_level = -1) {
        enable_logging = enable;

        if (!log_dir.empty()) {
            log_directory = log_dir;
        }

        if (verbose_level >= 0) {
            log_level = verbose_level;
        }

        if (enable_logging) {
            createLogDirectory();
        }
    }

    // Rest of your existing methods...
    // Reset assignment to initial state
    void resetAssignment() {
        // Reset link volumes
        for (auto& link_pair : network.links) {
            link_pair.second.volume = 0.0;
            link_pair.second.current_time = link_pair.second.free_flow_time;
        }

        // Reset route flows and costs
        for (auto& route : network.routes) {
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

    // Update route flows using logit model with MSA
    void updateRouteFlows(int iteration) {
        // Store current flows as previous flows
        std::vector<double> previous_flows;
        previous_flows.reserve(network.routes.size());

        for (size_t i = 0; i < network.routes.size(); i++) {
            previous_flows.push_back(network.routes[i].flow);
        }

        // Calculate auxiliary flows using logit model
        std::vector<double> auxiliary_flows(network.routes.size(), 0.0);

        for (auto& od : network.od_pairs) {
            double demand = od.demand;

            // Skip if no demand
            if (demand <= 0) continue;

            // Calculate denominator for logit model
            double exp_sum = 0.0;
            for (int route_idx : od.route_indices) {
                double cost = network.routes[route_idx].cost;
                exp_sum += exp(-logit_theta * cost);
            }

            // Calculate auxiliary flow for each route
            for (int route_idx : od.route_indices) {
                double cost = network.routes[route_idx].cost;

                if (exp_sum > 0) {
                    auxiliary_flows[route_idx] = demand * exp(-logit_theta * cost) / exp_sum;
                }
                else {
                    // Distribute evenly if all costs are too high
                    auxiliary_flows[route_idx] = demand / od.route_indices.size();
                }
            }
        }

        // Apply MSA to update flows
        double step_size = 1.0 / (iteration + 1.0); // Classic MSA step size

        for (size_t i = 0; i < network.routes.size(); i++) {
            // New flow = (1-step_size) * previous flow + step_size * auxiliary flow
            network.routes[i].flow = (1.0 - step_size) * previous_flows[i] + step_size * auxiliary_flows[i];
        }
    }

    // Update route flows using logit model with smart line search
    void updateRouteFlows_LS(int iteration) {
        // Store current flows as previous flows
        std::vector<double> previous_flows;
        previous_flows.reserve(network.routes.size());
        for (size_t i = 0; i < network.routes.size(); i++) {
            previous_flows.push_back(network.routes[i].flow);
        }

        // Calculate auxiliary flows using logit model
        std::vector<double> auxiliary_flows(network.routes.size(), 0.0);
        for (auto& od : network.od_pairs) {
            double demand = od.demand;
            // Skip if no demand
            if (demand <= 0) continue;

            // Calculate denominator for logit model
            double exp_sum = 0.0;
            for (int route_idx : od.route_indices) {
                double cost = network.routes[route_idx].cost;
                exp_sum += exp(-logit_theta * cost);
            }

            // Calculate auxiliary flow for each route
            for (int route_idx : od.route_indices) {
                double cost = network.routes[route_idx].cost;
                if (exp_sum > 0) {
                    auxiliary_flows[route_idx] = demand * exp(-logit_theta * cost) / exp_sum;
                }
                else {
                    // Distribute evenly if all costs are too high
                    auxiliary_flows[route_idx] = demand / od.route_indices.size();
                }
            }
        }

        // Smart line search instead of fixed MSA
        // First, try different step sizes and find the one that minimizes the gap
        std::vector<double> candidate_step_sizes;
        double base_step = 1.0 / (iteration + 1.0); // Classical MSA step as base

        // Generate candidate step sizes
        candidate_step_sizes.push_back(base_step);                  // Classical MSA
        candidate_step_sizes.push_back(std::min(1.0, 2.0 * base_step)); // More aggressive
        candidate_step_sizes.push_back(std::min(1.0, 0.5 * base_step)); // More conservative

        // For later iterations, add larger steps to escape local minima
        if (iteration > 5) {
            candidate_step_sizes.push_back(std::min(1.0, 0.1 + base_step)); // Jump to escape local minimum
        }

        // If we're making good progress, be more aggressive
        if (iteration > 1 && last_gap_improvement > 0.1) {
            candidate_step_sizes.push_back(std::min(1.0, 0.8)); // Very aggressive step
        }

        // Try each step size and keep the best one
        double best_gap = std::numeric_limits<double>::max();
        double best_step_size = base_step; // Default to classical MSA

        for (double step_size : candidate_step_sizes) {
            // Temporarily apply this step size
            for (size_t i = 0; i < network.routes.size(); i++) {
                network.routes[i].flow = (1.0 - step_size) * previous_flows[i] + step_size * auxiliary_flows[i];
            }

            // Update link volumes based on new route flows
            updateLinkVolumes();

            // Update link travel times
            updateLinkTravelTimes();

            // Update route costs
            updateRouteCosts();

            // Calculate gap for this step size
            double current_gap = calculateEquilibriumGap();

            // If this is the best gap so far, remember this step size
            if (current_gap < best_gap) {
                best_gap = current_gap;
                best_step_size = step_size;
            }

            // Restore previous flows for next trial
            for (size_t i = 0; i < network.routes.size(); i++) {
                network.routes[i].flow = previous_flows[i];
            }

            // Restore original network state
            updateLinkVolumes();
            updateLinkTravelTimes();
            updateRouteCosts();
        }

        // Now apply the best step size found
        std::cout << "  Using step size: " << std::fixed << std::setprecision(4) << best_step_size << std::endl;

        for (size_t i = 0; i < network.routes.size(); i++) {
            network.routes[i].flow = (1.0 - best_step_size) * previous_flows[i] + best_step_size * auxiliary_flows[i];
        }

        // Store gap improvement for next iteration
        last_gap_improvement = previous_gap - best_gap;
        previous_gap = best_gap;
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
        flows.reserve(network.routes.size()); // Pre-allocate for efficiency

        for (size_t i = 0; i < network.routes.size(); i++) {
            flows.push_back(network.routes[i].flow);
        }

        return flows;
    }

    // Calculate gap between iterations
    double calculateGap(const std::vector<double>& old_flows) {
        double gap = 0.0;
        double total_flow = 0.0;

        for (size_t i = 0; i < old_flows.size(); i++) {
            gap += std::abs(network.routes[i].flow - old_flows[i]);

            if (gap > 1000)
            {
                int idebug = 1;
            }
            total_flow += network.routes[i].flow;
        }

        return total_flow > 0 ? gap / total_flow : gap;
    }

    // Calculate relative gap for user equilibrium
    double calculateEquilibriumGap() {
        double current_total_cost = 0.0;   // Current total system travel time
        double minimum_total_cost = 0.0;   // Minimum possible total system travel time
        double total_demand = 0.0;         // Total system demand

        // For each OD pair
        for (const auto& od : network.od_pairs) {
            if (od.demand <= 0) continue;

            // Get minimum cost route for this OD pair
            double min_route_cost = std::numeric_limits<double>::max();
            for (int route_idx : od.route_indices) {
                min_route_cost = std::min(min_route_cost, network.routes[route_idx].cost);
            }

            // Calculate costs
            for (int route_idx : od.route_indices) {
                // Current cost: flow � actual route cost
                current_total_cost += network.routes[route_idx].flow * network.routes[route_idx].cost;

                // Minimum cost: flow � minimum route cost (as if all flow used best route)
                minimum_total_cost += network.routes[route_idx].flow * min_route_cost;
            }

            total_demand += od.demand;
        }

        // Prevent division by zero
        if (current_total_cost == 0 || total_demand == 0) {
            return 1.0;  // Return maximum gap
        }

        // Calculate relative gap
        double relative_gap = (current_total_cost - minimum_total_cost) / current_total_cost;

        // Add some debugging output
        std::cout << "  Current system cost: " << current_total_cost
            << ", Minimum possible cost: " << minimum_total_cost
            << ", Gap: " << relative_gap << std::endl;

        return relative_gap;
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

private:
    // Logging methods

    // Create log directory if it doesn't exist
    void createLogDirectory() {
        // This is a simple implementation; you might need to enhance it based on your system
        // For Windows: system(("mkdir " + log_directory).c_str());
        // For Linux/Mac: system(("mkdir -p " + log_directory).c_str());

        // For cross-platform approach, you'd typically use filesystem library (C++17)
        // or a third-party library
    }

    // Get current timestamp as string
    std::string getTimestamp() {
        std::time_t now = std::time(nullptr);
        std::tm localTime;

        // Use localtime_s (safer version)
#ifdef _WIN32
        localtime_s(&localTime, &now);
#else
        localtime_r(&now, &localTime);
#endif

        std::ostringstream oss;
        oss << std::put_time(&localTime, "%Y%m%d_%H%M%S");
        return oss.str();
    }

    // Log network state (initial or final)
    void logNetworkState(const std::string& state) {
        std::string filename = log_directory + "network_state_" + state + "_" + getTimestamp() + ".csv";
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error: Could not open log file " << filename << std::endl;
            return;
        }

        // Log settings
        file << "# NetFlow Traffic Assignment System\n";
        file << "# State: " << state << "\n";
        file << "# Time: " << getTimestamp() << "\n";
        file << "# Settings:\n";
        file << "# - Logit theta: " << logit_theta << "\n";
        file << "# - Max iterations: " << max_iterations << "\n";
        file << "# - Convergence threshold: " << convergence_threshold << "\n";
        file << "# - Number of nodes: " << network.num_nodes << "\n";
        file << "# - Number of links: " << network.num_links << "\n";
        file << "# - Number of zones: " << network.num_zones << "\n";
        file << "# - Number of modes: " << network.num_modes << "\n";
        file << "# - Number of routes: " << network.routes.size() << "\n";
        file << "# - Number of OD pairs: " << network.od_pairs.size() << "\n\n";

        // Node section
        file << "[Nodes]\n";
        file << "node_id,zone_id,x_coord,y_coord\n";
        for (const auto& pair : network.nodes) {
            const Node& node = pair.second;
            file << node.node_id << "," << node.zone_id << "," << node.x_coord << "," << node.y_coord << "\n";
        }
        file << "\n";

        // Link section
        file << "[Links]\n";
        file << "link_id,from_node_id,to_node_id,free_flow_time,lanes,capacity,vdf_alpha,vdf_beta,volume,current_time\n";
        for (const auto& pair : network.links) {
            const Link& link = pair.second;
            file << link.link_id << "," << link.from_node_id << "," << link.to_node_id << ","
                << link.free_flow_time << "," << link.lanes << "," << link.capacity << ","
                << link.vdf_alpha << "," << link.vdf_beta << "," << link.volume << "," << link.current_time << "\n";
        }
        file << "\n";

        // Route section
        file << "[Routes]\n";
        file << "mode,origin_zone,dest_zone,route_id,flow,cost,link_sequence\n";
        for (const auto& route : network.routes) {
            file << route.mode << "," << route.origin_zone << "," << route.dest_zone << ","
                << route.route_id << "," << route.flow << "," << route.cost << ",";

            // Link sequence
            for (size_t i = 0; i < route.link_ids.size(); i++) {
                file << route.link_ids[i];
                if (i < route.link_ids.size() - 1) {
                    file << ";";
                }
            }
            file << "\n";
        }
        file << "\n";

        // OD pairs section
        file << "[OD_Pairs]\n";
        file << "mode,origin,destination,demand\n";
        for (const auto& od : network.od_pairs) {
            file << od.mode << "," << od.origin << "," << od.destination << "," << od.demand << "\n";
        }

        file.close();
    }

    // Log link travel times for current iteration
    void logLinkTravelTimes(int iteration) {
        std::string filename = log_directory + "link_travel_times_iter" + std::to_string(iteration) + ".csv";
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error: Could not open log file " << filename << std::endl;
            return;
        }

        file << "link_id,from_node_id,to_node_id,free_flow_time,current_time,ratio\n";
        for (const auto& pair : network.links) {
            const Link& link = pair.second;
            double ratio = link.current_time / link.free_flow_time;
            file << link.link_id << "," << link.from_node_id << "," << link.to_node_id << ","
                << link.free_flow_time << "," << link.current_time << "," << ratio << "\n";
        }

        file.close();
    }

    // Log link volumes for current iteration
    void logLinkVolumes(int iteration) {
        std::string filename = log_directory + "link_volumes_iter" + std::to_string(iteration) + ".csv";
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error: Could not open log file " << filename << std::endl;
            return;
        }

        file << "link_id,from_node_id,to_node_id,capacity,volume,vc_ratio\n";
        for (const auto& pair : network.links) {
            const Link& link = pair.second;
            double vc_ratio = link.volume / (link.capacity * link.lanes);
            file << link.link_id << "," << link.from_node_id << "," << link.to_node_id << ","
                << (link.capacity * link.lanes) << "," << link.volume << "," << vc_ratio << "\n";
        }

        file.close();
    }

    // Log route costs for current iteration
    void logRouteCosts(int iteration) {
        std::string filename = log_directory + "route_costs_iter" + std::to_string(iteration) + ".csv";
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error: Could not open log file " << filename << std::endl;
            return;
        }

        file << "mode,origin_zone,dest_zone,route_id,flow,cost,flow_to_cost_ratio\n";

        // First, sort routes by flow (volume) for easy identification of major routes
        std::vector<std::pair<int, double>> sorted_routes;
        for (size_t i = 0; i < network.routes.size(); i++) {
            sorted_routes.push_back(std::make_pair(i, network.routes[i].flow));
        }

        //// Sort in descending order of flow
        //std::sort(sorted_routes.begin(), sorted_routes.end(),
        //    [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
        //        return a.second > b.second;
        //    });

        // Counter for routes reported
        int count = 0;

        // Log major routes (flow > 100) and at least the first 5 routes
        for (const auto& pair : sorted_routes) {
            int route_index = pair.first;
            const Route& route = network.routes[route_index];

            // Log if it's a major route (flow > 100) or one of the first 5 routes
            if (route.flow > 100 || count < 10) {
                double flow_to_cost_ratio = (route.cost > 0) ? route.flow / route.cost : 0;

                file << route.mode << ","
                    << route.origin_zone << ","
                    << route.dest_zone << ","
                    << route.route_id << ","
                    << route.flow << ","
                    << route.cost << ","
                    << flow_to_cost_ratio << "\n";

                count++;
            }
        }

        // Add a summary line at the end
        file << "# Summary: Logged " << count << " routes out of " << network.routes.size()
            << " total routes, focusing on major routes (flow > 100) and top 5 routes by flow\n";

        file.close();

        std::cout << "Logged " << count << " major routes for iteration " << iteration << std::endl;
    }

    // Log route flows for current iteration
    void logRouteFlows(int iteration) {
        std::string filename = log_directory + "route_flows_iter" + std::to_string(iteration) + ".csv";
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error: Could not open log file " << filename << std::endl;
            return;
        }

        file << "mode,origin_zone,dest_zone,route_id,flow,cost,percent_of_od_demand\n";

        // First, sort routes by flow for easy identification of major routes
        std::vector<std::pair<int, double>> sorted_routes;
        for (size_t i = 0; i < network.routes.size(); i++) {
            sorted_routes.push_back(std::make_pair(i, network.routes[i].flow));
        }

        // Sort in descending order of flow
        std::sort(sorted_routes.begin(), sorted_routes.end(),
            [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                return a.second > b.second;
            });

        // Counter for routes reported
        int count = 0;

        // Log major routes (flow > 100) and at least the first 5 routes
        for (const auto& pair : sorted_routes) {
            int route_index = pair.first;
            const Route& route = network.routes[route_index];

            // Find corresponding OD pair to calculate percentage
            double od_demand = 0.0;
            for (const auto& od : network.od_pairs) {
                if (od.mode == route.mode && od.origin == route.origin_zone && od.destination == route.dest_zone) {
                    od_demand = od.demand;
                    break;
                }
            }

            double percentage = (od_demand > 0) ? (route.flow / od_demand) * 100.0 : 0.0;

            // Log if it's a major route (flow > 100) or one of the first 5 routes
            if (route.flow > 100 || count < 5) {
                file << route.mode << ","
                    << route.origin_zone << ","
                    << route.dest_zone << ","
                    << route.route_id << ","
                    << route.flow << ","
                    << route.cost << ","
                    << percentage << "\n";

                count++;
            }
        }

        // Add a summary line at the end
        file << "# Summary: Logged " << count << " routes out of " << network.routes.size()
            << " total routes, focusing on major routes (flow > 100) and top 5 routes by flow\n";

        file.close();

        std::cout << "Logged " << count << " major route flows for iteration " << iteration << std::endl;
    }

    // Log OD flows (aggregated route flows) for current iteration
    void logODFlows(int iteration) {
        std::string filename = log_directory + "od_flows_iter" + std::to_string(iteration) + ".csv";
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error: Could not open log file " << filename << std::endl;
            return;
        }

        file << "mode,origin,destination,demand,assigned_flow,percent_assigned,num_routes,min_route_cost,max_route_cost\n";

        // Create a vector of OD pairs with their total flow for sorting
        std::vector<std::tuple<int, double, double>> od_with_flows; // index, demand, flow

        for (size_t i = 0; i < network.od_pairs.size(); i++) {
            const ODPair& od = network.od_pairs[i];
            double total_flow = 0.0;
            for (int route_idx : od.route_indices) {
                total_flow += network.routes[route_idx].flow;
            }
            od_with_flows.push_back(std::make_tuple(i, od.demand, total_flow));
        }

        // Sort by flow in descending order
        std::sort(od_with_flows.begin(), od_with_flows.end(),
            [](const auto& a, const auto& b) {
                return std::get<2>(a) > std::get<2>(b);
            });

        // Counter for OD pairs reported
        int count = 0;

        // Log major OD pairs (flow > 100) and at least the first 5 OD pairs
        for (const auto& tuple : od_with_flows) {
            int od_index = std::get<0>(tuple);
            double demand = std::get<1>(tuple);
            double total_flow = std::get<2>(tuple);

            const ODPair& od = network.od_pairs[od_index];

            // Calculate percentage assigned
            double percent_assigned = (demand > 0) ? (total_flow / demand) * 100.0 : 0.0;

            // Find min and max route costs for this OD pair
            double min_cost = std::numeric_limits<double>::max();
            double max_cost = 0.0;

            for (int route_idx : od.route_indices) {
                double cost = network.routes[route_idx].cost;
                min_cost = std::min(min_cost, cost);
                max_cost = std::max(max_cost, cost);
            }

            if (od.route_indices.empty()) {
                min_cost = 0.0; // No routes
            }

            // Log if it's a major OD pair (flow > 100) or one of the first 5 OD pairs
            if (total_flow > 100 || count < 5) {
                file << od.mode << ","
                    << od.origin << ","
                    << od.destination << ","
                    << demand << ","
                    << total_flow << ","
                    << percent_assigned << ","
                    << od.route_indices.size() << ","
                    << min_cost << ","
                    << max_cost << "\n";

                count++;
            }
        }

        // Add a summary line
        file << "# Summary: Logged " << count << " OD pairs out of " << network.od_pairs.size()
            << " total OD pairs, focusing on major OD pairs (flow > 100) and top 5 OD pairs by flow\n";

        file.close();

        std::cout << "Logged " << count << " major OD flows for iteration " << iteration << std::endl;
    }

    // Log iteration summary
    void logIterationSummary(int iteration, double gap) {
        std::string filename = log_directory + "iteration_summary.csv";
        bool file_exists = std::ifstream(filename).good();

        std::ofstream file(filename, std::ios::app);

        if (!file.is_open()) {
            std::cerr << "Error: Could not open log file " << filename << std::endl;
            return;
        }

        // Write header if file is new
        if (!file_exists) {
            file << "iteration,gap,total_vht,total_vkt,max_vc_ratio,avg_vc_ratio\n";
        }

        // Calculate metrics
        double total_vht = 0.0; // Vehicle Hours Traveled
        double total_vkt = 0.0; // Vehicle Kilometers Traveled
        double max_vc_ratio = 0.0;
        double sum_vc_ratio = 0.0;
        int count_links = 0;

        for (const auto& pair : network.links) {
            const Link& link = pair.second;
            total_vht += link.volume * link.current_time; // Assuming time unit is hours
            // For VKT, you would need link length - add if available
            // total_vkt += link.volume * link.length;

            double vc_ratio = link.volume / (link.capacity * link.lanes);
            max_vc_ratio = std::max(max_vc_ratio, vc_ratio);
            sum_vc_ratio += vc_ratio;
            count_links++;
        }

        double avg_vc_ratio = count_links > 0 ? sum_vc_ratio / count_links : 0.0;

        file << iteration << "," << gap << "," << total_vht << "," << total_vkt << ","
            << max_vc_ratio << "," << avg_vc_ratio << "\n";

        file.close();
    }

    // Log comprehensive summary report at the end
    void logSummaryReport(int total_iterations) {
        std::string filename = log_directory + "summary_report_" + getTimestamp() + ".txt";
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error: Could not open log file " << filename << std::endl;
            return;
        }

        file << "NetFlow Traffic Assignment Summary Report\n";
        file << "=========================================\n\n";

        file << "Run Information:\n";
        file << "  Time: " << getTimestamp() << "\n";
        file << "  Total Iterations: " << total_iterations << "\n";
        file << "  Convergence Threshold: " << convergence_threshold << "\n";
        file << "  Logit Theta: " << logit_theta << "\n\n";

        file << "Network Statistics:\n";
        file << "  Number of Nodes: " << network.num_nodes << "\n";
        file << "  Number of Links: " << network.num_links << "\n";
        file << "  Number of Zones: " << network.num_zones << "\n";
        file << "  Number of Modes: " << network.num_modes << "\n";
        file << "  Number of Routes: " << network.routes.size() << "\n";
        file << "  Number of OD Pairs: " << network.od_pairs.size() << "\n\n";

        // Calculate summary statistics
        double total_demand = 0.0;
        double total_assigned = 0.0;
        for (const auto& od : network.od_pairs) {
            total_demand += od.demand;
            double od_flow = 0.0;
            for (int route_idx : od.route_indices) {
                od_flow += network.routes[route_idx].flow;
            }
            total_assigned += od_flow;
        }

        double total_vht = 0.0;
        double total_vkt = 0.0;
        double max_vc_ratio = 0.0;
        int congested_links = 0;

        for (const auto& pair : network.links) {
            const Link& link = pair.second;
            total_vht += link.volume * link.current_time;
            // total_vkt += link.volume * link.length; // Add if length is available

            double vc_ratio = link.volume / (link.capacity * link.lanes);
            max_vc_ratio = std::max(max_vc_ratio, vc_ratio);

            if (vc_ratio > 0.9) {
                congested_links++;
            }
        }

        file << "Assignment Results:\n";
        file << "  Total Demand: " << total_demand << "\n";
        file << "  Total Assigned Flow: " << total_assigned << "\n";
        file << "  Total VHT: " << total_vht << "\n";
        if (total_vkt > 0) {
            file << "  Total VKT: " << total_vkt << "\n";
        }
        file << "  Maximum V/C Ratio: " << max_vc_ratio << "\n";
        file << "  Number of Congested Links (V/C > 0.9): " << congested_links << "\n\n";

        // Top congested links
        file << "Top 10 Most Congested Links:\n";
        file << "  Link ID, From Node, To Node, Volume, Capacity, V/C Ratio, Travel Time\n";

        std::vector<std::pair<int, double>> link_congestion;
        for (const auto& pair : network.links) {
            const Link& link = pair.second;
            double vc_ratio = link.volume / (link.capacity * link.lanes);
            link_congestion.emplace_back(link.link_id, vc_ratio);
        }

        // Sort by V/C ratio in descending order
        std::sort(link_congestion.begin(), link_congestion.end(),
            [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                return a.second > b.second;
            });

        // Print top 10 or less if fewer links
        int top_count = std::min(10, static_cast<int>(link_congestion.size()));
        for (int i = 0; i < top_count; i++) {
            int link_id = link_congestion[i].first;
            const Link& link = network.links.at(link_id);
            double vc_ratio = link.volume / (link.capacity * link.lanes);

            file << "  " << link.link_id << ", " << link.from_node_id << ", " << link.to_node_id
                << ", " << link.volume << ", " << (link.capacity * link.lanes) << ", "
                << vc_ratio << ", " << link.current_time << "\n";
        }

        file << "\nNote: Detailed logs are available in the log directory.\n";

        file.close();
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
    TrafficAssignmentLog assignment(network, 0.1, 100, 0.001, false, "./logs/", 0);

    // Run assignment
    assignment.assign();

    // Output results
    auto link_volumes = assignment.getLinkVolumes();
    auto route_costs = assignment.getRouteCosts();

    // Calculate system-wide metrics
    double total_vmt = 0.0;  // Vehicle Miles Traveled (or distance units)
    double total_vht = 0.0;  // Vehicle Hours Traveled
    double total_volume = 0.0;

    // Calculate link-based metrics
    for (const auto& link_pair : network.links) {
        const Link& link = link_pair.second;
        // Assuming you have link.length available. If not, you can approximate or omit VMT
        // double link_vmt = link.volume * link.length;
        // total_vmt += link_vmt;

        // Calculate VHT
        double link_vht = link.volume * link.current_time;
        total_vht += link_vht;

        // Accumulate volume
        total_volume += link.volume;
    }

    // Calculate average travel time (hours)
    double avg_travel_time = (total_volume > 0) ? total_vht / total_volume : 0.0;

    // Print system-wide metrics
    std::cout << "\nSystem-wide Metrics:" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    // std::cout << "Total Vehicle Miles Traveled (VMT): " << total_vmt << std::endl;
    std::cout << "Total Vehicle Hours Traveled (VHT): " << total_vht << std::endl;
    std::cout << "Total Volume: " << total_volume << " vehicles" << std::endl;
    std::cout << "Average Travel Time: " << (avg_travel_time * 60) << " minutes" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;

    // Sort routes by volume for better output
    std::vector<std::tuple<int, int, int, int, double, double>> sorted_routes;
    for (const auto& tuple : route_costs) {
        int mode = std::get<0>(tuple);
        int origin = std::get<1>(tuple);
        int dest = std::get<2>(tuple);
        int route_id = std::get<3>(tuple);
        double cost = std::get<4>(tuple);

        // Find the corresponding route to get its volume
        double volume = 0.0;
        for (const auto& route : network.routes) {
            if (route.mode == mode && route.origin_zone == origin &&
                route.dest_zone == dest && route.route_id == route_id) {
                volume = route.flow;
                break;
            }
        }

        sorted_routes.push_back(std::make_tuple(mode, origin, dest, route_id, cost, volume));
    }

    // Sort by volume in descending order
    std::sort(sorted_routes.begin(), sorted_routes.end(),
        [](const auto& a, const auto& b) {
            return std::get<5>(a) > std::get<5>(b);
        });

    // Print top routes by volume
    std::cout << "\nTop Routes by Volume:" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "Mode | Origin | Dest | Route ID | Volume | Cost (min)" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;

    // Print at most 10 routes
    int count = 0;
    for (const auto& tuple : sorted_routes) {
        if (count >= 10) break;

        int mode = std::get<0>(tuple);
        int origin = std::get<1>(tuple);
        int dest = std::get<2>(tuple);
        int route_id = std::get<3>(tuple);
        double cost = std::get<4>(tuple);
        double volume = std::get<5>(tuple);

        // Only print routes with volume > 0
        if (volume > 0) {
            std::cout << std::setw(4) << mode << " | "
                << std::setw(6) << origin << " | "
                << std::setw(4) << dest << " | "
                << std::setw(8) << route_id << " | "
                << std::setw(6) << std::fixed << std::setprecision(1) << volume << " | "
                << std::setw(9) << std::fixed << std::setprecision(2) << cost << std::endl;
            count++;
        }
    }

    if (count == 0) {
        std::cout << "No routes with positive volume found." << std::endl;
    }

    std::cout << "------------------------------------------------------" << std::endl;

    // Export results
    //  exportNetworkToJSON(network, "network_results.json");

    return 0;
}
