#include <memory>
#include "../stdafx.h"
#include "LoadingFunctions.h"

/**
 * @file LoadingFunctions.cpp
 * @brief Implementation of all loading Functions used during the Project.
 */

std::vector<Vertex> Vertexes;
std::vector<Edge> Edges;


constexpr double kEarthRadiusKm = 6371.0;

double ToRadians(double degrees) {
    return degrees * M_PI / 180.0;
}

/**
 * @brief Removes termination characters (e.g., '\r', '\n', '\0') from the given string.
 * @param str The string from which to remove termination characters.
 * @complexity O(n), where n is the length of the string.
 */
void Remove_terminations(std::string& str)
{
    size_t pos = str.find('\r');
    if (pos != std::string::npos) {
        str.erase(pos);
    }
    pos = str.find('\n');
    if (pos != std::string::npos) {
        str.erase(pos);
    }
    pos = str.find('\0');
    if (pos != std::string::npos) {
        str.erase(pos);
    }
}

std::string extractLetters(const std::string filename) {
    size_t start = filename.find('_');
    size_t end = filename.find('.');
    if (start != std::string::npos && end != std::string::npos && start < end) {
        return filename.substr(start + 1, end - start - 1);
    } else {
        return "";
    }
}


void LoadToyGraphs(Graph * g , const std::string& path , const int& graph){

    std::string file_name;

    int edgeId = 0;

    switch (graph) {
        case 0:
            file_name = "shipping.csv";
            break;
        case 1:
            file_name = "stadiums.csv";
            break;
        case 2:
            file_name = "tourism.csv";
            break;
        default:
            std::cerr << "Choose a valid graph";
            return;
    }

    std::string full_path = path + '/' + file_name;

    std::ifstream file(full_path);
    if (!file.is_open()) {
        std::cerr << "Failed to open the CSV file." << std::endl;
    }

    std::string line;
    getline(file, line);

    while (getline(file, line)) {

        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        std::istringstream lineStream(line);

        std::string id_v1 , id_v2 , distance;

        getline(lineStream , id_v1 , ',');
        getline(lineStream , id_v2 , ',');
        getline(lineStream , distance , ',');

        int id_V1 = stoi(id_v1);
        int id_V2 = stoi(id_v2);
        double v_distances = stod(distance);

        g->addVertex(id_V1 , 0 , 0);
        g->addVertex(id_V2 , 0 , 0);
        g->addEdge(id_V1 , id_V2 ,edgeId, v_distances);

        edgeId++;
    }

    file.close();
}


void LoadRealWorldGraphs(Graph * g , const std::string& path , const int& graph){
    std::string file_name;
    int edgeId = 0;

    switch (graph) {
        case 0:
            file_name = "/graph1/nodes.csv";
            break;
        case 1:
            file_name = "/graph2/nodes.csv";
            break;
        case 2:
            file_name = "/graph3/nodes.csv";
            break;
        default:
            std::cerr << "Choose a valid graph";
            return;
    }

    std::string full_path = path+ '/' + file_name;

    std::ifstream file(full_path);
    if (!file.is_open()) {
        std::cerr << "Failed to open the CSV file." << std::endl;
    }

    std::string line;
    getline(file, line);

    while (getline(file, line)) {

        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        std::istringstream lineStream(line);

        std::string id , intitude , latitude;

        getline(lineStream , id , ',');
        getline(lineStream , intitude , ',');
        getline(lineStream , latitude , ',');

        int id_V1 = stoi(id);
        double intitude_ = stod(intitude);
        double latitude_ = stod(latitude);

        g->addVertex(id_V1 , intitude_ , latitude_);
    }

    file.close();

    switch (graph) {
        case 0:
            file_name = "/graph1/edges.csv";
            break;
        case 1:
            file_name = "/graph2/edges.csv";
            break;
        case 2:
            file_name = "/graph3/edges.csv";
            break;
        default:
            std::cerr << "Choose a valid graph";
            return;
    }

    full_path = path+ '/' + file_name;

    std::ifstream file2(full_path);
    if (!file2.is_open()) {
        std::cerr << "Failed to open the CSV file." << std::endl;
    }

    std::string line2;
    getline(file2, line2);

    while (getline(file2, line2)) {

        line2.erase(std::remove(line2.begin(), line2.end(), '\r'), line2.end());

        std::istringstream lineStream(line2);

        std::string id_v1 , id_v2 , distance;

        std::getline(lineStream , id_v1 , ',');
        std::getline(lineStream , id_v2 , ',');
        std::getline(lineStream , distance , ',');

        int id_V1 = stoi(id_v1);
        int id_V2 = stoi(id_v2);
        double v_distances = stod(distance);

        Vertex v1 = Vertex(id_V1 , 0 , 0);
        Vertex v2 = Vertex(id_V2 , 0 , 0);

        Edge edge = Edge(&v1 , &v2 ,edgeId, v_distances);

        g->addBidirectionalEdge(v1.getId() , v2.getId() , edge.getId() , edge.getWeight());

        edgeId += 2;
    }

    file2.close();

}

void LoadMediumGraphs(Graph * g , const std::string& path , const int& graph){
    std::string file_name2 = "/nodes.csv";

    std::string full_path = path + file_name2;

    std::string  file_name;

    int edgeId = 0;

    switch (graph) {
        case 0:
            file_name = "edges_25.csv";
            break;
        case 1:
            file_name = "edges_50.csv";
            break;
        case 2:
            file_name = "edges_75.csv";
            break;
        case 3:
            file_name = "edges_100.csv";
            break;
        case 4:
            file_name = "edges_200.csv";
            break;
        case 5:
            file_name = "edges_300.csv";
            break;
        case 6:
            file_name = "edges_400.csv";
            break;
        case 7:
            file_name = "edges_500.csv";
            break;
        case 8:
            file_name = "edges_600.csv";
            break;
        case 9:
            file_name = "edges_700.csv";
            break;
        case 10:
            file_name = "edges_800.csv";
            break;
        case 11:
            file_name = "edges_900.csv";
            break;
        default:
            std::cerr << "Choose a valid graph";
            return;
    }

    std::string numbers = extractLetters(file_name);
    int i = std::stoi(numbers);

    std::ifstream file(full_path);
    if (!file.is_open()) {
        std::cerr << "Failed to open the CSV file." << std::endl;
    }

    std::string line;
    getline(file, line);

    while (getline(file, line) && i > 0) {

        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        std::istringstream lineStream(line);
        std::vector<std::string> tokens;
        std::string token;

        while (getline(lineStream, token, ',')) {
            tokens.push_back(token);
        }

        std::string id = tokens[0];
        std::string longitude_ = tokens[1];
        std::string latitude = tokens[2];

        Remove_terminations(latitude);

        int id_V1 = stoi(id);
        double longitude = stod(longitude_);
        double latitude_ = stod(latitude);

        g->addVertex(id_V1 ,longitude, latitude_ );
        i--;
    }

    file.close();


    full_path = path+ '/' + file_name;

    std::ifstream file2(full_path);
    if (!file2.is_open()) {
        std::cerr << "Failed to open the CSV file." << std::endl;
    }

    std::string line2;

    while (getline(file2, line2)) {

        line2.erase(std::remove(line2.begin(), line2.end(), '\r'), line2.end());

        std::istringstream lineStream(line2);

        std::string id_v1 , id_v2 , distance;

        std::getline(lineStream , id_v1 , ',');
        std::getline(lineStream , id_v2 , ',');
        std::getline(lineStream , distance , ',');

        int id_V1 = stoi(id_v1);
        int id_V2 = stoi(id_v2);
        double v_distances = stod(distance);

        Vertex v1 = Vertex(id_V1 , 0 , 0);
        Vertex v2 = Vertex(id_V2 , 0 , 0);

        Edge edge = Edge(&v1 , &v2 , edgeId ,  v_distances);

        g->addBidirectionalEdge(v1.getId() ,v2.getId() ,edge.getId(), edge.getWeight());

        edgeId += 2;
    }

    file2.close();
}