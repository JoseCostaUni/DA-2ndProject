#include "UI.h"
#include "../Logic/Logic.h"
#include "../Logic/Algorithms.h"
#include <thread>
#include <chrono>
#include <iomanip>
#include "map"
#include <functional>

UI::UI() {}

/**
 * @brief Clears the console screen by outputting multiple newline characters.
 *
 * @note Time complexity: O(1)
 */


void UI::clear_screen() {
    int i = 0;
    while(i != 100) {
        std::cout << std::endl;
        i++;
    }
}

/**
 * @brief Validates user input to ensure it is within a specified range.
 *
 * @param op Reference to the character variable where the input will be stored.
 * @param lower_bound Lower bound of the valid input range.
 * @param upper_bound Upper bound of the valid input range.
 * @return true if input is valid, false otherwise.
 *
 * @note Time complexity: O(1)
 */
bool UI::validate_input(char &op, const char lower_bound, const char upper_bound) {
    std::string tempValue;
    while(true){
        std::cin >> tempValue;
        std::cout << "\n";
        op = tempValue[0];
        op = std::toupper(op);
        if (std::cin.fail() || tempValue.length() != 1) {
            std::cout << "Introduce a valid option (" << lower_bound << "-" << upper_bound << "): ";
        }else{
            break;
        }
    }
    while (op < lower_bound || op > upper_bound) {
        std::cout << "Introduce a valid option (" << lower_bound << "-" << upper_bound << "): ";
        std::cin.clear();
        std::cin >> op;
    }
    return true;
}


/**
 * @brief Validates user input to ensure it contains only digits.
 *
 * @param index Reference to the integer variable where the input will be stored.
 * @return true if input is valid, false otherwise.
 *
 * @note Time complexity: O(n), where n is the number of characters in the input.
 */
bool UI::validate_int_input(int &index) {
    std::string tempValue;
    char op;

    while (true) {
        std::cin >> tempValue;
        std::cout << "\n";

        // Check if the input contains only digits
        bool isNumber = true;
        for (char c : tempValue) {
            if (!std::isdigit(c)) {
                isNumber = false;
                break;
            }
        }

        if (!isNumber) {
            std::cout << "Input should contain only numbers. Please try again: ";
        } else {
            // Convert the string to an integer
            index = std::stoi(tempValue);
            break;
        }
    }

    return true;
}
/**
 * @brief Loads the required dataset for analysis and starts the program.
 *
 * This function prompts the user to choose a dataset for analysis and loads the corresponding data.
 * It then waits for the user to press 'A' to start the program.
 *
 * @param ui The UI object to perform operations on.
 *
 * @note Time complexity: O(1)
 */
void UI::loading_stuff(UI &ui) {
    char op , secondOp;
    std::string path;
    std::cout << "Which data set do you wish to use during the analysis?" << std::endl
      << "A. Toy Data Sets" << std::endl
      << "B. Medium Data Sets" << std::endl
      << "C. Real World Data Sets" << std::endl
      << "Insert the letter: " ;
    validate_input(op,'A','C');

    switch (op) {
        case 'A':
            path = "../Graphs/Toy-Graphs/Toy-Graphs";
            std::cout << "Which Data Set do you wish to use?" << std::endl
                << "A. shipping.csv" << std::endl
                << "B. stadiums.csv" << std::endl
                << "C. tourism.csv" << std::endl
                << "D. myGraph.csv" << std::endl
                << "Insert the letter: " ;
            validate_input(secondOp,'A','D');
            switch (secondOp) {
                case 'A':
                    LoadToyGraphs(&g , path , 0);
                    file_path = path + "/shipping.csv";
                    break;
                case 'B':
                    LoadToyGraphs(&g , path , 1);
                    file_path = path + "/stadiums.csv";
                    break;
                case 'C':
                    LoadToyGraphs(&g , path , 2);
                    file_path = path + "/tourism.csv";
                    break;
                case 'D':
                    LoadToyGraphs(&g , path , 3);
                    file_path = path + "/myGraph.csv";
                    break;
                default:
                    std::cerr << "Invalid Input";
                    return;
            }
            break;
        case 'B':
            path = "../Graphs/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs";
            std::cout << "Which Data Set do you wish to use?" << std::endl
                      << "A. edges25.csv" << std::endl
                      << "B. edges50.csv" << std::endl
                      << "C. edges75.csv" << std::endl
                      << "D. edges100.csv" << std::endl
                      << "E. edges200.csv" << std::endl
                      << "F. edges300.csv" << std::endl
                      << "G. edges400.csv" << std::endl
                      << "H. edges500.csv" << std::endl
                      << "I. edges600.csv" << std::endl
                      << "J. edges700.csv" << std::endl
                      << "K. edges800.csv" << std::endl
                      << "L. edges900.csv" << std::endl
                      << "Insert the letter: " ;
            validate_input(secondOp,'A','L');
            switch (secondOp) {
                case 'A':
                    LoadMediumGraphs(&g , path , 0);
                    break;
                case 'B':
                    LoadMediumGraphs(&g , path , 1);
                    break;
                case 'C':
                    LoadMediumGraphs(&g , path , 2);
                    break;
                case 'D':
                    LoadMediumGraphs(&g , path , 3);
                    break;
                case 'E':
                    LoadMediumGraphs(&g , path , 4);
                    break;
                case 'F':
                    LoadMediumGraphs(&g , path , 5);
                    break;
                case 'G':
                    LoadMediumGraphs(&g , path , 6);
                    break;
                case 'H':
                    LoadMediumGraphs(&g , path , 7);
                    break;
                case 'I':
                    LoadMediumGraphs(&g , path , 8);
                    break;
                case 'J':
                    LoadMediumGraphs(&g , path , 9);
                    break;
                case 'K':
                    LoadMediumGraphs(&g , path , 10);
                    break;
                case 'L':
                    LoadMediumGraphs(&g , path , 11);
                    break;
                default:
                    std::cerr << "Invalid Input";
                    return;
            }

            break;
        case 'C':
            path = "../Graphs/Real-world Graphs/Real-world Graphs";
            std::cout << "Which Data Set do you wish to use?" << std::endl
                      << "A. Graph1" << std::endl
                      << "B. Graph2" << std::endl
                      << "C. Graph3" << std::endl
                      << "Insert the letter: " ;
            validate_input(secondOp,'A','C');
            switch (secondOp) {
                case 'A':
                    LoadRealWorldGraphs(&g , path , 0);
                    file_path = path + "/Graph1";
                    break;
                case 'B':
                    LoadRealWorldGraphs(&g , path , 1);
                    file_path = path + "/Graph2";
                    break;
                case 'C':
                    LoadRealWorldGraphs(&g , path , 2);
                    file_path = path + "/Graph3";
                    break;
                default:
                    std::cerr << "Invalid Input";
                    return;
            }
            break;
        default:
            std::cerr << "Invalid input";
            return;
    }


    std::cout << "Load Finished" << std::endl;
    std::cout << "Press A to start the program: ";
    op;
    validate_input(op,'A','A');
    menu_start();
}

/**
 * @brief Allows the user to change the dataset being used.
 *
 * This function clears the current graph data and prompts the user to choose a new dataset
 * for analysis. It then loads the corresponding data and starts the program again.
 *
 * @note Time complexity: O(1)
 */
void UI::changeDataSet(){

    this->g.clear();

    char op , secondOp;
    std::string path;
    std::cout << "Which data set do you wish to switch to?" << std::endl
              << "A. Toy Data Sets" << std::endl
              << "B. Medium Data Sets" << std::endl
              << "C. Real World Data Sets" << std::endl
              << "Insert the letter: " ;
    validate_input(op,'A','C');

    switch (op) {
        case 'A':
            path = "../Graphs/Toy-Graphs/Toy-Graphs";
            std::cout << "Which Data Set do you wish to use?" << std::endl
                      << "A. shipping.csv" << std::endl
                      << "B. stadiums.csv" << std::endl
                      << "C. tourism.csv" << std::endl
                      << "D. myGraph.csv" << std::endl
                      << "Insert the letter: " ;
            validate_input(secondOp,'A','D');
            switch (secondOp) {
                case 'A':
                    LoadToyGraphs(&g , path , 0);
                    file_path = path + "/shipping.csv";
                    break;
                case 'B':
                    LoadToyGraphs(&g , path , 1);
                    file_path = path + "/stadiums.csv";
                    break;
                case 'C':
                    LoadToyGraphs(&g , path , 2);
                    file_path = path + "/tourism.csv";
                    break;
                case 'D':
                    LoadToyGraphs(&g , path , 3);
                    file_path = path + "/myGraph.csv";
                    break;
                default:
                    std::cerr << "Invalid Input";
                    return;
            }
            break;
        case 'B':
            path = "../Graphs/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs";
            std::cout << "Which Data Set do you wish to use?" << std::endl
                      << "A. edges25.csv" << std::endl
                      << "B. edges50.csv" << std::endl
                      << "C. edges75.csv" << std::endl
                      << "D. edges100.csv" << std::endl
                      << "E. edges200.csv" << std::endl
                      << "F. edges300.csv" << std::endl
                      << "G. edges400.csv" << std::endl
                      << "H. edges500.csv" << std::endl
                      << "I. edges600.csv" << std::endl
                      << "J. edges700.csv" << std::endl
                      << "K. edges800.csv" << std::endl
                      << "L. edges900.csv" << std::endl
                      << "Insert the letter: " ;
            validate_input(secondOp,'A','L');
            switch (secondOp) {
                case 'A':
                    LoadMediumGraphs(&g , path , 0);
                    file_path = path + "/edges25.csv";
                    break;
                case 'B':
                    LoadMediumGraphs(&g , path , 1);
                    file_path = path + "/edges50.csv";
                    break;
                case 'C':
                    LoadMediumGraphs(&g , path , 2);
                    file_path = path + "/edges75.csv";
                    break;
                case 'D':
                    LoadMediumGraphs(&g , path , 3);
                    file_path = path + "/edges100.csv";
                    break;
                case 'E':
                    LoadMediumGraphs(&g , path , 4);
                    file_path = path + "/edges200.csv";
                    break;
                case 'F':
                    LoadMediumGraphs(&g , path , 5);
                    file_path = path + "/edges300.csv";
                    break;
                case 'G':
                    LoadMediumGraphs(&g , path , 6);
                    file_path = path + "/edges400.csv";
                    break;
                case 'H':
                    LoadMediumGraphs(&g , path , 7);
                    file_path = path + "/edges500.csv";
                    break;
                case 'I':
                    LoadMediumGraphs(&g , path , 8);
                    file_path = path + "/edges600.csv";
                    break;
                case 'J':
                    LoadMediumGraphs(&g , path , 9);
                    file_path = path + "/edges700.csv";
                    break;
                case 'K':
                    LoadMediumGraphs(&g , path , 10);
                    file_path = path + "/edges800.csv";
                    break;
                case 'L':
                    LoadMediumGraphs(&g , path , 11);
                    file_path = path + "/edges900.csv";
                    break;
                default:
                    std::cerr << "Invalid Input";
                    return;
            }

            break;
        case 'C':
            path = "../Graphs/Real-world Graphs/Real-world Graphs";
            std::cout << "Which Data Set do you wish to use?" << std::endl
                      << "A. Graph1" << std::endl
                      << "B. Graph2" << std::endl
                      << "C. Graph3" << std::endl
                      << "Insert the letter: " ;
            validate_input(secondOp,'A','C');
            switch (secondOp) {
                case 'A':
                    LoadRealWorldGraphs(&g , path , 0);
                    file_path = path + "/Graph1";
                    break;
                case 'B':
                    LoadRealWorldGraphs(&g , path , 1);
                    file_path = path + "/Graph2";
                    break;
                case 'C':
                    LoadRealWorldGraphs(&g , path , 2);
                    file_path = path + "/Graph3";
                    break;
                default:
                    std::cerr << "Invalid Input";
                    return;
            }
            break;
        default:
            std::cerr << "Invalid input";
            return;
    }


    std::cout << "Load Finished" << std::endl;
    std::cout << "Press A to start the program: ";
    op;
    validate_input(op,'A','A');
    main_menu();
}

/**
 * @brief Displays the startup menu of the program.
 *
 * This function displays the startup menu of the program, prompting the user to either proceed
 * to the application or close it.
 *
 * @note Time complexity: O(1)
 */
void UI::menu_start() {
    char op;
    std::cout << "##############################################################################" << std::endl
              << "#  @@@@@   @@@@@   @@@@@       @@@@@   @@@@@   @       @   @  @@@@@   @@@@   #" << std::endl
              << "#    @     @       @   @       @       @   @   @       @   @  @       @   @  #" << std::endl
              << "#    @     @@@@@   @@@@@       @@@@@   @   @   @       @   @  @@@@@   @@@@   #" << std::endl
              << "#    @         @   @               @   @   @   @        @ @   @       @  @   #" << std::endl
              << "#    @     @@@@@   @           @@@@@   @@@@@   @@@@@     @    @@@@@   @   @  #" << std::endl
              << "##############################################################################" << std::endl << '\n'
              << "Welcome to the Routing Algorithms App for Ocean Shipping and Urban Deliveries (TSP), what would you like to do?" << std::endl
              << "A. Proceed to the application" << std::endl
              << "C. Close the application" << std::endl
              << "Insert the letter: " ;
    validate_input(op,'A','B');
    switch(op){
        case 'A':
            main_menu();
            break;
        case 'B':
            std::cout << "Thanks for using our Routing Algorithms App!" << std::endl << "\n"
                 << "Made by: " << std::endl
                 << "Ângelo Oliveira || 202207798" << std::endl
                 << "José Costa      || 202207871" << std::endl
                 << "Bruno Fortes    || 202209730" << std::endl;
            std::this_thread::sleep_for(std::chrono::seconds(2)); // Espera 2 segundos antes de fechar o terminal
            exit(0);

    }
}

/**
 * @brief Displays options related to graph information and handles user input.
 *
 * This function displays options related to graph information, such as printing node information,
 * graph basic information, calculating Harverstein distance, or returning to the main menu.
 *
 * @note Time complexity: O(1)
 */
void UI::GraphOptionsMenu() {
    char op;
    std::cout << "What would you like to know?" << std::endl;
    std::cout << "A. Print all Nodes information (Not Recommended for large graphs)" << std::endl
              << "B. Print Graph basic information (Number of nodes and edges)" << std::endl
              << "C. Calculate Harverstein Distance between two nodes" << std::endl
              << "D. Go back to the main menu" << std::endl
              << "Insert your choice:";

    validate_input(op, 'A', 'D');
    int index = 0, index2 = 0;
    Vertex *testV1;
    Vertex *testV2;
    switch (op) {
        case 'A':
            g.printNodesContente();
            back_menu_GraphOptions();
            break;
        case 'B':
            g.printGraphInfo();
            back_menu_GraphOptions();
            break;
        case 'C':

            std::cout << "Please introduce the index of the source vertex [0 - " << g.getVertexSet().size() - 1 << "]: ";

            while (true){
                validate_int_input(index);
                if(index < 0 || index >= g.getVertexSet().size()){
                    std::cout << "Please introduce a valid index [0 - " << g.getVertexSet().size() - 1 << "]: ";
                }else{
                    break;
                }
            }

            std::cout << "Please introduce the index of the destination vertex [0 - " << g.getVertexSet().size() - 1 << "]: ";

            while (true){
                validate_int_input(index2);
                if(index < 0 || index >= g.getVertexSet().size()){
                    std::cout << "Please introduce a valid index [0 - " << g.getVertexSet().size() - 1 << "]: ";
                }else{
                    break;
                }
            }

            testV1 = g.findVertex(index);
            testV2 = g.findVertex(index2);

            std::cout << "Distance = "
                      << Harverstein(testV1->getLongitude(), testV1->getLatitude(), testV2->getLongitude(),
                                     testV2->getLatitude()) << std::endl;

            back_menu_GraphOptions();
            break;
        case 'D':
            main_menu();
            break;
        default:
            std::cerr << "Error";
    }
}

/**
 * @brief Displays the main menu options and handles user input.
 *
 * This function displays the main menu options and waits for user input to select an action.
 * Depending on the user's choice, it performs various operations such as changing the dataset,
 * accessing graph information options, calculating the Traveling Salesman Problem (TSP) using different algorithms,
 * or exiting the program.
 *
 * Time Complexity: O(N + E), where N is the number of vertices and E is the number of edges in the graph.
 */
void UI::main_menu(){
    clear_screen();
    char op;
    std::cout << "What would you like to know?" <<std::endl;
    std::cout << "A. Switch Data Set" <<std::endl
              << "B. Graph Information Options" << std::endl
              << "C. Calculate TSP with Backtracking Brute Force" << std::endl
              << "D. Calculate TSP with Triangular Approximation (Beware that this will make the graph fully connected which may take some time)" << std::endl
              << "E. Calculate TSP with Nearest Neighbour algorithm (Beware that this will make the graph fully connected which may take some time)" << std::endl
              << "F. Calculate TSP for incomplete Graphs" << std::endl
             << "G Exit the program" << std::endl
             << "Insert your choice:";

    validate_input(op, 'A', 'G');

    double weight = 0 , lastWeight = DBL_MAX;
    std::vector<Vertex *> articulationPoints;
    std::vector<std::vector<Vertex *>> scc;
    Vertex * source;
    Vertex * dest;
    bool isFullyConnected = false;
    switch(op){
        case 'A':
            changeDataSet();
            main_menu();
            break;
        case 'B':
            GraphOptionsMenu();
            break;
        case 'C':

            if(fullyConnected){
                pathSelector();
                fullyConnected = false;
            }

            tspBruteForce(&g);
            back_menu();
            break;
        case 'D':
            source = g.findVertex(0);

            if(g.makeFullyConnected()){
                isFullyConnected = true;
            }else {
                std::cout << "Not possible to make fully connected\n";
                back_menu();
            }

            TSP = TriangularApproximationHeuristic(&g , source , nullptr);
            std::cout << source->getId() ;
            std::cout << " -------> ";
            for(Edge * e : TSP){
                weight += e->getWeight();
                std::cout << e->getDestination()->getId() ;
                if(e->getDestination()->getId() != 0){
                    std::cout << "-------> ";
                }
            }
            std::cout << std::endl;
            std::cout << std::fixed;
            std::cout << std::setprecision(2);
            std::cout << "Weight :" << weight << std::endl;

            back_menu();
            break;
        case 'E':

            if(!fullyConnected){
                if(g.makeFullyConnected()){
                    fullyConnected = true;
                }else {
                    std::cout << "Not possible to make fully connected\n";
                    back_menu();
                }
            }

            source = g.findVertex(0);

            TSP = NearestNeighbour(&g , source);

            std::cout << source->getId() ;
            std::cout << " -------> ";
            for(Edge * e : TSP){
                weight += e->getWeight();
                std::cout << e->getDestination()->getId() ;
                if(e->getDestination()->getId() != 0){
                    std::cout << "-------> ";
                }
            }
            std::cout << std::endl;
            std::cout << std::fixed;
            std::cout << std::setprecision(2);
            std::cout << "Weight :" << weight << std::endl;

            back_menu();
            break;
        case 'F':

            BackTrackMenu();
            break;
        case 'G':
            std::cout << "Thanks for using our TSP solver tool!" <<std::endl << "\n"
                      << "Made by: " <<std::endl
                      << "Ângelo Oliveira || 202207798" <<std::endl
                      << "José Costa      || 202207871" <<std::endl
                      << "Bruno Fortes    || 202209730" << std::endl;
            std::this_thread::sleep_for(std::chrono::seconds(2)); // Espera 2 segundos antes de fechar o terminal
            exit(0);
        default:
            std::cerr << "Error";
    }
}

/**
 * @brief Retrieves the current graph stored in the UI.
 *
 * This function returns the current graph stored in the UI.
 *
 * @return The current graph stored in the UI.
 *
 * Time Complexity: O(1).
 */
Graph UI::getGraph() const {
    return g;
}

/**
 * @brief Displays an option to go back to the main menu.
 *
 * This function prompts the user to press 'A' to go back to the main menu.
 * Upon user input, it redirects to the main menu.
 *
 * Time Complexity: O(1).
 */
void UI::back_menu(){
    char op;
    std::cout << "Press A to go back to the menu: ";
    validate_input(op,'A','A');
    main_menu();
}

/**
 * @brief Displays an option to go back to the Graph Options menu.
 *
 * This function prompts the user to press 'A' to go back to the Graph Options menu.
 * Upon user input, it redirects to the Graph Options menu.
 *
 * Time Complexity: O(1).
 */
void UI::back_menu_GraphOptions() {
    char op;
    std::cout << "Press A to go back to the menu: ";
    validate_input(op,'A','A');
    GraphOptionsMenu();
}

/**
 * @brief Displays the BackTrack menu and prompts the user to choose an option.
 *
 * This function presents the user with options to perform TSP calculations using different algorithms.
 * The user can select either:
 * - A. Calculate TSP with Nearest Neighbour with Backtrack algorithm for incomplete Graphs without Two Opt Optimization.
 * - B. Calculate TSP with Nearest Neighbour with Backtrack algorithm for incomplete Graphs with Two Opt Optimization (Beware that this will take 10x more than the previous).
 * - C. Go back to the main menu.
 *
 * Upon user selection, the function performs the corresponding action.
 *
 * Time Complexity: O(n), where n is the size of the vertex set in the graph.
 */
void UI::BackTrackMenu() {
    char op;
    std::cout << "Which one would you like to use?" << std::endl;
    std::cout << "A. Calculate TSP with Nearest Neighbour with Backtrack algorithm for incomplete Graphs without Two Opt Optimization " << std::endl
              << "B. Calculate TSP with Nearest Neighbour with Backtrack algorithm for incomplete Graphs with Two Opt Optimization(Beware that this will take 10x more than the previous)" << std::endl
              << "C. Go back to the main menu" << std::endl
              << "Insert your choice:";

    validate_input(op, 'A', 'C');
    int index = 0, index2 = 0;
    Vertex *source;
    std::vector<Vertex * > path;
    switch (op) {
        case 'A':

            if(fullyConnected){
                pathSelector();
                fullyConnected = false;
            }

            std::cout << "Please introduce the index of the source vertex [0 - " << g.getVertexSet().size() - 1 << "]: ";

            while (true){
                validate_int_input(index);
                if(index < 0 || index >= g.getVertexSet().size()){
                    std::cout << "Please introduce a valid index [0 - " << g.getVertexSet().size() - 1 << "]: ";
                }else{
                    break;
                }
            }

            source = g.findVertex(index);

            path.clear();
            nn_with_backtracking(&g , source , path);

            back_menu_BacktrackingOptions();
            break;
        case 'B':

            if(fullyConnected){
                pathSelector();
                fullyConnected = false;
            }

            std::cout << "Please introduce the index of the source vertex [0 - " << g.getVertexSet().size() - 1 << "]: ";

            while (true){
                validate_int_input(index);
                if(index < 0 || index >= g.getVertexSet().size()){
                    std::cout << "Please introduce a valid index [0 - " << g.getVertexSet().size() - 1 << "]: ";
                }else{
                    break;
                }
            }

            source = g.findVertex(index);

            path.clear();
            nn_with_backtrackingAndTwoOpt(&g , source , path);

            back_menu_BacktrackingOptions();
            break;
        case 'C':
            main_menu();
            break;
        default:
            std::cerr << "Error";
    }
}

/**
 * @brief Navigates back to the BackTrack menu upon user input.
 *
 * This function prompts the user to press 'A' to return to the BackTrack menu.
 * It validates the input and calls the BackTrackMenu function accordingly.
 *
 * Time Complexity: O(1)
 */
void UI::back_menu_BacktrackingOptions() {
    char op;
    std::cout << "Press A to go back to the BackTrack menu: ";
    validate_input(op,'A','A');
    BackTrackMenu();
}

/**
 * @brief Selects the appropriate action based on the file path loaded.
 *
 * If the file path is empty, it prints an error message and returns.
 * Otherwise, it extracts the last word from the file path and determines the action
 * based on this word.
 *
 * If the last word corresponds to a specific command (e.g., "Graph1", "shipping"),
 * it clears the graph, loads the corresponding data set, and performs the associated action.
 *
 * @note This function assumes that the file path has been previously loaded.
 *
 * @note Time complexity: O(n), where n is the length of the file path.
 */
void UI::pathSelector() {

    if(file_path.empty()){
        std::cerr << "No file path was loaded" << std::endl;
        return;
    }

    std::size_t lastSlashPos = file_path.find_last_of('/');
    std::string command;
    if (lastSlashPos != std::string::npos) {
        command = file_path.substr(lastSlashPos + 1);
        std::cout << "Last word: " << command << std::endl;
    } else {
        command = "error";
        std::cout << "No '/' found in the path" << std::endl;
    }

    std::string path;
    std::map<std::string, std::function<void()>> commandMap;

    if(command.empty()){
        std::cerr << "Invalid file path" << std::endl;
        return;
    }else if(command == "Graph1"){
        g.clear();
        path = "../Graphs/Real-world Graphs/Real-world Graphs";
        LoadRealWorldGraphs(&g , path , 0);
    }else if(command == "Graph2") {
        g.clear();
        path = "../Graphs/Real-world Graphs/Real-world Graphs";
        LoadRealWorldGraphs(&g, path, 1);
    }else if(command == "Graph3") {
        g.clear();
        path = "../Graphs/Real-world Graphs/Real-world Graphs";
        LoadRealWorldGraphs(&g, path, 2);
    }else if(command == "shipping") {
        g.clear();
        path = "../Graphs/Toy-Graphs/Toy-Graphs";
        LoadToyGraphs(&g , path , 0);
    }else if(command == "stadiums") {
        g.clear();
        path = "../Graphs/Toy-Graphs/Toy-Graphs";
        LoadToyGraphs(&g, path, 1);
    }else if(command == "tourism") {
        g.clear();
        path = "../Graphs/Toy-Graphs/Toy-Graphs";
        LoadToyGraphs(&g, path, 2);
    }else{
        std::cerr << "Invalid file path" << std::endl;
    }

}



