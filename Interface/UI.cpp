#include "UI.h"
#include "../Logic/Logic.h"
#include "../Logic/Algorithms.h"
#include <thread>
#include <chrono>
#include <iomanip>

UI::UI() {}

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
 * @brief Loads data and initializes the program.
 *
 * @param ui Reference to the UI object.
 *
 * @note Time complexity: Depends on the loading process, including file I/O and graph creation.
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
                    break;
                case 'B':
                    LoadToyGraphs(&g , path , 1);
                    break;
                case 'C':
                    LoadToyGraphs(&g , path , 2);
                    break;
                case 'D':
                    LoadToyGraphs(&g , path , 3);
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
                    break;
                case 'B':
                    LoadRealWorldGraphs(&g , path , 1);
                    break;
                case 'C':
                    LoadRealWorldGraphs(&g , path , 2);
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
                    break;
                case 'B':
                    LoadToyGraphs(&g , path , 1);
                    break;
                case 'C':
                    LoadToyGraphs(&g , path , 2);
                    break;
                case 'D':
                    LoadToyGraphs(&g , path , 3);
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
                    break;
                case 'B':
                    LoadRealWorldGraphs(&g , path , 1);
                    break;
                case 'C':
                    LoadRealWorldGraphs(&g , path , 2);
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
 * @brief Displays the main menu of the program.
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
 * @brief Displays the main menu of the program after clearing the screen.
 *
 * @note Time complexity: O(1)
 */
void UI::main_menu(){
    clear_screen();
    char op;
    std::cout << "What would you like to know?" <<std::endl;
    std::cout << "A. Switch Data Set" <<std::endl
              << "B. Graph Information Options" << std::endl
              << "C. Calculate TSP with Backtracking Brute Force" << std::endl
              << "D. Calculate TSP with Triangular Aprox" << std::endl
              << "E. Calculate TSP with Nearest Neighbour algorithm" << std::endl
              << "F. Calculate TSP with Nearest Neighbour with Backtrack algorithm for incomplete Graphs" << std::endl
             << "G Exit the program" << std::endl
             << "Insert your choice:";

    validate_input(op, 'A', 'G');

    int index = 0 , index2 = 0;
    Vertex * testV1 ;
    Vertex * testeV2;
    double weight = 0 , lastWeight = DBL_MAX;
    std::vector<Vertex *> articulationPoints;
    std::vector<std::vector<Vertex *>> scc;
    double bestValue1 = 0.0 , bestValue2 = 0.0;
    Vertex * source;
    Vertex * dest;
    switch(op){
        case 'A':
            changeDataSet();
            main_menu();
            break;
        case 'B':
            GraphOptionsMenu();
            break;
        case 'C':
            tspBruteForce(&g);
            back_menu();
            break;
        case 'D':
            source = g.findVertex(0);

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

            articulationPoints.clear();
            nn_with_backtracking(&g , source , articulationPoints);

            back_menu();
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
 * @brief Retrieves the graph associated with the UI object.
 *
 * @return The graph object.
 *
 * @note Time complexity: O(1)
 */
Graph UI::getGraph() const {
    return g;
}

void UI::back_menu(){
    char op;
    std::cout << "Press A to go back to the menu: ";
    validate_input(op,'A','A');
    main_menu();
}

void UI::back_menu_GraphOptions() {
    char op;
    std::cout << "Press A to go back to the menu: ";
    validate_input(op,'A','A');
    GraphOptionsMenu();
}


