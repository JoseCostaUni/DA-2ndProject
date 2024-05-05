#include "UI.h"
#include "../Logic/Logic.h"
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

/**
 * @brief Loads data and initializes the program.
 *
 * @param ui Reference to the UI object.
 *
 * @note Time complexity: Depends on the loading process, including file I/O and graph creation.
 */
void UI::loading_stuff(UI &ui) {
    char op;
    std::string path;
    std::cout << "Which data set do you wish to use during the analysis?" << std::endl
      << "A. Small Data Set" << std::endl
      << "B. Large Data Set" << std::endl
      << "Insert the letter: " ;
    validate_input(op,'A','B');

    if(op == 'A'){
        path = "SmallDataSet";
    }else{
        path = "LargeDataSet";
    }

    LoadFireStations(path);
    LoadWaterReservoirs(path);
    LoadPipes(path);
    LoadCities(path);
    createGraph(&g);

    DeliverySite supersource = DeliverySite("SuperSource");
    DeliverySite supersink = DeliverySite("SuperSink");
    DeliverySite dummy = DeliverySite("Empty");
    createSuperSourceSink(&g,supersource,supersink);
    inital_max_flow = edmondsKarp(&g,supersource,supersink, dummy);
    removeSuperSourceSink(&g,supersource,supersink);

    for(Vertex<DeliverySite>* v: g.getVertexSet()){

        int sumFlow = calculate_incoming_flow(v);
        v->setIncomingFlow(sumFlow);

        if(v->getInfo().getNodeType() == CITY){
            codeToFlow[v->getInfo().getCode()] = sumFlow;
        }
    }


    std::cout << "Load Finished" << std::endl;
    std::cout << "Press A to start the program: ";
    op;
    validate_input(op,'A','A');
    menu_start();
}

/**
 * @brief Displays the main menu of the program.
 *
 * @note Time complexity: O(1)
 */
void UI::menu_start() {
    char op;
    std::cout << "######################################################################" << std::endl
         << "@ @ @  @@@@@  @@@@@  @@@@@  @@@@@           @@@@@   @@@    @@@   @    " << std::endl
         << "@ @ @  @   @    @    @@     @  @@             @    @   @  @   @  @    " << std::endl
         << "@ @ @  @@@@@    @    @@@@@  @@@@              @    @   @  @   @  @    " << std::endl
         << "@ @ @  @   @    @    @@     @@ @@             @    @   @  @   @  @    " << std::endl
         << " @@@   @   @    @    @@@@@  @@ @@             @     @@@    @@@   @@@@@" << std::endl
         << "######################################################################" << std::endl << '\n'
         << "Welcome to the Analysis Tool for Water Supply Management, what would you like to do?" << std::endl
         << "A. Proceed to the application" << std::endl
         << "B. Close the application" << std::endl
         << "Insert the letter: " ;
    validate_input(op,'A','B');
    switch(op){
        case 'A':
            main_menu();
            break;
        case 'B':
            std::cout << "Thanks for using our analysis tool for water supply management!" << std::endl << "\n"
                 << "Made by: " << std::endl
                 << "Ângelo Oliveira || 202207798" << std::endl
                 << "José Costa      || 202207871" << std::endl
                 << "Bruno Fortes    || 202209730" << std::endl;
            std::this_thread::sleep_for(std::chrono::seconds(2)); // Espera 2 segundos antes de fechar o terminal
            exit(0);

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
   std::cout << "A. Run Max Flow algorithm" <<std::endl
             << "B. Check if every city meets it's water demand" <<std::endl
             << "C. Check heuristic stats of the max flow" <<std::endl
             << "D. Evaluate network's resiliency" <<std::endl
             << "E. Exit the program" << std::endl
             << "Insert your choice:";

    validate_input(op, 'A', 'E');
    switch(op){
        case 'A':
            max_flow();
            break;
        case 'B':
            check_demand();
            break;
        case 'C':
            check_heuristic();
            break;
        case 'D':
            evaluate_resiliency();
            break;
        case 'E':
           std::cout << "Thanks for using our water management tool!" <<std::endl << "\n"
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
Graph<DeliverySite> UI::getGraph() const {
    return g;
}

/**
 * @brief Executes the max flow algorithm based on user input.
 *
 * @note Time complexity: Depends on the execution of the max flow algorithm, including graph traversal and flow calculations.
 */
void UI::max_flow(){
    char op;
    std::cout << "What would you like to do?" << std::endl
              << "A. Calculate max flow for the entire network" << std::endl
              << "B. Calculate max flow for a specific city" << std::endl
              << "Insert your choice:";

    validate_input(op, 'A', 'B');
    switch (op) {
        case 'A': {
            char op1;
            std:: cout << "Would you like to only get the max flow result or the listing of max flows of the cities" << std::endl
                       << "A. Only the result" << std::endl
                       << "B. Listing for all cities" << std::endl
                       << "Insert your choice:" << std::endl;
            validate_input(op1, 'A', 'B');

            DeliverySite supersource = DeliverySite("SuperSource");
            DeliverySite supersink = DeliverySite("SuperSink");
            DeliverySite dummy = DeliverySite("Empty");
            createSuperSourceSink(&g,supersource,supersink);
            double max_flow = edmondsKarp(&g,supersource,supersink, dummy);
            removeSuperSourceSink(&g,supersource,supersink);

            switch (op1) {
                case 'A': {
                    std::cout << "The max flow for the entire network is: " << max_flow << std::endl;
                    break;
                }
                case 'B': {
                    std::cout << std::left << std::setw(20) << "City Name"
                              << std::setw(20) << "City Code"
                              << std::setw(20) << "Flow" << std::endl;

                    for(auto v : g.getVertexSet()){
                        if(v->getInfo().getNodeType() == CITY){
                            std::cout << std::left << std::setw(20) << v->getInfo().getName()
                                      << std::setw(20) << v->getInfo().getCode()
                                      << std::setw(20) << v->calculateIncomingFlow() << std::endl;
                        }

                    }
                    std::cout << "The max flow for the entire network is: " << max_flow << std::endl;

                    break;
                }
            }
            break;
        }
        case 'B':{
            std::string code;
            bool foundVertex = false;
            DeliverySite supersource = DeliverySite("SuperSource");

            while (!foundVertex) {
                std::cout << "Insert the code of the city: " << std::endl;
                std::cin >> code;

                DeliverySite sink = DeliverySite(code);

                if (!g.findVertex(sink)) {
                    std::cout << "Error: City with code '" << code << "' not found. Please try again." << std::endl;
                } else {
                    foundVertex = true;
                }
            }
            DeliverySite sink = DeliverySite(code);

            createSuperSource(&g,supersource);
            DeliverySite dummy = DeliverySite("Empty");
            double max_flow = edmondsKarp(&g,supersource,sink,dummy);
            removeSuperSource(&g,supersource);

            std::cout << "The max flow for the city " << code << " is: " << max_flow << std::endl;

            break;
        }
    }
    back_menu();
};

/**
 * @brief Checks if every city meets its water demand.
 *
 * @note Time complexity: Depends on the graph traversal and demand checking process.
 */
void UI::check_demand(){

    DeliverySite supersource = DeliverySite("SuperSource");
    DeliverySite supersink = DeliverySite("SuperSink");
    DeliverySite dummy = DeliverySite("Empty");
    createSuperSourceSink(&g,supersource,supersink);
    edmondsKarp(&g,supersource,supersink,dummy);
    removeSuperSourceSink(&g,supersource,supersink);

    std::cout << "The following cities don't receive enough water : " << std::endl;

    std::cout << std::left << std::left << std::setw(20) << "City Name"
              << std::setw(20) << "City Code"
              << std::setw(20) << "Demand"
              << std::setw(20) << "Flow"
              << std::setw(20) << "Defecit" << std::endl;

    for(Vertex<DeliverySite>* ds: g.getVertexSet()){
        int sumFlow = calculate_incoming_flow(ds);
        ds->setIncomingFlow(sumFlow);

        int difference = ds->getInfo().getDemand() - ds->getIncomingFlow();

        if(ds->getInfo().getNodeType() == CITY && difference > 0 ){
            std::cout << std::left << std::setw(20) << ds->getInfo().getName()
                      << std::setw(20) << ds->getInfo().getCode()
                      << std::setw(20) << ds->getInfo().getDemand()
                      << std::setw(20) << ds->getIncomingFlow()
                      << std::setw(20) << abs(difference) << std::endl;
        }
    }
    std::cout << std::endl;
    back_menu();
}

/**
 * @brief Calculates and displays heuristic stats of the max flow.
 *
 * @note Time complexity: Depends on the execution of the heuristic algorithm.
 */
void UI::check_heuristic(){
    DeliverySite supersource = DeliverySite("SuperSource");
    DeliverySite supersink = DeliverySite("SuperSink");
    DeliverySite dummy = DeliverySite("Empty");

    createSuperSourceSink(&g,supersource,supersink);
    edmondsKarp(&g,supersource,supersink,dummy);
    removeSuperSourceSink(&g,supersource,supersink);

    heuristic(&g);
    std::cout << std::endl;
    back_menu();
}

/**
 * @brief Redirects the user back to the main menu.
 *
 * @note Time complexity: O(1)
 */
void UI::back_menu(){
    char op;
    std::cout << "Press A to go back to the menu: ";
    validate_input(op,'A','A');
    main_menu();
}

/**
 * @brief Displays the menu to evaluate the network's resiliency and performs the chosen evaluation.
 *
 * @note Time complexity: Depends on the chosen evaluation method and associated computations.
 */
void UI::evaluate_resiliency() {
    char op;
    std:: cout << "How would you like to evaluate the resiliency?" << std::endl
               << "A. Water Reservoir out of comission" << std::endl
               << "B. Pumping Station out of comission /in maintenance" << std::endl
               << "C. Pipeline/s out of comission / ruptured" << std::endl
               << "Insert your choice:" << std::endl;
    validate_input(op, 'A', 'C');
    switch(op){
        case 'A':{
            std::string code;
            bool foundVertex = false;

            while (!foundVertex) {
                std::cout << "Insert the code of the water reservoir: " << std::endl;
                std::cin >> code;

                DeliverySite water_reservoir = DeliverySite(code);

                if (!g.findVertex(water_reservoir)) {
                    std::cout << "Error: Water Reservoir with code '" << code << "' not found. Please try again." << std::endl;
                } else {
                    foundVertex = true;
                }
            }

            std::cout << "Would you like to try the balancing algorithm without the whole network?" << std::endl
                      << "A. Yes" << std::endl
                      << "B. No" << std::endl
                      << "Insert your choice:" << std::endl;
            validate_input(op, 'A', 'B');

            switch(op){
                case 'A':
                    redistributeWithoutMaxFlowVersion2(code);
                    break;
                case 'B':
                    redistributeWithoutMaxFlow(code);
                    break;

            }
            break;
        }
        case 'B':{
            std::string code;
            bool foundVertex = false;

            while (!foundVertex) {
                std::cout << "Insert the code of the pumping station: " << std::endl;
                std::cin >> code;

                DeliverySite pump_station = DeliverySite(code);

                if (!g.findVertex(pump_station)) {
                    std::cout << "Error: Pumping Station with code '" << code << "' not found. Please try again." << std::endl;
                } else {
                    foundVertex = true;
                }
            }

            DeliverySite supersource = DeliverySite("SuperSource");
            DeliverySite supersink = DeliverySite("SuperSink");
            DeliverySite pump_station = DeliverySite(code);
            createSuperSourceSink(&g,supersource,supersink);
            double max_flow = edmondsKarp(&g,supersource,supersink,pump_station);
            removeSuperSourceSink(&g,supersource,supersink);

            std::cout << "The max flow of the network removing " << code << " is: " << max_flow << std::endl;

            std::cout << std::left << std::setw(20) << "City Name"
                      << std::setw(20) << "City Code"
                      << std::setw(20) << "New Units"
                      << std::setw(20) << "New Flow"
                      << std::setw(20) << "Old Flow" << std::endl;

            for(Vertex<DeliverySite>* v: g.getVertexSet()){
                if(v->getInfo().getNodeType() == CITY){

                    int sumFlow = calculate_incoming_flow(v);

                    v->setIncomingFlow(sumFlow);
                    int result = v->getIncomingFlow() - codeToFlow[v->getInfo().getCode()];
                    if(result < 0){
                        std::cout << std::left << std::setw(20) << v->getInfo().getName()
                                  << std::setw(20) << v->getInfo().getCode()
                                  << std::setw(20) << abs(result)
                                  << std::setw(20) << v->getIncomingFlow()
                                  << std::setw(20) << codeToFlow[v->getInfo().getCode()] << std::endl;
                    }
                }
            }

            break;
        }
        case 'C':{
            std::string code1;
            std::string code2;
            bool foundEdge = false;
            std::vector<Edge<DeliverySite>*> pointerVector;

            while (!foundEdge) {
                Edge<DeliverySite>* edgeFound = nullptr;
                std::cout << "Insert the code of the first delivery site: " << std::endl;
                std::cin >> code1;

                std::cout << "Insert the code of the second delivery site: " << std::endl;
                std::cin >> code2;

                DeliverySite ds1 = DeliverySite(code1);

                DeliverySite ds2 = DeliverySite(code2);

                if (!g.findVertex(ds1) && !g.findVertex(ds2)) {
                    std::cout << "Error: Delivery Sites not found. Please try again." << std::endl;
                }

                for(auto edge: g.getEdges()){
                    if(edge->getOrig()->getInfo().getCode() == ds1.getCode() && edge->getDest()->getInfo().getCode() == ds2.getCode()){
                        foundEdge = true;
                        edgeFound = edge;
                        break;
                    }
                }

                if(!foundEdge){
                    std::cout << "Error: Pipe not found. Please try again." << std::endl;
                }else{
                    if(edgeFound->getReverse() != nullptr){
                        pointerVector.push_back(edgeFound->getReverse());
                    }
                    pointerVector.push_back(edgeFound);
                    std::cout << "Do you want to add another pipeline to the out of comission/ruptured list ?" << std::endl
                              << "A. Yes" << std::endl
                              << "B. No" << std::endl;
                    char op1;
                    validate_input(op1,'A','B');
                    if(op1 == 'A'){
                        foundEdge = false;
                    }
                }
            }

            DeliverySite supersource = DeliverySite("SuperSource");
            DeliverySite supersink = DeliverySite("SuperSink");
            createSuperSourceSink(&g,supersource,supersink);
            double max_flow = edmondsKarpPipe(&g,supersource,supersink,pointerVector);
            removeSuperSourceSink(&g,supersource,supersink);

            std::cout << "The max flow of the network is " << max_flow << " removing the pipes: " << std::endl;
            for(auto pipe: pointerVector){
                std::cout << pipe->getOrig()->getInfo().getCode() << "-" << pipe->getDest()->getInfo().getCode() << std::endl;
            }
            std::cout << std::endl;

            std::cout << std::left << std::setw(20) << "City Name"
                      << std::setw(20) << "City Code"
                      << std::setw(20) << "Required Units"
                      << std::setw(20) << "New Flow"
                      << std::setw(20) << "Old Flow" << std::endl;


            for(Vertex<DeliverySite>* v: g.getVertexSet()){
                if(v->getInfo().getNodeType() == CITY){

                    int sumFlow = calculate_incoming_flow(v);

                    v->setIncomingFlow(sumFlow);
                    int result = v->getIncomingFlow() - codeToFlow[v->getInfo().getCode()];
                    if(result < 0){
                        std::cout << std::left << std::setw(20) << v->getInfo().getName()
                                  << std::setw(20) << v->getInfo().getCode()
                                  << std::setw(20) << abs(result)
                                  << std::setw(20) << v->getIncomingFlow()
                                  << std::setw(20) << codeToFlow[v->getInfo().getCode()] << std::endl;
                    }
                }
            }
            std::cout << std::endl;
        }

        break;
    }

    back_menu();
}


/**
 * @brief Redistributes water without calculating the max flow.
 *
 * @param wr Code of the water reservoir to be removed.
 *
 * @note Time complexity: Depends on the graph traversal and redistribution process.
 */
void UI::redistributeWithoutMaxFlow(std::string wr) {

    DeliverySite reservoir_ds = DeliverySite(wr);
    Vertex<DeliverySite>* ds = g.findVertex(reservoir_ds);
    DeliverySite supersource = DeliverySite("SuperSource");
    DeliverySite supersink = DeliverySite("SuperSink");

    createSuperSourceSink(&g,supersource,supersink);
    inital_max_flow = edmondsKarp(&g,supersource,supersink, reservoir_ds);
    removeSuperSourceSink(&g,supersource,supersink);

    std::cout << "The max flow of the network removing " << ds->getInfo().getCode() << " is: " << inital_max_flow << std::endl;

    std::cout << std::left << std::setw(20) << "City Name"
              << std::setw(20) << "City Code"
              << std::setw(20) << "Required Units"
              << std::setw(20) << "New Flow"
              << std::setw(20) << "Old Flow" << std::endl;

    for (Vertex<DeliverySite> *v: g.getVertexSet()) {
        if (v->getInfo().getNodeType() == CITY) {
            int sumFlow = 0;
            for(auto e : v->getIncoming()){
                sumFlow += e->getFlow();
            }
            v->setIncomingFlow(sumFlow);

            int result = v->getIncomingFlow() - codeToFlow[v->getInfo().getCode()];
            if (result < 0) {
                std::cout << std::left << std::setw(20) << v->getInfo().getName()
                          << std::setw(20) << v->getInfo().getCode()
                          << std::setw(20) << abs(result)
                          << std::setw(20) << v->getIncomingFlow()
                          << std::setw(20) << codeToFlow[v->getInfo().getCode()] << std::endl;
            }
        }
    }

}

/**
 * @brief Redistributes water without calculating the max flow (version 2).
 *
 * @param wr_code Code of the water reservoir to be removed.
 *
 * @note Time complexity: Depends on the graph traversal and redistribution process.
 */
void UI::redistributeWithoutMaxFlowVersion2(std::string& wr_code){
    DeliverySite reservoir_ds = DeliverySite(wr_code);

    std::vector<std::vector<Edge<DeliverySite>*>> paths;

    Vertex<DeliverySite>* ds = g.findVertex(reservoir_ds);
    std::vector<Edge<DeliverySite>*> path;

    for(auto ver: g.getVertexSet()){
        ver->setVisited(false);
    }

    findAllPathsRedistribute(&g,ds,path,paths);

    double maxflow = redistributeWaterWithoutMaxFlow2(&g,paths);

    std::cout << "The max flow of the network removing " << ds->getInfo().getCode() << " is: " << maxflow << std::endl;

    std::cout << std::left << std::setw(20) << "City Name"
              << std::setw(20) << "City Code"
              << std::setw(20) << "Required Units"
              << std::setw(20) << "New Flow"
              << std::setw(20) << "Old Flow" << std::endl;

    for (Vertex<DeliverySite> *v: g.getVertexSet()) {
        if (v->getInfo().getNodeType() == CITY) {
            int sumFlow = 0;
            for(auto e : v->getIncoming()){
                sumFlow += e->getFlow();
            }
            v->setIncomingFlow(sumFlow);

            int result = v->getIncomingFlow() - codeToFlow[v->getInfo().getCode()];
            if (result < 0) {
                std::cout << std::left << std::setw(20) << v->getInfo().getName()
                          << std::setw(20) << v->getInfo().getCode()
                          << std::setw(20) << abs(result)
                          << std::setw(20) << v->getIncomingFlow()
                          << std::setw(20) << codeToFlow[v->getInfo().getCode()] << std::endl;
            }
        }
    }

}

/**
 * @brief Calculates the total incoming flow for a given vertex.
 *
 * @param v Pointer to the vertex for which incoming flow needs to be calculated.
 * @return The total incoming flow for the vertex.
 *
 * @note Time complexity: O(E), where E is the number of edges incident to the vertex.
 */
int UI::calculate_incoming_flow(Vertex<DeliverySite>* v){
    int sumFlow = 0;
    for(Edge<DeliverySite>* p : v->getIncoming()){
        if(!p->isSelected()) {
            sumFlow += p->getFlow();
        }
    }
    return sumFlow;
}
