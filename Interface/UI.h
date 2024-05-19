#ifndef PROJETO_UI_H
#define PROJETO_UI_H

#include "../Logic/LoadingFunctions.h"
#include "../stdafx.h"
#include "../Graph.h"

/**
 * @brief User Interface class for managing interaction with the water distribution system.
 */
 class UI {
public:
    /**
   * @brief Default constructor for the UI class.
   */
    UI();
     /**
     * @brief Displays the main menu and handles user input.
     */
    void menu_start();
     /**
      * @brief Clears the console screen.
      */
    static void clear_screen();
     /**
     * @brief Loads necessary data and settings for the UI.
     *
     * @param ui Reference to the UI object.
     */
    void loading_stuff(UI &ui);


    /**
     * @brief Displays the options to change the data set.
     * @return none.
     * @param none.
     * @complexity O(1).
     */
    void changeDataSet();
     /**
    * @brief Validates user input within a specified range.
    *
    * @param op Reference to the input character to be validated.
    * @param lower_bound The lower bound of the valid range.
    * @param upper_bound The upper bound of the valid range.
    * @return True if input is valid, False otherwise.
    */
    static bool validate_input(char &op, const char lower_bound, const char upper_bound);

    void GraphOptionsMenu();
    void BackTrackMenu();

    /**
     * @brief Validates user input as an integer.
     *
     * @param index Reference to the input integer to be validated.
     * @return True if input is valid, False otherwise.
     */

    static bool validate_int_input(int &index);
     /**
     * @brief Displays the main menu options.
     * @return none.
     * @param none.
     * @complexity O(1).
     */
    void main_menu();
     /**
     * @brief Displays the back menu options.
     * @param none.
     * @return none.
     * @complexity O(1).
     */
    void back_menu();

    /**
    * @brief Displays the Graph Options menu options.
    * @param none.
    * @return none.
    * @complexity O(1).
    */
    void back_menu_GraphOptions();

     /**
     * @brief Displays the BackTracking Options menu options.
     * @param none.
     * @return none.
     * @complexity O(1).
     */
    void back_menu_BacktrackingOptions();

    /**
     * @brief returns the graph.
     * @param none.
     * @return the graph.
     * @complexity O(1).
     */
    Graph getGraph() const;
    void pathSelector();
private:
    Graph g;
    std::vector<Edge *> TSP;
    std::string file_path;
    bool fullyConnected = false;
};


#endif //PROJETO_UI_H
