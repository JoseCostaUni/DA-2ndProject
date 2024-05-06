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

    static bool validate_int_input(int &index);
     /**
     * @brief Displays the main menu options.
     */
    void main_menu();
     /**
     * @brief Executes the max flow algorithm.
     */
    void back_menu();

    Graph getGraph() const;
private:
    Graph g;
    std::vector<Edge *> TSP;
};


#endif //PROJETO_UI_H
