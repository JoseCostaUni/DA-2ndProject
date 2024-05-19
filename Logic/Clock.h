


#ifndef PROJETO2_CLOCK_H
#define PROJETO2_CLOCK_H

#include "../stdafx.h"

/**
 * @class Clock
 * @file Clock.h
 * @brief A simple clock class to measure elapsed time.
 *
 * The Clock class provides functionality to start a timer and measure the elapsed time
 * since the timer was started.
 */

class Clock {
public:

    /**
     * @brief Default constructor for the Clock class.
     */
    Clock() = default;

    /**
    * @brief Starts the clock.
    *
    * This function records the current time as the start time.
    */
    void start() {
        startTime = std::chrono::high_resolution_clock::now();
    }

    /**
     * @brief Prints the elapsed time since the clock was started.
     *
     * This function calculates the elapsed time since the start time and prints it to the standard output.
     */
    void elapsed() {
        auto endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsedSeconds = endTime - startTime;
        std::cout << "Elapsed time: " << elapsedSeconds.count() << " seconds" << std::endl;
    }

private:
    std::chrono::high_resolution_clock::time_point startTime;

};


#endif //PROJETO2_CLOCK_H
