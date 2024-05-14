


#ifndef PROJETO2_CLOCK_H
#define PROJETO2_CLOCK_H

#include "../stdafx.h"

class Clock {
public:
    Clock() = default;

    void start() {
        startTime = std::chrono::high_resolution_clock::now();
    }

    void elapsed() {
        auto endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsedSeconds = endTime - startTime;

        std::cout << "Elapsed time: " << elapsedSeconds.count() << "s" << std::endl;
    }

private:
    std::chrono::high_resolution_clock::time_point startTime;

};


#endif //PROJETO2_CLOCK_H
