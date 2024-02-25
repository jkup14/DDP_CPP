#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <vector>

class Timer {
    public:
        ~Timer () ;

        void Start();

        long long Stop();

        long long avg_Time();

        long long total_Time();

    private:
        std::chrono::time_point<std::chrono::high_resolution_clock> m_startTime;
        std::vector<long long> times;
};

#endif