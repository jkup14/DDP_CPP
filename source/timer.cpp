#include <chrono>

class Timer {
    public:
        Timer() {
            m_startTime = std::chrono::high_resolution_clock::now();
        }

        ~Timer () {
            Stop();
        }

        long long Stop() {
            auto endTime = std::chrono::high_resolution_clock::now();
            auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_startTime).time_since_epoch().count();
            auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTime).time_since_epoch().count();
            auto duration = end - start;
            return duration * 0.001;
        }

    private:
        std::chrono::time_point<std::chrono::high_resolution_clock> m_startTime;
};