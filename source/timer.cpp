#include <chrono>
#include <numeric>
#include "../include/timer.hpp"

Timer::~Timer () {
    Stop();
}

void Timer::Start() {
    m_startTime = std::chrono::high_resolution_clock::now();
}

long long Timer::Stop() {
    auto endTime = std::chrono::high_resolution_clock::now();
    auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_startTime).time_since_epoch().count();
    auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTime).time_since_epoch().count();
    auto duration = end - start;
    times.push_back(duration);
    return duration;
}

long long Timer::avg_Time() {
    return std::accumulate(times.begin(), times.end(), 0)/times.size();
}

long long Timer::total_Time() {
    return std::accumulate(times.begin(), times.end(), 0);
}
