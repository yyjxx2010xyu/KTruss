#ifndef TIMER
#define TIMER

#include <iostream>
#include <chrono>
#include <iomanip>
#include "tools.h"

class Timer
{
private:
    uint64_t _startTime;
    uint64_t _endTime;
    uint64_t _lastTime;
    bool _timerState;
public:
    explicit Timer()
    {
        _startTime = _endTime = _lastTime = 0;
        _timerState = false;
    }
    uint64_t GetSysTime()
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    }
    void StartTimer()
    {
        _timerState = true;
        _startTime = _lastTime = GetSysTime();
    }
    void EndTimer()
    {
        _endTime = GetSysTime();
        _timerState = false;
    }
    void UpdateLastTimer()
    {
        _lastTime = GetSysTime();
    }
    void GapTimer(const char * str)
    {
        uint64_t curTime = GetSysTime();
        std::cout << CYAN << "Timer Info" << RESET;
        std::cout << std::setfill('-') << std::setw(15) << str << " : "<< \
        std::fixed << std::setprecision(4) <<  static_cast<double>(curTime - _lastTime)/ 1e6 << "s";
        std::cout << " Total Time Usage: " << static_cast<double>(curTime - _startTime) / 1e6 << "s" << std::endl;
        _lastTime = curTime;
        std::cout << std::setfill(' ');
    }
    double QueryTimer()
    {
        return static_cast<double>(_endTime - _startTime) / 1e6;
    }
};

#endif