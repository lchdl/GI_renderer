#pragma once
#include <cstdint>
#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <unistd.h>
#endif

class Timer {
protected:
    int64_t last, now, freq;

protected:
    int64_t _getTicks() {
        LARGE_INTEGER ticks;
        QueryPerformanceCounter(&ticks);
        return ticks.QuadPart;
    }
    int64_t _getFreq() {
        LARGE_INTEGER freq;
        QueryPerformanceFrequency(&freq);
        return freq.QuadPart;
    }

public:
    Timer() {
        last = _getTicks();
        now = last;
        freq = _getFreq();
    }
    double tick() {
        /* get elapsed time since last tick() call */
        now = _getTicks();
        int64_t diff = now - last;
        double dt = double(diff) / double(freq);
        last = now;
        return dt;
    }
    double query() {
        /* get elapsed time since last tick() call but it is not recorded */
        now = _getTicks();
        int64_t diff = now - last;
        double dt = double(diff) / double(freq);
        return dt;
    }
};
