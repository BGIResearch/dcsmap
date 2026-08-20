#include "simp_logger.h"

#include <cstdarg>
#include <ctime>
#include <cstdio>
#include <mutex>

// Singleton mutex behind ConsoleIoMutex(): shared by safe_fprintf_stderr and
// the task-log replay in process.cc so console output blocks stay atomic.
std::mutex& ConsoleIoMutex() {
    static std::mutex mu;
    return mu;
}

static char* get_formatted_time(char *buffer, int buffer_size) {
    time_t now;
    struct tm *local_time;

    // Get current time
    time(&now);

    // Convert to local time
    local_time = localtime(&now);

    // Format time as "YYYY-MM-DD HH:MM:SS"
    size_t ret = strftime(buffer, buffer_size, "%Y-%m-%d %H:%M:%S", local_time);
    if (ret == 0) {
        fprintf(stderr, "Error! get_formatted_time buffer size small.\n");
    }

    return buffer;
}


int safe_fprintf_stderr(const char *format, ...) {
    va_list args;
    int result;

    // Lock the shared console mutex, perform fprintf, unlock.
    std::lock_guard<std::mutex> guard(ConsoleIoMutex());
    const int buffer_size = 4096;
    char buffer[buffer_size];
    fprintf(stderr, "[%s] ", get_formatted_time(buffer, buffer_size));
    va_start(args, format);
    result = vfprintf(stderr, format, args);
    va_end(args);

    // Flush to ensure complete output before unlocking
    fflush(stderr);

    return result;
}
