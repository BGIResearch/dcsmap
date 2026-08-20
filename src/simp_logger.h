#ifndef CWORKFLOW_ENGINE_SIMP_LOGGER_H_
#define CWORKFLOW_ENGINE_SIMP_LOGGER_H_

#include <mutex>

// Thread-safe, timestamped stderr logging. safe_fprintf_stderr prepends
// "[YYYY-MM-DD HH:MM:SS] "; the SIMP_LOG_* macros add "[level] [func:line] ".

// Shared mutex serializing all writes that go to the process console (stderr
// log lines, and the task-log replay written by process.cc). Holding it makes
// a multi-line block -- e.g. one task's replayed stdout/stderr -- atomic
// w.r.t. other tasks and w.r.t. SIMP_LOG_* lines, so blocks never interleave.
std::mutex& ConsoleIoMutex();

int safe_fprintf_stderr(const char *format, ...);

#define SIMP_LOG_INFO(msg, ...) \
    do { \
        safe_fprintf_stderr("[info] [%s:%d] " msg "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    } while (0)

#define SIMP_LOG_WARN(msg, ...) \
    do { \
        safe_fprintf_stderr("[warning] [%s:%d] " msg "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    } while (0)

#define SIMP_LOG_ERROR(msg, ...) \
    do { \
        safe_fprintf_stderr("[error] [%s:%d] " msg "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    } while (0)


#endif  // CWORKFLOW_ENGINE_SIMP_LOGGER_H_
