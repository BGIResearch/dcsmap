#pragma once

#include <atomic>
#include <string>
#include <thread>
#include <vector>

// Result of running a Task: the bash exit code and the list of expected output
// file paths (absolute). `outputs` is only filled when `exit_code == 0` and all
// expected files exist on disk.
struct TaskResult {
    int exit_code = -1;
    std::vector<std::string> outputs;
};

// A Task wraps a single bash script that runs in its own job subdirectory
// `<work_root>/task-<name>`. It writes `command.sh`, runs it with bash, and
// redirects stdout/stderr to `command.sh.log.o` / `command.sh.log.e`.
//
// Threading model: start() spawns a std::thread running run(); join() waits for
// it. This is intentionally minimal (no scheduler) - the caller chains tasks
// with start()/join() to express dependencies and parallelism.
class Task {
public:
    Task(std::string name,
         std::string work_root,
         std::string script,
         std::vector<std::string> expected_outputs = {});

    // Spawn the worker thread. Must be called at most once.
    void start();

    // Wait for the worker thread to finish.
    void join();

    // Read-only accessors (valid after join()).
    const std::string& name() const { return name_; }
    const std::string& jobDir() const { return job_dir_; }
    const TaskResult& result() const { return result_; }

private:
    void run();  // executed on thread_

    std::string name_;
    std::string work_root_;
    std::string job_dir_;          // <work_root>/task-<name>
    std::string script_;
    std::vector<std::string> expected_outputs_;

    std::thread thread_;
    TaskResult result_;
};
