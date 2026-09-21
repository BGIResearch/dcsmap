#include "Task.h"
#include "simp_logger.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace {

bool dirExists(const std::string& p) {
    struct stat st;
    return stat(p.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

bool makeDir(const std::string& p) {
    // mkdir -p
    std::string acc;
    for (size_t i = 0; i < p.size(); ++i) {
        acc += p[i];
        if (p[i] == '/' && !acc.empty()) {
            mkdir(acc.c_str(), 0755);
        }
    }
    mkdir(p.c_str(), 0755);
    return dirExists(p);
}

bool fileExists(const std::string& p) {
    struct stat st;
    return stat(p.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

std::string joinPath(const std::string& a, const std::string& b) {
    if (a.empty()) return b;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

// Write the task return code to <jobdir>/command.sh.rc.
void writeRc(const std::string& jobdir, int code) {
    std::ofstream f(joinPath(jobdir, "command.sh.rc"));
    if (f) f << code << "\n";
}

// Print the last N lines of a task's command.sh.log.e to stderr, to help
// diagnose failures.
void dumpLogTail(const std::string& loge, const std::string& task_name, int n = 100) {
    std::ifstream f(loge);
    if (!f) {
        SIMP_LOG_ERROR("cannot open task log.e: %s", loge.c_str());
        return;
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(f, line)) lines.push_back(line);
    size_t total = lines.size();
    size_t start = total > static_cast<size_t>(n) ? total - n : 0;
    SIMP_LOG_ERROR("task '%s': last %zu lines of %s (total %zu)",
                   task_name.c_str(), (total - start), loge.c_str(), total);
    for (size_t i = start; i < total; ++i) {
        std::fprintf(stderr, "%s\n", lines[i].c_str());
    }
    SIMP_LOG_ERROR("task '%s': end of %s", task_name.c_str(), loge.c_str());
}

}  // namespace

Task::Task(std::string name,
           std::string work_root,
           std::string script,
           std::vector<std::string> expected_outputs)
    : name_(std::move(name)),
      work_root_(std::move(work_root)),
      script_(std::move(script)),
      expected_outputs_(std::move(expected_outputs)) {
    job_dir_ = joinPath(work_root_, "task-" + name_);
}

void Task::start() {
    thread_ = std::thread(&Task::run, this);
}

void Task::join() {
    if (thread_.joinable()) thread_.join();
}

void Task::run() {
    // 1. create job subdir
    if (!makeDir(job_dir_)) {
        SIMP_LOG_ERROR("Task '%s': cannot create job dir %s", name_.c_str(), job_dir_.c_str());
        result_.exit_code = -1;
        writeRc(job_dir_, result_.exit_code);
        return;
    }

    // 2. write command.sh
    // Inject a per-task temp directory preamble (WDL-style): create a tmp dir
    // inside the job dir, export TMPDIR and _JAVA_OPTIONS. Cleanup is handled by
    // cleanWorkDir (no EXIT trap).
    const std::string marker = "set -euo pipefail\n";
    std::string preamble =
        "#!/bin/bash\n"
        "set -euo pipefail\n"
        "cd '" + job_dir_ + "'\n"
        "_tmp_dir=\"$(mktemp -d '" + job_dir_ + "/tmp.XXXXXX')\"\n"
        "export TMPDIR=\"${_tmp_dir}\"\n"
        "export _JAVA_OPTIONS=\"-Djava.io.tmpdir=${_tmp_dir}\"\n"
        "# resolve an awk implementation: prefer mawk, fall back to awk\n"
        "if command -v mawk >/dev/null 2>&1; then\n"
        "    AWK=mawk\n"
        "elif command -v awk >/dev/null 2>&1; then\n"
        "    AWK=awk\n"
        "else\n"
        "    echo \"[ERROR] neither mawk nor awk found in PATH\" >&2\n"
        "    exit 1\n"
        "fi\n"
        "export AWK\n"
        "# resolve a GNU-compatible timer: it must exist *and* accept -f, since\n"
        "# /usr/bin/time may be absent (slim images) or a BSD/busybox variant\n"
        "TIME_CMD=\"\"\n"
        "for _t in /usr/bin/time \"$(command -v gtime || true)\" \"$(command -v time || true)\"; do\n"
        "    [ -n \"${_t}\" ] && [ -x \"${_t}\" ] || continue\n"
        "    if \"${_t}\" -f \"%e\" true >/dev/null 2>&1; then TIME_CMD=\"${_t}\"; break; fi\n"
        "done\n"
        "unset _t\n"
        "if [ -z \"${TIME_CMD}\" ]; then\n"
        "    echo \"[WARN] no GNU-compatible 'time' found; running without resource timing\" >&2\n"
        "fi\n"
        "# run a command under the timer, or unmeasured when no timer is available\n"
        "timed() {\n"
        "    local _label=\"$1\"; shift\n"
        "    if [ -n \"${TIME_CMD}\" ]; then\n"
        "        \"${TIME_CMD}\" -f \"[${_label}] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M\" \"$@\"\n"
        "    else\n"
        "        \"$@\"\n"
        "    fi\n"
        "}\n";
    std::string full_script;
    if (script_.rfind(marker, 0) == 0) {
        full_script = preamble + script_.substr(marker.size());
    } else {
        full_script = preamble + script_;
    }

    const std::string cmd_sh = joinPath(job_dir_, "command.sh");
    {
        std::ofstream f(cmd_sh, std::ios::binary);
        if (!f) {
            SIMP_LOG_ERROR("Task '%s': cannot write %s", name_.c_str(), cmd_sh.c_str());
            result_.exit_code = -1;
            writeRc(job_dir_, result_.exit_code);
            return;
        }
        f << full_script;
    }

    // 3. run bash command.sh > log.o 2> log.e
    // The script itself cds into its job dir (see preamble), so each bash
    // subprocess has its own cwd. We must NOT chdir() here because chdir() is
    // process-global and would race between parallel task threads.
    const std::string log_o = joinPath(job_dir_, "command.sh.log.o");
    const std::string log_e = joinPath(job_dir_, "command.sh.log.e");
    std::ostringstream shell;
    // Use bash to run the script; redirect stdout/stderr to the log files.
    shell << "bash " << " '" << cmd_sh << "' > '" << log_o << "' 2> '" << log_e
          << "'";
    SIMP_LOG_INFO("task '%s' start (dir=%s)", name_.c_str(), job_dir_.c_str());
    int rc = std::system(shell.str().c_str());
    int code = WIFEXITED(rc) ? WEXITSTATUS(rc) : -1;
    result_.exit_code = code;
    writeRc(job_dir_, code);
    SIMP_LOG_INFO("task '%s' end (exit=%d)", name_.c_str(), code);

    if (code != 0) {
        SIMP_LOG_ERROR("Task '%s' failed (exit %d). See logs: %s | %s",
                       name_.c_str(), code, log_o.c_str(), log_e.c_str());
        dumpLogTail(log_e, name_);
        return;
    }

    // 4. verify expected outputs exist
    for (const auto& out : expected_outputs_) {
        if (!fileExists(out)) {
            SIMP_LOG_ERROR("Task '%s': expected output missing: %s", name_.c_str(), out.c_str());
            result_.exit_code = -1;
            return;
        }
        result_.outputs.push_back(out);
    }
}
