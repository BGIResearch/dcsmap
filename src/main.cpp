#include "Scripts.h"
#include "Task.h"
#include "simp_logger.h"

#include <cstdio>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <sstream>
#include <algorithm>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#ifndef DCSMAP_VERSION
#define DCSMAP_VERSION "unknown"
#endif
#ifndef DCSMAP_GIT_COMMIT
#define DCSMAP_GIT_COMMIT "unknown"
#endif
static const char* kVersion = DCSMAP_VERSION "-" DCSMAP_GIT_COMMIT;

namespace {

bool fileExists(const std::string& p) {
    struct stat st;
    return stat(p.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

bool makeDirP(const std::string& p) {
    std::string acc;
    for (size_t i = 0; i < p.size(); ++i) {
        acc += p[i];
        if (p[i] == '/' && !acc.empty()) mkdir(acc.c_str(), 0755);
    }
    mkdir(p.c_str(), 0755);
    struct stat st;
    return stat(p.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

std::string joinPath(const std::string& a, const std::string& b) {
    if (a.empty()) return b;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

std::string canonical(const std::string& p) {
    char buf[PATH_MAX];
    if (realpath(p.c_str(), buf)) return std::string(buf);
    return p;
}

std::string absolute(const std::string& p) {
    if (!p.empty() && p[0] == '/') return canonical(p);
    char cwd[PATH_MAX];
    if (!getcwd(cwd, sizeof(cwd))) return p;
    return canonical(std::string(cwd) + "/" + p);
}

// Make a path absolute without resolving symlinks. Needed for --ref-fasta:
// the linear aligner (bwa-mem2) looks up its index (.0123/.bwt.2bit.64/...)
// next to the ref path, and the index may live in a directory reached
// through a symlink. realpath() would collapse that symlink and break index
// lookup. We only normalize "." and ".." lexically and prepend cwd for
// relative paths.
std::string absoluteNoSymlink(const std::string& p) {
    std::string base = p;
    if (p.empty() || p[0] != '/') {
        char cwd[PATH_MAX];
        if (!getcwd(cwd, sizeof(cwd))) return p;
        base = std::string(cwd) + "/" + p;
    }
    std::vector<std::string> parts;
    std::string seg;
    for (size_t i = 0; i <= base.size(); ++i) {
        if (i == base.size() || base[i] == '/') {
            if (!seg.empty()) {
                if (seg == ".") { /* skip */ }
                else if (seg == "..") { if (!parts.empty()) parts.pop_back(); }
                else parts.push_back(seg);
            }
            seg.clear();
        } else {
            seg += base[i];
        }
    }
    std::string out;
    for (const auto& p2 : parts) out += "/" + p2;
    return out.empty() ? std::string("/") : out;
}

std::string dirName(const std::string& p) {
    size_t s = p.find_last_of('/');
    return (s == std::string::npos) ? "." : p.substr(0, s);
}

std::string firstLine(const std::string& path) {
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty()) return line;
    }
    return line;
}

int runCapture(const std::string& cmd, std::string& out) {
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return -1;
    char tmp[1024];
    while (fgets(tmp, sizeof(tmp), pipe)) out += tmp;
    int rc = pclose(pipe);
    return WIFEXITED(rc) ? WEXITSTATUS(rc) : -1;
}

std::string mktempWork(const std::string& dir) {
    makeDirP(dir);
    std::ostringstream cmd;
    cmd << "mktemp -d -p '" << dir << "' work.XXXXXX 2>/dev/null";
    std::string out;
    if (runCapture(cmd.str(), out) == 0) {
        while (!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
        if (!out.empty()) return out;
    }
    std::string fb;
    runCapture("mktemp -d -t work.XXXXXX 2>/dev/null", fb);
    while (!fb.empty() && (fb.back() == '\n' || fb.back() == '\r')) fb.pop_back();
    return fb;
}

// Remove task outputs from <work_root>/task-<name>/, keeping only the per-task
// command.sh, command.sh.log.o, command.sh.log.e and command.sh.rc files.
// Only known task directory names are considered.
bool cleanWorkDir(const std::string& work_root) {
    namespace fs = std::filesystem;
    const std::vector<std::string> keep = {
        "command.sh", "command.sh.log.o", "command.sh.log.e", "command.sh.rc"
    };
    const std::vector<std::string> task_names = {
        "linear_align_extract", "vg_haplotype", "vg_giraffe", "merge_bams", "realign"
    };
    std::error_code ec;

    if (work_root.empty()) {
        SIMP_LOG_ERROR("refusing to clean: work root is empty");
        return false;
    }

    fs::path work = fs::canonical(fs::path(work_root), ec);
    if (ec || work.empty() || work == work.root_path()) {
        SIMP_LOG_ERROR("refusing to clean unsafe work root: %s", work_root.c_str());
        return false;
    }

    bool success = true;
    for (const auto& task_name : task_names) {
        fs::path task = work / ("task-" + task_name);
        fs::file_status task_status = fs::symlink_status(task, ec);
        if (ec) {
            SIMP_LOG_ERROR("cannot inspect task directory: %s", task.string().c_str());
            ec.clear();
            success = false;
            continue;
        }
        if (!fs::exists(task_status)) continue;
        if (fs::is_symlink(task_status) || !fs::is_directory(task_status)) {
            SIMP_LOG_WARN("skip unsafe task directory: %s", task.string().c_str());
            success = false;
            continue;
        }

        fs::path canonical_task = fs::canonical(task, ec);
        if (ec || canonical_task.parent_path() != work) {
            SIMP_LOG_ERROR("refusing to clean task directory outside work root: %s",
                           task.string().c_str());
            ec.clear();
            success = false;
            continue;
        }

        fs::directory_iterator entries(canonical_task, ec);
        if (ec) {
            SIMP_LOG_ERROR("cannot list task directory: %s",
                           canonical_task.string().c_str());
            ec.clear();
            success = false;
            continue;
        }
        for (const auto& entry : entries) {
            const std::string name = entry.path().filename().string();
            bool preserve = false;
            for (const auto& kept_name : keep) {
                if (name == kept_name) {
                    preserve = true;
                    break;
                }
            }
            if (preserve) continue;

            // remove_all removes a symlink itself and does not follow its target.
            fs::remove_all(entry.path(), ec);
            if (ec) {
                SIMP_LOG_ERROR("failed to remove work entry: %s (%s)",
                               entry.path().string().c_str(), ec.message().c_str());
                ec.clear();
                success = false;
            }
        }
    }
    return success;
}

struct Args {
    std::string fq1, fq2;
    std::string ref_fasta;
    std::string gbz, hapl, graph_ref_contigs;
    std::string extract_model;
    std::string sample_name = "SAMPLE";
    std::string platform = "DNBSEQ";
    int threads = 32;
    std::string tools_root;
    std::string java_home;
    std::string mode = "dcsmap";
    std::string out_bam;
    std::string work_dir;
    bool parallel = true;
    bool clean = true;
};

// Print one option line: "<opt><pad>desc" with the description column at 46.
// Continuation lines in `desc` are indented to column 46.
void printOpt(const std::string& opt, const std::vector<std::string>& desc) {
    const size_t kCol = 40;
    std::string first = desc.empty() ? std::string() : desc[0];
    std::string line = opt;
    if (line.size() < kCol) line.append(kCol - line.size(), ' ');
    else line += "  ";
    line += first;
    std::cerr << line << "\n";
    for (size_t i = 1; i < desc.size(); ++i) {
        std::cerr << std::string(kCol, ' ') << desc[i] << "\n";
    }
}

void usage() {
    std::cerr <<
"Program: dcsmap\n"
"version: " << kVersion << "\n"
"\n"
"Usage: dcsmap [-option]\n"
"\n"
"Required:\n";
    printOpt("  --fq1 <file>",                {"input read1 fastq path"});
    printOpt("  --fq2 <file>",                {"input read2 fastq path"});
    printOpt("  --out-bam <file>",            {"output bam path (work dir defaults to its dirname)"});
    printOpt("  --ref-fasta <file>",          {"reference fasta path (with .fai + .dict + aligner index)"});
    printOpt("  --gbz <file>",                {"vg gbz graph path"});
    printOpt("  --hapl <file>",               {"vg hapl index path"});
    printOpt("  --graph-ref-contigs <file>",  {"vg ref-paths file"});
    printOpt("  --extract-model <file>",      {"model file for extract-bam"});
    std::cerr << "\nOptions:\n";
    printOpt("  --sample-name <str>",         {"sample name (default: SAMPLE)"});
    printOpt("  --platform <str>",            {"sequencing platform (default: DNBSEQ)"});
    printOpt("  --threads <int>",             {"number of threads to use (default: 32)"});
    printOpt("  --mode <str>",                {"workflow mode (default: dcsmap)",
                                              "available options: {dcsmap, dcsmap-m1}"});
    printOpt("  --tools-root <dir>",          {"tools root dir (libexec + jar)",
                                              "default: DCSTOOLS_HOME env or parent-of-exe-dir"});
    printOpt("  --java-home <dir>",           {"JAVA home (must be Java 8)",
                                              "default: JAVA_HOME env"});
    printOpt("  --work-dir <dir>",            {"work dir (default: out-bam-dir/work.XXXXXX)"});
    printOpt("  --parallel <bool>",           {"run linear_align_extract and vg_haplotype in parallel (default: true)",
                                              "set to false to run vg_haplotype first, then linear_align_extract"});
    printOpt("  --clean <bool>",              {"clean work dir on success, keeping command.sh/logs/rc (default: true)"});
    printOpt("  -h, --help",                  {"display help message"});
    printOpt("  --version",                   {"display version message"});
}

bool parseBool(const std::string& v, bool& out, std::string& err, const char* key) {
    std::string s;
    s.reserve(v.size());
    for (char c : v) s.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    if (s == "true")  { out = true;  return true; }
    if (s == "false") { out = false; return true; }
    err = std::string("invalid value for ") + key + " (expected true|false): " + v;
    return false;
}

bool parseArgs(int argc, char** argv, Args& a, std::string& err) {
    auto need = [&](int i, const char* key) -> bool {
        if (i + 1 >= argc) { err = std::string("missing value for ") + key; return false; }
        return true;
    };
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto next = [&](std::string& dst) { dst = argv[++i]; };
        if (k == "--fq1") { if (!need(i,k.c_str())) return false; next(a.fq1); }
        else if (k == "--fq2") { if (!need(i,k.c_str())) return false; next(a.fq2); }
        else if (k == "--ref-fasta") { if (!need(i,k.c_str())) return false; next(a.ref_fasta); }
        else if (k == "--gbz") { if (!need(i,k.c_str())) return false; next(a.gbz); }
        else if (k == "--hapl") { if (!need(i,k.c_str())) return false; next(a.hapl); }
        else if (k == "--graph-ref-contigs") { if (!need(i,k.c_str())) return false; next(a.graph_ref_contigs); }
        else if (k == "--extract-model") { if (!need(i,k.c_str())) return false; next(a.extract_model); }
        else if (k == "--sample-name") { if (!need(i,k.c_str())) return false; next(a.sample_name); }
        else if (k == "--platform") { if (!need(i,k.c_str())) return false; next(a.platform); }
        else if (k == "--threads") { if (!need(i,k.c_str())) return false; a.threads = std::stoi(argv[++i]); }
        else if (k == "--tools-root") { if (!need(i,k.c_str())) return false; next(a.tools_root); }
        else if (k == "--java-home") { if (!need(i,k.c_str())) return false; next(a.java_home); }
        else if (k == "--mode") { if (!need(i,k.c_str())) return false; next(a.mode); }
        else if (k == "--out-bam") { if (!need(i,k.c_str())) return false; next(a.out_bam); }
        else if (k == "--work-dir") { if (!need(i,k.c_str())) return false; next(a.work_dir); }
        else if (k == "--parallel") { if (!need(i,k.c_str())) return false; if (!parseBool(argv[++i], a.parallel, err, k.c_str())) return false; }
        else if (k == "--clean") { if (!need(i,k.c_str())) return false; if (!parseBool(argv[++i], a.clean, err, k.c_str())) return false; }
        else if (k == "--help" || k == "-h") { usage(); std::exit(0); }
        else if (k == "--version") { std::cerr << "version: " << kVersion << "\n"; std::exit(0); }
        else { err = "unknown option: " + k; return false; }
    }
    return true;
}

std::string resolveToolsRoot(const Args& a) {
    if (!a.tools_root.empty()) return absolute(a.tools_root);
    const char* env = std::getenv("DCSTOOLS_HOME");
    if (env && env[0]) return absolute(env);
    char exe[PATH_MAX];
    ssize_t n = readlink("/proc/self/exe", exe, sizeof(exe) - 1);
    if (n > 0) { exe[n] = '\0'; return dirName(dirName(canonical(exe))); }
    return std::string("./tools");
}

bool checkJava8(const std::string& java_home, std::string& err) {
    if (java_home.empty()) { err = "java_home is empty (set --java-home or JAVA_HOME)"; return false; }
    std::string javaBin = joinPath(joinPath(java_home, "bin"), "java");
    std::string out;
    std::ostringstream cmd;
    cmd << "'" << javaBin << "' -version 2>&1";
    int rc = runCapture(cmd.str(), out);
    if (rc != 0) { err = "failed to run " + javaBin + " -version (exit " + std::to_string(rc) + ")"; return false; }
    if (out.find("1.8.") != std::string::npos || out.find("\"1.8") != std::string::npos) return true;
    err = "java at " + java_home + " is not Java 8. Output:\n" + out;
    return false;
}

std::string resolveJavaHome(const Args& a, std::string& err) {
    std::string jh = a.java_home;
    if (jh.empty()) {
        const char* env = std::getenv("JAVA_HOME");
        if (env && env[0]) jh = env;
    }
    if (jh.empty()) { err = "JAVA_HOME not set and --java-home not given"; return jh; }
    jh = absolute(jh);
    if (!checkJava8(jh, err)) return std::string();
    return jh;
}

void deriveRefPrefix(const std::string& graph_ref_contigs,
                     std::string& set_reference,
                     std::string& ref_path_prefix) {
    std::string line = firstLine(graph_ref_contigs);
    std::string derived_ref;
    size_t pos = line.find("#0#");
    if (pos != std::string::npos) derived_ref = line.substr(0, pos);
    set_reference = derived_ref;
    ref_path_prefix = set_reference.empty() ? std::string() : (set_reference + "#0#");
}

}  // namespace

// Run a single task synchronously (start + join). Returns true on success.
// On failure the Task itself logs the error and dumps its log.e tail.
bool runSync(Task& t) {
    t.start();
    t.join();
    return t.result().exit_code == 0;
}

int orchestrate(Config& c, const Args& a) {
    SIMP_LOG_INFO("==== workflow start (mode=%s parallel=%s threads=%d) ====",
                  a.mode.c_str(), a.parallel ? "true" : "false", c.threads);

    // ---- T1: linear_align_extract ; T2: vg_haplotype ----
    std::vector<std::string> t1_out, t2_out;
    Task t1("linear_align_extract", c.work_root, buildLinearAlignExtract(c, t1_out), t1_out);
    Task t2("vg_haplotype", c.work_root, buildVgHaplotype(c, t2_out), t2_out);

    if (a.parallel) {
        t1.start();
        t2.start();
        t1.join();
        t2.join();
        if (t1.result().exit_code != 0) return 1;
        if (t2.result().exit_code != 0) return 1;
    } else {
        if (!runSync(t2)) return 1;  // vg_haplotype first
        if (!runSync(t1)) return 1;  // then linear_align_extract
    }

    // T1 outputs
    std::string extract_fq1 = t1_out[0];
    std::string extract_fq2 = t1_out[1];
    std::string linear_bam = t1_out[2];
    // T2 outputs
    std::string hap_gbz = t2_out[0];
    std::string hap_dist = t2_out[1];
    std::string hap_min = t2_out[2];
    std::string hap_zipcodes = t2_out[3];

    // ---- T3: vg_giraffe ----
    std::vector<std::string> t3_out;
    Task t3("vg_giraffe", c.work_root,
            buildVgGiraffe(c, extract_fq1, extract_fq2, hap_gbz, hap_dist, hap_min, hap_zipcodes, t3_out),
            t3_out);
    if (!runSync(t3)) return 1;
    std::string giraffe_bam = t3_out[0];
    std::string giraffe_bai = t3_out[1];

    std::string final_bam = c.out_bam;
    if (a.mode == "dcsmap") {
        // realign (intermediate) on giraffe bam -> abra2 bam in job dir
        std::vector<std::string> ra_out;
        Task ra("realign", c.work_root,
                buildRealign(c, giraffe_bam, giraffe_bai, "abra2", "", ra_out), ra_out);
        if (!runSync(ra)) return 1;
        std::string abra2_bam = ra_out[0];
        std::string abra2_bai = ra_out[1];

        // merge_bams (FINAL) -> out-bam
        std::vector<std::string> mg_out;
        Task mg("merge_bams", c.work_root,
                buildMergeBams(c, linear_bam, abra2_bam, final_bam, mg_out), mg_out);
        if (!runSync(mg)) return 1;
    } else {
        // dcsmap-m1: merge_bams (intermediate) on linear + giraffe bam
        std::vector<std::string> mg_out;
        Task mg("merge_bams", c.work_root,
                buildMergeBams(c, linear_bam, giraffe_bam, "", mg_out), mg_out);
        if (!runSync(mg)) return 1;
        std::string merged_bam = mg_out[0];
        std::string merged_bai = mg_out[1];

        // realign (FINAL) on merged bam -> out-bam
        std::vector<std::string> ra_out;
        Task ra("realign", c.work_root,
                buildRealign(c, merged_bam, merged_bai, "abra2", final_bam, ra_out), ra_out);
        if (!runSync(ra)) return 1;
    }

    SIMP_LOG_INFO("DONE. final bam: %s", final_bam.c_str());

    if (a.clean) {
        SIMP_LOG_INFO("cleaning work dir (keeping command.sh/logs/rc): %s", c.work_root.c_str());
        if (!cleanWorkDir(c.work_root)) {
            SIMP_LOG_WARN("work directory cleanup was incomplete: %s", c.work_root.c_str());
        }
    }
    return 0;
}

int main(int argc, char** argv) {
    Args a;
    std::string err;
    if (!parseArgs(argc, argv, a, err)) { SIMP_LOG_ERROR("%s", err.c_str()); usage(); return 1; }

    auto require = [&](const std::string& v, const char* name) {
        if (v.empty()) { SIMP_LOG_ERROR("missing required --%s", name); usage(); std::exit(1); }
    };
    require(a.fq1, "fq1");
    require(a.fq2, "fq2");
    require(a.out_bam, "out-bam");
    require(a.ref_fasta, "ref-fasta");
    require(a.gbz, "gbz");
    require(a.hapl, "hapl");
    require(a.graph_ref_contigs, "graph-ref-contigs");
    require(a.extract_model, "extract-model");

    if (a.mode != "dcsmap" && a.mode != "dcsmap-m1") {
        SIMP_LOG_ERROR("--mode must be dcsmap or dcsmap-m1 (got: %s)", a.mode.c_str());
        usage();
        return 1;
    }

    std::string tools_root = resolveToolsRoot(a);
    if (!makeDirP(tools_root)) { SIMP_LOG_ERROR("tools_root not found: %s", tools_root.c_str()); usage(); return 1; }

    std::string java_home = resolveJavaHome(a, err);
    if (java_home.empty()) { SIMP_LOG_ERROR("%s", err.c_str()); usage(); return 1; }

    Config c;
    c.fq1 = absolute(a.fq1);
    c.fq2 = absolute(a.fq2);
    c.ref_fasta = absoluteNoSymlink(a.ref_fasta);
    c.ref_fasta_fai = c.ref_fasta + ".fai";
    c.ref_fasta_dict = c.ref_fasta + ".dict";
    c.gbz = absolute(a.gbz);
    c.hapl = absolute(a.hapl);
    c.graph_ref_contigs = absolute(a.graph_ref_contigs);
    c.extract_model = absolute(a.extract_model);
    c.sample_name = a.sample_name;
    c.platform = a.platform;
    c.threads = a.threads;
    c.samtools_sort_threads = a.parallel ? 4 : a.threads;
    c.vg_haplotype_threads = a.parallel ? std::max(a.threads / 2, 4) : a.threads;
    c.tools_root = tools_root;
    c.java_home = java_home;
    c.out_bam = absolute(a.out_bam);

    deriveRefPrefix(c.graph_ref_contigs, c.set_reference, c.ref_path_prefix);
    SIMP_LOG_INFO("set_reference='%s' ref_path_prefix='%s'",
                  c.set_reference.c_str(), c.ref_path_prefix.c_str());

    auto checkFile = [&](const std::string& p, const char* what) {
        if (!fileExists(p)) { SIMP_LOG_ERROR("%s not found: %s", what, p.c_str()); usage(); std::exit(1); }
    };
    checkFile(c.fq1, "fq1");
    checkFile(c.fq2, "fq2");
    checkFile(c.ref_fasta, "ref-fasta");
    checkFile(c.gbz, "gbz");
    checkFile(c.hapl, "hapl");
    checkFile(c.graph_ref_contigs, "graph-ref-contigs");
    checkFile(c.extract_model, "extract-model");

    makeDirP(c.out_bam.empty() ? std::string(".") : dirName(c.out_bam));
    if (!a.work_dir.empty()) {
        c.work_root = absolute(a.work_dir);
        makeDirP(c.work_root);
    } else {
        c.work_root = mktempWork(dirName(c.out_bam));
        if (c.work_root.empty()) { SIMP_LOG_ERROR("failed to create work dir under %s", dirName(c.out_bam).c_str()); usage(); return 1; }
    }
    SIMP_LOG_INFO("work_root=%s", c.work_root.c_str());
    SIMP_LOG_INFO("out_bam=%s", c.out_bam.c_str());

    int rc = orchestrate(c, a);
    SIMP_LOG_INFO("==== workflow end, return=%d ====", rc);
    return rc;
}
