// pcr.cpp
//
// In-silico PCR: find every primer-pair amplicon, on both strands, in a
// single pass per record. C++ port of numba_pcr.py, with both fixes that
// were verified against it kept in place:
//   - an 'N' in the SUBJECT sequence costs a real mismatch (mask = 0),
//     it never matches a primer base for free.
//   - only the NEAREST valid end is paired with each start, avoiding the
//     spurious "merged" amplicon bug found in the original pairing logic.
//
// Build:
//   g++ -O3 -std=c++17 -pthread -o pcr pcr.cpp
//
// Usage matches numba_pcr.py's CLI:
//   ./pcr -i refs.fasta -f FWDPRIMER -r REVPRIMER -m 2 \
//         --min-length 1 --max-length 1000 -t 8 -o amplicons.fasta

#include <algorithm>
#include <array>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <deque>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// ---------------------------------------------------------------------
// IUPAC bit masks: A=1, C=2, G=4, T=8. Combinations for degenerate codes.
// N is intentionally 0 -- see note above.
// ---------------------------------------------------------------------
static std::array<uint8_t, 256> build_mask_table() {
    std::array<uint8_t, 256> table{};
    table.fill(0);
    auto setmask = [&](const std::string &chars, uint8_t mask) {
        for (char c : chars) table[(unsigned char)c] = mask;
    };
    setmask("Aa", 1);
    setmask("Cc", 2);
    setmask("Gg", 4);
    setmask("TtUu", 8);
    setmask("Rr", 1 | 4);
    setmask("Yy", 2 | 8);
    setmask("Ss", 2 | 4);
    setmask("Ww", 1 | 8);
    setmask("Kk", 4 | 8);
    setmask("Mm", 1 | 2);
    setmask("Bb", 2 | 4 | 8);
    setmask("Dd", 1 | 4 | 8);
    setmask("Hh", 1 | 2 | 8);
    setmask("Vv", 1 | 2 | 4);
    // "Nn" left at 0 on purpose.
    return table;
}
static const std::array<uint8_t, 256> MASK_TABLE = build_mask_table();

// ---------------------------------------------------------------------
// Reverse complement, IUPAC-aware.
// ---------------------------------------------------------------------
static std::array<char, 256> build_comp_table() {
    std::array<char, 256> table{};
    for (int i = 0; i < 256; i++) table[i] = (char)i;
    auto setcomp = [&](char a, char b) {
        table[(unsigned char)a] = b;
        table[(unsigned char)tolower((unsigned char)a)] = (char)tolower((unsigned char)b);
    };
    setcomp('A', 'T'); setcomp('T', 'A'); setcomp('C', 'G'); setcomp('G', 'C');
    setcomp('U', 'A');
    setcomp('R', 'Y'); setcomp('Y', 'R'); setcomp('S', 'S'); setcomp('W', 'W');
    setcomp('K', 'M'); setcomp('M', 'K');
    setcomp('B', 'V'); setcomp('V', 'B'); setcomp('D', 'H'); setcomp('H', 'D');
    setcomp('N', 'N');
    return table;
}
static const std::array<char, 256> COMP_TABLE = build_comp_table();

static std::string reverse_complement(const std::string &seq) {
    std::string out(seq.size(), 'N');
    size_t n = seq.size();
    for (size_t i = 0; i < n; i++) out[i] = COMP_TABLE[(unsigned char)seq[n - 1 - i]];
    return out;
}

static std::string to_upper(const std::string &s) {
    std::string out = s;
    for (auto &c : out) c = (char)toupper((unsigned char)c);
    return out;
}

static std::vector<uint8_t> to_mask_array(const std::string &seq) {
    std::vector<uint8_t> masks(seq.size());
    for (size_t i = 0; i < seq.size(); i++) masks[i] = MASK_TABLE[(unsigned char)seq[i]];
    return masks;
}

// ---------------------------------------------------------------------
// Sliding-window search. Returns hit start positions (ascending, since we
// scan left to right) and the mismatch count at each.
// ---------------------------------------------------------------------
struct Hits {
    std::vector<int64_t> starts;
    std::vector<int16_t> errors;
};

static Hits find_hits(const std::vector<uint8_t> &seq_masks,
                       const std::vector<uint8_t> &primer_masks,
                       int max_mismatches) {
    Hits result;
    int64_t n = (int64_t)seq_masks.size();
    int64_t k = (int64_t)primer_masks.size();
    if (n < k) return result;

    for (int64_t i = 0; i <= n - k; i++) {
        int mismatches = 0;
        for (int64_t j = 0; j < k; j++) {
            if ((seq_masks[(size_t)(i + j)] & primer_masks[(size_t)j]) == 0) {
                mismatches++;
                if (mismatches > max_mismatches) break;
            }
        }
        if (mismatches <= max_mismatches) {
            result.starts.push_back(i);
            result.errors.push_back((int16_t)mismatches);
        }
    }
    return result;
}

// ---------------------------------------------------------------------
// FASTA record + streaming reader
// ---------------------------------------------------------------------
struct FastaRecord {
    std::string name;
    std::string seq;
};

class FastaReader {
public:
    explicit FastaReader(const std::string &path) : in(path) {
        if (!in) throw std::runtime_error("cannot open input file: " + path);
        std::string line;
        while (std::getline(in, line)) {
            strip_cr(line);
            if (!line.empty() && line[0] == '>') {
                next_header = line;
                break;
            }
        }
    }

    bool next(FastaRecord &out) {
        if (next_header.empty()) return false;
        std::string header = next_header.substr(1);
        size_t sp = header.find_first_of(" \t");
        out.name = (sp == std::string::npos) ? header : header.substr(0, sp);
        out.seq.clear();
        next_header.clear();

        std::string line;
        while (std::getline(in, line)) {
            strip_cr(line);
            if (line.empty()) continue;
            if (line[0] == '>') {
                next_header = line;
                return true;
            }
            out.seq += line;
        }
        return true;  // last record; next_header stays empty
    }

private:
    std::ifstream in;
    std::string next_header;
    static void strip_cr(std::string &s) {
        while (!s.empty() && (s.back() == '\r' || s.back() == '\n')) s.pop_back();
    }
};

// ---------------------------------------------------------------------
// Core algorithm: find every amplicon in one sequence, both strands.
// ---------------------------------------------------------------------
struct Amplicon {
    std::string header;
    std::string seq;
};

static std::vector<Amplicon> extract_from_sequence(
    const std::string &seq_id,
    const std::string &seq_upper,
    const std::vector<uint8_t> &fwd_masks,
    const std::vector<uint8_t> &rev_rc_masks,
    int64_t fwd_len,
    int64_t rev_len,
    int mismatches,
    int64_t min_len,
    int64_t max_len,
    bool include_primers) {
    std::vector<Amplicon> results;

    std::string rc_seq = reverse_complement(seq_upper);
    const std::pair<const char *, const std::string *> directions[2] = {
        {"forward", &seq_upper}, {"reverse", &rc_seq}};

    for (auto &dir_pair : directions) {
        const char *direction = dir_pair.first;
        const std::string &strand_seq = *dir_pair.second;

        std::vector<uint8_t> strand_masks = to_mask_array(strand_seq);

        Hits f_hits = find_hits(strand_masks, fwd_masks, mismatches);
        Hits r_hits = find_hits(strand_masks, rev_rc_masks, mismatches);

        if (f_hits.starts.empty() || r_hits.starts.empty()) continue;

        for (size_t f_idx = 0; f_idx < f_hits.starts.size(); f_idx++) {
            int64_t f_start = f_hits.starts[f_idx];
            int64_t f_end = f_start + fwd_len;

            int64_t min_r_start = f_end + min_len;
            int64_t max_r_start = f_end + max_len;

            // nearest valid end only -- mirrors the fixed Python version
            auto lo_it = std::lower_bound(r_hits.starts.begin(), r_hits.starts.end(), min_r_start);
            if (lo_it == r_hits.starts.end()) continue;
            int64_t r_start = *lo_it;
            if (r_start > max_r_start) continue;

            size_t r_idx = (size_t)(lo_it - r_hits.starts.begin());
            int64_t r_end = r_start + rev_len;

            int64_t amp_start = include_primers ? f_start : f_end;
            int64_t amp_end = include_primers ? r_end : r_start;

            std::string amplicon = strand_seq.substr((size_t)amp_start, (size_t)(amp_end - amp_start));

            std::ostringstream header;
            header << seq_id << "_sub[" << (amp_start + 1) << ".." << amp_end << "]"
                   << " direction=" << direction
                   << " forward_error=" << f_hits.errors[f_idx]
                   << " reverse_error=" << r_hits.errors[r_idx]
                   << " forward_match=" << strand_seq.substr((size_t)f_start, (size_t)fwd_len)
                   << " reverse_match=" << strand_seq.substr((size_t)r_start, (size_t)rev_len);

            results.push_back({header.str(), amplicon});
        }
    }
    return results;
}

// ---------------------------------------------------------------------
// Bounded thread-safe queue for streaming records to worker threads
// without loading the whole input file into memory at once.
// ---------------------------------------------------------------------
class RecordQueue {
public:
    explicit RecordQueue(size_t max_size) : max_size_(max_size) {}

    void push(FastaRecord rec) {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_not_full_.wait(lock, [&] { return queue_.size() < max_size_ || done_pushing_; });
        queue_.push_back(std::move(rec));
        lock.unlock();
        cv_not_empty_.notify_one();
    }

    bool pop(FastaRecord &out) {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_not_empty_.wait(lock, [&] { return !queue_.empty() || done_pushing_; });
        if (queue_.empty()) return false;
        out = std::move(queue_.front());
        queue_.pop_front();
        lock.unlock();
        cv_not_full_.notify_one();
        return true;
    }

    void close() {
        std::unique_lock<std::mutex> lock(mtx_);
        done_pushing_ = true;
        lock.unlock();
        cv_not_empty_.notify_all();
        cv_not_full_.notify_all();
    }

private:
    std::deque<FastaRecord> queue_;
    size_t max_size_;
    bool done_pushing_ = false;
    std::mutex mtx_;
    std::condition_variable cv_not_full_;
    std::condition_variable cv_not_empty_;
};

// ---------------------------------------------------------------------
// CLI argument parsing (minimal, no external dependency)
// ---------------------------------------------------------------------
struct Args {
    std::string input;
    std::string forward;
    std::string reverse_primer;
    int mismatches = 0;
    int64_t min_length = 0;
    int64_t max_length = -1;
    bool include_primers = false;
    int width = 0;
    std::string output;
    int threads = 1;
};

static void print_usage_and_exit() {
    std::cerr << "Usage: pcr -i INPUT.fasta -f FWD_PRIMER -r REV_PRIMER -m MISMATCHES "
                 "--max-length N [--min-length N] [-o OUTPUT.fasta] [-t THREADS] "
                 "[--include-primers] [-w WIDTH]\n";
    std::exit(1);
}

static Args parse_args(int argc, char **argv) {
    Args args;
    bool have_max_length = false;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        auto next_val = [&]() -> std::string {
            if (i + 1 >= argc) print_usage_and_exit();
            return argv[++i];
        };
        if (a == "-i" || a == "--input") args.input = next_val();
        else if (a == "-f" || a == "--forward") args.forward = next_val();
        else if (a == "-r" || a == "--reverse") args.reverse_primer = next_val();
        else if (a == "-m" || a == "--mismatches") args.mismatches = std::stoi(next_val());
        else if (a == "--min-length") args.min_length = std::stoll(next_val());
        else if (a == "--max-length") { args.max_length = std::stoll(next_val()); have_max_length = true; }
        else if (a == "--include-primers") args.include_primers = true;
        else if (a == "-w" || a == "--width") args.width = std::stoi(next_val());
        else if (a == "-o" || a == "--output") args.output = next_val();
        else if (a == "-t" || a == "--threads") args.threads = std::stoi(next_val());
        else if (a == "-h" || a == "--help") print_usage_and_exit();
        else { std::cerr << "Unknown argument: " << a << "\n"; print_usage_and_exit(); }
    }
    if (args.input.empty() || args.forward.empty() || args.reverse_primer.empty() || !have_max_length) {
        print_usage_and_exit();
    }
    return args;
}

static void write_fasta(std::ostream &out, const std::string &header, const std::string &seq, int width) {
    out << ">" << header << "\n";
    if (width <= 0) {
        out << seq << "\n";
    } else {
        for (size_t i = 0; i < seq.size(); i += (size_t)width) {
            out << seq.substr(i, (size_t)width) << "\n";
        }
    }
}

// ---------------------------------------------------------------------
// main
// ---------------------------------------------------------------------
int main(int argc, char **argv) {
    Args args = parse_args(argc, argv);

    std::string fwd = to_upper(args.forward);
    std::string rev = to_upper(args.reverse_primer);
    std::string rev_rc = reverse_complement(rev);

    std::vector<uint8_t> fwd_masks = to_mask_array(fwd);
    std::vector<uint8_t> rev_rc_masks = to_mask_array(rev_rc);
    int64_t fwd_len = (int64_t)fwd.size();
    int64_t rev_len = (int64_t)rev_rc.size();

    std::ostream *out_stream = &std::cout;
    std::ofstream out_file;
    if (!args.output.empty()) {
        out_file.open(args.output);
        if (!out_file) {
            std::cerr << "cannot open output file: " << args.output << "\n";
            return 1;
        }
        out_stream = &out_file;
    }

    std::mutex output_mutex;
    int threads = std::max(1, args.threads);

    if (threads == 1) {
        // single-threaded path -- no queue/thread overhead
        FastaReader reader(args.input);
        FastaRecord rec;
        while (reader.next(rec)) {
            std::string seq_upper = to_upper(rec.seq);
            auto amplicons = extract_from_sequence(
                rec.name, seq_upper, fwd_masks, rev_rc_masks, fwd_len, rev_len,
                args.mismatches, args.min_length, args.max_length, args.include_primers);
            for (auto &amp : amplicons) write_fasta(*out_stream, amp.header, amp.seq, args.width);
        }
    } else {
        RecordQueue queue((size_t)(threads * 4));
        std::atomic<bool> reading_done{false};

        std::thread reader_thread([&] {
            FastaReader reader(args.input);
            FastaRecord rec;
            while (reader.next(rec)) {
                queue.push(std::move(rec));
            }
            queue.close();
        });

        std::vector<std::thread> workers;
        for (int t = 0; t < threads; t++) {
            workers.emplace_back([&] {
                FastaRecord rec;
                while (queue.pop(rec)) {
                    std::string seq_upper = to_upper(rec.seq);
                    auto amplicons = extract_from_sequence(
                        rec.name, seq_upper, fwd_masks, rev_rc_masks, fwd_len, rev_len,
                        args.mismatches, args.min_length, args.max_length, args.include_primers);
                    if (!amplicons.empty()) {
                        std::lock_guard<std::mutex> lock(output_mutex);
                        for (auto &amp : amplicons) write_fasta(*out_stream, amp.header, amp.seq, args.width);
                    }
                }
            });
        }

        reader_thread.join();
        for (auto &w : workers) w.join();
    }

    if (out_file.is_open()) out_file.close();
    return 0;
}