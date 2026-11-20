#ifndef RABBITSKETCH_NUMA_MEMORY_H
#define RABBITSKETCH_NUMA_MEMORY_H

#include "api/OpenMPCompat.h"
#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstring>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#ifdef __linux__
#include <linux/mempolicy.h>
#include <sched.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif

namespace Sketch {
struct NumaIndexOptions {
    bool replicate = false;
    // Sum of the node-local index mappings, not whole-process memory.
    size_t max_replica_bytes = size_t(64) << 30;
};

namespace Numa {
inline std::vector<int> parseList(const std::string& text) {
    std::vector<int> result;
    std::istringstream input(text);
    std::string part;
    while (std::getline(input, part, ',')) {
        if (part.empty()) continue;
        const auto dash = part.find('-');
        int lo = std::stoi(part), hi = dash == std::string::npos ? lo : std::stoi(part.substr(dash + 1));
        if (lo < 0 || hi < lo || hi > 1048575) throw std::invalid_argument("invalid CPU/node list");
        for (int i = lo; i <= hi; ++i) result.push_back(i);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}
inline std::string readFile(const std::string& path) {
    std::ifstream stream(path);
    return std::string(std::istreambuf_iterator<char>(stream), {});
}

struct Mask {
    std::vector<unsigned long> words;
    static constexpr size_t width = sizeof(unsigned long) * 8;
    explicit Mask(size_t bits = 1024) : words((bits + width - 1) / width, 0) {}
    size_t bits() const { return words.size() * width; }
    size_t bytes() const { return words.size() * sizeof(unsigned long); }
    void add(size_t bit) { if (bit >= bits()) words.resize(bit / width + 1); words[bit / width] |= 1UL << (bit % width); }
    bool has(size_t bit) const { return bit < bits() && (words[bit / width] & (1UL << (bit % width))); }
    bool empty() const { for (auto word : words) if (word) return false; return true; }
    Mask intersection(const Mask& other) const {
        Mask out(bits());
        for (size_t i = 0; i < std::min(words.size(), other.words.size()); ++i) out.words[i] = words[i] & other.words[i];
        return out;
    }
    static Mask affinity() {
#ifdef __linux__
        for (size_t bits = 1024; bits <= 1048576; bits *= 2) {
            Mask out(bits);
            if (syscall(SYS_sched_getaffinity, 0, out.bytes(), out.words.data()) >= 0) return out;
            if (errno != EINVAL) break;
        }
#endif
        throw std::runtime_error("cannot read thread CPU affinity");
    }
    bool applyAffinity() const noexcept {
#ifdef __linux__
        return syscall(SYS_sched_setaffinity, 0, bytes(), words.data()) == 0;
#else
        return false;
#endif
    }
};

struct Node { int id; Mask cpus; };
struct Topology {
    std::vector<Node> nodes;
    size_t policy_bits = 1024;
    static Topology detect(int threads) {
        Topology out;
#ifdef __linux__
        const auto online = parseList(readFile("/sys/devices/system/node/online"));
        auto memory = parseList(readFile("/sys/devices/system/node/has_memory"));
        auto possible = parseList(readFile("/sys/devices/system/node/possible"));
        if (!possible.empty()) out.policy_bits = std::max<size_t>(1024, possible.back() + 1);
        std::string status = readFile("/proc/self/status");
        const auto start = status.find("Mems_allowed_list:");
        if (start != std::string::npos) {
            const auto end = status.find('\n', start);
            auto allowed = parseList(status.substr(start + 18, end - start - 18));
            std::vector<int> common;
            std::set_intersection(memory.begin(), memory.end(), allowed.begin(), allowed.end(), std::back_inserter(common));
            memory.swap(common);
        }
        std::vector<Node> all;
        for (int id : online) if (std::binary_search(memory.begin(), memory.end(), id)) {
            Mask cpus;
            for (int cpu : parseList(readFile("/sys/devices/system/node/node" + std::to_string(id) + "/cpulist"))) cpus.add(cpu);
            if (!cpus.empty()) all.push_back({id, std::move(cpus)});
        }
        std::vector<unsigned char> used(all.size(), 0);
        // Query the actual OpenMP team: its workers can have narrower affinity
        // masks than the process launcher, or fewer workers than requested.
        #pragma omp parallel num_threads(std::max(1, threads))
        {
            try {
                Mask affinity = Mask::affinity();
                std::vector<size_t> eligible;
                for (size_t i = 0; i < all.size(); ++i)
                    if (!affinity.intersection(all[i].cpus).empty()) eligible.push_back(i);
                if (!eligible.empty()) {
                    const size_t selected = eligible[size_t(omp_get_thread_num()) % eligible.size()];
                    #pragma omp critical(rabbitsketch_numa_discovery)
                    { used[selected] = 1; }
                }
            } catch (...) {}
        }
        for (size_t i = 0; i < all.size(); ++i) if (used[i]) out.nodes.push_back(std::move(all[i]));
#else
        (void)threads;
#endif
        return out;
    }
};

// Restrict, never expand, the caller's affinity. Restore both affinity and
// per-thread memory policy when leaving the query region.
class ThreadBinding {
    Mask previous_cpu_, previous_memory_;
    int previous_policy_ = 0;
    bool cpu_changed_ = false, policy_changed_ = false;
public:
    bool active = false;
    ThreadBinding(const Node& node, size_t policy_bits)
        : previous_cpu_(Mask::affinity()), previous_memory_(policy_bits) {
#ifdef __linux__
        const Mask local = previous_cpu_.intersection(node.cpus);
        Mask target(policy_bits); target.add(node.id);
        if (local.empty() || !local.applyAffinity()) return;
        cpu_changed_ = true;
        if (syscall(SYS_get_mempolicy, &previous_policy_, previous_memory_.words.data(),
                    previous_memory_.bits(), nullptr, 0) != 0) return;
        if (syscall(SYS_set_mempolicy, MPOL_BIND | MPOL_F_STATIC_NODES,
                    target.words.data(), target.bits()) != 0) return;
        policy_changed_ = true;
        active = true;
#else
        (void)node;
#endif
    }
    ~ThreadBinding() {
#ifdef __linux__
        if (policy_changed_)
            syscall(SYS_set_mempolicy, previous_policy_, previous_policy_ == MPOL_DEFAULT ? nullptr : previous_memory_.words.data(),
                    previous_policy_ == MPOL_DEFAULT ? 0 : previous_memory_.bits());
        if (cpu_changed_) previous_cpu_.applyAffinity();
#endif
    }
    ThreadBinding(const ThreadBinding&) = delete;
    ThreadBinding& operator=(const ThreadBinding&) = delete;
};

struct Budget {
    size_t limit, mapped = 0;
    explicit Budget(size_t bytes) : limit(bytes) {}
};

// Dedicated mmap regions prevent allocator reuse from leaving a nominally
// local copy in old interleaved heap pages. mbind applies before first touch.
class Arena {
    struct Block { uint8_t* data; size_t bytes, used; };
    std::vector<Block> blocks_;
    std::shared_ptr<Budget> budget_;
    int node_;
    size_t policy_bits_;
public:
    Arena(int node, size_t policy_bits, std::shared_ptr<Budget> budget)
        : budget_(std::move(budget)), node_(node), policy_bits_(policy_bits) {}
    ~Arena() {
#ifdef __linux__
        for (const auto& block : blocks_) { munmap(block.data, block.bytes); budget_->mapped -= block.bytes; }
#endif
    }
    void* allocate(size_t bytes, size_t alignment) {
        if (!bytes) return nullptr;
        if (!alignment || (alignment & (alignment - 1))) throw std::bad_alloc();
        if (!blocks_.empty()) {
            auto& block = blocks_.back();
            if (block.used <= SIZE_MAX - (alignment - 1)) {
                const size_t begin = (block.used + alignment - 1) & ~(alignment - 1);
                if (begin <= block.bytes && bytes <= block.bytes - begin) {
                    block.used = begin + bytes; return block.data + begin;
                }
            }
        }
#ifdef __linux__
        const size_t page = static_cast<size_t>(sysconf(_SC_PAGESIZE));
        if (alignment > page || bytes > SIZE_MAX - page) throw std::bad_alloc();
        size_t chunk = std::max<size_t>(65536, std::min<size_t>(size_t(64) << 20,
            blocks_.empty() ? 65536 : blocks_.back().bytes * (blocks_.back().bytes < (size_t(32) << 20) ? 2 : 1)));
        chunk = std::max(chunk, (bytes + page - 1) / page * page);
        if (budget_->mapped > budget_->limit || chunk > budget_->limit - budget_->mapped)
            throw std::length_error("NUMA index replica memory budget exceeded");
        void* region = mmap(nullptr, chunk, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (region == MAP_FAILED) throw std::bad_alloc();
        Mask target(policy_bits_); target.add(node_);
        if (syscall(SYS_mbind, region, chunk, MPOL_BIND | MPOL_F_STATIC_NODES,
                    target.words.data(), target.bits(), 0) != 0) {
            const int error = errno; munmap(region, chunk);
            throw std::runtime_error("NUMA index mbind failed: " + std::string(std::strerror(error)));
        }
        try { blocks_.push_back({static_cast<uint8_t*>(region), chunk, bytes}); }
        catch (...) { munmap(region, chunk); throw; }
        budget_->mapped += chunk;
        return region;
#else
        throw std::runtime_error("NUMA-local allocation requires Linux");
#endif
    }
    template<class T>
    T* copy(const T* values, size_t count) {
        if (!count) return nullptr;
        if (count > SIZE_MAX / sizeof(T)) throw std::bad_alloc();
        T* out = static_cast<T*>(allocate(count * sizeof(T), alignof(T)));
        std::copy_n(values, count, out);
        return out;
    }
    int node() const { return node_; }
    size_t mappedBytes() const { size_t total = 0; for (const auto& b : blocks_) total += b.bytes; return total; }
    // Diagnostic only, never in the query loop. Query the actual backing node
    // of touched addresses, in addition to the VMA allocation policy.
    bool placementMatches() const {
#ifdef __linux__
        for (const auto& block : blocks_) if (block.used) {
            for (size_t offset : {size_t(0), block.used / 2, block.used - 1}) {
                int node = -1;
                if (syscall(SYS_get_mempolicy, &node, nullptr, 0, block.data + offset,
                            MPOL_F_NODE | MPOL_F_ADDR) != 0 || node != node_) return false;
            }
        }
        return true;
#else
        return false;
#endif
    }
};

template<class T>
struct Span {
    T* pointer = nullptr;
    size_t length = 0;
    T* data() const { return pointer; }
    size_t size() const { return length; }
    bool empty() const { return !length; }
    T* begin() const { return pointer; }
    T* end() const { return length ? pointer + length : pointer; }
    T& operator[](size_t i) const { return pointer[i]; }
};

template<class T>
struct Allocator {
    using value_type = T;
    std::shared_ptr<Arena> arena;
    explicit Allocator(std::shared_ptr<Arena> owner = {}) : arena(std::move(owner)) {}
    template<class U> Allocator(const Allocator<U>& other) : arena(other.arena) {}
    T* allocate(size_t count) {
        if (!arena || count > SIZE_MAX / sizeof(T)) throw std::bad_alloc();
        return static_cast<T*>(arena->allocate(count * sizeof(T), alignof(T)));
    }
    void deallocate(T*, size_t) noexcept {}
    template<class U> bool operator==(const Allocator<U>& other) const { return arena == other.arena; }
    template<class U> bool operator!=(const Allocator<U>& other) const { return !(*this == other); }
};

}} // namespace Sketch::Numa
#endif
