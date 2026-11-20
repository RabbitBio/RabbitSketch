#ifndef RABBITSKETCH_NUMA_INDEX_H
#define RABBITSKETCH_NUMA_INDEX_H

#include "NumaMemory.h"
#include "PostingCodec.h"
#include "phmap.h"
#include <iostream>
#include <chrono>

namespace Sketch {

// Immutable hash-index clone. Its hash table and encoded lists both use
// dedicated mappings bound to the same NUMA node.
template<class Key, class Range>
struct NumaHashIndexReplica {
    using StandardMap = phmap::flat_hash_map<Key, Range>;
    using Map = phmap::flat_hash_map<Key, Range, typename StandardMap::hasher,
        typename StandardMap::key_equal, Numa::Allocator<typename StandardMap::value_type>>;
    std::shared_ptr<Numa::Arena> arena;
    Map postIdx;
    Numa::Span<uint8_t> encodedPosts;
    template<class Source>
    NumaHashIndexReplica(const Source& source, std::shared_ptr<Numa::Arena> owner)
        : arena(std::move(owner)), postIdx(0, source.postIdx.hash_function(), source.postIdx.key_eq(),
            Numa::Allocator<typename StandardMap::value_type>(arena)) {
        postIdx.reserve(source.postIdx.size());
        postIdx.insert(source.postIdx.begin(), source.postIdx.end());
        encodedPosts = {arena->copy(source.encodedPosts.data(), source.encodedPosts.size()), source.encodedPosts.size()};
    }
    bool findPostingRange(const Key& key, size_t& off, uint32_t& count) const {
        const auto it = postIdx.find(key);
        if (it == postIdx.end()) return false;
        off = it->second.off; count = it->second.cnt; return true;
    }
    template<class Visit>
    void forEachPosting(size_t off, uint32_t count, Visit visit) const {
        if (count) PostingCodec::forEachCount(encodedPosts.data() + off, count, visit);
    }
};

template<class Index>
class NumaIndex {
    using Replica = typename Index::NumaReplica;
    std::unique_ptr<Index> shared_;
    Numa::Topology topology_;
    std::vector<std::unique_ptr<Replica>> replicas_;
    std::shared_ptr<Numa::Budget> budget_;
    size_t keys_ = 0;
    mutable std::atomic<size_t> binding_failures_{0};
public:
    size_t totalPostings = 0;
    class Local {
        const Index* shared_ = nullptr;
        const Replica* replica_ = nullptr;
        std::unique_ptr<Numa::ThreadBinding> binding_;
        friend class NumaIndex;
    public:
        int node = -1;
        bool locally_bound = false;
        Local() = default;
        Local(Local&&) = default;
        template<class Key>
        bool findPostingRange(const Key& key, size_t& off, uint32_t& count) const {
            return replica_ ? replica_->findPostingRange(key, off, count)
                            : shared_->findPostingRange(key, off, count);
        }
        template<class Visit>
        void forEachPosting(size_t off, uint32_t count, Visit visit) const {
            if (replica_) replica_->forEachPosting(off, count, visit);
            else shared_->forEachPosting(off, count, visit);
        }
    };
    NumaIndex(Index&& source, int threads, NumaIndexOptions options = {})
        : shared_(new Index(std::move(source))) {
        keys_ = shared_->numaKeyCount();
        totalPostings = shared_->totalPostings;
        if (!options.replicate || !totalPostings) return;
        const auto started = std::chrono::steady_clock::now();
        topology_ = Numa::Topology::detect(threads);
        if (topology_.nodes.size() < 2) {
            std::cerr << "NUMA index: shared (fewer than two usable query nodes)\n";
            return;
        }
        const size_t estimate = shared_->numaStorageBytes();
        if (estimate > options.max_replica_bytes / topology_.nodes.size())
            throw std::length_error("NUMA index replicas exceed --numa-index-max-gib");
        budget_ = std::make_shared<Numa::Budget>(options.max_replica_bytes);
        replicas_.reserve(topology_.nodes.size());
        for (const auto& node : topology_.nodes) {
            auto arena = std::make_shared<Numa::Arena>(node.id, topology_.policy_bits, budget_);
            if (replicas_.empty()) {
                replicas_.emplace_back(new Replica(*shared_, arena));
                // Release the interleaved index before creating further copies:
                // never retain an unnecessary original plus all node replicas.
                shared_.reset();
            } else replicas_.emplace_back(new Replica(*replicas_.front(), arena));
        }
        std::cerr << "NUMA index: " << replicas_.size() << " replicas, nodes=";
        for (size_t i = 0; i < topology_.nodes.size(); ++i)
            std::cerr << (i ? "," : "") << topology_.nodes[i].id;
        std::cerr << ", mapped_bytes=" << budget_->mapped
                  << ", copy_seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count()
                  << ", original index released, local query binding enabled\n";
    }
    ~NumaIndex() {
        if (binding_failures_.load())
            std::cerr << "NUMA index: " << binding_failures_.load()
                      << " query workers could not bind locally; those workers used a shared replica\n";
    }
    size_t keyCount() const { return keys_; }
    size_t replicaCount() const { return replicas_.size(); }
    size_t replicaBytes() const { return budget_ ? budget_->mapped : 0; }
    size_t bindingFailures() const { return binding_failures_.load(); }
    std::vector<int> replicaNodes() const {
        std::vector<int> ids;
        for (const auto& node : topology_.nodes) if (!replicas_.empty()) ids.push_back(node.id);
        return ids;
    }
    bool placementMatches() const {
        for (const auto& replica : replicas_) if (!replica->arena->placementMatches()) return false;
        return true;
    }
    Local local() const {
        Local result;
        if (shared_) { result.shared_ = shared_.get(); return result; }
        result.replica_ = replicas_.front().get();
        try {
            const auto affinity = Numa::Mask::affinity();
            std::vector<size_t> eligible;
            for (size_t i = 0; i < topology_.nodes.size(); ++i)
                if (!affinity.intersection(topology_.nodes[i].cpus).empty()) eligible.push_back(i);
            if (!eligible.empty()) {
                const size_t selected = eligible[size_t(omp_get_thread_num()) % eligible.size()];
                result.binding_.reset(new Numa::ThreadBinding(topology_.nodes[selected], topology_.policy_bits));
                result.locally_bound = result.binding_->active;
                result.replica_ = replicas_[selected].get();
                result.node = topology_.nodes[selected].id;
                if (result.locally_bound) return result;
                result.binding_.reset();
            }
        } catch (...) {}
        if (binding_failures_.fetch_add(1) == 0)
            std::cerr << "NUMA index: local query binding failed; using a remote replica for affected workers\n";
        return result;
    }
};

template<class Index>
NumaIndex<Index> makeNumaIndex(Index&& index, int threads, NumaIndexOptions options = {}) {
    return NumaIndex<Index>(std::move(index), threads, options);
}

template<class Index>
const Index& localIndex(const Index& index) { return index; }
template<class Index>
typename NumaIndex<Index>::Local localIndex(const NumaIndex<Index>& index) { return index.local(); }

} // namespace Sketch
#endif
