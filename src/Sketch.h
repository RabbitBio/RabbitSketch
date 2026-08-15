#ifndef Sketch_H
#define Sketch_H
#include "shuffle.h"
#include "phmap.h"
#include <map>
#include <vector>
#include <string>
//#include <string.h>
#include <stdint.h>
#include <float.h>
#include <unordered_set>
#include "MinHash.h"
#include "histoSketch.h"
#include "HyperLogLog.h"
#include "SetSketch.h"
#include "rank/RankMetadata.h"
//#include "Kssd.h"
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif  

#define COMPONENT_SZ 7
//#define MIN_SUBCTX_DIM_SMP_SZ 4096
#define _64MASK 0xffffffffffffffffLLU
#define CTX_SPC_USE_L 8
#define LD_FCTR 0.6
#define DEFAULT_CHAR_KSSD -1



//typedef struct fileInfo
//{
//	string fileName;
//	uint64_t fileSize;
//} fileInfo_t;

//typedef struct sketch
//{
//	string fileName;
//	string seqName;
//	string comment;
//	int id;
//	vector<uint32_t> hashSet;
//	vector<uint64_t> hashSet64;
//} sketch_t;


//for converting Kssd sketch into RabbitKSSD sketch format
typedef struct co_dirstat
{
  unsigned int shuf_id;
  bool koc;
  int kmerlen;
  int dim_rd_len;
  int comp_num;
  int infile_num;
  uint64_t all_ctx_ct;
} co_dstat_t;

typedef struct setResult
{
  int common;
  int size0;
  int size1;
  double jaccard;
} setResult_t;

struct DistInfo
{
  string refName;
  int common;
  int refSize;
  double jorc;
  double dist;
};
struct cmpDistInfo
{
  bool operator()(DistInfo d1, DistInfo d2){
    return d1.dist < d2.dist;
  }
};
/// \brief Sketch namespace
namespace Sketch{
  struct fileInfo_t {
    std::string fileName;
    uint64_t fileSize;
  };
  struct sketch_t {
    std::string fileName;
    std::vector<uint32_t> hashSet;
    std::vector<uint64_t> hashSet64;
    int id;
  };

  struct sketchInfo_t
  {
    int id;
    int half_k;
    int half_subk;
    int drlevel;
    int genomeNumber;
  };
  typedef uint64_t hash_t;

  struct WMHParameters
  {
    int kmerSize;
    int sketchSize;
    int windowSize;
    double * r; 
    double * c; 
    double * b; 
  };

  struct Reference
  {
    // no sequence for now
    std::string name;
    std::string comment;
    uint64_t length;
    HashList hashesSorted;
    std::vector<uint32_t> counts;
  };

  struct MashLite {
    std::string fileName;               
    int id;                            
    //std::vector<uint32_t> hashList;    
    std::vector<uint64_t> hashList64;  

    MashLite() : id(0) {}
  };

  void transMinHashes(vector<MashLite>& sketches, sketchInfo_t& info, string dictFile, string indexFile, int numThreads);
  void index_tridist_MinHash(vector<MashLite>& sketches, sketchInfo_t& info, string refSketchOut, string outputFile, int kmer_size, double maxDist, int isContainment, int numThreads);
  void saveMinHashes(vector<MashLite>& sketches, sketchInfo_t& info, string outputFile);
  // Sketching seqeunces using minhash method
  class MinHash
  {

    public:
      /// minhash init with parameters
      MinHash(int k = 21, int size = 1000, uint32_t seed = 42, bool rc = true):
        kmerSpace(0.0), minHashHeap(nullptr), kmerSize(k),
        alphabetSize(4), seed(seed), sketchSize(size),
        noncanonical(!rc), use64(true)
    {
      if (kmerSize < 1 || kmerSize > 32)
        throw std::invalid_argument("MinHash k-mer size must be in [1, 32]");
      if (sketchSize < 2)
        throw std::invalid_argument("MinHash sketch size must be at least 2");
      minHashHeap = new MinHashHeap(use64, sketchSize);

      this->kmerSpace = pow(alphabetSize, kmerSize);

      this->totalLength = 0;

      this->needToList = true;
    }

      // Explicit destructor: minHashHeap is a raw `new`-allocated pointer
      // member.  The compiler-generated default destructor would leak it,
      // costing ~16-64 KiB per sketch (a real bug for large-N runs).
      ~MinHash() {
          delete minHashHeap;
          minHashHeap = nullptr;
      }
      MinHash(const MinHash&)            = delete;
      MinHash& operator=(const MinHash&) = delete;

      static std::unique_ptr<MinHash> fromHashes(
          int kmer_size, uint32_t maximum_sketch_size, uint32_t seed,
          bool reverse_complement, const std::vector<uint64_t>& hashes,
          double estimated_cardinality, uint64_t total_length);

      /* for the containment of sequences(genomes).
       * the size of minHashHeap(as sketchSize) is proportatd with the sequence(genome) length.
       * addbyxxm 2021/9/18
       */
      int id;
      MashLite toLite(vector<uint64_t> hashL) const {
        MashLite lite;
        lite.fileName = fileName;
        lite.id = id;
        //lite.hashList = hashList;
        lite.hashList64 = hashL;
        return lite;
      }
      /// minhash is updatable with multiple sequences
      void update(char * seq);

      /// merge two minhashes
      void merge(MinHash& msh);

      /// return the jaccard index
      double jaccard(MinHash * msh);			

      double containJaccard(MinHash * msh);

      double containDistance(MinHash * msh);
      /// return mutation distance defined in Mash instead of jaccard distance
      double distance(MinHash * msh);

      //print hash values for debug
      void printMinHashes();

      //get hash values for saving
      vector<uint64_t> storeMinHashes();

      //load hash valued from files
      void loadMinHashes(vector<uint64_t> hashArr);

      /// return totalSeqence length, including multiple updates
      uint64_t getTotalLength() const {return totalLength;}

      ///Estimate the cardinality count
      int count(){
        return static_cast<int>(cardinality());
      }
      double cardinality() {
        if (minHashHeap && needToList)
          estimatedCardinality_ = minHashHeap->estimateSetSize();
        return estimatedCardinality_;
      }
      Rank::RankMetadata metadata() const;
      // parameters

      /// return kmerSize
      int getKmerSize() const { return kmerSize; }

      // return alphabet size
      //uint32_t getAlphabetSize() { return alphabetSize; }

      // return whether to preserve case
      //bool isPreserveCase() { return preserveCase; }

      // return whether to use 64bit hash
      //bool isUse64() { return use64; }
      //save file name
      string fileName;
      /// return hash seed
      uint32_t getSeed() const {return seed; }

      /// return sketch size
      uint32_t getMaxSketchSize() const {return sketchSize; }

      /// return whether to use reverse complement
      bool isReverseComplement() const { return !noncanonical; }
      bool isSealed() const noexcept { return sealed_; }

      /// Materialize the hash list AND release the build-only scratch state.
      ///
      /// Must be called once after the last update() and before pairwise
      /// distance queries.  Two effects:
      ///   1. heapToList(): heap → sorted unique bottom-k vector (compact).
      ///   2. delete minHashHeap: frees ~16-64 KiB of scratch hash table /
      ///      priority queue per sketch.  For 200k sketches this is ~3-12 GiB
      ///      reclaimed before the O(N²) distance loop starts.
      ///
      /// One-way operation: calling update() after finalize() is UB.
      /// Mirrors Kssd::finalize() and BinDash::finalize().
      void finalize() {
          if (sealed_) return;
          if (minHashHeap && needToList)
              estimatedCardinality_ = minHashHeap->estimateSetSize();
          ensureHeapToListed();
          if (minHashHeap) {
              delete minHashHeap;
              minHashHeap = nullptr;
          }
          // Compact the sorted hash vectors to their exact final size.
          // heapToList() already does std::move(reserve(sketchSize)),
          // so capacity should be ≤ sketchSize already; shrink_to_fit
          // is a no-op when size == capacity but cheap regardless.
          reference.hashesSorted.hashes64.shrink_to_fit();
          reference.hashesSorted.hashes32.shrink_to_fit();
          sealed_ = true;
      }

      /// Returns the finalized bottom-k sorted hash list for inverted-index use.
      /// Calls finalize() internally; result is valid until the sketch is modified.
      /// The returned reference points into the internal Reference struct — zero copy.
      const std::vector<uint64_t>& getHashesSorted() {
          finalize();
          return reference.hashesSorted.hashes64;
      }

      /// test whether this minhash is empty
      bool isEmpty() { 
        ensureHeapToListed();
        if(this->reference.hashesSorted.size() <= 0)
          return true;
        else
          return false;
      }

      /// get sketch size, it should be less than max sketch size
      int getSketchSize() {
        ensureHeapToListed();
        return this->reference.hashesSorted.size();
      }

    private:
      bool needToList = true;
      double kmerSpace;
      MinHashHeap * minHashHeap;
      Reference reference;
      uint64_t totalLength = 0;
      double estimatedCardinality_ = 0.0;
      bool sealed_ = false;

      double pValue(uint64_t x, uint64_t lengthRef, uint64_t lengthQuery, double kmerSpace, uint64_t sketchSize);
      void heapToList();
      void ensureHeapToListed();

      //parameters
      int kmerSize;
      uint32_t alphabetSize = 4; //nuc sequences
      uint32_t seed;
      uint64_t sketchSize; //minHashesPerWindow
      bool noncanonical;
      bool use64; //always true (remove support for hash32)

      //FIXME: perserveCase is not included in constructor
      bool preserveCase = false;
  };

  std::tuple<int, int, int, int*> read_shuffled_file(std::string filepath);
  struct kssd_parameter_t {
    int half_k;
    int half_subk;
    int drlevel;
    int rev_add_move;
    int half_outctx_len;
    int * shuffled_dim;
    // Owns the malloc-allocated table while preserving the historical raw
    // pointer field used by existing APIs. Copies share the immutable table.
    std::shared_ptr<int> shuffled_dim_owner;
    int dim_start;
    int dim_end;
    unsigned int kmer_size;
    int hashSize;
    int hashLimit;
    uint64_t domask;
    uint64_t tupmask;
    uint64_t undomask0;
    uint64_t undomask1;
    uint64_t shuffle_fingerprint;
    phmap::flat_hash_map<uint32_t, int> shuffled_map;

  private:
    static int _checked_dim_end(int half_k_, int half_subk_, int drlevel_) {
      if (half_k_ < half_subk_ || half_k_ > 16 ||
          half_subk_ < 3 || half_subk_ > 7 || drlevel_ < 0)
        throw std::invalid_argument(
          "kssd_parameter_t requires half_subk<=half_k<=16, "
          "half_subk in [3,7], and drlevel>=0");
      if (half_subk_ - drlevel_ < 3)
        throw std::invalid_argument(
          "kssd_parameter_t: half_subk - drlevel must be >= 3");
      return 1 << (4 * (half_subk_ - drlevel_));
    }

    // shared init: compute bit-masks and populate shuffled_map from shuffled_dim
    void _init_masks_and_map() {
      int comp_bittl = 64 - 4 * half_k;
      tupmask   = _64MASK >> comp_bittl;
      domask    = (tupmask >> (4 * half_outctx_len)) << (2 * half_outctx_len);
      uint64_t undomask = (tupmask ^ domask) & tupmask;
      undomask1 = undomask & (tupmask >> ((half_k + half_subk) * 2));
      undomask0 = undomask ^ undomask1;
      int dim_size  = 1 << (4 * half_subk);
      int dim_limit = 1 << (4 * (half_subk - drlevel));
      shuffled_map.reserve(dim_size);
      shuffle_fingerprint = UINT64_C(1469598103934665603);
      for (int t = 0; t < dim_size; t++) {
        shuffle_fingerprint ^= static_cast<uint32_t>(shuffled_dim[t]);
        shuffle_fingerprint *= UINT64_C(1099511628211);
        if (shuffled_dim[t] >= 0 && shuffled_dim[t] < dim_limit)
          shuffled_map.emplace(t, shuffled_dim[t]);
      }
    }

  public:
    // ── Constructor 1: no shuffle file required (recommended) ──────────────
    // The shuffle dictionary is generated in memory with the fixed canonical
    // seed (348842630), guaranteeing reproducible results without extra files.
    kssd_parameter_t(int half_k_ = 10, int half_subk_ = 6, int drlevel_ = 3)
      : half_k(half_k_), half_subk(half_subk_), drlevel(drlevel_),
        rev_add_move(4 * half_k_ - 2),
        half_outctx_len(half_k_ - half_subk_), shuffled_dim(nullptr),
        dim_start(0),
        dim_end(_checked_dim_end(half_k_, half_subk_, drlevel_)),
        kmer_size(2 * half_k_),
        hashSize(2000), hashLimit(static_cast<int>(2000 * LD_FCTR))
    {
      shuffled_dim = generate_shuffle_dim(half_subk_);
      shuffled_dim_owner = std::shared_ptr<int>(
        shuffled_dim, [](int* pointer) { std::free(pointer); });
      _init_masks_and_map();
    }

    // ── Constructor 2: load shuffle dictionary from file (legacy) ──────────
    [[deprecated("Use the 3-argument constructor; shuffle file is no longer needed.")]]
    kssd_parameter_t(int half_k_, int half_subk_, int drlevel_,
                     const string& shuffle_file)
      : half_k(half_k_), half_subk(half_subk_), drlevel(drlevel_),
        rev_add_move(4 * half_k_ - 2),
        half_outctx_len(half_k_ - half_subk_), shuffled_dim(nullptr),
        dim_start(0),
        dim_end(_checked_dim_end(half_k_, half_subk_, drlevel_)),
        kmer_size(2 * half_k_),
        hashSize(2000), hashLimit(static_cast<int>(2000 * LD_FCTR))
    {
      auto result  = Sketch::read_shuffled_file(shuffle_file);
      shuffled_dim = std::get<3>(result);
      if (shuffled_dim == nullptr)
        throw std::invalid_argument(
          "kssd_parameter_t: failed to load shuffle dictionary");
      shuffled_dim_owner = std::shared_ptr<int>(
        shuffled_dim, [](int* pointer) { std::free(pointer); });
      _init_masks_and_map();
    }
  };

  struct KssdLite {
    std::string fileName;               
    int id;                            
    std::vector<uint32_t> hashList;    
    std::vector<uint64_t> hashList64;  

    KssdLite() : id(0) {}
  };



  class Kssd{
    public:


      //Kssd(int half_k, int half_subk, int drlevel, std::vector<int> shuffled_dim)
      //	:params_(half_k, half_subk, drlevel, std::move(shuffled_dim)),
      //	dim_size_(1 << (4 * half_subk)), use64_((half_k - drlevel) > 8)
      // Takes a shared_ptr so all Kssd objects built from the same parameters
      // share one copy of the two large tables (shuffled_dim + shuffled_map).
      // Previously a value-copy was made per object, wasting ~192 MB each.
      Kssd(std::shared_ptr<const kssd_parameter_t> params,
           uint64_t seed = 42)
        : id(0),
        params_ptr_(requireParams(std::move(params))),
        half_k_(params_ptr_->half_k),
        half_subk_(params_ptr_->half_subk),
        drlevel_(params_ptr_->drlevel),
        shuffled_dim_(params_ptr_->shuffled_dim),
        dim_size_(1 << (4 * params_ptr_->half_subk)),
        use64((params_ptr_->half_k - params_ptr_->drlevel) > 8),
        half_outctx_len(params_ptr_->half_outctx_len),
        rev_add_move(params_ptr_->rev_add_move),
        kmer_size(params_ptr_->kmer_size),
        dim_start(params_ptr_->dim_start),
        dim_end(params_ptr_->dim_end),
        hashSize(params_ptr_->hashSize),
        hashLimit(params_ptr_->hashLimit),
        component_num(0),
        tupmask(params_ptr_->tupmask),
        domask(params_ptr_->domask),
        undomask0(params_ptr_->undomask0),
        undomask1(params_ptr_->undomask1),
        shuffled_map_(&params_ptr_->shuffled_map),
        seed_(seed)
    {
    }


      Kssd(const Kssd&) = delete;
      Kssd& operator=(const Kssd&) = delete;

      KssdLite toLite() const {
        KssdLite lite;
        lite.fileName = fileName;
        lite.id = id;
        lite.hashList = hashList;
        lite.hashList64 = hashList64;
        return lite;
      }


      Kssd(Kssd&&) = default;
      Kssd& operator=(Kssd&&) = default;

      ~Kssd() = default;
      phmap::flat_hash_set<uint32_t> hashSet;
      phmap::flat_hash_set<uint64_t> hashSet64;
      std::vector<uint64_t> hashList64;
      std::vector<uint32_t> hashList;

      std::string fileName;
      vector<uint32_t> storeHashes();
      vector<uint64_t> storeHashes64();
      const vector<uint32_t>& getHashes() const { return hashList; }
      const vector<uint64_t>& getHashes64() const { return hashList64; }
      void update(const char* seq);
      uint64_t getShuffleFingerprint() const {
        return params_ptr_->shuffle_fingerprint;
      }
      double admissionProbability() const;
      uint32_t admissionDimensionCount() const {
        return static_cast<uint32_t>(dim_end - dim_start);
      }
      uint32_t dimensionUniverseSize() const {
        return static_cast<uint32_t>(dim_size_);
      }
      Rank::RankMetadata metadata() const;
      /**
       * Exact state-only projection to a sparser (larger) drlevel.
       * The native reduced tuples retain the shuffled dimension rank needed
       * for filtering/remapping.  A denser state cannot be reconstructed.
       */
      Kssd project(int target_drlevel) const;
      // Compact sketch storage after all update() calls: shrinks the sorted
      // hashList vector to its exact size (dropping any capacity overshoot
      // accumulated by multi-contig merges) and releases the hashSet's heap.
      // Single-threaded; call once per sketch before pairwise distance loop.
      void finalize();
      bool isSealed() const noexcept { return sealed_; }
      bool existFile(string fileName);
      bool isFastaList(string inputList);
      bool isFastqList(string inputList);
      bool isFastaGZList(string inputList);
      bool isFastqGZList(string inputList);
      int id;
      // Different drlevels compare at max(drlevel_a, drlevel_b), the common
      // exactly-projectable sparse resolution.
      double jaccard(Kssd* kssd);
      double distance(Kssd* kssd);

      void loadHashes64(vector<uint64_t> hashArr);
      void loadHashes(vector<uint32_t> hashArr);
      //std::tuple<int, int, int, std::unique_ptr<int[]>> read_shuffled_file(std::string filepath);
      bool cmpSketch(sketch_t s1, sketch_t s2);

      //for result accuracy testing
      bool cmpSketchName(sketch_t s1, sketch_t s2);
      bool isSketchFile(string inputFile);
      //	void saveSketches(vector<Kssd*>& sketches, sketchInfo_t& info, string outputFile);
      void readSketches(vector<Kssd*>& sketches, sketchInfo_t& info, string inputFile);
      void printSketches(vector<Kssd*>& sketches, string outputFile);
      void printInfos(vector<Kssd*>& sketches, string outputFile);
      void convertSketch(vector<Kssd*>& sketches, sketchInfo_t& info, string inputDir, int numThreads);
      void tri_dist(vector<Kssd*>& sketches, string outputFile, int kmer_size, double maxDist, int numThreads);
      void dist(vector<Kssd*>& ref_sketches, vector<sketch_t>& query_sketches, string outputFile, int kmer_size, double maxDist, int numThreads);
      int get_half_subk() const;
      int get_drlevel() const;
      int get_half_k() const;

    private:
      static std::shared_ptr<const kssd_parameter_t> requireParams(
          std::shared_ptr<const kssd_parameter_t> params) {
        if (!params)
          throw std::invalid_argument("Kssd parameters must not be null");
        return params;
      }

      // Shared ownership of the parameter block (shuffled_dim + shuffled_map).
      // All Kssd objects built from the same parameters point to one instance.
      std::shared_ptr<const kssd_parameter_t> params_ptr_;
      int half_k_;
      int half_subk_;
      int drlevel_;
      int * shuffled_dim_;          // raw pointer into params_ptr_->shuffled_dim
      int dim_size_;
      bool use64=  (half_k_ - drlevel_) > 8;
      int half_outctx_len;
      int rev_add_move;
      int kmer_size;
      int dim_start;
      int dim_end;
      int hashSize;
      int hashLimit;
      int component_num;
      int comp_bittl;

      uint64_t tupmask;
      uint64_t domask;
      uint64_t undomask;
      uint64_t undomask0;
      uint64_t undomask1;
      int get_hashSize(int half_k, int drlevel);
      void SetToList64();
      void SetToList();

      // Pointer into params_ptr_->shuffled_map — no per-object copy.
      const phmap::flat_hash_map<uint32_t, int>* shuffled_map_;

      uint64_t seed_;
      bool sealed_ = false;

  };

  void index_tridist(vector<KssdLite>& sketches, sketchInfo_t& info, string refSketchOut, string outputFile, int kmer_size, double maxDist, int isContainment, int numThreads);
  void saveSketches(vector<KssdLite>& sketches, Sketch::sketchInfo_t& info, std::string outputFile);
  void transSketches(vector<KssdLite>& sketches, sketchInfo_t& info, string dictFile, string indexFile, int numThreads);
  //	struct KSSDParameters
  //	{
  //		int half_k;
  //		int half_subk;
  //		int drlevel;
  //		int* shuffled_dim;
  //	
  //		int* shuffleN(int n, int base);
  //		int* shuffle(int arr[], int length);
  //		int get_hashSize(int half_k, int drlevel);
  //		int hashSize;
  //	
  //		KSSDParameters(int half_k_=10, int half_subk_=6, int drlevel_=3):
  //			half_k(half_k_), half_subk(half_subk_), drlevel(drlevel_)
  //		{
  //			int dim_size = 1 << 4 * half_subk;
  //			shuffled_dim = (int*)malloc(sizeof(int) * dim_size);
  //			shuffled_dim = shuffleN(dim_size, 0);
  //			shuffled_dim = shuffle(shuffled_dim, dim_size);
  //			hashSize = get_hashSize(half_k, drlevel);
  //		}
  //	
  //	};
  //
  //	class KSSD
  //	{
  //		public:
  //			KSSD(KSSDParameters parameter)
  //			{
  //				for(int i = 0; i< 128; i++)
  //				{
  //					BaseMap[i] = DEFAULT_CHAR_KSSD;
  //				}
  //				BaseMap['a'] = 0;
  //				BaseMap['c'] = 1;
  //				BaseMap['g'] = 2;
  //				BaseMap['t'] = 3;
  //				BaseMap['A'] = 0;
  //				BaseMap['C'] = 1;
  //				BaseMap['G'] = 2;
  //				BaseMap['T'] = 3;
  //				shuffled_dim = parameter.shuffled_dim;
  //				half_k = parameter.half_k;
  //				half_subk = parameter.half_subk;
  //				drlevel = parameter.drlevel;
  //				hashSize = parameter.hashSize;
  //	
  //	
  //				half_outctx_len = half_k - half_subk;
  //				rev_addmove = 4 * half_k - 2;
  //				kmer_size = 2 * half_k;
  //				dim_start = 0;
  //				dim_end = MIN_SUBCTX_DIM_SMP_SZ;
  //				hashLimit = hashSize * LD_FCTR;
  //				component_num = half_k - drlevel > COMPONENT_SZ ? 1LU << 4 * (half_k - drlevel - COMPONENT_SZ) : 1;
  //				comp_bittl = 64 - 4 * half_k;
  //				tupmask = _64MASK >> comp_bittl;
  //				domask = (tupmask >> (4 * half_outctx_len)) << (2 * half_outctx_len);
  //				undomask = (tupmask ^ domask) & tupmask;
  //				undomask1 = undomask &	(tupmask >> ((half_k + half_subk) * 2));
  //				undomask0 = undomask ^ undomask1;
  //	
  //			}
  //	
  //			void update(char* seq);
  //			double jaccard(KSSD* kssd);
  //			double distance(KSSD* kssd);
  //			void printHashes();
  //			vector<uint64_t> storeHashes();
  //			void loadHashes(vector<uint64_t> hashArr);
  //			int get_half_k();
  //			int get_half_subk();
  //			int get_drlevel();
  //	
  //		private:
  //			int half_k;
  //			int half_subk;
  //			int drlevel;
  //			int half_outctx_len;
  //			int rev_addmove;
  //			int kmer_size;
  //			int dim_start;
  //			int dim_end;
  //			int hashSize;
  //			int hashLimit;
  //			Reference reference;
  //			int component_num;
  //			int comp_bittl;
  //			int BaseMap[128];
  //			int* shuffled_dim;
  //	
  //			uint64_t tupmask;
  //			uint64_t domask;
  //			uint64_t undomask;
  //			uint64_t undomask0;
  //			uint64_t undomask1;
  //	
  //			vector<uint64_t> hashList;
  //			unordered_set<uint64_t> hashSet;
  //	
  //			int get_hashSize(int half_k, int drlevel);
  //			void SetToList();
  //	
  //	};


  class WMinHash{
    public:
      //WMinHash(int k = 21, int size = 50, int windowSize = 20, double paraDWeight = 0.0):
      WMinHash(WMHParameters parameters, double paraDWeight = 0.0):
        kmerSize(parameters.kmerSize), histoSketchSize(parameters.sketchSize), minimizerWindowSize(parameters.windowSize), paraDecayWeight(paraDWeight)
    {	
      //numBins = pow(kmerSize, alphabetSize); //need to be confirmed
      //histoDimension = pow(kmerSize, alphabetSize); //need to be confirmed
      numBins = pow(kmerSize, 4); //consult from HULK
      histoDimension = numBins; //need to be confirmed

      binsArr = (double *)malloc(numBins * sizeof(double));
      memset(binsArr, 0, numBins*sizeof(double));//init binsArr

      int g = ceil(2 / EPSILON);
      int d = ceil(log(1 - DELTA) / log(0.5));
      countMinSketch = (double *)malloc(d * g *sizeof(double));
      memset(countMinSketch, 0, d*g*sizeof(double));//init countMinSketch 

      r = parameters.r;
      c = parameters.c;
      b = parameters.b;

      //	//the r, c, b and getCWS need to be outClass
      //	r = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
      //	c = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
      //	b = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
      //	getCWS(r, c, b, histoSketchSize, histoDimension);

      ////for getCWS test debug
      //	FILE * fp;
      //	fp = fopen("cws.txt", "w");
      //	for(int i = 0; i < histoSketchSize*histoDimension; i++){
      //		fprintf(fp, "%ld\t%lf\t%lf\t%lf\n", i, r[i], c[i], b[i]);
      //	}

      histoSketches = (uint32_t *) malloc (histoSketchSize * sizeof(uint32_t));
      histoWeight = (double *) malloc (histoSketchSize * sizeof(double));
      for(int i = 0; i < histoSketchSize; i++){
        histoWeight[i] = DBL_MAX;
      }
      //memset(histoWeight, 0, histoSketchSize * sizeof(double));

      //add the applyConceptDrift and decayWeight.
      if(paraDecayWeight < 0.0 || paraDecayWeight > 1.0){
        cerr << "the paraDecayWeight must between 0.0 and 1.0 " << endl;
        exit(1);
      }
      else{
        applyConceptDrift = true;
      }
      if(paraDecayWeight == 1.0){
        applyConceptDrift = false;
      }
      decayWeight = 1.0;
      if(applyConceptDrift){
        decayWeight = exp(-paraDecayWeight);
      }

      needToCompute = true;

    }


      ~WMinHash();

      /// weightedMinHash is updatable with multiple sequences
      void update(char * seq);

      void computeHistoSketch();

      /// return the jaccard index
      double wJaccard(WMinHash * wmh);

      /// return the distance which is negative correlation with jaccard index
      double distance(WMinHash * wmh);

      /// print hash value fo debug
      void getWMinHash();

      //void setKmerSize(int kmerSizeNew) { kmerSize = kmerSizeNew; }
      //void setAlphabetSize(int alphabetSizeNew) { alphabetSize = alphabetSizeNew; }
      //void setNumBins(int numBinsNew) { numBins = numBinsNew; }
      //void setMinimizerWindowSize(int minimizerWindowSizeNew) { minimizerWindowSize = minimizerWindowSizeNew; }
      //void setHistoSketchSize(int histoSketchSizeNew); //{ histoSketchSize = histoSketchSizeNew; }
      //void setHistoDimension(int histoDimensionNew); //{ histoDimension = histoDimensionNew; }
      //void setParaDecayWeight(double paraDecayWeightNew) { paraDecayWeight = paraDecayWeightNew; }
      //void setApplyConceptDrift(bool applyConceptDriftNew) { applyConceptDrift = applyConceptDriftNew; }

      /// return kmerSize
      int getKmerSize() { return kmerSize; }

      /// return the minimizer window size
      int getMinimizerWindowSize() { return minimizerWindowSize; }

      //int getAlphabetSize() { return alphabetSize; }
      //int getNumBins() { return numBins; }
      //int getMinimizerWindowSize() { return minimizerWindowSize; }
      //int getHistoSketchSize() { return histoSketchSize; }
      //int getHistoDimension() { return histoDimension; }
      //double getParaDecayWeight() { return paraDecayWeight; }
      //bool isApplyComceptDrift() { return applyConceptDrift; }




    private:
      WMHParameters parameters;
      int kmerSize;
      int alphabetSize = 4;
      int numBins;
      int minimizerWindowSize;
      int histoSketchSize;
      int histoDimension;
      double paraDecayWeight;
      bool applyConceptDrift;
      double decayWeight;

      bool needToCompute = true;
      double * binsArr;
      double * countMinSketch; 
      std::vector<uint64_t> sketches;
      std::vector<Bin> kmerSpectrums;
      uint32_t * histoSketches;
      double * histoWeight;


      double * r;
      double * c;
      double * b;



  };

  //OMinHash
  struct OSketch {
    //std::string       name;
    int               k = 0, l = 0, m = 0;
    std::vector<char> data;
    std::vector<char> rcdata;

    bool operator==(const OSketch& rhs) const {
      return k == rhs.k && l == rhs.l && m == rhs.m && data == rhs.data && rcdata == rhs.rcdata;
    }

  };

  /// Sketching and compare sequences or strings using Order MinHash algorithm.
  class OrderMinHash{

    public:
      /// OrderMinHash constructor
      OrderMinHash() = default;
      /// OrderMinHash constructor for sketching sequences using default parameters
      OrderMinHash(char * seqNew);
      explicit OrderMinHash(const std::string& sequence);
      ~OrderMinHash() = default;

      /** Restore a validated immutable native sketch. */
      static OrderMinHash fromSketch(const OSketch& sketch,
          uint64_t seed, bool reverse_complement);

      /// return sketch result in `OSketch` type
      OSketch getSektch() const { return sk; }
      const OSketch& getSketch() const noexcept { return sk; }

      /** \rst
        Build a `OrderMinHash` sketch.
        `seqNew` is NULL pointer in default.
        If seqNew is NULL pointer, buildSketch() will rebuild sketh using old data.
        This is useful when chaning parameters and build a new sketch.
        \endrst
        */
      void buildSketch(char * seqNew);
      void buildSketch(const std::string& sequence);

      /** 
        Return similarity between two `OrderMinHash` sketches. In `OrderMinHash` class, there is no jaccard function provied. Because `OrderMinHash` is a proxy of edit distance instead of jaccard index.
        */
      double similarity(const OrderMinHash & omh2) const;

      /// Return distance between two `OrderMinHash` sketches
      double distance(const OrderMinHash & omh2) const
      {
        return (double)1.0 - similarity(omh2);
      }

      /// Set parameter `kmerSize`: default 21.
      void setK(int k){ m_k = k; sealed_ = false; }

      /// Set parameter `l`: default 2 (normally 2 - 5).
      void setL(int l){ m_l = l; sealed_ = false; }

      /// Set parameter `m`: default 500.
      void setM(int m){ m_m = m; sealed_ = false; }

      /// Set seed value for random generator: default 32.
      void setSeed(uint64_t seedNew) { mtSeed = seedNew; sealed_ = false; }

      /** 
        Choose whether to deal with reverse complement sequences: default false.
        Reverse complement is normally used in biological sequences such as DNA or protein sequences.
        */
      void setReverseComplement(bool isRC){rc = isRC; sealed_ = false;}

      /// Return parameter `kmerSize`.
      int getK() const noexcept {return m_k;}

      /// Return parameter `l`.
      int getL() const noexcept {return m_l;}

      /// Return parameter `m`.
      int getM() const noexcept {return m_m;}

      /// Return random generator seed value.
      uint64_t getSeed() const noexcept { return mtSeed; }

      /// Test whether to deal with reverse complement kmers.
      bool isReverseComplement() const noexcept {return rc;}
      bool isSealed() const noexcept { return sealed_; }
      Rank::RankMetadata metadata() const;

    private:
      std::string sequence_;
      std::string rc_sequence_;
      //Parameters parameters;
      int m_k = 21, m_l = 2, m_m = 500;
      //reverse complement
      bool rc = false;
      bool sealed_ = false;
      OSketch sk;
      uint64_t mtSeed = 32; //default value

      void sketch();

      inline void compute_sketch(char * ptr, const char * seq);

      double compare_sketches(const OSketch& sk1, const OSketch& sk2, 
          ssize_t m = -1, bool circular = false) const;
      double compare_sketch_pair(const char* p1, const char* p2,
          unsigned m, unsigned k, unsigned l, bool circular) const;

  };

  class HyperLogLog{

    public:
      explicit HyperLogLog(int np, int kmerlen = 32);

      // Witness-tracking constructor: when true, every register update also
      // records the 64-bit hash of the k-mer that "won" the register.
      // Used by the inverted-index --hll path (mirrors SetSketch witnesses).
      // Memory cost: 8 * (1 << np) bytes.
      HyperLogLog(int np, bool track_witnesses, int kmerlen = 32);

      /// Seed-aware constructor. The four-argument order avoids ambiguity with
      /// the historical (precision, track_witnesses, kmer_length) overload.
      HyperLogLog(int np, int kmerlen, uint64_t seed, bool track_witnesses);

      /** Restore a validated native register state without original sequence. */
      static HyperLogLog fromRegisters(
          uint32_t precision,
          int kmer_size,
          uint64_t seed,
          const std::vector<uint8_t>& registers);

      ~HyperLogLog(){};
      void update(char* seq);
      /**
       * Add an already coordinated 64-bit RankStream fingerprint.
       * This bypasses sequence canonicalization and hashing; callers are
       * responsible for following the contract reported by metadata().
       */
      void addFingerprint(uint64_t fingerprint);
      /// Exact HLL precision down-projection. Projected lower-precision states
      /// intentionally drop witnesses because native high-p winners do not
      /// contain enough information to reconstruct low-p winner identities.
      HyperLogLog project(uint32_t target_precision) const;
      HyperLogLog merge(const HyperLogLog &other) const;
      void printSketch();
      double distance(const HyperLogLog &h2) const {return 1.0 - jaccard_index(h2);}
      //double jaccard_index(HyperLogLog &h2); 
      double jaccard_index(const HyperLogLog &h2) const;
      /// Cardinality estimate (k-mer count). Used for size-based pre-filtering.
      double cardinality() const { return creport(); }
      /// Raw register array. Used by external LSH banding.
      const std::vector<uint8_t>& getCore() const { return core_; }
      uint32_t getPrecision() const { return np_; }
      int getKmerSize() const { return kmerLen_; }
      uint64_t getSeed() const { return seed_; }
      Rank::RankMetadata metadata() const;
      /// Per-register "winning" k-mer hash – only populated when constructed
      /// with track_witnesses=true. Used by inverted-index --index mode as
      /// high-entropy candidate keys (mirrors SetSketch witnesses).
      /// Returns an empty vector when witnesses are not being tracked.
      const std::vector<uint64_t>& getWitnesses() const { return witnesses_; }
      bool tracksWitnesses() const { return track_witnesses_; }
      /// Release the witness array (e.g. after KMV-key extraction in
      /// --index mode) to reclaim ~8 bytes/register before the verify pass.
      void clearWitnesses() {
          std::vector<uint64_t>().swap(witnesses_);
          track_witnesses_ = false;
      }

      /// Fraction of registers whose values are identical in both sketches.
      /// Computed with SIMD (AVX-512BW / AVX2 / scalar fallback).
      /// Cheap O(m) pre-filter: for Jaccard > J the fraction is empirically
      /// well above J*0.3, so this call is ~30-60x faster than distance().
      double equalRegisterFraction(const HyperLogLog& other) const;

      /// Fast-path distance: applies a cheap equal-register pre-filter before
      /// calling the full Ertl joint MLE.
      /// Returns -1.0 if the pair is provably dissimilar
      ///   (equalRegisterFraction < min_jaccard * prefilter_factor).
      /// Returns the Jaccard distance (1 - Jaccard) otherwise.
      /// prefilter_factor = 0.3 is very conservative (near-zero false negatives).
      double distanceFiltered(const HyperLogLog& other,
                              double min_jaccard,
                              double prefilter_factor = 0.3) const;

      /// Containment of *this in other: |A ∩ B| / |A|.
      ///   C(A⊆B) = (|A| + |B| - |AUB|) / |A|
      /// Returns 0 if cardinality of this sketch is zero.
      double containment(const HyperLogLog& other) const;

      /// Average Nucleotide Identity estimated from Jaccard similarity.
      ///   ANI = (2J / (1+J))^(1/kmer_size)
      /// @param kmer_size  k-mer length used during sketching (default 32).
      double ani(const HyperLogLog& other, int kmer_size = 32) const;

    protected:
      std::vector<uint8_t>  core_;//sketchInfo; 
      std::vector<uint64_t> witnesses_;          // populated only when track_witnesses_
      mutable double value_; //cardinality
      uint32_t np_; // 10-20
      int      kmerLen_;    // k-mer length used in update() (default 32)
      uint64_t seed_;       // coordinated RankStream seed
      mutable uint8_t is_calculated_;
      EstimationMethod                        estim_;
      JointEstimationMethod                  jestim_;
      bool                                   track_witnesses_;
      //HashStruct                                 hf_;

    private:
      uint32_t p() const {return np_;}//verification
      uint32_t q() const {return (sizeof(uint64_t) * CHAR_BIT) - np_;}
      uint64_t m() const {return static_cast<uint64_t>(1) << np_;}
      size_t size() const {return size_t(m());}
      bool get_is_ready() const {return is_calculated_;}
      const auto &core()    const {return core_;}
      EstimationMethod get_estim()       const {return  estim_;}
      JointEstimationMethod get_jestim() const {return jestim_;}
      void set_estim(EstimationMethod val) { estim_ = std::max(val, ERTL_MLE);}
      void set_jestim(JointEstimationMethod val) { jestim_ = val;}
      void set_jestim(uint16_t val) {set_jestim(static_cast<JointEstimationMethod>(val));}
      void set_estim(uint16_t val)  {estim_  = static_cast<EstimationMethod>(val);}

      // Returns cardinality estimate. Sums if not calculated yet.
      double creport() const {
        csum();
        return value_;
      }
      double report() noexcept {
        csum();
        return creport();
      }


      //private:
      void add(uint64_t hashval);
      void addh(const std::string &element);
      double alpha()          const {return make_alpha(m());}
      static double small_range_correction_threshold(uint64_t m) {return 2.5 * m;}
      double union_size(const HyperLogLog &other) const;
      // Call sum to recalculate if you have changed contents.
      void csum() const { if(!is_calculated_) sum(); }
      void sum() const {
        const auto counts(sum_counts(core_)); // std::array<uint32_t, 64>  // K->C
        value_ = calculate_estimate(counts, estim_, m(), np_, alpha(), 1e-2);
        is_calculated_ = 1;
      }
      std::array<uint32_t,64> sum_counts(const std::vector<uint8_t> &sketchInfo) const;
      double calculate_estimate(const std::array<uint32_t,64> &counts, EstimationMethod estim, uint64_t m, uint32_t p, double alpha, double relerr) const; 
      template<typename T>
        void compTwoSketch(const std::vector<uint8_t> &sketch1, const std::vector<uint8_t> &sketch2, T &c1, T &c2, T &cu, T &cg1, T &cg2, T &ceq) const;
      template<typename T>
        double ertl_ml_estimate(const T& c, unsigned p, unsigned q, double relerr=1e-2) const; 
      template<typename HllType>
        std::array<double, 3> ertl_joint(const HllType &h1, const HllType &h2) const; 



  };


  //std::tuple<int, int, int, std::unique_ptr<int[]>> read_shuffled_file(std::string filepath);
  //std::tuple<int, int, int, std::unique_ptr<int[]>> read_shuffled_file_wrapper(std::string filepath) {
  //    return read_shuffled_file(filepath);
  //}
}//namespace sketch

#endif //Sketch_h
