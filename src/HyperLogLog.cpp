//#include <zlib.h>
#include "HyperLogLog.h"
#include "Sketch.h"
#include "api/RuntimeInfo.h"
#include "MurmurHash3.h"
#include "rank/RankStream.h"
#include "rank/CanonicalKmer.h"
#include <immintrin.h>
#include <stdexcept>

using namespace Sketch;

namespace {

static size_t hll_register_count(int precision) {
	if (precision < 4 || precision > 20) {
		throw std::invalid_argument("HLL precision must be in [4, 20]");
	}
	return size_t(1) << precision;
}

static void validate_hll_kmer_size(int kmer_size) {
	if (kmer_size < 1 || kmer_size > 32) {
		throw std::invalid_argument("HLL k-mer size must be in [1, 32]");
	}
}

} // namespace

HyperLogLog::HyperLogLog(int np, int kmerlen)
	: HyperLogLog(np, kmerlen, 42, false) {}

HyperLogLog::HyperLogLog(int np, bool track_witnesses, int kmerlen)
	: HyperLogLog(np, kmerlen, 42, track_witnesses) {}

HyperLogLog::HyperLogLog(int np, int kmerlen, uint64_t seed,
	                     bool track_witnesses)
	: core_(hll_register_count(np), 0),
	  witnesses_(track_witnesses ? hll_register_count(np) : 0, 0),
	  value_(0.0), np_(static_cast<uint32_t>(np)), kmerLen_(kmerlen),
	  seed_(seed), is_calculated_(0),
	  estim_(EstimationMethod::ERTL_MLE),
	  jestim_(JointEstimationMethod::ERTL_JOINT_MLE),
	  track_witnesses_(track_witnesses) {
	validate_hll_kmer_size(kmerlen);
}

HyperLogLog HyperLogLog::fromRegisters(
	uint32_t precision,
	int kmer_size,
	uint64_t seed,
	const std::vector<uint8_t>& registers) {
	HyperLogLog result(static_cast<int>(precision), kmer_size, seed, false);
	if (registers.size() != result.core_.size())
		throw std::invalid_argument(
			"HyperLogLog::fromRegisters register count does not match precision");
	const uint8_t maximum_rank = static_cast<uint8_t>(65 - precision);
	if (std::any_of(registers.begin(), registers.end(),
	                [maximum_rank](uint8_t rank) {
		                return rank > maximum_rank;
	                }))
		throw std::invalid_argument(
			"HyperLogLog::fromRegisters found an invalid register rank");
	result.core_ = registers;
	result.value_ = 0.0;
	result.is_calculated_ = 0;
	return result;
}

// ── SIMD equal-register counter ───────────────────────────────────────────────
// Counts positions where a[i] == b[i] for arrays of n uint8_t registers.
// Compile-time dispatch: AVX-512BW (64/cycle) → AVX2 (32/cycle) → scalar.
static int hll_count_equal_regs(const uint8_t* __restrict__ a,
                                 const uint8_t* __restrict__ b,
                                 int n)
{
	return static_cast<int>(Runtime::countEqualBytes(
		a, b, static_cast<size_t>(n)));
}


std::array<uint32_t,64> HyperLogLog::sum_counts(const std::vector<uint8_t> &sketchInfo) const {
	std::array<uint32_t,64> sum_count {0};//default 64
	for(uint64_t i=0; i<sketchInfo.size(); ++i){
		// lzt can be 1..(65-np_), so with small p value can be 64; clamp to avoid OOB
		uint8_t b = sketchInfo[i];
		sum_count[b >= 64 ? 63u : b]++;
	}
	return sum_count;
}





template<typename T>
void HyperLogLog::compTwoSketch(const std::vector<uint8_t> &sketch1, const std::vector<uint8_t> &sketch2, T &c1, T &c2, T &cu, T &cg1, T &cg2, T &ceq) const {
	assert(sketch1.size() == sketch2.size());
	std::array<uint32_t, 64> c1l{0}, c2l{0}, c1g{0}, c2g{0};

	const uint64_t sz = sketch1.size();
	constexpr uint64_t BLOCK = 512;
	const uint8_t* __restrict__ s1 = sketch1.data();
	const uint8_t* __restrict__ s2 = sketch2.data();

	// Block-wise local histograms to reduce random writes into c1l/c2l/c1g/c2g/ceq
	// (better cache locality, fewer store-forwarding stalls).
	for (uint64_t start = 0; start < sz; start += BLOCK) {
		const uint64_t end = (start + BLOCK < sz) ? (start + BLOCK) : sz;
		std::array<uint32_t, 64> l1l{0}, l2l{0}, l1g{0}, l2g{0}, leq{0};

		uint64_t i = start;
		const uint64_t end4 = start + ((end - start) & ~uint64_t(3));
		for (; i < end4; i += 4) {
			const uint8_t a0 = s1[i],     b0 = s2[i];
			const uint8_t a1 = s1[i+1],   b1 = s2[i+1];
			const uint8_t a2 = s1[i+2],   b2 = s2[i+2];
			const uint8_t a3 = s1[i+3],   b3 = s2[i+3];

			const int lt0 = (a0 < b0), eq0 = (a0 == b0), gt0 = 1 - lt0 - eq0;
			const int lt1 = (a1 < b1), eq1 = (a1 == b1), gt1 = 1 - lt1 - eq1;
			const int lt2 = (a2 < b2), eq2 = (a2 == b2), gt2 = 1 - lt2 - eq2;
			const int lt3 = (a3 < b3), eq3 = (a3 == b3), gt3 = 1 - lt3 - eq3;

			l1l[a0] += lt0;  l2g[b0] += lt0;  l1g[a0] += gt0;  l2l[b0] += gt0;  leq[a0] += eq0;
			l1l[a1] += lt1;  l2g[b1] += lt1;  l1g[a1] += gt1;  l2l[b1] += gt1;  leq[a1] += eq1;
			l1l[a2] += lt2;  l2g[b2] += lt2;  l1g[a2] += gt2;  l2l[b2] += gt2;  leq[a2] += eq2;
			l1l[a3] += lt3;  l2g[b3] += lt3;  l1g[a3] += gt3;  l2l[b3] += gt3;  leq[a3] += eq3;
		}
		for (; i < end; ++i) {
			const uint8_t a = s1[i], b = s2[i];
			const int lt = (a < b), eq = (a == b), gt = 1 - lt - eq;
			l1l[a] += lt;  l2g[b] += lt;  l1g[a] += gt;  l2l[b] += gt;  leq[a] += eq;
		}

		for (int k = 0; k < 64; ++k) {
			c1l[k] += l1l[k];  c2l[k] += l2l[k];
			c1g[k] += l1g[k];  c2g[k] += l2g[k];
			ceq[k] += leq[k];
		}
	}

	for (int i = 0; i < 64; ++i) {
		c1[i]  = c1l[i] + ceq[i] + c1g[i];
		c2[i]  = c2l[i] + ceq[i] + c2g[i];
		cu[i]  = c1g[i] + ceq[i] + c2g[i];
		cg1[i] = c1g[i];
		cg2[i] = c2g[i];
	}
}


// Rolling-hash optimized update():
//   - No seqRev buffer: rev_enc is computed from seq by complement-on-read (saves
//     full reverse-complement pass and cache pressure).
	//   - No upfront to-upper: the shared LUT encodes A/a->0, C/c->1, G/g->2, T/t->3
//     so rolling uses one table lookup per base.
//   - Window invalid_count: only when invalid_count==0 do we add the k-mer (N etc.).
//   - Canonical = min(fwd_enc, rev_enc) with encoding A=0,C=1,G=2,T=3 (lex order).
//   - Optional future: per-thread or block-local core_ buffer, merge at end, to
//     reduce random writes when doing multi-threaded batch updates (p=12..16).
	#define ENC(c)   (Rank::encodeBase((uint8_t)(c)))
	#define COMP(e)  (Rank::complementBase((e)))
	#define VALID(e) (Rank::validBase((e)))

	void HyperLogLog::update(char* seq) {
	if (seq == nullptr) {
		throw std::invalid_argument("HyperLogLog::update sequence must not be null");
	}
	const uint64_t LENGTH = strlen(seq);
	const int KMERLEN = kmerLen_;
	// For k<32 the forward encoding only uses 2*k bits; mask off stale upper bits.
	const uint64_t kmer_mask = (KMERLEN >= 32) ? UINT64_MAX : ((1ULL << (2 * KMERLEN)) - 1);
	if (LENGTH < (uint64_t)KMERLEN) return;
	is_calculated_ = 0;

	uint32_t qq = q();

	// Initialize rolling encodings for the k-mer at position 0.
	// fwd_enc: MSB = seq[0], LSB = seq[KMERLEN-1]. rev_enc: from complement of
	// seq[KMERLEN-1..0], so MSB = comp(seq[KMERLEN-1]), LSB = comp(seq[0]).
	uint64_t fwd_enc = 0, rev_enc = 0;
	int invalid_count = 0;
	for (int k = 0; k < KMERLEN; k++) {
		uint8_t ef = ENC(seq[k]);
		if (!VALID(ef)) invalid_count++;
		fwd_enc = ((fwd_enc << 2) | (VALID(ef) ? (ef & 3u) : 0u)) & kmer_mask;
		uint8_t er = VALID(ef) ? (COMP(ef) & 3u) : 0u;
		rev_enc = (rev_enc >> 2) | (static_cast<uint64_t>(er) << (2 * (KMERLEN - 1)));
	}

#if defined __AVX512F__  && defined __AVX512DQ__
	__m512i vconst0 = _mm512_set1_epi64(0xff51afd7ed558ccd);
	__m512i vconst1 = _mm512_set1_epi64(0xc4ceb9fe1a85ec53);
#endif
#if defined __AVX512F__  && defined __AVX512CD__
	__m512i v1 = _mm512_set1_epi64(1);
#endif

	const int lanes = 8;
	const uint64_t N = ((LENGTH - KMERLEN) / lanes) * lanes;

	// Main 8-lane loop.  Rolling: fwd/rev from seq only; rev = complement by LUT.
	// Only add k-mer when invalid_count == 0 (no N etc. in window).
	for (uint64_t i = 0; i < N; i += lanes)
	{
		uint64_t resv[8];
		bool lane_valid[8];
		for (int j = 0; j < lanes; j++)
		{
			lane_valid[j] = (invalid_count == 0);
			resv[j] = lane_valid[j] ? ((fwd_enc <= rev_enc) ? fwd_enc : rev_enc) : 0;

			// Roll: leaving = seq[i+j], entering = seq[i+j+KMERLEN].
			uint8_t ef_out = ENC(seq[i + j]);
			uint8_t ef_in  = ENC(seq[i + j + KMERLEN]);
			if (!VALID(ef_out)) invalid_count--;
			if (!VALID(ef_in))  invalid_count++;
		fwd_enc = ((fwd_enc << 2) | (VALID(ef_in) ? (ef_in & 3u) : 0u)) & kmer_mask;
		uint8_t er_in = VALID(ef_in) ? (COMP(ef_in) & 3u) : 0u;
		rev_enc = (rev_enc >> 2) | ((uint64_t)er_in << (2 * (KMERLEN - 1)));
		}

		uint64_t hashvalv[8];
		#if defined __AVX512F__ && defined __AVX512DQ__
		__m512i vb = _mm512_loadu_si512((void*)resv);
			__m512i vseed = _mm512_set1_epi64((int64_t)seed_);
		__m512i va = _mm512_xor_epi64(vb, vseed);
		__m512i vtmp = _mm512_srli_epi64(va, 33);
		vb = _mm512_xor_epi64(va, vtmp);
		va = _mm512_mullo_epi64(vb, vconst0);
		vtmp = _mm512_srli_epi64(va, 33);
		vb = _mm512_xor_epi64(va, vtmp);
		va = _mm512_mullo_epi64(vb, vconst1);
		vtmp = _mm512_srli_epi64(va, 33);
		vb = _mm512_xor_epi64(va, vtmp);
		_mm512_storeu_si512(hashvalv, vb);
		#else
		for (int j = 0; j < lanes; j++)
				hashvalv[j] = Rank::RankStream::fmix64(resv[j], seed_);
		#endif

		uint64_t indexv[8], lztv[8];
		#if defined __AVX512CD__ && defined __AVX512F__
		__m512i vhash = _mm512_loadu_si512((void*)hashvalv);
		__m512i vindex = _mm512_srli_epi64(vhash, (uint8_t)qq);
		_mm512_storeu_si512(indexv, vindex);
		__m512i vlzhash = _mm512_slli_epi64(vhash, 1);
		vhash = _mm512_or_epi64(vlzhash, v1);
		vlzhash = _mm512_slli_epi64(vhash, (uint8_t)(np_ - 1));
		vhash = _mm512_lzcnt_epi64(vlzhash);
		vlzhash = _mm512_add_epi64(vhash, v1);
		_mm512_storeu_si512(lztv, vlzhash);
		#else
		for (int j = 0; j < lanes; j++) {
			indexv[j] = hashvalv[j] >> qq;
			lztv[j] = clz(((hashvalv[j] << 1) | 1) << (np_ - 1)) + 1;
		}
		#endif

		// Witness-tracking branch hoisted out of the lane loop so the non-tracking
		// path keeps its original `std::max` (no compare-then-store penalty).
		if (track_witnesses_) {
			for (int j = 0; j < lanes; j++) {
				if (!lane_valid[j]) continue;
				const uint32_t idx = indexv[j];
				const uint8_t  lzt = (uint8_t)lztv[j];
				if (lzt > core_[idx]) {
					core_[idx] = lzt;
					witnesses_[idx] = hashvalv[j];
				}
#if LZ_COUNTER
				++clz_counts_[clz(((hashvalv[j] << 1) | 1) << (np_ - 1)) + 1];
#endif
			}
		} else {
			for (int j = 0; j < lanes; j++) {
				if (!lane_valid[j]) continue;
				core_[indexv[j]] = std::max(core_[indexv[j]], (uint8_t)lztv[j]);
#if LZ_COUNTER
				++clz_counts_[clz(((hashvalv[j] << 1) | 1) << (np_ - 1)) + 1];
#endif
			}
		}
	}

	// Remainder: continue rolling from position N; only add when invalid_count==0.
		for (uint64_t i = N; i <= LENGTH - KMERLEN; ++i)
		{
			if (invalid_count == 0) {
				uint64_t res = (fwd_enc <= rev_enc) ? fwd_enc : rev_enc;
				uint64_t hashval = Rank::RankStream::fmix64(res, seed_);
			const uint32_t index = hashval >> qq;
			const uint8_t lzt = clz(((hashval << 1) | 1) << (np_ - 1)) + 1;
			if (track_witnesses_) {
				if (lzt > core_[index]) {
					core_[index] = lzt;
					witnesses_[index] = hashval;
				}
			} else {
				core_[index] = std::max(core_[index], lzt);
			}
#if LZ_COUNTER
			++clz_counts_[clz(((hashval << 1) | 1) << (np_ - 1)) + 1];
#endif
		}
			if (i < LENGTH - KMERLEN) {
				uint8_t ef_out = ENC(seq[i]);
				uint8_t ef_in  = ENC(seq[i + KMERLEN]);
				if (!VALID(ef_out)) invalid_count--;
				if (!VALID(ef_in))  invalid_count++;
				fwd_enc = ((fwd_enc << 2) | (VALID(ef_in) ? (ef_in & 3u) : 0u)) & kmer_mask;
				uint8_t er_in = VALID(ef_in) ? (COMP(ef_in) & 3u) : 0u;
				rev_enc = (rev_enc >> 2) | ((uint64_t)er_in << (2 * (KMERLEN - 1)));
			}
	}
	#undef ENC
	#undef COMP
	#undef VALID
}
 
//void HyperLogLog::update(char* seq) {
//    const uint64_t LENGTH = strlen(seq);
//    #pragma omp parallel for
//    for(uint64_t i = 0; i < LENGTH; i++){
//        if(seq[i] > 96 && seq[i] < 123){
//            seq[i] -= 32;
//        }
//    }
//
//    char* seqRev;
//    seqRev = new char[LENGTH];
//    char table[4] = {'T','G','A','C'};
//    #pragma omp parallel for
//    for ( uint64_t i = 0; i < LENGTH; i++ )
//    {
//        char base = seq[i];
//        base >>= 1;
//        base &= 0x03;
//        seqRev[LENGTH - i - 1] = table[base];
//    }
//
//    const int KMERLEN = 32;
//    if(LENGTH < KMERLEN) {
//        delete [] seqRev;
//        return;
//    }
//
//    #pragma omp parallel for
//    for(uint64_t i=0; i<LENGTH-KMERLEN; ++i) {
//        char kmer_fwd[KMERLEN+1];
//        char kmer_rev[KMERLEN+1];
//        memcpy(kmer_fwd, seq+i, KMERLEN);
//        memcpy(kmer_rev, seqRev+LENGTH-i-KMERLEN, KMERLEN);
//        kmer_fwd[KMERLEN] = '\0';
//        kmer_rev[KMERLEN] = '\0';
//        if(memcmp(kmer_fwd, kmer_rev, KMERLEN) <= 0) {
//            addh(kmer_fwd);
//        } else {
//            addh(kmer_rev);
//        }
//    }
//
//    delete [] seqRev;
//}
//

HyperLogLog HyperLogLog::project(uint32_t target_precision) const {
	if (target_precision < 4 || target_precision > np_) {
		throw std::invalid_argument(
			"HyperLogLog::project target precision must be in [4, source precision]");
	}
	if (target_precision == np_) return *this;

	HyperLogLog ret(static_cast<int>(target_precision), kmerLen_, seed_, false);
	const uint32_t delta = np_ - target_precision;
	const uint32_t suffix_mask = (uint32_t(1) << delta) - 1;

	for (uint32_t high_index = 0; high_index < core_.size(); ++high_index) {
		const uint8_t high_rank = core_[high_index];
		if (high_rank == 0) continue;

		const uint32_t low_index = high_index >> delta;
		const uint32_t suffix = high_index & suffix_mask;
		uint8_t candidate = 0;
		if (suffix == 0) {
			candidate = static_cast<uint8_t>(delta + high_rank);
		} else {
#if defined(__GNUC__) || defined(__clang__)
			const uint32_t leading = static_cast<uint32_t>(__builtin_clz(suffix))
				- (32u - delta);
#else
			uint32_t leading = 0;
			uint32_t bit = uint32_t(1) << (delta - 1);
			while ((suffix & bit) == 0) { ++leading; bit >>= 1; }
#endif
			candidate = static_cast<uint8_t>(leading + 1);
		}
		ret.core_[low_index] = std::max(ret.core_[low_index], candidate);
	}
	ret.is_calculated_ = 0;
	return ret;
}

HyperLogLog HyperLogLog::merge(const HyperLogLog &other) const {
	metadata().requireSameFamily(other.metadata(), "HyperLogLog::merge",
	                             /*allow_resolution_change=*/true);
	const uint32_t common_precision = std::min(np_, other.np_);
	HyperLogLog left = (np_ == common_precision) ? *this : project(common_precision);
	HyperLogLog right = (other.np_ == common_precision)
		? other : other.project(common_precision);

	const bool keep_witnesses = left.track_witnesses_ && right.track_witnesses_;
	HyperLogLog ret(static_cast<int>(common_precision), kmerLen_, seed_,
	                keep_witnesses);
	for (uint64_t i = 0; i < ret.m(); ++i) {
		ret.core_[i] = std::max(left.core_[i], right.core_[i]);
		if (keep_witnesses) {
			if (left.core_[i] > right.core_[i])
				ret.witnesses_[i] = left.witnesses_[i];
			else if (left.core_[i] < right.core_[i])
				ret.witnesses_[i] = right.witnesses_[i];
			else
				ret.witnesses_[i] = std::min(
					left.witnesses_[i], right.witnesses_[i]);
		}
	}
	ret.is_calculated_ = 0;
	return ret;
}


//TODO: int hash
//HyperLogLog::inline void addh(uint64_t element) {
//	element = hash(element); //TODO: which hf_
//	add(element);
//}
inline void HyperLogLog::addh(const std::string &element) {
	//add(std::hash<std::string>{}(element));
	uint64_t res[2];
	int len = element.length();
	const uint32_t seed = static_cast<uint32_t>(seed_);
	MurmurHash3_x64_128(element.c_str(), len, seed, res);
	add(res[0]);
	
}

//TODO: different hash function
//hash() {
//}

//TODO: clz is a function in clz.h
inline void HyperLogLog::add(uint64_t hashval) {
	is_calculated_ = 0;
	const uint32_t index(hashval >> q());
	const uint8_t lzt(clz(((hashval << 1)|1) << (np_ - 1)) + 1);
	if (track_witnesses_) {
		if (lzt > core_[index]) {
			core_[index] = lzt;
			witnesses_[index] = hashval;
		}
	} else {
		core_[index] = std::max(core_[index], lzt);
	}
#if LZ_COUNTER
	++clz_counts_[clz(((hashval << 1)|1) << (np_ - 1)) + 1];
#endif
}

void HyperLogLog::addFingerprint(uint64_t fingerprint) {
	add(fingerprint);
}
//Added by liumy to show sketch for testing. 
void HyperLogLog::printSketch(){
	fprintf(stdout,"Sketch info: [");
	int vecSize = core_.size();
	for(int i=0; i<vecSize-1; ++i)
		fprintf( stdout," %u,",  core_[i] );
	fprintf(stdout," %u ]\n", core_[vecSize-1]);
}

double HyperLogLog::union_size(const HyperLogLog &other) const {
	metadata().requireSameFamily(other.metadata(), "HyperLogLog::union_size");
	if(jestim_ != JointEstimationMethod::ERTL_JOINT_MLE) {
		assert(m() == other.m()|| !std::fprintf(stderr, "sizes don't match! Size1: %zu. Size2: %zu\n", m(), other.m()));
		std::array<uint32_t,64> counts{0};
		std::vector<uint8_t> unionCore(m(),0);
		for(uint64_t i=0; i<m(); ++i){
			unionCore[i] = std::max(core_[i],other.core_[i]); 
		}
		counts = sum_counts(unionCore);
		return calculate_estimate(counts, get_estim(), m(), p(), alpha(), 1e-2);
	}
	//std::fprintf(stderr, "jestim is ERTL_JOINT_MLE: %s\n", JESTIM_STRINGS[jestim_]);
	const auto full_counts = ertl_joint(*this, other);
	return full_counts[0] + full_counts[1] + full_counts[2];
}


double HyperLogLog::jaccard_index(const HyperLogLog &h2) const {
	metadata().requireSameFamily(h2.metadata(), "HyperLogLog::jaccard_index",
	                             /*allow_resolution_change=*/true);
	if (np_ != h2.np_) {
		const uint32_t common_precision = std::min(np_, h2.np_);
		HyperLogLog left = (np_ == common_precision) ? *this : project(common_precision);
		HyperLogLog right = (h2.np_ == common_precision)
			? h2 : h2.project(common_precision);
		return left.jaccard_index(right);
	}
	if(jestim_ == JointEstimationMethod::ERTL_JOINT_MLE) {
		auto full_cmps = ertl_joint(*this, h2);
		const double denom = full_cmps[0] + full_cmps[1] + full_cmps[2];
		if(denom <= 0.) return 0.;
		return full_cmps[2] / denom;
	}
	const double us = union_size(h2);
	if(us <= 0.) return 0.;
	const double ret = (creport() + h2.creport() - us) / us;
	return std::max(0., ret);
}

template<typename T>
double HyperLogLog::ertl_ml_estimate(const T& c, unsigned p, unsigned q, double relerr) const {
	/*
	   Note --
	   Putting all these optimizations together finally gives the new cardinality estimation
	   algorithm presented as Algorithm 8. The algorithm requires mainly only elementary
	   operations. For very large cardinalities it makes sense to use the strong (46) instead
	   of the weak lower bound (47) as second starting point for the secant method. The
	   stronger bound is a much better approximation especially for large cardinalities, where
	   the extra logarithm evaluation is amortized by savings in the number of iteration cycles.
	   -Ertl paper.
TODO:  Consider adding this change to the method. This could improve our performance for other
*/

#if DEBUG
	fprintf(stdout,"[W:%s:%d] Counts: ",__PRETTY_FUNCTION__, __LINE__);
	for(int i=0; i<64; i++)
		fprintf(stdout,"%d, ", c[i]);
	fprintf(stdout,"\n");
#endif

	const uint64_t m = 1ull << p;
	if (c[q+1] == m) return std::numeric_limits<double>::infinity();

	int kMin, kMax;
	for(kMin=0; c[kMin]==0; ++kMin);
	int kMinPrime = std::max(1, kMin);
	for(kMax=q+1; kMax && c[kMax]==0; --kMax);
	int kMaxPrime = std::min(static_cast<int>(q), kMax);
	double z = 0.;
	for(int k = kMaxPrime; k >= kMinPrime; z = 0.5*z + c[k--]);
	z = ldexp(z, -kMinPrime);
	unsigned cPrime = c[q+1];
	if(q) cPrime += c[kMaxPrime];
	double gprev;
	double x;
	double a = z + c[0];
	int mPrime = m - c[0];
	gprev = z + ldexp(c[q+1], -q); // Reuse gprev, setting to 0 after.
	x = gprev <= 1.5*a ? mPrime/(0.5*gprev+a): (mPrime/gprev)*std::log1p(gprev/a);
	gprev = 0;
	double deltaX = x;
	relerr /= std::sqrt(m);
	while(deltaX > x*relerr) {
		int kappaMinus1;
		frexp(x, &kappaMinus1);
		double xPrime = ldexp(x, -std::max(static_cast<int>(kMaxPrime+1), kappaMinus1+2));
		double xPrime2 = xPrime*xPrime;
		double h = xPrime - xPrime2/3 + (xPrime2*xPrime2)*(1./45. - xPrime2/472.5);
		for(int k = kappaMinus1; k >= kMaxPrime; --k) {
			double hPrime = 1. - h;
			h = (xPrime + h*hPrime)/(xPrime+hPrime);
			xPrime += xPrime;
		}
		double g = cPrime*h;
		for(int k = kMaxPrime-1; k >= kMinPrime; --k) {
			double hPrime = 1. - h;
			h = (xPrime + h*hPrime)/(xPrime+hPrime);
			xPrime += xPrime;
			g += c[k] * h;
		}
		g += x*a;
		if(gprev < g && g <= mPrime) deltaX *= (g-mPrime)/(gprev-g);
		else                         deltaX  = 0;
		x += deltaX;
		gprev = g;
	}
	return x*m;
}
template<typename HllType>
std::array<double, 3> HyperLogLog::ertl_joint(const HllType &h1, const HllType &h2) const {
	assert(h1.m() == h2.m() || !std::fprintf(stderr, "sizes don't match! Size1: %zu. Size2: %zu\n", h1.size(), h2.size()));
	std::array<double, 3> ret;
	if(h1.get_jestim() != JointEstimationMethod::ERTL_JOINT_MLE) {
		// intersection & union
		ret[2] = h1.union_size(h2);
		ret[0] = h1.creport();
		ret[1] = h2.creport();
		ret[2] = ret[0] + ret[1] - ret[2];
		ret[0] -= ret[2];
		ret[1] -= ret[2];
		ret[2] = std::max(ret[2], 0.);
		return ret;
	}
	//    using ertl_ml_estimate;
	auto p = h1.p();
	auto q = h1.q();
	std::array<uint32_t, 64> c1{0}, c2{0}, cu{0}, ceq{0}, cg1{0}, cg2{0};
	//TODO: K->C5
	//joint_unroller ju;
	//ju.sum_arrays(h1.core(), h2.core(), c1, c2, cu, cg1, cg2, ceq);
	compTwoSketch(h1.core(), h2.core(), c1, c2, cu, cg1, cg2, ceq);
	const double cAX = h1.get_is_ready() ? h1.creport() : ertl_ml_estimate(c1, h1.p(), h1.q(), 1e-2);
	const double cBX = h2.get_is_ready() ? h2.creport() : ertl_ml_estimate(c2, h2.p(), h2.q(), 1e-2);
	const double cABX = ertl_ml_estimate(cu, h1.p(), h1.q(), 1e-2);
	// std::fprintf(stderr, "Made initials: %lf, %lf, %lf\n", cAX, cBX, cABX);
	std::array<uint32_t, 64> countsAXBhalf;
	std::array<uint32_t, 64> countsBXAhalf;
	countsAXBhalf[q] = h1.m();
	countsBXAhalf[q] = h1.m();
	for(unsigned _q = 0; _q < q; ++_q) {
		// Handle AXBhalf
		countsAXBhalf[_q] = cg1[_q] + ceq[_q] + cg2[_q + 1];
		assert(countsAXBhalf[q] >= countsAXBhalf[_q]);
		countsAXBhalf[q] -= countsAXBhalf[_q];

		// Handle BXAhalf
		countsBXAhalf[_q] = cg2[_q] + ceq[_q] + cg1[_q + 1];
		assert(countsBXAhalf[q] >= countsBXAhalf[_q]);
		countsBXAhalf[q] -= countsBXAhalf[_q];
	}
	double cAXBhalf = ertl_ml_estimate(countsAXBhalf, p, q - 1, 1e-2);
	double cBXAhalf = ertl_ml_estimate(countsBXAhalf, p, q - 1, 1e-2);
	//std::fprintf(stderr, "Made halves: %lf, %lf\n", cAXBhalf, cBXAhalf);
	ret[0] = cABX - cBX;
	ret[1] = cABX - cAX;
	double cX1 = (1.5 * cBX + 1.5*cAX - cBXAhalf - cAXBhalf);
	double cX2 = 2.*(cBXAhalf + cAXBhalf) - 3.*cABX;
	ret[2] = std::max(0., 0.5 * (cX1 + cX2));
	return ret;
}




//template<typename CountArrType>
double HyperLogLog::calculate_estimate(const std::array<uint32_t,64> &counts,
		EstimationMethod estim, uint64_t m, uint32_t p, double alpha, double relerr) const {
	assert(estim <= 3);

#if DEBUG
	fprintf(stdout,"[W:%s:%d] Counts: ",__PRETTY_FUNCTION__, __LINE__);
	fprintf(stdout,"Counts: ");
	for(int i=0; i<64; i++)
		fprintf(stdout,"%d, ", counts[i]);
	fprintf(stdout,"\n");
#endif

#if ENABLE_COMPUTED_GOTO
	static constexpr void *arr [] {&&ORREST, &&ERTL_IMPROVED_EST, &&ERTL_MLE_EST};
	goto *arr[estim];
ORREST: {
#else
			switch(estim) {
				case ORIGINAL: {
#endif
								   assert(estim != ERTL_MLE);
								   double sum = counts[0];
								   for(unsigned i = 1; i < 64; ++i) if(counts[i]) sum += std::ldexp(counts[i], -i); // 64 - p because we can't have more than that many leading 0s. This is just a speed thing.
								   //for(unsigned i = 1; i < 64 - p + 1; ++i) sum += std::ldexp(counts[i], -i); // 64 - p because we can't have more than that many leading 0s. This is just a speed thing.
								   double value(alpha * m * m / sum);
								   if(value < small_range_correction_threshold(m)) {
									   if(counts[0]) {
#if DEBUG
										   std::fprintf(stderr, "[W:%s:%d] Small value correction. Original estimate %lf. New estimate %lf.\n",
												   __PRETTY_FUNCTION__, __LINE__, value, m * std::log(static_cast<double>(m) / counts[0]));
#endif
										   value = m * std::log(static_cast<double>(m) / counts[0]);
									   }
								   } else if(value > LARGE_RANGE_CORRECTION_THRESHOLD) {
									   // Reuse sum variable to hold correction.
									   // I do think I've seen worse accuracy with the large range correction, but I would need to rerun experiments to be sure.
									   sum = -std::pow(2.0L, 32) * std::log1p(-std::ldexp(value, -32));
									   if(!std::isnan(sum)) value = sum;
#if DEBUG
									   else std::fprintf(stderr, "[W:%s:%d] Large range correction returned nan. Defaulting to regular calculation.\n", __PRETTY_FUNCTION__, __LINE__);
#endif
								   }
								   return value;
							   }
#if ENABLE_COMPUTED_GOTO
ERTL_IMPROVED_EST: {
#else
					   case ERTL_IMPROVED: {
#endif
											   static const double divinv = 1. / (2.L*std::log(2.L));
											   double z = m * gen_tau(static_cast<double>((m-counts[64 - p + 1]))/static_cast<double>(m));
											   for(unsigned i = 64-p; i; z += counts[i--], z *= 0.5); // Reuse value variable to avoid an additional allocation.
											   z += m * gen_sigma(static_cast<double>(counts[0])/static_cast<double>(m));
											   return m * divinv * m / z;
										   }
#if ENABLE_COMPUTED_GOTO
ERTL_MLE_EST: return ertl_ml_estimate(counts, p, 64 - p, relerr);
#else
				   case ERTL_MLE: return ertl_ml_estimate(counts, p, 64 - p, relerr);
				   default: return 0.0;
			   }
#endif
		}


// ── equalRegisterFraction ─────────────────────────────────────────────────────
double HyperLogLog::equalRegisterFraction(const HyperLogLog& other) const
{
	metadata().requireSameFamily(other.metadata(),
	                             "HyperLogLog::equalRegisterFraction",
	                             /*allow_resolution_change=*/true);
	if (np_ != other.np_) {
		const uint32_t common_precision = std::min(np_, other.np_);
		HyperLogLog left = (np_ == common_precision) ? *this : project(common_precision);
		HyperLogLog right = (other.np_ == common_precision)
			? other : other.project(common_precision);
		return left.equalRegisterFraction(right);
	}
	const int n = (int)core_.size();
	if(n == 0) return 0.0;
	return (double)hll_count_equal_regs(core_.data(), other.core_.data(), n) / n;
}

// ── distanceFiltered ──────────────────────────────────────────────────────────
// Applies a cheap SIMD equal-register pre-check before the expensive Ertl MLE.
// Returns -1.0 when the pair is provably below min_jaccard; exact distance otherwise.
double HyperLogLog::distanceFiltered(const HyperLogLog& other,
                                     double min_jaccard,
                                     double prefilter_factor) const
{
	if(equalRegisterFraction(other) < min_jaccard * prefilter_factor)
		return -1.0;
	return distance(other);
}

// ── containment ─────────────────────────────────────────────────────────────
// C(this ⊆ other) = (|A| + |B| - |AUB|) / |A|
double HyperLogLog::containment(const HyperLogLog& other) const
{
	const double card_a = creport();
	if (card_a <= 0.0) return 0.0;
	const double card_b = other.creport();
	const double j = jaccard_index(other);
	if (j <= 0.0) return 0.0;
	const double estimate = j * (card_a + card_b) / (card_a * (1.0 + j));
	return std::min(1.0, std::max(0.0, estimate));
}

// ── ani ──────────────────────────────────────────────────────────────────────
// ANI = (2J / (1+J))^(1/kmer_size)   (Mash / Ondov et al. 2016)
double HyperLogLog::ani(const HyperLogLog& other, int kmer_size) const
{
	const double j = jaccard_index(other);
	if (j <= 0.0) return 0.0;
	if (j >= 1.0) return 1.0;
	return std::pow(2.0 * j / (1.0 + j), 1.0 / static_cast<double>(kmer_size));
}

Rank::RankMetadata HyperLogLog::metadata() const {
	Rank::RankMetadata meta;
	meta.backend = Rank::Backend::HLL;
	meta.hash_profile = Rank::HashProfile::Murmur3Fmix64V1;
	meta.weight_semantics = Rank::WeightSemantics::UnweightedSet;
	meta.resolution_kind = Rank::ResolutionKind::RegisterPrecisionBits;
	meta.fingerprint_bits = 64;
	meta.kmer_size = static_cast<uint16_t>(kmerLen_);
	meta.resolution = np_;
	meta.seed = seed_;
	return meta;
}
