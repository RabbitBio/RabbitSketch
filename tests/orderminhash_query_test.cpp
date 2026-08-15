#include "Sketch.h"
#include "api/SketchBuilder.h"
#include "estimators/UnifiedQuery.h"

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Exception, typename Function>
void expectThrows(Function function, const std::string& message) {
    try {
        function();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error(message);
}

Sketch::OrderMinHash build(const std::string& sequence, int samples = 64) {
    Sketch::OrderMinHash result;
    result.setK(5);
    result.setL(2);
    result.setM(samples);
    result.setSeed(42);
    result.setReverseComplement(true);
    result.buildSketch(sequence);
    return result;
}

void testLifecycleAndNativeQuery() {
    const std::string sequence =
        "ACGTACGTTGCAACGTACGTTGCAACGTACGTTGCA";
    Sketch::OrderMinHash left = build(sequence, 64);
    Sketch::OrderMinHash right = build(sequence, 32);
    check(left.isSealed() &&
              left.metadata().backend == Sketch::Rank::Backend::OrderMinHash &&
              left.metadata().resolution == 64,
          "OrderMinHash lifecycle metadata is incomplete");

    Sketch::Query::Request request;
    request.metric = Sketch::Query::Metric::OrderSimilarity;
    const auto result = Sketch::Query::query(left, right, request);
    check(result.estimate_available && result.value == 1.0 &&
              result.order_distance == 0.0 &&
              result.effective_samples == 32 && result.plan.exact_projection,
          "OrderMinHash unified prefix query is incorrect");

    request.metric = Sketch::Query::Metric::Jaccard;
    check(!Sketch::Query::query(left, right, request).estimate_available,
          "OrderMinHash incorrectly exposed a set Jaccard estimand");

    left.setM(16);
    check(!left.isSealed(),
          "changing OrderMinHash parameters did not invalidate lifecycle state");
    expectThrows<std::logic_error>(
        [&]() { (void)left.similarity(right); },
        "unbuilt OrderMinHash parameters were silently compared");
}

void testOwnershipAndFastxBuilder() {
    const std::string sequence =
        "ACGTACGTTGCAACGTACGTTGCAACGTACGTTGCA";
    char mutable_sequence[] =
        "ACGTACGTTGCAACGTACGTTGCAACGTACGTTGCA";
    Sketch::OrderMinHash owned(mutable_sequence);
    mutable_sequence[0] = 'T';
    Sketch::OrderMinHash expected(sequence);
    check(owned.similarity(expected) == 1.0,
          "OrderMinHash retained a borrowed sequence pointer");

    Sketch::API::SketchConfig config;
    config.algorithm = Sketch::API::Algorithm::OrderMinHash;
    config.sampling = Sketch::API::SamplingMode::AlgorithmNative;
    config.aggregation = Sketch::API::AggregationMode::OneSketchPerRecord;
    config.ambiguous_policy = Sketch::API::AmbiguousPolicy::RejectRecord;
    config.kmer_size = 5;
    config.order_l = 2;
    config.order_m = 32;
    Sketch::IO::FastxRecord record;
    record.name = "ordered";
    record.sequence = sequence;
    Sketch::API::MultiSketchBuilder builder({config});
    builder.update(record);
    const auto built = builder.finish("ordered", "memory", record.name);
    Sketch::Query::Request request;
    request.metric = Sketch::Query::Metric::OrderSimilarity;
    check(built.size() == 1 &&
              built[0].stats.accepted_kmers ==
                  sequence.size() - config.kmer_size + 1 &&
              built[0].sketch.metadata().backend ==
                  Sketch::Rank::Backend::OrderMinHash &&
              built[0].sketch.query(built[0].sketch, request).value == 1.0,
          "unified record builder did not produce a queryable OrderMinHash");
}

} // namespace

int main() {
    try {
        testLifecycleAndNativeQuery();
        testOwnershipAndFastxBuilder();
        std::cout << "orderminhash_query_test: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "orderminhash_query_test: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
