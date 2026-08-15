#ifndef __PYBIND_H__
#define __PYBIND_H__

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <limits>
#include "api/SketchBuilder.h"
#include "api/SketchConfig.h"
#include "api/RuntimeInfo.h"
#include "api/Version.h"
#include "BinDash.h"
#include "fastkmv.h"
#include "fracminhash.h"
#include "probmh.h"
#include "estimators/UnifiedQuery.h"
#include "rank/RankMetadata.h"
#include "rank/RankStream.h"
#include "io/FastxReader.h"
namespace py = pybind11;



PYBIND11_MODULE(rabbitsketch, m) {
  m.doc() = "rabbitsketch pybind";
  m.attr("__version__") = Sketch::API::VERSION;

  py::enum_<Sketch::Rank::Backend>(m, "RankBackend")
    .value("FastKMV", Sketch::Rank::Backend::FastKMV)
    .value("KSSD", Sketch::Rank::Backend::KSSD)
    .value("ProbMinHash", Sketch::Rank::Backend::ProbMinHash)
    .value("SetSketch", Sketch::Rank::Backend::SetSketch)
    .value("HLL", Sketch::Rank::Backend::HLL)
    .value("BinDash", Sketch::Rank::Backend::BinDash)
    .value("LegacyMinHash", Sketch::Rank::Backend::LegacyMinHash)
    .value("FracMinHash", Sketch::Rank::Backend::FracMinHash)
    .value("OrderMinHash", Sketch::Rank::Backend::OrderMinHash);

  py::enum_<Sketch::Rank::Canonicalization>(m, "RankCanonicalization")
    .value("Lexicographic2BitReverseComplementV1",
           Sketch::Rank::Canonicalization::Lexicographic2BitReverseComplementV1)
    .value("ForwardOnlyV1",
           Sketch::Rank::Canonicalization::ForwardOnlyV1);

  py::enum_<Sketch::Rank::HashProfile>(m, "RankHashProfile")
    .value("CanonicalBits", Sketch::Rank::HashProfile::CanonicalBits)
    .value("Murmur3Fmix64V1", Sketch::Rank::HashProfile::Murmur3Fmix64V1)
    .value("Murmur3Fmix64Twice",
           Sketch::Rank::HashProfile::Murmur3Fmix64Twice)
    .value("LegacyMurmur3X64V1",
           Sketch::Rank::HashProfile::LegacyMurmur3X64V1)
    .value("BinDashDoubleFmixV1",
           Sketch::Rank::HashProfile::BinDashDoubleFmixV1)
    .value("ProbMinHashRaceV1",
           Sketch::Rank::HashProfile::ProbMinHashRaceV1)
    .value("FracRankStreamV1",
           Sketch::Rank::HashProfile::FracRankStreamV1)
    .value("OrderMinHashV1", Sketch::Rank::HashProfile::OrderMinHashV1);

  py::enum_<Sketch::Rank::WeightSemantics>(m, "RankWeightSemantics")
    .value("UnweightedSet", Sketch::Rank::WeightSemantics::UnweightedSet)
    .value("Frequency", Sketch::Rank::WeightSemantics::Frequency)
    .value("UserSupplied", Sketch::Rank::WeightSemantics::UserSupplied);

  py::enum_<Sketch::API::Algorithm>(m, "Algorithm")
    .value("LegacyMinHash", Sketch::API::Algorithm::LegacyMinHash)
    .value("FastKMV", Sketch::API::Algorithm::FastKMV)
    .value("FracMinHash", Sketch::API::Algorithm::FracMinHash)
    .value("HyperLogLog", Sketch::API::Algorithm::HyperLogLog)
    .value("SetSketch", Sketch::API::Algorithm::SetSketch)
    .value("KSSD", Sketch::API::Algorithm::KSSD)
    .value("ProbMinHash", Sketch::API::Algorithm::ProbMinHash)
    .value("BinDash", Sketch::API::Algorithm::BinDash)
    .value("OrderMinHash", Sketch::API::Algorithm::OrderMinHash);

  py::enum_<Sketch::API::MoleculeType>(m, "MoleculeType")
    .value("DNA", Sketch::API::MoleculeType::DNA)
    .value("RNA", Sketch::API::MoleculeType::RNA)
    .value("Protein", Sketch::API::MoleculeType::Protein)
    .value("Dayhoff", Sketch::API::MoleculeType::Dayhoff)
    .value("HydrophobicPolar",
           Sketch::API::MoleculeType::HydrophobicPolar);

  py::enum_<Sketch::API::AmbiguousPolicy>(m, "AmbiguousPolicy")
    .value("SkipKmer", Sketch::API::AmbiguousPolicy::SkipKmer)
    .value("RejectRecord", Sketch::API::AmbiguousPolicy::RejectRecord)
    .value("MapToUnknown", Sketch::API::AmbiguousPolicy::MapToUnknown);

  py::enum_<Sketch::API::AggregationMode>(m, "AggregationMode")
    .value("OneSketchPerFile",
           Sketch::API::AggregationMode::OneSketchPerFile)
    .value("OneSketchPerRecord",
           Sketch::API::AggregationMode::OneSketchPerRecord)
    .value("OneSketchPerCollection",
           Sketch::API::AggregationMode::OneSketchPerCollection);

  py::enum_<Sketch::API::SamplingMode>(m, "SamplingMode")
    .value("BottomK", Sketch::API::SamplingMode::BottomK)
    .value("Scaled", Sketch::API::SamplingMode::Scaled)
    .value("RegisterPrecision",
           Sketch::API::SamplingMode::RegisterPrecision)
    .value("KssdReduction", Sketch::API::SamplingMode::KssdReduction)
    .value("AlgorithmNative", Sketch::API::SamplingMode::AlgorithmNative);

  py::enum_<Sketch::API::Capability>(m, "Capability")
    .value("SequenceUpdate", Sketch::API::Capability::SequenceUpdate)
    .value("RecordUpdate", Sketch::API::Capability::RecordUpdate)
    .value("HashUpdate", Sketch::API::Capability::HashUpdate)
    .value("WeightedUpdate", Sketch::API::Capability::WeightedUpdate)
    .value("Seal", Sketch::API::Capability::Seal)
    .value("Cardinality", Sketch::API::Capability::Cardinality)
    .value("Jaccard", Sketch::API::Capability::Jaccard)
    .value("JaccardDistance", Sketch::API::Capability::JaccardDistance)
    .value("MashDistance", Sketch::API::Capability::MashDistance)
    .value("Containment", Sketch::API::Capability::Containment)
    .value("Intersection", Sketch::API::Capability::Intersection)
    .value("Union", Sketch::API::Capability::Union)
    .value("ANI", Sketch::API::Capability::ANI)
    .value("WeightedJaccard", Sketch::API::Capability::WeightedJaccard)
    .value("WeightedContainment",
           Sketch::API::Capability::WeightedContainment)
    .value("OrderSimilarity", Sketch::API::Capability::OrderSimilarity)
    .value("OrderDistance", Sketch::API::Capability::OrderDistance)
    .value("Merge", Sketch::API::Capability::Merge)
    .value("Project", Sketch::API::Capability::Project)
    .value("SearchKeys", Sketch::API::Capability::SearchKeys);

  py::class_<Sketch::API::SketchConfig>(m, "SketchConfig")
    .def(py::init<>())
    .def_readwrite("algorithm", &Sketch::API::SketchConfig::algorithm)
    .def_readwrite("molecule", &Sketch::API::SketchConfig::molecule)
    .def_readwrite("ambiguous_policy",
                   &Sketch::API::SketchConfig::ambiguous_policy)
    .def_readwrite("aggregation", &Sketch::API::SketchConfig::aggregation)
    .def_readwrite("sampling", &Sketch::API::SketchConfig::sampling)
    .def_readwrite("weight_semantics",
                   &Sketch::API::SketchConfig::weight_semantics)
    .def_readwrite("kmer_size", &Sketch::API::SketchConfig::kmer_size)
    .def_readwrite("seed", &Sketch::API::SketchConfig::seed)
    .def_readwrite("resolution", &Sketch::API::SketchConfig::resolution)
    .def_readwrite("scaled", &Sketch::API::SketchConfig::scaled)
    .def_readwrite("canonical", &Sketch::API::SketchConfig::canonical)
    .def_readwrite("track_abundance",
                   &Sketch::API::SketchConfig::track_abundance)
    .def_readwrite("homopolymer_compressed",
                   &Sketch::API::SketchConfig::homopolymer_compressed)
    .def_readwrite("min_abundance",
                   &Sketch::API::SketchConfig::min_abundance)
    .def_readwrite("max_abundance",
                   &Sketch::API::SketchConfig::max_abundance)
    .def_readwrite("minimum_base_quality",
                   &Sketch::API::SketchConfig::minimum_base_quality)
    .def_readwrite("bindash_bits",
                   &Sketch::API::SketchConfig::bindash_bits)
    .def_readwrite("probminhash_max_l",
                   &Sketch::API::SketchConfig::probminhash_max_l)
    .def_readwrite("setsketch_base",
                   &Sketch::API::SketchConfig::setsketch_base)
    .def_readwrite("setsketch_a",
                   &Sketch::API::SketchConfig::setsketch_a)
    .def_readwrite("kssd_half_subk",
                   &Sketch::API::SketchConfig::kssd_half_subk)
    .def_readwrite("order_l", &Sketch::API::SketchConfig::order_l)
    .def_readwrite("order_m", &Sketch::API::SketchConfig::order_m)
    .def_readwrite("memory_budget_bytes",
                   &Sketch::API::SketchConfig::memory_budget_bytes)
    .def_readwrite("target_relative_error",
                   &Sketch::API::SketchConfig::target_relative_error)
    .def("validate", &Sketch::API::SketchConfig::validate)
    .def("fingerprint", &Sketch::API::SketchConfig::fingerprint);

  py::class_<Sketch::API::Capabilities>(m, "Capabilities")
    .def_readonly("algorithm", &Sketch::API::Capabilities::algorithm)
    .def_readonly("flags", &Sketch::API::Capabilities::flags)
    .def("supports", &Sketch::API::Capabilities::supports);

  m.def("capabilities_for", &Sketch::API::capabilitiesFor,
        py::arg("algorithm"));

  py::class_<Sketch::Runtime::RuntimeInfo>(m, "RuntimeInfo")
    .def_readonly("architecture",
                  &Sketch::Runtime::RuntimeInfo::architecture)
    .def_readonly("compiler", &Sketch::Runtime::RuntimeInfo::compiler)
    .def_readonly("portable_baseline",
                  &Sketch::Runtime::RuntimeInfo::portable_baseline)
    .def_readonly("runtime_dispatch_available",
                  &Sketch::Runtime::RuntimeInfo::runtime_dispatch_available)
    .def_readonly("cpu_sse2", &Sketch::Runtime::RuntimeInfo::cpu_sse2)
    .def_readonly("cpu_sse41", &Sketch::Runtime::RuntimeInfo::cpu_sse41)
    .def_readonly("cpu_avx2", &Sketch::Runtime::RuntimeInfo::cpu_avx2)
    .def_readonly("cpu_avx512f",
                  &Sketch::Runtime::RuntimeInfo::cpu_avx512f)
    .def_readonly("cpu_avx512bw",
                  &Sketch::Runtime::RuntimeInfo::cpu_avx512bw)
    .def_readonly("selected_byte_path",
                  &Sketch::Runtime::RuntimeInfo::selected_byte_path)
    .def_readonly("selected_u64_path",
                  &Sketch::Runtime::RuntimeInfo::selected_u64_path)
    .def_readonly("environment_override",
                  &Sketch::Runtime::RuntimeInfo::environment_override)
    .def_readonly("environment_override_honored",
                  &Sketch::Runtime::RuntimeInfo::environment_override_honored);

  m.def("runtime_info", []() {
    return Sketch::Runtime::runtimeInfo();
  });

  py::class_<Sketch::API::BuildStats>(m, "BuildStats")
    .def_readonly("records", &Sketch::API::BuildStats::records)
    .def_readonly("input_bases", &Sketch::API::BuildStats::input_bases)
    .def_readonly("candidate_kmers",
                  &Sketch::API::BuildStats::candidate_kmers)
    .def_readonly("accepted_kmers",
                  &Sketch::API::BuildStats::accepted_kmers)
    .def_readonly("skipped_ambiguous_kmers",
                  &Sketch::API::BuildStats::skipped_ambiguous_kmers)
    .def_readonly("skipped_low_quality_kmers",
                  &Sketch::API::BuildStats::skipped_low_quality_kmers)
    .def_readonly("distinct_kmers",
                  &Sketch::API::BuildStats::distinct_kmers)
    .def_readonly("accepted_distinct_kmers",
                  &Sketch::API::BuildStats::accepted_distinct_kmers)
    .def_readonly("skipped_abundance_kmers",
                  &Sketch::API::BuildStats::skipped_abundance_kmers)
    .def_readonly("sealed", &Sketch::API::BuildStats::sealed);

  py::enum_<Sketch::IO::FastxFormat>(m, "FastxFormat")
    .value("Unknown", Sketch::IO::FastxFormat::Unknown)
    .value("Fasta", Sketch::IO::FastxFormat::Fasta)
    .value("Fastq", Sketch::IO::FastxFormat::Fastq);

  py::class_<Sketch::IO::FastxRecord>(m, "FastxRecord")
    .def(py::init<>())
    .def_readwrite("index", &Sketch::IO::FastxRecord::index)
    .def_readwrite("name", &Sketch::IO::FastxRecord::name)
    .def_readwrite("comment", &Sketch::IO::FastxRecord::comment)
    .def_readwrite("sequence", &Sketch::IO::FastxRecord::sequence)
    .def_readwrite("quality", &Sketch::IO::FastxRecord::quality)
    .def_property_readonly("has_quality",
                           &Sketch::IO::FastxRecord::hasQuality);

  py::class_<Sketch::IO::FastxReadStats>(m, "FastxReadStats")
    .def_readonly("records", &Sketch::IO::FastxReadStats::records)
    .def_readonly("bases", &Sketch::IO::FastxReadStats::bases)
    .def_readonly("quality_bases",
                  &Sketch::IO::FastxReadStats::quality_bases);

  py::class_<Sketch::IO::FastxReader>(m, "FastxReader")
    .def(py::init<std::string>(), py::arg("path"))
    .def("read", [](Sketch::IO::FastxReader& reader) -> py::object {
      Sketch::IO::FastxRecord record;
      bool available = false;
      {
        py::gil_scoped_release release;
        available = reader.next(record);
      }
      if (!available) return py::none();
      return py::cast(std::move(record));
    })
    .def("__iter__", [](Sketch::IO::FastxReader& reader)
      -> Sketch::IO::FastxReader& { return reader; },
      py::return_value_policy::reference_internal)
    .def("__next__", [](Sketch::IO::FastxReader& reader) {
      Sketch::IO::FastxRecord record;
      bool available = false;
      {
        py::gil_scoped_release release;
        available = reader.next(record);
      }
      if (!available) throw py::stop_iteration();
      return record;
    })
    .def_property_readonly("path", &Sketch::IO::FastxReader::path)
    .def_property_readonly("format", &Sketch::IO::FastxReader::format)
    .def_property_readonly("stats", &Sketch::IO::FastxReader::stats,
                           py::return_value_policy::reference_internal);

  py::enum_<Sketch::Rank::ResolutionKind>(m, "RankResolutionKind")
    .value("BottomK", Sketch::Rank::ResolutionKind::BottomK)
    .value("RegisterPrecisionBits",
           Sketch::Rank::ResolutionKind::RegisterPrecisionBits)
    .value("KssdDrLevel", Sketch::Rank::ResolutionKind::KssdDrLevel)
    .value("RegisterCount", Sketch::Rank::ResolutionKind::RegisterCount)
    .value("Scaled", Sketch::Rank::ResolutionKind::Scaled)
    .value("OrderSampleCount",
           Sketch::Rank::ResolutionKind::OrderSampleCount);

  py::enum_<Sketch::Rank::RankStream::Domain>(m, "RankDomain")
    .value("Admission", Sketch::Rank::RankStream::Domain::Admission)
    .value("Bucket", Sketch::Rank::RankStream::Domain::Bucket)
    .value("Register", Sketch::Rank::RankStream::Domain::Register)
    .value("WeightedRace", Sketch::Rank::RankStream::Domain::WeightedRace);

  py::class_<Sketch::Rank::RankStream>(m, "RankStream")
    .def(py::init<uint64_t>(), py::arg("seed") = 42)
    .def_property_readonly("seed", &Sketch::Rank::RankStream::seed)
    .def("fingerprint", &Sketch::Rank::RankStream::fingerprint,
         py::arg("canonical_code"))
    .def("derive", &Sketch::Rank::RankStream::derive,
         py::arg("fingerprint"), py::arg("domain"), py::arg("counter") = 0)
    .def("exponential_rank", &Sketch::Rank::RankStream::exponentialRank,
         py::arg("fingerprint"), py::arg("counter"), py::arg("weight"));

  py::class_<Sketch::Rank::RankMetadata>(m, "RankMetadata")
    .def_property_readonly("schema_version",
      [](const Sketch::Rank::RankMetadata& meta) { return meta.schema_version; })
    .def_property_readonly("rank_stream_version",
      [](const Sketch::Rank::RankMetadata& meta) { return meta.rank_stream_version; })
    .def_property_readonly("backend",
      [](const Sketch::Rank::RankMetadata& meta) { return meta.backend; })
    .def_property_readonly("canonicalization",
      [](const Sketch::Rank::RankMetadata& meta) {
        return meta.canonicalization;
      })
    .def_property_readonly("hash_profile",
      [](const Sketch::Rank::RankMetadata& meta) { return meta.hash_profile; })
    .def_property_readonly("weight_semantics",
      [](const Sketch::Rank::RankMetadata& meta) {
        return meta.weight_semantics;
      })
    .def_property_readonly("resolution_kind",
      [](const Sketch::Rank::RankMetadata& meta) {
        return meta.resolution_kind;
      })
    .def_property_readonly("fingerprint_bits",
      [](const Sketch::Rank::RankMetadata& meta) { return meta.fingerprint_bits; })
    .def_property_readonly("kmer_size",
      [](const Sketch::Rank::RankMetadata& meta) { return meta.kmer_size; })
    .def_property_readonly("resolution",
      [](const Sketch::Rank::RankMetadata& meta) { return meta.resolution; })
    .def_property_readonly("seed",
      [](const Sketch::Rank::RankMetadata& meta) { return meta.seed; })
    .def_property_readonly("parameter_fingerprint",
      [](const Sketch::Rank::RankMetadata& meta) {
        return meta.parameter_fingerprint;
      })
    .def("coordinated_with", &Sketch::Rank::RankMetadata::coordinatedWith)
    .def("coordination_mismatch",
         &Sketch::Rank::RankMetadata::coordinationMismatch)
    .def("same_family_mismatch",
         &Sketch::Rank::RankMetadata::sameFamilyMismatch,
         py::arg("other"), py::arg("allow_resolution_change") = false)
    .def("to_bytes", [](const Sketch::Rank::RankMetadata& meta) {
      const auto bytes = meta.toBytes();
      return py::bytes(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    })
    .def_static("from_bytes", [](py::bytes raw) {
      const std::string bytes = raw;
      return Sketch::Rank::RankMetadata::fromBytes(
        reinterpret_cast<const uint8_t*>(bytes.data()), bytes.size());
    });
  m.def("save_sketches", &Sketch::saveSketches, 
      py::arg("sketches"), py::arg("info"), py::arg("filename"),
      "Save sketches to a file");

  m.def("trans_sketches", &Sketch::transSketches, 
      py::arg("sketches"), py::arg("info"), py::arg("dict_file"), py::arg("index_file"), py::arg("num_threads"),
      "Transform sketches with additional parameters");

  m.def("index_dict", &Sketch::index_tridist, 
      py::arg("sketches"), py::arg("info"), py::arg("ref_sketch_out"), py::arg("output_file"), 
      py::arg("kmer_size"), py::arg("max_dist"), py::arg("is_containment"), py::arg("num_threads"),
      "Compute index dictionary for sketches");


  py::class_<Sketch::MinHash>(m, "MinHash")
    //.def(py::init<>())
    //.def(py::init<int>(), py::arg("k"))
    //.def(py::init<int, int>(), py::arg("k"), py::kwonly(), py::arg("size"))
    .def(py::init<int, int, uint32_t>(), py::arg("kmer") = 21, py::arg("size")=1000, py::arg("seed")=42)
    //.def(py::init<int, int, uint32_t>())
    .def("update", &Sketch::MinHash::update,
         py::call_guard<py::gil_scoped_release>())
    .def("merge", &Sketch::MinHash::merge,
         py::call_guard<py::gil_scoped_release>())
    .def("finalize", &Sketch::MinHash::finalize,
         py::call_guard<py::gil_scoped_release>())
    .def("jaccard", &Sketch::MinHash::jaccard,
         py::call_guard<py::gil_scoped_release>())
    .def("distance", &Sketch::MinHash::distance,
         py::call_guard<py::gil_scoped_release>())
    .def("get_total_length", &Sketch::MinHash::getTotalLength)
    .def("print_min_hashes", &Sketch::MinHash::printMinHashes)
    .def("count", &Sketch::MinHash::count)
    .def("cardinality", &Sketch::MinHash::cardinality)
    .def("metadata", &Sketch::MinHash::metadata)

    //parameters
    //.def("setKmerSize", &Sketch::MinHash::setKmerSize)
    //.def("setAlphabetSize", &Sketch::MinHash::setAlphabetSize)
    //.def("setPreserveCase", &Sketch::MinHash::setPreserveCase)
    //.def("setUse64", &Sketch::MinHash::setUse64)
    //.def("setSeed", &Sketch::MinHash::setSeed)
    //.def("setSketchSize", &Sketch::MinHash::setSketchSize)
    //.def("setNoncanonical", &Sketch::MinHash::setNoncanonical)
    .def("get_kmer_size", &Sketch::MinHash::getKmerSize)
    .def("get_seed", &Sketch::MinHash::getSeed)
    .def("get_max_sketch_size", &Sketch::MinHash::getMaxSketchSize)
    .def("get_sketch_size", &Sketch::MinHash::getSketchSize)
    .def("is_empty", &Sketch::MinHash::isEmpty)
    .def("is_reverse_complement", &Sketch::MinHash::isReverseComplement)
    .def_property_readonly("sealed", &Sketch::MinHash::isSealed)
    .def_property_readonly("hashes", [](Sketch::MinHash& sketch) {
      return sketch.getHashesSorted();
    })
    ;
  py::class_<Sketch::OSketch>(m, "OSketch")
    .def(py::init<>())
    .def_readwrite("k", &Sketch::OSketch::k)
    .def_readwrite("l", &Sketch::OSketch::l)
    .def_readwrite("m", &Sketch::OSketch::m)
    .def_readwrite("data", &Sketch::OSketch::data)
    .def_readwrite("rcdata", &Sketch::OSketch::rcdata)
    .def("__eq__", &Sketch::OSketch::operator==);

  py::class_<Sketch::OrderMinHash>(m, "OrderMinHash")
    .def(py::init<>())
    .def(py::init<const std::string&>(), py::arg("sequence"))
    .def("build_sketch", [](Sketch::OrderMinHash& sketch,
                              const std::string& sequence) {
      sketch.buildSketch(sequence);
    }, py::arg("sequence"), py::call_guard<py::gil_scoped_release>())
    .def("similarity", &Sketch::OrderMinHash::similarity,
         py::call_guard<py::gil_scoped_release>())
    .def("distance", &Sketch::OrderMinHash::distance,
         py::call_guard<py::gil_scoped_release>())
    .def("set_k", &Sketch::OrderMinHash::setK)
    .def("set_l", &Sketch::OrderMinHash::setL)
    .def("set_m", &Sketch::OrderMinHash::setM)
    .def("set_seed", &Sketch::OrderMinHash::setSeed)
    .def("set_reverse_complement", &Sketch::OrderMinHash::setReverseComplement)
    .def("metadata", &Sketch::OrderMinHash::metadata)
    .def_property_readonly("kmer_size", &Sketch::OrderMinHash::getK)
    .def_property_readonly("l", &Sketch::OrderMinHash::getL)
    .def_property_readonly("sample_count", &Sketch::OrderMinHash::getM)
    .def_property_readonly("seed", &Sketch::OrderMinHash::getSeed)
    .def_property_readonly("reverse_complement",
                           &Sketch::OrderMinHash::isReverseComplement)
    .def_property_readonly("sealed", &Sketch::OrderMinHash::isSealed)
    .def_property_readonly("sketch", &Sketch::OrderMinHash::getSektch)
    ;

  py::class_<Sketch::HyperLogLog>(m, "HyperLogLog")
    .def(py::init<int, int>(), py::arg("precision") = 10,
         py::arg("kmer_size") = 32)
    .def(py::init([](int precision, int kmer_size, uint64_t seed) {
             return Sketch::HyperLogLog(
                 precision, kmer_size, seed, false);
         }),
         py::arg("precision"), py::arg("kmer_size"), py::arg("seed"))
    .def("jaccard", &Sketch::HyperLogLog::jaccard_index,
         py::call_guard<py::gil_scoped_release>())
    .def("update", &Sketch::HyperLogLog::update,
         py::call_guard<py::gil_scoped_release>())
    .def("add_fingerprint", &Sketch::HyperLogLog::addFingerprint)
    .def("project", &Sketch::HyperLogLog::project)
    .def("merge", &Sketch::HyperLogLog::merge,
         py::call_guard<py::gil_scoped_release>())
    .def("distance", &Sketch::HyperLogLog::distance,
         py::call_guard<py::gil_scoped_release>())
    .def("cardinality", &Sketch::HyperLogLog::cardinality)
    .def("containment", &Sketch::HyperLogLog::containment)
    .def("metadata", &Sketch::HyperLogLog::metadata)
    .def_property_readonly("precision", &Sketch::HyperLogLog::getPrecision)
    .def_property_readonly("kmer_size", &Sketch::HyperLogLog::getKmerSize)
    .def_property_readonly("seed", &Sketch::HyperLogLog::getSeed)
    .def_property_readonly("core", &Sketch::HyperLogLog::getCore)
    ;

  py::class_<Sketch::FastKMV>(m, "FastKMV")
    .def(py::init<uint32_t, int, uint64_t>(),
         py::arg("size") = 1024, py::arg("kmer_size") = 21,
         py::arg("seed") = 42)
    .def("update", [](Sketch::FastKMV& sketch, const std::string& sequence) {
      sketch.update(sequence.data(), sequence.size());
    }, py::call_guard<py::gil_scoped_release>())
    .def("finalize", &Sketch::FastKMV::finalize,
         py::call_guard<py::gil_scoped_release>())
    .def("project", &Sketch::FastKMV::project)
    .def("merge", &Sketch::FastKMV::merge,
         py::call_guard<py::gil_scoped_release>())
    .def("jaccard", &Sketch::FastKMV::jaccard,
         py::arg("other"), py::arg("min_jaccard") = 0.0,
         py::call_guard<py::gil_scoped_release>())
    .def("distance", &Sketch::FastKMV::distance,
         py::arg("other"),
         py::arg("max_distance") = std::numeric_limits<double>::infinity(),
         py::call_guard<py::gil_scoped_release>())
    .def("cardinality", &Sketch::FastKMV::cardinality)
    .def("containment", &Sketch::FastKMV::containment)
    .def("ani", &Sketch::FastKMV::ani)
    .def("metadata", &Sketch::FastKMV::metadata)
    .def_property_readonly("size", &Sketch::FastKMV::getK)
    .def_property_readonly("kmer_size", &Sketch::FastKMV::getKmerSize)
    .def_property_readonly("seed", &Sketch::FastKMV::getSeed)
    .def_property_readonly("sealed", &Sketch::FastKMV::isSealed);

  py::class_<Sketch::BinDash>(m, "BinDash")
    .def(py::init<uint32_t, int, uint32_t, uint64_t>(),
         py::arg("sketchsize64") = 32, py::arg("kmer_size") = 21,
         py::arg("bbits") = 16, py::arg("seed") = 42)
    .def("update", [](Sketch::BinDash& sketch,
                       const std::string& sequence) {
      sketch.update(sequence.data(), sequence.size());
    }, py::call_guard<py::gil_scoped_release>())
    .def("finalize", &Sketch::BinDash::finalize,
         py::call_guard<py::gil_scoped_release>())
    .def("jaccard", &Sketch::BinDash::jaccard,
         py::call_guard<py::gil_scoped_release>())
    .def("distance", &Sketch::BinDash::distance)
    .def("cardinality", &Sketch::BinDash::cardinality)
    .def("containment", &Sketch::BinDash::containment)
    .def("ani", &Sketch::BinDash::ani)
    .def("metadata", &Sketch::BinDash::metadata)
    .def_property_readonly("num_bins", &Sketch::BinDash::getNumBins)
    .def_property_readonly("bbits", &Sketch::BinDash::getBbits)
    .def_property_readonly("kmer_size", &Sketch::BinDash::getKmerSize)
    .def_property_readonly("seed", &Sketch::BinDash::getSeed)
    .def_property_readonly("sealed", &Sketch::BinDash::isSealed)
    .def_property_readonly("packed_signatures",
      [](const Sketch::BinDash& sketch) {
        return sketch.packedSignatures();
      });

  py::class_<Sketch::FracMinHash>(m, "FracMinHash")
    .def(py::init<uint32_t, int, uint64_t>(),
         py::arg("scaled") = 1000, py::arg("kmer_size") = 21,
         py::arg("seed") = 42)
    .def("update", [](Sketch::FracMinHash& sketch,
                       const std::string& sequence) {
      sketch.update(sequence.data(), sequence.size());
    }, py::call_guard<py::gil_scoped_release>())
    .def("add_hash", &Sketch::FracMinHash::addHash,
         py::arg("coordinated_fingerprint"))
    .def("finalize", &Sketch::FracMinHash::finalize,
         py::call_guard<py::gil_scoped_release>())
    .def("project", &Sketch::FracMinHash::project,
         py::arg("target_scaled"))
    .def("merge", &Sketch::FracMinHash::merge,
         py::call_guard<py::gil_scoped_release>())
    .def("jaccard", &Sketch::FracMinHash::jaccard,
         py::call_guard<py::gil_scoped_release>())
    .def("distance", &Sketch::FracMinHash::distance)
    .def("cardinality", &Sketch::FracMinHash::cardinality)
    .def("containment", &Sketch::FracMinHash::containment)
    .def("ani", &Sketch::FracMinHash::ani)
    .def("metadata", &Sketch::FracMinHash::metadata)
    .def_property_readonly("scaled", &Sketch::FracMinHash::getScaled)
    .def_property_readonly("kmer_size", &Sketch::FracMinHash::getKmerSize)
    .def_property_readonly("seed", &Sketch::FracMinHash::getSeed)
    .def_property_readonly("threshold", &Sketch::FracMinHash::threshold)
    .def_property_readonly("sampling_probability",
                           &Sketch::FracMinHash::samplingProbability)
    .def_property_readonly("size", &Sketch::FracMinHash::size)
    .def_property_readonly("hashes", &Sketch::FracMinHash::hashes)
    .def_property_readonly("sealed", &Sketch::FracMinHash::isSealed);

  py::class_<Sketch::ProbMinHash4>(m, "ProbMinHash")
    .def(py::init<uint32_t, int, uint64_t, uint32_t,
                  Sketch::Rank::WeightSemantics>(),
         py::arg("registers") = 1024, py::arg("kmer_size") = 21,
         py::arg("seed") = 42, py::arg("max_l") = 0,
         py::arg("weight_semantics") =
             Sketch::Rank::WeightSemantics::UnweightedSet)
    .def("update", [](Sketch::ProbMinHash4& sketch,
                       const std::string& sequence) {
      sketch.update(sequence.data(), sequence.size());
    }, py::call_guard<py::gil_scoped_release>())
    .def("update_weighted", [](Sketch::ProbMinHash4& sketch,
                                const std::string& sequence,
                                double weight_each) {
      sketch.updateWeighted(sequence.data(), sequence.size(), weight_each);
    }, py::arg("sequence"), py::arg("weight_each"),
       py::call_guard<py::gil_scoped_release>())
    .def("update_entropy", [](Sketch::ProbMinHash4& sketch,
                               const std::string& sequence, double w_min) {
      sketch.updateEntropy(sequence.data(), sequence.size(), w_min);
    }, py::arg("sequence"), py::arg("w_min") = 0.1,
       py::call_guard<py::gil_scoped_release>())
    .def("add_hash", &Sketch::ProbMinHash4::addHash,
         py::arg("hash"), py::arg("weight") = 1.0)
    .def("finalize", &Sketch::ProbMinHash4::finalize,
         py::call_guard<py::gil_scoped_release>())
    .def("jaccard", &Sketch::ProbMinHash4::jaccard,
         py::call_guard<py::gil_scoped_release>())
    .def("weighted_jaccard", &Sketch::ProbMinHash4::jaccard_weighted,
         py::call_guard<py::gil_scoped_release>())
    .def("distance", &Sketch::ProbMinHash4::distance)
    .def("containment", &Sketch::ProbMinHash4::containment)
    .def("ani", &Sketch::ProbMinHash4::ani)
    .def("merge", &Sketch::ProbMinHash4::merge,
         py::call_guard<py::gil_scoped_release>())
    .def("metadata", &Sketch::ProbMinHash4::metadata)
    .def_property_readonly("register_count", &Sketch::ProbMinHash4::getM)
    .def_property_readonly("kmer_size", &Sketch::ProbMinHash4::getKmerSize)
    .def_property_readonly("seed", &Sketch::ProbMinHash4::getSeed)
    .def_property_readonly("max_l", &Sketch::ProbMinHash4::getMaxL)
    .def_property_readonly("total_weight", &Sketch::ProbMinHash4::total_weight)
    .def_property_readonly("weight_semantics",
                           &Sketch::ProbMinHash4::weightSemantics)
    .def_property_readonly("sealed", &Sketch::ProbMinHash4::isSealed)
    .def_property_readonly("register_values",
      [](const Sketch::ProbMinHash4& sketch) {
        return std::vector<double>(sketch.getRegisters(),
                                   sketch.getRegisters() + sketch.getM());
      })
    .def_property_readonly("winner_hashes",
      [](const Sketch::ProbMinHash4& sketch) {
        return std::vector<uint64_t>(sketch.getWinners(),
                                     sketch.getWinners() + sketch.getM());
      });

  py::class_<Sketch::SetSketch>(m, "SetSketch")
    .def(py::init([](int precision, double base, double a,
                     int kmer_size, uint64_t seed) {
      return std::unique_ptr<Sketch::SetSketch>(new Sketch::SetSketch(
          precision, base, a, kmer_size, false, seed));
    }),
         py::arg("precision") = 14, py::arg("base") = 2.0,
         py::arg("a") = 5.0, py::arg("kmer_size") = 32,
         py::arg("seed") = 42)
    .def("update", [](Sketch::SetSketch& sketch, const std::string& sequence) {
      std::string mutable_sequence = sequence;
      sketch.update(mutable_sequence.data(), mutable_sequence.size());
    }, py::call_guard<py::gil_scoped_release>())
    .def("add_fingerprint", &Sketch::SetSketch::addFingerprint)
    .def("project", &Sketch::SetSketch::project)
    .def("merge", &Sketch::SetSketch::merge,
         py::call_guard<py::gil_scoped_release>())
    .def("jaccard", &Sketch::SetSketch::jaccard_index,
         py::call_guard<py::gil_scoped_release>())
    .def("distance", &Sketch::SetSketch::distance)
    .def("cardinality", &Sketch::SetSketch::cardinality)
    .def("containment", &Sketch::SetSketch::containment)
    .def("metadata", &Sketch::SetSketch::metadata)
    .def_property_readonly("precision", &Sketch::SetSketch::getPrecision)
    .def_property_readonly("kmer_size", &Sketch::SetSketch::getKmerSize)
    .def_property_readonly("seed", &Sketch::SetSketch::getSeed)
    .def_property_readonly("core", &Sketch::SetSketch::getCore);

  py::class_<Sketch::kssd_parameter_t,
             std::shared_ptr<Sketch::kssd_parameter_t>>(m, "kssd_parameter_t")
    .def(py::init<int, int, int>(), py::arg("half_k") = 10,
         py::arg("half_subk") = 6, py::arg("drlevel") = 3)
    .def(py::init<int, int, int, string>())
    .def_readwrite("half_k", &Sketch::kssd_parameter_t::half_k)
    .def_readwrite("half_subk", &Sketch::kssd_parameter_t::half_subk)
    .def_readwrite("drlevel", &Sketch::kssd_parameter_t::drlevel)
    .def_readwrite("shuffled_dim", &Sketch::kssd_parameter_t::shuffled_dim)
    .def_readwrite("hash_size", &Sketch::kssd_parameter_t::hashSize);


  //py::class_<Sketch::KssdLite>(m, "KssdLite")
  //  .def(py::init<>())  
  //  .def_readwrite("fileName", &Sketch::KssdLite::fileName)
  //  .def_readwrite("id", &Sketch::KssdLite::id)
  //  .def_readwrite("hashList", &Sketch::KssdLite::hashList)
  //  .def_readwrite("hashList64", &Sketch::KssdLite::hashList64);



  py::class_<Sketch::Kssd>(m, "Kssd")
    .def(py::init([](std::shared_ptr<Sketch::kssd_parameter_t> params,
                     uint64_t seed) {
      return std::make_unique<Sketch::Kssd>(
        std::move(params), seed);
    }), py::arg("parameters"),
        py::arg("seed") = 42)
    .def("update", [](Sketch::Kssd& sketch, const std::string& sequence) {
      sketch.update(sequence.c_str());
    }, py::call_guard<py::gil_scoped_release>())
    .def("jaccard", &Sketch::Kssd::jaccard,
         py::call_guard<py::gil_scoped_release>())
    .def("distance", &Sketch::Kssd::distance)
    .def("project", &Sketch::Kssd::project,
         py::arg("target_drlevel"))
    .def("finalize", &Sketch::Kssd::finalize,
         py::call_guard<py::gil_scoped_release>())
    .def_property_readonly("sealed", &Sketch::Kssd::isSealed)
    .def("store_hashes", &Sketch::Kssd::storeHashes)
    .def("store_hashes64", &Sketch::Kssd::storeHashes64)
    .def("metadata", &Sketch::Kssd::metadata)
    .def_property_readonly("admission_probability",
      &Sketch::Kssd::admissionProbability)
    .def_property_readonly("admission_dimension_count",
      &Sketch::Kssd::admissionDimensionCount)
    .def_property_readonly("dimension_universe_size",
      &Sketch::Kssd::dimensionUniverseSize)
    .def_property_readonly("shuffle_fingerprint",
      &Sketch::Kssd::getShuffleFingerprint)
    .def("get_halfk", &Sketch::Kssd::get_half_k)
    .def("get_half_subk", &Sketch::Kssd::get_half_subk)
    .def("toLite", &Sketch::Kssd::toLite)
    .def("get_drlevel", &Sketch::Kssd::get_drlevel)
    .def_readwrite("fileName", &Sketch::Kssd::fileName);

  py::enum_<Sketch::Query::Metric>(m, "QueryMetric")
    .value("Jaccard", Sketch::Query::Metric::Jaccard)
    .value("Intersection", Sketch::Query::Metric::Intersection)
    .value("LeftContainment", Sketch::Query::Metric::LeftContainment)
    .value("RightContainment", Sketch::Query::Metric::RightContainment)
    .value("ANI", Sketch::Query::Metric::ANI)
    .value("JaccardDistance", Sketch::Query::Metric::JaccardDistance)
    .value("MashDistance", Sketch::Query::Metric::MashDistance)
    .value("Union", Sketch::Query::Metric::Union)
    .value("MaxContainment", Sketch::Query::Metric::MaxContainment)
    .value("AverageContainment", Sketch::Query::Metric::AverageContainment)
    .value("LeftContainmentANI", Sketch::Query::Metric::LeftContainmentANI)
    .value("RightContainmentANI", Sketch::Query::Metric::RightContainmentANI)
    .value("MaxContainmentANI", Sketch::Query::Metric::MaxContainmentANI)
    .value("WeightedJaccard", Sketch::Query::Metric::WeightedJaccard)
    .value("WeightedContainment", Sketch::Query::Metric::WeightedContainment)
    .value("OrderSimilarity", Sketch::Query::Metric::OrderSimilarity)
    .value("OrderDistance", Sketch::Query::Metric::OrderDistance);

  py::enum_<Sketch::Query::Support>(m, "QuerySupport")
    .value("Supported", Sketch::Query::Support::Supported)
    .value("Unsupported", Sketch::Query::Support::Unsupported)
    .value("Incompatible", Sketch::Query::Support::Incompatible);

  py::enum_<Sketch::Query::Method>(m, "QueryMethod")
    .value("None", Sketch::Query::Method::None)
    .value("FastKmvNative", Sketch::Query::Method::FastKmvNative)
    .value("HllNative", Sketch::Query::Method::HllNative)
    .value("SetSketchNative", Sketch::Query::Method::SetSketchNative)
    .value("KssdNative", Sketch::Query::Method::KssdNative)
    .value("ProbMinHashNative", Sketch::Query::Method::ProbMinHashNative)
    .value("BinDashNative", Sketch::Query::Method::BinDashNative)
    .value("LegacyMinHashNative", Sketch::Query::Method::LegacyMinHashNative)
    .value("FracMinHashNative", Sketch::Query::Method::FracMinHashNative)
    .value("OrderMinHashNative", Sketch::Query::Method::OrderMinHashNative);

  py::class_<Sketch::Query::Request>(m, "QueryRequest")
    .def(py::init<>())
    .def_readwrite("metric", &Sketch::Query::Request::metric)
    .def_readwrite("confidence_level",
                   &Sketch::Query::Request::confidence_level);

  py::class_<Sketch::Query::SketchDescriptor>(m, "SketchDescriptor")
    .def_readonly("metadata", &Sketch::Query::SketchDescriptor::metadata);

  py::class_<Sketch::Query::Plan>(m, "QueryPlan")
    .def_readonly("metric", &Sketch::Query::Plan::metric)
    .def_readonly("support", &Sketch::Query::Plan::support)
    .def_readonly("method", &Sketch::Query::Plan::method)
    .def_readonly("left_backend", &Sketch::Query::Plan::left_backend)
    .def_readonly("right_backend", &Sketch::Query::Plan::right_backend)
    .def_readonly("common_resolution",
                  &Sketch::Query::Plan::common_resolution)
    .def_readonly("exact_projection",
                  &Sketch::Query::Plan::exact_projection)
    .def_readonly("reason", &Sketch::Query::Plan::reason)
    .def_property_readonly("executable", &Sketch::Query::Plan::executable)
    .def_property_readonly("metric_name", [](const Sketch::Query::Plan& plan) {
      return std::string(Sketch::Query::metricName(plan.metric));
    })
    .def_property_readonly("support_name", [](const Sketch::Query::Plan& plan) {
      return std::string(Sketch::Query::supportName(plan.support));
    })
    .def_property_readonly("method_name", [](const Sketch::Query::Plan& plan) {
      return std::string(Sketch::Query::methodName(plan.method));
    });

  py::class_<Sketch::Query::Result>(m, "QueryResult")
    .def_readonly("plan", &Sketch::Query::Result::plan)
    .def_readonly("estimate_available",
                  &Sketch::Query::Result::estimate_available)
    .def_readonly("value", &Sketch::Query::Result::value)
    .def_readonly("jaccard", &Sketch::Query::Result::jaccard)
    .def_readonly("jaccard_distance",
                  &Sketch::Query::Result::jaccard_distance)
    .def_readonly("mash_distance",
                  &Sketch::Query::Result::mash_distance)
    .def_readonly("ani", &Sketch::Query::Result::ani)
    .def_readonly("intersection", &Sketch::Query::Result::intersection)
    .def_readonly("union_size", &Sketch::Query::Result::union_size)
    .def_readonly("left_containment",
                  &Sketch::Query::Result::left_containment)
    .def_readonly("right_containment",
                  &Sketch::Query::Result::right_containment)
    .def_readonly("max_containment",
                  &Sketch::Query::Result::max_containment)
    .def_readonly("average_containment",
                  &Sketch::Query::Result::average_containment)
    .def_readonly("left_containment_ani",
                  &Sketch::Query::Result::left_containment_ani)
    .def_readonly("right_containment_ani",
                  &Sketch::Query::Result::right_containment_ani)
    .def_readonly("max_containment_ani",
                  &Sketch::Query::Result::max_containment_ani)
    .def_readonly("weighted_jaccard",
                  &Sketch::Query::Result::weighted_jaccard)
    .def_readonly("weighted_left_containment",
                  &Sketch::Query::Result::weighted_left_containment)
    .def_readonly("weighted_right_containment",
                  &Sketch::Query::Result::weighted_right_containment)
    .def_readonly("weighted_intersection",
                  &Sketch::Query::Result::weighted_intersection)
    .def_readonly("weighted_union",
                  &Sketch::Query::Result::weighted_union)
    .def_readonly("order_similarity",
                  &Sketch::Query::Result::order_similarity)
    .def_readonly("order_distance",
                  &Sketch::Query::Result::order_distance)
    .def_readonly("left_weight", &Sketch::Query::Result::left_weight)
    .def_readonly("right_weight", &Sketch::Query::Result::right_weight)
    .def_readonly("left_cardinality",
                  &Sketch::Query::Result::left_cardinality)
    .def_readonly("right_cardinality",
                  &Sketch::Query::Result::right_cardinality)
    .def_readonly("effective_samples",
                  &Sketch::Query::Result::effective_samples)
    .def_readonly("confidence_level",
                  &Sketch::Query::Result::confidence_level)
    .def_readonly("plugin_standard_error",
                  &Sketch::Query::Result::plugin_standard_error)
    .def_readonly("approximate_ci95_lower",
                  &Sketch::Query::Result::approximate_ci95_lower)
    .def_readonly("approximate_ci95_upper",
                  &Sketch::Query::Result::approximate_ci95_upper)
    .def_readonly("normal_ci95_eligible",
                  &Sketch::Query::Result::normal_ci95_eligible)
    .def_readonly("warnings", &Sketch::Query::Result::warnings);

  py::class_<Sketch::API::BuiltSketch>(m, "BuiltSketch")
    .def_property_readonly("algorithm",
                           &Sketch::API::BuiltSketch::algorithm)
    .def_property_readonly(
        "metadata", &Sketch::API::BuiltSketch::metadata,
        py::return_value_policy::reference_internal)
    .def_property_readonly("descriptor",
                           &Sketch::API::BuiltSketch::descriptor)
    .def("query", &Sketch::API::BuiltSketch::query,
         py::arg("other"),
         py::arg("request") = Sketch::Query::Request(),
         py::call_guard<py::gil_scoped_release>());

  py::class_<Sketch::API::BuildResult>(m, "BuildResult")
    .def_readonly("sketch", &Sketch::API::BuildResult::sketch)
    .def_readonly("label", &Sketch::API::BuildResult::label)
    .def_readonly("config", &Sketch::API::BuildResult::config)
    .def_readonly("stats", &Sketch::API::BuildResult::stats)
    .def_readonly("source", &Sketch::API::BuildResult::source)
    .def_readonly("record_name", &Sketch::API::BuildResult::record_name);

  py::class_<Sketch::API::MultiSketchBuilder>(m, "MultiSketchBuilder")
    .def(py::init<std::vector<Sketch::API::SketchConfig>>(),
         py::arg("configs"))
    .def("update", &Sketch::API::MultiSketchBuilder::update,
         py::arg("record"), py::call_guard<py::gil_scoped_release>())
    .def("finish", &Sketch::API::MultiSketchBuilder::finish,
         py::arg("label_prefix"), py::arg("source") = std::string(),
         py::arg("record_name") = std::string(),
         py::call_guard<py::gil_scoped_release>())
    .def_property_readonly("lanes", &Sketch::API::MultiSketchBuilder::lanes)
    .def_property_readonly("finished",
                           &Sketch::API::MultiSketchBuilder::finished);

  m.def("build_fastx", &Sketch::API::buildFastx,
        py::arg("path"), py::arg("configs"),
        py::call_guard<py::gil_scoped_release>());
  m.def("build_fastx_files", &Sketch::API::buildFastxFiles,
        py::arg("paths"), py::arg("configs"),
        py::call_guard<py::gil_scoped_release>());

  m.def("describe_fastkmv", [](const Sketch::FastKMV& sketch) {
    return Sketch::Query::describe(sketch);
  });
  m.def("describe_bindash", [](const Sketch::BinDash& sketch) {
    return Sketch::Query::describe(sketch);
  });
  m.def("describe_probminhash", [](const Sketch::ProbMinHash4& sketch) {
    return Sketch::Query::describe(sketch);
  });
  m.def("describe_legacy_minhash", [](const Sketch::MinHash& sketch) {
    return Sketch::Query::describe(sketch);
  });
  m.def("describe_fracminhash", [](const Sketch::FracMinHash& sketch) {
    return Sketch::Query::describe(sketch);
  });
  m.def("describe_orderminhash", [](const Sketch::OrderMinHash& sketch) {
    return Sketch::Query::describe(sketch);
  });
  m.def("describe_hll", [](const Sketch::HyperLogLog& sketch) {
    return Sketch::Query::describe(sketch);
  });
  m.def("describe_setsketch", [](const Sketch::SetSketch& sketch) {
    return Sketch::Query::describe(sketch);
  });
  m.def("describe_kssd", [](const Sketch::Kssd& sketch) {
    return Sketch::Query::describe(sketch);
  });
  m.def("plan_query", &Sketch::Query::planQuery,
        py::arg("left"), py::arg("right"),
        py::arg("request") = Sketch::Query::Request());

  m.def("query", [](const Sketch::FastKMV& left,
                     const Sketch::FastKMV& right,
                     const Sketch::Query::Request& request) {
    py::gil_scoped_release release;
    return Sketch::Query::query(left, right, request);
  }, py::arg("left"), py::arg("right"),
     py::arg("request") = Sketch::Query::Request());
  m.def("query", [](const Sketch::BinDash& left,
                     const Sketch::BinDash& right,
                     const Sketch::Query::Request& request) {
    py::gil_scoped_release release;
    return Sketch::Query::query(left, right, request);
  }, py::arg("left"), py::arg("right"),
     py::arg("request") = Sketch::Query::Request());
  m.def("query", [](const Sketch::ProbMinHash4& left,
                     const Sketch::ProbMinHash4& right,
                     const Sketch::Query::Request& request) {
    py::gil_scoped_release release;
    return Sketch::Query::query(left, right, request);
  }, py::arg("left"), py::arg("right"),
     py::arg("request") = Sketch::Query::Request());
  m.def("query", [](Sketch::MinHash& left,
                     Sketch::MinHash& right,
                     const Sketch::Query::Request& request) {
    py::gil_scoped_release release;
    return Sketch::Query::query(left, right, request);
  }, py::arg("left"), py::arg("right"),
     py::arg("request") = Sketch::Query::Request());
  m.def("query", [](const Sketch::FracMinHash& left,
                     const Sketch::FracMinHash& right,
                     const Sketch::Query::Request& request) {
    py::gil_scoped_release release;
    return Sketch::Query::query(left, right, request);
  }, py::arg("left"), py::arg("right"),
     py::arg("request") = Sketch::Query::Request());
  m.def("query", [](const Sketch::HyperLogLog& left,
                     const Sketch::HyperLogLog& right,
                     const Sketch::Query::Request& request) {
    py::gil_scoped_release release;
    return Sketch::Query::query(left, right, request);
  }, py::arg("left"), py::arg("right"),
     py::arg("request") = Sketch::Query::Request());
  m.def("query", [](const Sketch::SetSketch& left,
                     const Sketch::SetSketch& right,
                     const Sketch::Query::Request& request) {
    py::gil_scoped_release release;
    return Sketch::Query::query(left, right, request);
  }, py::arg("left"), py::arg("right"),
     py::arg("request") = Sketch::Query::Request());
  m.def("query", [](Sketch::Kssd& left, Sketch::Kssd& right,
                     const Sketch::Query::Request& request) {
    py::gil_scoped_release release;
    return Sketch::Query::query(left, right, request);
  }, py::arg("left"), py::arg("right"),
     py::arg("request") = Sketch::Query::Request());
  m.def("query", [](const Sketch::OrderMinHash& left,
                     const Sketch::OrderMinHash& right,
                     const Sketch::Query::Request& request) {
    py::gil_scoped_release release;
    return Sketch::Query::query(left, right, request);
  }, py::arg("left"), py::arg("right"),
     py::arg("request") = Sketch::Query::Request());
  py::class_<Sketch::sketchInfo_t>(m, "SketchInfo")
    .def(py::init<>())  
    .def_readwrite("id", &Sketch::sketchInfo_t::id)
    .def_readwrite("half_k", &Sketch::sketchInfo_t::half_k)
    .def_readwrite("half_subk", &Sketch::sketchInfo_t::half_subk)
    .def_readwrite("drlevel", &Sketch::sketchInfo_t::drlevel)
    .def_readwrite("genomeNumber", &Sketch::sketchInfo_t::genomeNumber)
    .def("__repr__", [](const Sketch::sketchInfo_t &info) {
        return "<SketchInfo id=" + std::to_string(info.id) +
        ", half_k=" + std::to_string(info.half_k) +
        ", half_subk=" + std::to_string(info.half_subk) +
        ", drlevel=" + std::to_string(info.drlevel) +
        ", genomeNumber=" + std::to_string(info.genomeNumber) + ">";
        });



    py::class_<Sketch::KssdLite>(m, "KssdLite")
        .def(py::init<>())
        .def_readwrite("fileName", &Sketch::KssdLite::fileName)
        .def_readwrite("id", &Sketch::KssdLite::id)
        .def_readwrite("hashList", &Sketch::KssdLite::hashList)
        .def_readwrite("hashList64", &Sketch::KssdLite::hashList64)
        .def("__getstate__", [](const Sketch::KssdLite &self) {
            return py::make_tuple(self.fileName, self.id, self.hashList, self.hashList64);
        })
        .def("__setstate__", [](Sketch::KssdLite &self, py::tuple t) {
            if (t.size() != 4) {
                throw std::runtime_error("Invalid state for KssdLite");
            }
            new (&self) Sketch::KssdLite(); 
            self.fileName = t[0].cast<std::string>();
            self.id = t[1].cast<int>();
            self.hashList = t[2].cast<std::vector<uint32_t>>();
            self.hashList64 = t[3].cast<std::vector<uint64_t>>();
        });


}

#endif // __PYBIND_H__
