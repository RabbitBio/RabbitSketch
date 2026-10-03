// Model identity, public-query consistency and sparse-state bias regression.
#include "SetSketch.h"
#include "estimators/PartitionedSetSketch.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

static void require(bool yes, const char* message) {
    if (!yes) throw std::runtime_error(message);
}
static uint64_t fingerprint(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
static void populate(Sketch::SetSketch& a, Sketch::SetSketch& b,
                     uint64_t na,uint64_t nb,uint64_t common,uint64_t seed) {
    for (uint64_t i=0;i<na+nb-common;++i) {
        uint64_t h=fingerprint(i+seed*100000000ULL);
        if (i<common) { a.addFingerprint(h); b.addFingerprint(h); }
        else if (i<na) a.addFingerprint(h);
        else b.addFingerprint(h);
    }
}
int main() {
    try {
        // Independently compute the joint PMF by four CDF corners. Include
        // empty registers, fully shared sets, and the censored upper tail.
        unsigned model_cases=0, api_cases=0;
        for (double base : {1.5,2.0,3.0}) for (double scale : {.7,5.0}) {
            Sketch::detail::PartitionedSetSketchModel model(base,scale,5);
            for (double la : {.001,.1,1.,10.,100.}) for (double lb : {.001,.1,1.,10.,100.})
            for (double fraction : {0.,.2,1.}) {
                const double shared=fraction*std::min(la,lb);
                auto cdf = [&](int ra,int rb) -> long double {
                    if (ra<0 || rb<0) return 0.0L;
                    auto exceed = [&](int r) -> long double {
                        return r==6 ? 0.0L : -std::expm1(-static_cast<long double>(scale)*
                            std::pow(static_cast<long double>(base),-r));
                    };
                    const long double ca=exceed(ra),cb=exceed(rb);
                    return std::exp(-la*ca-lb*cb+shared*std::min(ca,cb));
                };
                double total=0.0;
                for (unsigned ra=0;ra<=6;++ra) for (unsigned rb=0;rb<=6;++rb) {
                    const long double expected=cdf(ra,rb)-cdf(int(ra)-1,rb)-
                        cdf(ra,int(rb)-1)+cdf(int(ra)-1,int(rb)-1);
                    const double measured=std::exp(model.logJoint(ra,rb,la,lb,shared));
                    require(std::isfinite(measured) && measured>=0.0,"invalid joint probability");
                    require(std::abs(measured-expected)<2e-12L,"PMF disagrees with independent joint CDF");
                    total+=measured;
                }
                require(std::abs(total-1.0)<1e-11,"joint PMF does not normalize");
                ++model_cases;
            }
        }
        for (int p : {4,8,12,16}) for (int nb : {100,1000,10000,100000})
        for (double overlap : {0.,.1,1.}) {
            const int na=nb/10,common=int(na*overlap);
            Sketch::SetSketch a(p,2.,5.,20,false,42),b(p,2.,5.,20,false,42);
            populate(a,b,na,nb,common,17);
            double j=a.jaccard_index(b), reverse=b.jaccard_index(a);
            require(std::isfinite(j) && j>=0.0 && j<=1.0,"default J out of bounds");
            require(std::abs(j-reverse)<1e-7,"default J asymmetric");
            require(std::abs(a.jaccard_index(a)-1.0)<1e-12,"identical-state J not one");
            double flat=Sketch::SetSketch::jaccardFromCores(a.getCore().data(),b.getCore().data(),a.getM(),
                a.getBaseInvPow(),a.getFactor(),a.cardinality(),b.cardinality());
            require(flat>=0.0 && flat<=1.0,"flat inclusion-exclusion J out of bounds");
            require(std::abs(flat-a.jaccard_index_inclexcl(b))<1e-8,"flat/OO union estimators disagree");
            double joint=Sketch::SetSketch::jaccardFromCoresMLE(a.getCore().data(),b.getCore().data(),a.getM(),
                a.getBaseInvPow(),a.getFactor(),a.cardinality(),b.cardinality(),a.getBase(),a.getQ());
            require(std::abs(joint-j)<1e-12,"flat/OO joint estimators disagree");
            for (double threshold : {0.,.001,.01,.1}) {
                double screened=Sketch::SetSketch::jaccardFromCoresBatch(a.getCore().data(),b.getCore().data(),
                    a.getM(),a.getBaseInvPow(),a.getFactor(),a.cardinality(),b.cardinality(),threshold,nullptr,nullptr,16);
                require(flat<threshold || std::abs(screened-flat)<1e-12,"screening false negative");
            }
            ++api_cases;
        }
        double error=0.0;
        for (int rep=0;rep<20;++rep) {
            Sketch::SetSketch a(16,2.,5.,20,false,42),b(16,2.,5.,20,false,42);
            populate(a,b,10000,100000,110,100+rep);
            error+=std::abs(a.jaccard_index(b)-110.0/109890.0);
        }
        error/=20;
        require(error<.003,"sparse-state bias regression");
        std::cout << "{\"model_cases\":" << model_cases << ",\"api_cases\":" << api_cases
                  << ",\"sparse_replicates\":20,\"sparse_mae\":" << error << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
