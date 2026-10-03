// Metric/workload adapter dispatched by the existing RabbitSketch main -t.
#include "Sketch.h"
#include "fastkmv.h"
#include "fracminhash.h"
#include "BinDash.h"
#include "probmh.h"
#include "api/RuntimeInfo.h"
#include "../screening/real_input.h"
#include <omp.h>
#include <zlib.h>
#include <sys/resource.h>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <unordered_map>

namespace {
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
uint64_t rss(){rusage r{};getrusage(RUSAGE_SELF,&r);return uint64_t(r.ru_maxrss)*1024;}
std::string quote(const std::string& s){std::string o="\"";for(char c:s){if(c=='"'||c=='\\')o+='\\';if(c=='\n')o+="\\n";else o+=c;}return o+'"';}
struct Row{
 std::map<std::string,std::string> f;
 void s(std::string k,std::string v){f[k]=quote(v);}
 void n(std::string k,double v){std::ostringstream o;o<<std::setprecision(17)<<v;f[k]=std::isfinite(v)?o.str():"null";}
 void u(std::string k,uint64_t v){f[k]=std::to_string(v);}
 void emit(){bool first=true;std::cout<<'{';for(auto& kv:f){if(!first)std::cout<<',';first=false;std::cout<<quote(kv.first)<<':'<<kv.second;}std::cout<<"}\n";}
};
std::vector<std::string> fields(const std::string& line){std::vector<std::string> v;std::istringstream in(line);std::string s;while(std::getline(in,s,'\t')){if(!s.empty()&&s.back()=='\r')s.pop_back();v.push_back(s);}return v;}
struct Object{std::string id,path,sequence;std::vector<std::pair<uint64_t,double>> weights;uint64_t digest=0;};
struct Pair{std::string id;size_t a,b;};
struct BinDashWithHll{
 Sketch::BinDash bins;Sketch::HyperLogLog cardinal;
 BinDashWithHll(int p,int k,uint64_t seed):bins(p/64,k,16,seed),cardinal(p==384?8:p==1536?10:12,k,seed,false){}
};
std::string fasta(const std::string& path){
 if(path.size()<3||path.substr(path.size()-3)!=".gz")return readScreeningFasta(path);
 gzFile f=gzopen(path.c_str(),"rb");if(!f)throw std::runtime_error("cannot open "+path);
 std::string raw;char buf[65536];int n;while((n=gzread(f,buf,sizeof(buf)))>0)raw.append(buf,n);
 if(n<0){gzclose(f);throw std::runtime_error("gzip read failed");}gzclose(f);
 std::string result;result.reserve(raw.size());bool start=true,header=false,record=false;
 for(unsigned char c:raw){if(start&&c=='>'){header=true;record=true;result+='N';}if(c=='\n'){start=true;header=false;continue;}start=false;if(header||std::isspace(c))continue;if(!record)throw std::runtime_error("FASTA header missing");c=std::toupper(c);result+=(c=='A'||c=='C'||c=='G'||c=='T')?c:'N';}
 if(!record)throw std::runtime_error("empty FASTA");return result;
}
uint64_t digest(const std::string& s){uint64_t h=14695981039346656037ULL;for(unsigned char c:s){h^=c;h*=1099511628211ULL;}return h;}
struct Config{
 std::string algorithm,mode;int p,k,threads,passes,half_subk;uint64_t seed;
 std::vector<Object> objects;std::vector<Pair> pairs;
 Row row(const char* type)const{Row r;r.s("type",type);r.s("algorithm",algorithm);r.s("mode",mode);r.u("parameter",p);r.u("k",k);r.u("seed",seed);return r;}
};
double derived(double j,double a,double b){if(!(a>0)||!std::isfinite(a)||!std::isfinite(b))return NAN;return j*(a+b)/(a*(1+j));}

template<class Make,class Update,class Finalize,class Bytes,class Empty,class Card,class Query,class Contain>
void measured(Config& c,Make make,Update update,Finalize finalize,Bytes bytes,Empty empty,Card card,Query query,Contain contain){
 using Ptr=decltype(make());
 const size_t n=c.objects.size();std::vector<Ptr> sketches(n);
 std::vector<double> build(n),cards(n),cardtime(n);std::vector<uint64_t> payload(n),empties(n);
 std::string error;auto start=Clock::now();
 #pragma omp parallel for schedule(dynamic,1) num_threads(c.threads)
 for(size_t i=0;i<n;++i){try{auto t=Clock::now();sketches[i]=make();update(*sketches[i],c.objects[i]);finalize(*sketches[i]);build[i]=elapsed(t);payload[i]=bytes(*sketches[i]);empties[i]=empty(*sketches[i]);}catch(const std::exception& e){
 #pragma omp critical
 {error=e.what();}
 }}
 const double build_wall=elapsed(start);if(!error.empty())throw std::runtime_error(error);
 // Materialize mutable cardinality caches before any concurrent reads.
 start=Clock::now();
 #pragma omp parallel for schedule(dynamic,1) num_threads(c.threads)
 for(size_t i=0;i<n;++i){auto t=Clock::now();cards[i]=card(*sketches[i]);cardtime[i]=elapsed(t);}
 const double prepare_wall=elapsed(start);
 for(size_t i=0;i<n;++i){auto r=c.row("object");r.s("object_id",c.objects[i].id);r.u("input_digest",c.objects[i].digest);r.u("input_bytes",c.objects[i].sequence.size());r.u("payload_bytes",payload[i]);r.u("empty",empties[i]);r.n("cardinality",cards[i]);r.n("build_seconds",build[i]);r.n("cardinality_seconds",cardtime[i]);r.emit();}
 auto r=c.row("build");r.u("objects",n);r.u("threads",c.threads);r.n("wall_seconds",build_wall);r.n("prepare_query_wall_seconds",prepare_wall);r.u("payload_bytes",std::accumulate(payload.begin(),payload.end(),uint64_t(0)));r.u("process_peak_rss_bytes",rss());r.emit();std::cout<<std::flush;
 if(c.mode=="allpairs"){
  for(int pass=0;pass<c.passes;++pass){double checksum=0;uint64_t invalid=0;auto t=Clock::now();
   #pragma omp parallel for schedule(dynamic,1) num_threads(c.threads) reduction(+:checksum,invalid)
   for(size_t i=0;i<n;++i)for(size_t j=i+1;j<n;++j){double v=query(*sketches[i],*sketches[j]);if(!std::isfinite(v)||v<0||v>1)++invalid;else checksum+=v;}
   auto q=c.row("allpairs");q.u("objects",n);q.u("pairs",n*(n-1)/2);q.u("pass",pass);q.n("wall_seconds",elapsed(t));q.n("checksum",checksum);q.u("invalid",invalid);q.u("threads",c.threads);q.emit();std::cout<<std::flush;
  }
 }else{
  struct Result{double j,cab,cba,nab,nba,t;};std::vector<Result> results(c.pairs.size());start=Clock::now();
  #pragma omp parallel for schedule(dynamic,1) num_threads(c.threads)
  for(size_t i=0;i<c.pairs.size();++i){auto p=c.pairs[i];auto& a=*sketches[p.a];auto& b=*sketches[p.b];auto t=Clock::now();double j=query(a,b);double qt=elapsed(t);results[i]={j,derived(j,cards[p.a],cards[p.b]),derived(j,cards[p.b],cards[p.a]),contain(a,b),contain(b,a),qt};}
  const double pairwall=elapsed(start);
  for(size_t i=0;i<c.pairs.size();++i){auto p=c.pairs[i];auto v=results[i];auto q=c.row("pair");q.s("pair_id",p.id);q.s("a",c.objects[p.a].id);q.s("b",c.objects[p.b].id);q.n("jaccard",v.j);q.n("containment_derived_ab",v.cab);q.n("containment_derived_ba",v.cba);q.n("containment_native_ab",v.nab);q.n("containment_native_ba",v.nba);q.n("query_seconds",v.t);q.u("payload_bytes_a",payload[p.a]);q.u("payload_bytes_b",payload[p.b]);q.u("empty_a",empties[p.a]);q.u("empty_b",empties[p.b]);q.emit();}
  auto q=c.row("pairs_summary");q.n("wall_seconds_including_containment",pairwall);q.u("pairs",c.pairs.size());q.emit();
 }
}
void exact(Config& c){
 std::vector<std::vector<uint64_t>> words(c.objects.size());auto t=Clock::now();
 #pragma omp parallel for schedule(dynamic,1) num_threads(c.threads)
 for(size_t i=0;i<words.size();++i)words[i]=exactScreeningWords(c.objects[i].sequence,c.k);
 for(size_t i=0;i<words.size();++i){auto r=c.row("object");r.s("object_id",c.objects[i].id);r.u("n",words[i].size());r.u("input_digest",c.objects[i].digest);r.emit();}
 std::vector<uint64_t> common(c.pairs.size());
 #pragma omp parallel for schedule(dynamic,1) num_threads(c.threads)
 for(size_t i=0;i<c.pairs.size();++i){auto p=c.pairs[i];auto& a=words[p.a];auto& b=words[p.b];if(a.size()*30<b.size()){uint64_t s=0;for(auto w:a)s+=std::binary_search(b.begin(),b.end(),w);common[i]=s;}else common[i]=exactScreeningIntersection(a,b);}
 for(size_t i=0;i<c.pairs.size();++i){auto p=c.pairs[i];auto r=c.row("pair");r.s("pair_id",p.id);r.u("n_a",words[p.a].size());r.u("n_b",words[p.b].size());r.u("intersection",common[i]);r.emit();}
 auto r=c.row("exact_summary");r.n("wall_seconds",elapsed(t));r.u("process_peak_rss_bytes",rss());r.emit();
}
}

int run_regimes_manifest(const std::string& path,int threads){try{
 std::ifstream in(path);std::string line;std::getline(in,line);if(line!="#regimes-v1")throw std::runtime_error("bad regimes marker");
 std::map<std::string,std::string> kv;while(std::getline(in,line)){if(line.empty()||line[0]=='#')continue;auto v=fields(line);if(v.size()!=2||!kv.emplace(v[0],v[1]).second)throw std::runtime_error("bad config row");}
 Config c;c.algorithm=kv.at("algorithm");c.mode=kv.at("mode");c.p=std::stoi(kv.at("parameter"));c.k=std::stoi(kv.at("k"));c.seed=std::stoull(kv.at("seed"));c.threads=threads;c.passes=kv.count("passes")?std::stoi(kv.at("passes")):3;c.half_subk=kv.count("half_subk")?std::stoi(kv.at("half_subk")):6;
 if(threads<1||c.k<1||c.k>32||c.p<1||c.passes<1)throw std::runtime_error("invalid configuration");
 std::ifstream objects(kv.at("objects"));if(!objects)throw std::runtime_error("objects missing");std::getline(objects,line);std::unordered_map<std::string,size_t> ids;
 while(std::getline(objects,line)){auto v=fields(line);if(v.size()<2||ids.count(v[0]))throw std::runtime_error("bad object row");ids[v[0]]=c.objects.size();c.objects.push_back({v[0],v[1],{}, {},0});}
 if(kv.count("limit")){auto n=std::stoull(kv.at("limit"));if(n>c.objects.size())throw std::runtime_error("limit exceeds inputs");c.objects.resize(n);}
 if(kv.at("pairs")!="-"){std::ifstream pairs(kv.at("pairs"));if(!pairs)throw std::runtime_error("pairs missing");std::getline(pairs,line);while(std::getline(pairs,line)){auto v=fields(line);if(v.size()<3)throw std::runtime_error("bad pair");auto a=ids.at(v[1]),b=ids.at(v[2]);if(a>=c.objects.size()||b>=c.objects.size())throw std::runtime_error("pair outside input limit");c.pairs.push_back({v[0],a,b});}}
 omp_set_dynamic(0);std::string error;auto start=Clock::now();
 #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
 for(size_t i=0;i<c.objects.size();++i){try{auto& o=c.objects[i];if(c.mode=="weighted"){
  std::ifstream w(o.path);if(!w)throw std::runtime_error("weights missing");uint64_t key;double weight;while(w>>key>>weight){if(!(weight>0)||!std::isfinite(weight))throw std::runtime_error("invalid weight");o.weights.emplace_back(key,weight);}std::sort(o.weights.begin(),o.weights.end());for(size_t j=1;j<o.weights.size();++j)if(o.weights[j-1].first==o.weights[j].first)throw std::runtime_error("duplicate weighted key");
 }else{o.sequence=fasta(o.path);o.digest=digest(o.sequence);}}catch(const std::exception& e){
 #pragma omp critical
 {error=e.what();}
 }}
 if(!error.empty())throw std::runtime_error(error);auto read=c.row("input");read.u("objects",c.objects.size());read.u("threads",threads);read.n("wall_seconds",elapsed(start));read.emit();
 if(c.algorithm=="exact"){exact(c);return 0;}
 auto noop=[](auto&){};auto finish=[](auto& s){s.finalize();};auto up=[](auto& s,Object& o){s.update(o.sequence.data(),o.sequence.size());};auto upstr=[](auto& s,Object& o){s.update(o.sequence.data());};
 auto card=[](auto& s){return s.cardinality();};auto query=[](auto& a,auto& b){return a.jaccard(b);};auto contain=[](auto& a,auto& b){return a.containment(b);};auto notempty=[](auto&){return false;};
 if(c.algorithm=="minhash")measured(c,[&]{return std::make_unique<Sketch::FastKMV>(c.p,c.k,c.seed);},up,finish,[](auto& s){return s.size()*8;},[](auto& s){return s.size()==0;},card,query,contain);
 else if(c.algorithm=="fracminhash")measured(c,[&]{return std::make_unique<Sketch::FracMinHash>(c.p,c.k,c.seed);},up,finish,[](auto& s){return s.size()*8;},[](auto& s){return s.size()==0;},card,query,contain);
 else if(c.algorithm=="bindash")measured(c,[&]{return std::make_unique<Sketch::BinDash>(c.p/64,c.k,16,c.seed);},up,finish,[](auto& s){return s.memoryBytes();},[](auto& s){return s.getRawNonempty()==0;},card,query,contain);
 else if(c.algorithm=="bindash_hll"){
  if(c.p!=384&&c.p!=1536&&c.p!=6144)throw std::runtime_error("invalid BinDash+HLL configuration");
  measured(c,[&]{return std::make_unique<BinDashWithHll>(c.p,c.k,c.seed);},[](auto& s,Object& o){s.bins.update(o.sequence.data(),o.sequence.size());s.cardinal.update(o.sequence.data());},[](auto& s){s.bins.finalize();},[](auto& s){return s.bins.memoryBytes()+s.cardinal.getCore().size();},[](auto& s){return s.bins.getRawNonempty()==0;},[](auto& s){return s.cardinal.cardinality();},[](auto& a,auto& b){return a.bins.jaccard(b.bins);},[](auto& a,auto& b){return derived(a.bins.jaccard(b.bins),a.cardinal.cardinality(),b.cardinal.cardinality());});
 }
 else if(c.algorithm=="hll")measured(c,[&]{return std::make_unique<Sketch::HyperLogLog>(c.p,c.k,c.seed,false);},upstr,noop,[](auto& s){return s.getCore().size();},notempty,card,[](auto& a,auto& b){return a.jaccard_index(b);},contain);
 else if(c.algorithm=="setsketch")measured(c,[&]{return std::make_unique<Sketch::SetSketch>(c.p,2.,5.,c.k,false,c.seed);},upstr,noop,[](auto& s){return s.getCore().size();},notempty,card,[](auto& a,auto& b){return a.jaccard_index(b);},contain);
 else if(c.algorithm=="kssd"){
  if(c.k%2||c.half_subk<3||c.half_subk>7||c.half_subk>c.k/2||c.p>c.half_subk-3)throw std::runtime_error("invalid KSSD resolution");start=Clock::now();auto params=std::make_shared<Sketch::kssd_parameter_t>(c.k/2,c.half_subk,c.p);auto setup=c.row("shared_setup");setup.n("wall_seconds",elapsed(start));setup.u("half_subk",c.half_subk);setup.emit();
  auto kc=[](auto& s){return (s.getHashes().size()+s.getHashes64().size())/s.admissionProbability();};
  measured(c,[&]{return std::make_unique<Sketch::Kssd>(params,c.seed);},upstr,finish,[](auto& s){return s.getHashes().size()*4+s.getHashes64().size()*8;},[](auto& s){return s.getHashes().empty()&&s.getHashes64().empty();},kc,[](auto& a,auto& b){return a.jaccard(&b);},[&](auto& a,auto& b){return derived(a.jaccard(&b),kc(a),kc(b));});
 }else if(c.algorithm=="probminhash"){
  measured(c,[&]{return std::make_unique<Sketch::ProbMinHash4>(c.p,c.k,c.seed,0,Sketch::Rank::WeightSemantics::Frequency);},[](auto& s,Object& o){for(auto v:o.weights)s.addHash(v.first,v.second);},finish,[](auto& s){return s.getM()*8;},notempty,[](auto& s){return s.total_weight();},[](auto& a,auto& b){return a.jaccard_weighted(b);},[](auto&,auto&){return NAN;});
 }else throw std::runtime_error("unknown algorithm");
 return 0;
 }catch(const std::exception& e){std::cerr<<"regimes: "<<e.what()<<'\n';return 2;}}
