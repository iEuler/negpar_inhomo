// Frozen-population Fourier resampling audit; no collision or moment projection.
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include "Resampler.h"
#include "RandomSampling.h"
using namespace coulomb;
using namespace coulomb::resampling;
using Values=std::array<double,6>;
Values estimate(const NeParticleGroup& group,double weight) {
    Values r{};
    for (auto kind : {ParticleKind::Positive,ParticleKind::Negative}) {
        double w=kind==ParticleKind::Positive?weight:-weight;
        for (const auto& p : group.list(kind)) {
            const auto& v=p.velocity(); double x=v[0],y=v[1],z=v[2];
            Values a{1.,x,x*x+y*y+z*z,x*x-.5*(y*y+z*z),x*x*x*x,std::cos(x)-std::cos(y)};
            for (int j=0;j<6;++j) r[j]+=w*a[j];
        }
    }
    return r;
}
int main(int argc,char** argv) {
    try {
        if (argc!=3 && argc!=4) throw std::invalid_argument("Usage: probe output.csv replicas [geometry]");
        const bool geometry=argc==4;
        if (geometry && std::string(argv[3])!="geometry") throw std::invalid_argument("Unknown probe mode");
        int replicas=std::stoi(argv[2]);
        if (replicas<8) throw std::invalid_argument("Use >=8 replicas");
        std::ifstream existing(argv[1]);
        if (existing.good()) throw std::invalid_argument("Output exists");
        std::ofstream out(argv[1]); if (!out) throw std::runtime_error("Cannot open output");
        out<<std::setprecision(17)<<"frequency,certified,wrapped,replica,seconds,attempts,envelope_increases,positive,negative";
        for(int j=0;j<6;++j) out<<",source_"<<j<<",sample_"<<j;
        out<<'\n';
        RandomContext initial; initial.reseed(610701); RandomSampling sample(initial);
        NeParticleGroup group;
        for(int i=0;i<512;++i) {
            group.pushBack(Particle1D3D({std::sqrt(2.)*sample.normal(),std::sqrt(.5)*sample.normal(),std::sqrt(.5)*sample.normal()}),ParticleKind::Positive);
            group.pushBack(Particle1D3D({sample.normal(),sample.normal(),sample.normal()}),ParticleKind::Negative);
        }
        const double inputWeight=.05/512.,outputWeight=.05/128.;
        auto source=estimate(group,inputWeight);
        for(int frequency : {4,8,12}) {
            for(int replica=0;replica<replicas;++replica) {
                // Alternate order to reduce systematic warm-cache timing differences.
                const int modes=geometry?3:2;
                for(int side=0;side<modes;++side) {
                    int mode=(replica+side)%modes;
                    const bool certified=mode!=0, wrapped=mode==2;
                    FourierResamplerConfig config;
                    config.frequencyCount=frequency;
                    config.effectiveParticleWeight=outputWeight;
                    config.sourceSignedParticleWeight=inputWeight;
                    config.envelope=certified?ResamplingEnvelope::CertifiedQuadratic:ResamplingEnvelope::LegacyAdaptive;
                    config.cellGeometry=wrapped?ResamplingCellGeometry::PeriodicWrapped:ResamplingCellGeometry::LegacyShifted;
                    RandomContext random; random.reseed(710000+replica);
                    FourierResamplerDiagnostics diagnostics;
                    auto start=std::chrono::steady_clock::now();
                    auto result=FourierResampler(group,config).resample(random,&diagnostics);
                    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                    auto measured=estimate(result,outputWeight);
                    out<<frequency<<','<<certified<<','<<wrapped<<','<<replica<<','<<seconds<<','<<diagnostics.attempts<<','<<diagnostics.envelopeIncreases<<','<<result.size(ParticleKind::Positive)<<','<<result.size(ParticleKind::Negative);
                    for(int j=0;j<6;++j) out<<','<<source[j]<<','<<measured[j];
                    out<<'\n';
                }
            }
            std::cout<<"frequency "<<frequency<<" complete\n"<<std::flush;
        }
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
