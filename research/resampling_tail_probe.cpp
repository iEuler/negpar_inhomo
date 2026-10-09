// Isolated core/tail audit. Optional bounded low-moment correction; no collisions.
#include <algorithm>
#include <array>
#include <cstdint>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include "ParticlePartition.h"
#include "BoundedMomentCorrection.h"
#include "ParticleGroupOperations.h"
#include "RandomSampling.h"
#include "Resampler.h"

using namespace coulomb;
using namespace coulomb::resampling;
constexpr std::size_t observableCount=20;
using Values=std::array<double,observableCount>;
constexpr std::array<std::array<double,3>,6> modes{{
    {{1,0,0}},{{0,1,0}},{{0,0,1}},{{1,1,0}},{{2,0,0}},{{0,2,0}}}};

Values estimate(const NeParticleGroup& group,double weight) {
    Values result{};
    for(auto kind:{ParticleKind::Positive,ParticleKind::Negative}) {
        const double w=kind==ParticleKind::Positive?weight:-weight;
        for(const auto& particle:group.list(kind)) {
            const auto& v=particle.velocity();
            const double x=v[0],y=v[1],z=v[2],r2=x*x+y*y+z*z;
            Values a{1.,x,y,z,r2,x*x-.5*(y*y+z*z),x*x*x*x,r2*r2};
            for(std::size_t k=0;k<modes.size();++k) {
                const double phase=modes[k][0]*x+modes[k][1]*y+modes[k][2]*z;
                a[8+2*k]=std::cos(phase); a[9+2*k]=std::sin(phase);
            }
            for(std::size_t j=0;j<result.size();++j) result[j]+=w*a[j];
        }
    }
    return result;
}

// Match the production partial-resampling thinning when weights increase.
void coarsenTail(NeParticleGroup& tail,double inputWeight,double outputWeight,RandomContext& random) {
    if(outputWeight==inputWeight) return;
    if(outputWeight<inputWeight) throw std::invalid_argument("This audit only increases weights");
    for(auto kind:{ParticleKind::Positive,ParticleKind::Negative}) {
        const int remove=RandomSampling(random).stochasticFloor(
            (1.-inputWeight/outputWeight)*tail.size(kind));
        for(int j=0;j<remove;++j) {
            const int index=static_cast<int>(RandomSampling(random).uniform()*tail.size(kind));
            tail.erase(index,kind);
        }
    }
}

int main(int argc,char** argv) {
    try {
        if(argc!=4 && argc!=5) throw std::invalid_argument("Usage: tail_probe output.csv replicas rounds [aligned|stratified|corrected]");
        const bool correctionStudy=argc==5 && std::string(argv[4])=="corrected";
        const bool stratifiedStudy=correctionStudy || (argc==5 && std::string(argv[4])=="stratified");
        const bool alignedStudy=argc==5;
        if(alignedStudy && !stratifiedStudy && std::string(argv[4])!="aligned") throw std::invalid_argument("Unknown study option");
        const int replicas=std::stoi(argv[2]),rounds=std::stoi(argv[3]);
        if(replicas<8 || rounds<1 || rounds>100) throw std::invalid_argument("Need >=8 replicas and 1..100 rounds for distinct seed schedules");
        const std::filesystem::path output(argv[1]);
        if(std::filesystem::exists(output)) throw std::invalid_argument("Output exists");
        std::ofstream out(output); if(!out) throw std::runtime_error("Cannot open output");
        out<<std::setprecision(17)<<"weight_ratio,frequency,cutoff,aligned,stratified,replica,round,seconds,cumulative_seconds,attempts,positive,negative,tail_positive,tail_negative,gate_accept,tail_moment_error,corrected,correction_status,correction_iterations,correction_removed,correction_rms,correction_max,correction_residual";
        for(std::size_t j=0;j<7;++j) out<<",call_target_"<<j<<",correction_moment_"<<j;
        for(std::size_t j=0;j<observableCount;++j) out<<",source_"<<j<<",sample_"<<j;
        out<<'\n';
        RandomContext initial; initial.reseed(610701); RandomSampling sample(initial);
        NeParticleGroup source; source.tprtM=source.t1M=source.t2M=source.t3M=1.;
        for(int i=0;i<512;++i) {
            source.pushBack(Particle1D3D({std::sqrt(2.)*sample.normal(),std::sqrt(.5)*sample.normal(),std::sqrt(.5)*sample.normal()}),ParticleKind::Positive);
            source.pushBack(Particle1D3D({sample.normal(),sample.normal(),sample.normal()}),ParticleKind::Negative);
        }
        std::ofstream sourceFile(output.parent_path()/"source_particles.csv");
        if(!sourceFile) throw std::runtime_error("Cannot archive source particles");
        sourceFile<<std::setprecision(17)<<"sign,vx,vy,vz\n";
        for(auto kind:{ParticleKind::Positive,ParticleKind::Negative})
            for(const auto& p:source.list(kind)) sourceFile<<(kind==ParticleKind::Positive?1:-1)<<','<<p.velocity(0)<<','<<p.velocity(1)<<','<<p.velocity(2)<<'\n';
        const double initialWeight=.05/512.;
        const auto target=estimate(source,initialWeight);
        constexpr std::array<double,4> cutoffs{{0.,2.5,3.,4.}}; // zero means full resampling
        for(int ratio:{1,4}) for(int frequency:{4,8,12}) {
            for(int replica=0;replica<replicas;++replica) {
                // Rotate mode order to limit systematic timing-order effects.
                const int methodCount=alignedStudy?7:4;
                for(int side=0;side<methodCount;++side) {
                    const int method=(replica+side)%methodCount;
                    const bool aligned=stratifiedStudy?method>0:method>=4;
                    const bool stratified=stratifiedStudy && (correctionStudy?method>0:method>=4);
                    const bool corrected=correctionStudy && method>=4;
                    const double cutoff=cutoffs[static_cast<std::size_t>(method>=4?method-3:method)];
                    auto current=source;
                    double cumulativeSeconds=0.;
                    for(int round=1;round<=rounds;++round) {
                        RandomContext random;
                        // Each call runs continuously on one context; modes share the initial seed.
                        random.reseed(static_cast<std::uint32_t>(910000+100*replica+round));
                        const double inputWeight=round==1?initialWeight:ratio*initialWeight;
                        const double outputWeight=ratio*initialWeight;
                        FourierResamplerConfig config;
                        config.frequencyCount=static_cast<std::size_t>(frequency);
                        config.sourceSignedParticleWeight=inputWeight;
                        config.effectiveParticleWeight=outputWeight;
                        config.envelope=ResamplingEnvelope::CertifiedQuadratic;
                        config.cellGeometry=ResamplingCellGeometry::PeriodicWrapped;
                        if(stratified) config.proposalAllocation=ResamplingProposalAllocation::Stratified;
                        if(aligned) {
                            const double radius=cutoff*std::sqrt(current.tprtM);
                            config.fixedVelocityBounds={current.u1M-radius,current.u1M+radius,
                                current.u2M-radius,current.u2M+radius,current.u3M-radius,current.u3M+radius};
                        }
                        const int oldPositive=current.size(ParticleKind::Positive),oldNegative=current.size(ParticleKind::Negative);
                        FourierResamplerDiagnostics diagnostics;
                        const auto start=std::chrono::steady_clock::now();
                        SignedLowMoments correctionTarget{};
                        if(corrected) correctionTarget=BoundedMomentCorrection::moments(current,inputWeight);
                        const auto wholeMomentTarget=correctionTarget;
                        auto partition=cutoff>0?ParticlePartitioning::split(current,cutoff):ParticlePartition{current,NeParticleGroup{}};
                        auto reconstructed=current;
                        for(auto kind:{ParticleKind::Positive,ParticleKind::Negative,ParticleKind::Full}) reconstructed.clear(kind);
                        if(partition.core.size(ParticleKind::Positive)+partition.core.size(ParticleKind::Negative)>0)
                            reconstructed=FourierResampler(partition.core,config).resample(random,&diagnostics);
                        const double reconstructionSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                        const auto tailBefore=estimate(partition.tail,inputWeight);
                        const auto mergeStart=std::chrono::steady_clock::now();
                        coarsenTail(partition.tail,inputWeight,outputWeight,random);
                        const int tailPositive=partition.tail.size(ParticleKind::Positive),tailNegative=partition.tail.size(ParticleKind::Negative);
                        MomentCorrectionResult correction;
                        if(corrected) {
                            const auto tailTarget=BoundedMomentCorrection::moments(partition.tail,outputWeight);
                            for(std::size_t j=0;j<correctionTarget.size();++j) correctionTarget[j]-=tailTarget[j];
                            correction=BoundedMomentCorrection{}.apply(reconstructed,correctionTarget,outputWeight,
                                {current.u1M,current.u2M,current.u3M},cutoff*std::sqrt(current.tprtM),random);
                        }
                        ParticleGroupOperations{}.mergeSigned(reconstructed,partition.tail);
                        reconstructed.tprtM=reconstructed.t1M=reconstructed.t2M=reconstructed.t3M=1.;
                        const double seconds=reconstructionSeconds+std::chrono::duration<double>(std::chrono::steady_clock::now()-mergeStart).count();
                        const auto tailAfter=estimate(partition.tail,outputWeight);
                        double tailError=0.;
                        for(std::size_t j=0;j<observableCount;++j) tailError=std::max(tailError,std::abs(tailAfter[j]-tailBefore[j]));
                        cumulativeSeconds+=seconds;
                        const auto measured=estimate(reconstructed,outputWeight);
                        SignedLowMoments correctedMoments{};
                        if(corrected) correctedMoments=BoundedMomentCorrection::moments(reconstructed,outputWeight);
                        for(double v:measured) if(!std::isfinite(v)) throw std::runtime_error("Nonfinite observable");
                        const bool accept=reconstructed.size(ParticleKind::Positive)<oldPositive && reconstructed.size(ParticleKind::Negative)<oldNegative;
                        out<<ratio<<','<<frequency<<','<<cutoff<<','<<aligned<<','<<stratified<<','<<replica<<','<<round<<','<<seconds<<','<<cumulativeSeconds<<','<<diagnostics.attempts<<','<<reconstructed.size(ParticleKind::Positive)<<','<<reconstructed.size(ParticleKind::Negative)<<','<<tailPositive<<','<<tailNegative<<','<<accept<<','<<tailError<<','<<corrected<<','<<(corrected?static_cast<int>(correction.status):-1)<<','<<correction.iterations<<','<<correction.removed<<','<<correction.rmsDisplacement<<','<<correction.maxDisplacement<<','<<correction.normalizedResidual;
                        for(std::size_t j=0;j<7;++j) out<<','<<wholeMomentTarget[j]<<','<<correctedMoments[j];
                        for(std::size_t j=0;j<observableCount;++j) out<<','<<target[j]<<','<<measured[j];
                        out<<'\n';
                        // Apply every candidate to expose accumulated reconstruction error.
                        // This is deliberately distinct from the production reduction gate.
                        current=std::move(reconstructed);
                    }
                }
            }
            out.flush();
            std::cout<<"ratio "<<ratio<<" frequency "<<frequency<<" complete\n"<<std::flush;
        }
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
