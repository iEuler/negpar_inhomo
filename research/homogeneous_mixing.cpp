// Collision-only homogeneous experiment. No advection, projection,
// synchronization, weight adaptation or resampling code is called.
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "Collisions.h"
#include "Grid.h"
#include "NegativeParticleCollisions.h"
#include "NegativeParticleSampling.h"
#include "ParticleGroup.h"
#include "RandomSampling.h"
#include "SimulationConfig.h"

using namespace coulomb;
using Clock = std::chrono::steady_clock;
constexpr int observableCount = 8;
using Values = std::array<double, observableCount>;

Values observable(const std::vector<double>& v) {
    const double x=v[0], y=v[1], z=v[2];
    return {x*x-.5*(y*y+z*z), x*x*x*x, std::cos(x)-std::cos(y),
            1., x, y, z, x*x+y*y+z*z};
}

Values particleEstimate(const NeParticleGroup& group, ParticleKind kind, double weight) {
    Values result{};
    for (const auto& particle : group.list(kind)) {
        const auto value=observable(particle.velocity());
        for (int j=0;j<observableCount;++j) result[j]+=weight*value[j];
    }
    return result;
}

int main(int argc, char** argv) {
    try {
        int replicas=1, steps=20, nf=2048, np=128;
        double epsilon=.3, dt=.01, strength=5.;
        std::uint32_t seed=91000;
        std::string mode="hdp", output;
        for (int i=1;i<argc;++i) {
            if (i+1>=argc) throw std::invalid_argument("Every option requires a value");
            const std::string key=argv[i], value=argv[++i];
            if (key=="--replicas") replicas=std::stoi(value);
            else if (key=="--steps") steps=std::stoi(value);
            else if (key=="--full-count") nf=std::stoi(value);
            else if (key=="--sign-count") np=std::stoi(value);
            else if (key=="--epsilon") epsilon=std::stod(value);
            else if (key=="--dt") dt=std::stod(value);
            else if (key=="--strength") strength=std::stod(value);
            else if (key=="--seed") seed=static_cast<std::uint32_t>(std::stoull(value));
            else if (key=="--mode") mode=value;
            else if (key=="--output") output=value;
            else throw std::invalid_argument("Unknown option: "+key);
        }
        if (replicas<1 || steps<0 || nf<2 || nf%2 || np<1 ||
            !std::isfinite(epsilon) || epsilon<=0 || epsilon>=1 ||
            !std::isfinite(dt) || dt<=0 || !std::isfinite(strength) || strength<0 ||
            (mode!="hdp" && mode!="pic") || output.empty() || (mode=="hdp" && nf<2*np))
            throw std::invalid_argument("Invalid homogeneous experiment arguments");
        if (std::filesystem::exists(output)) throw std::invalid_argument("Output already exists");
        std::ofstream file(output);
        if (!file) throw std::runtime_error("Cannot open output");
        file<<std::setprecision(17)<<"replica,step,time";
        for (const auto prefix : {"full", "signed", "mixed", "pic_cv", "mixed_cv"})
            for (int j=0;j<observableCount;++j) file<<','<<prefix<<'_'<<j;
        file<<",positive_count,negative_count,full_count,count_proxy_weight,initialization_seconds,collision_seconds,diagnostic_seconds\n";
        for (int replica=0;replica<replicas;++replica) {
            const auto start=Clock::now();
            RandomContext random;
            random.reseed(seed+static_cast<std::uint32_t>(replica)*7919u);
            RandomSampling sampling(random);
            NumericGridClass grid(1);
            grid.dx=1.; grid.neff=epsilon/np; grid.neffF=1./nf; grid.dt=dt;
            ParaClass parameters;
            parameters.dt=dt; parameters.coeffBinaryColl=strength;
            parameters.flagUseOpenMp=false;
            std::vector<NeParticleGroup> groups(1);
            auto& group=groups[0];
            group.rhoM=1.; group.rhoF=1.; group.rho=1.;
            group.tprtM=1.; group.t1M=group.t2M=group.t3M=1.;
            group.setXRange(0.,1.);
            const auto sample=[&](bool anisotropic) {
                const double sx=anisotropic?std::sqrt(2.):1.;
                const double sy=anisotropic?std::sqrt(.5):1.;
                return Particle1D3D({sx*sampling.normal(),sy*sampling.normal(),sy*sampling.normal()});
            };
            for (int j=0;j<nf;++j) group.pushBack(sample(sampling.uniform()<epsilon),ParticleKind::Full);
            if (mode=="hdp") {
                for (int j=0;j<np;++j) group.pushBack(sample(true),ParticleKind::Positive);
                for (int j=0;j<np;++j) group.pushBack(sample(false),ParticleKind::Negative);
                // M, full density, dt and collision coefficient stay fixed, so
                // the source bounds are constant throughout this experiment.
                if (strength>0) NegativeParticleSampling{}.updateBounds(groups,grid,parameters);
            }
            const double initialization=std::chrono::duration<double>(Clock::now()-start).count();
            double collisionTime=0., diagnosticTime=0.;
            Values initialFull{},initialMixed{};
            const Values exactInitial{1.5*epsilon,3.+9.*epsilon,
                epsilon*(std::exp(-1.)-std::exp(-.25)),1.,0.,0.,0.,3.};
            for (int step=0;step<=steps;++step) {
                auto diagStart=Clock::now();
                const auto full=particleEstimate(group,ParticleKind::Full,grid.neffF);
                auto signedValue=particleEstimate(group,ParticleKind::Positive,grid.neff);
                const auto negative=particleEstimate(group,ParticleKind::Negative,grid.neff);
                const Values maxwellian{0.,3.,0.,1.,0.,0.,0.,3.};
                for (int j=0;j<observableCount;++j) signedValue[j]+=maxwellian[j]-negative[j];
                const int nd=group.size(ParticleKind::Positive)+group.size(ParticleKind::Negative);
                const double vd=grid.neff*grid.neff*nd, vf=grid.neffF*grid.neffF*nf;
                const double weight=vd/(vd+vf);
                Values mixed{},picCv{},mixedCv{};
                for (int j=0;j<observableCount;++j)
                    mixed[j]=mode=="hdp"?weight*full[j]+(1.-weight)*signedValue[j]:full[j];
                if (step==0) {initialFull=full;initialMixed=mixed;}
                for (int j=0;j<observableCount;++j) {
                    picCv[j]=exactInitial[j]+full[j]-initialFull[j];
                    mixedCv[j]=exactInitial[j]+mixed[j]-initialMixed[j];
                }
                diagnosticTime+=std::chrono::duration<double>(Clock::now()-diagStart).count();
                file<<replica<<','<<step<<','<<step*dt;
                for (double value : full) { if (!std::isfinite(value)) throw std::runtime_error("Non-finite full estimate"); file<<','<<value; }
                for (double value : signedValue) { if (!std::isfinite(value)) throw std::runtime_error("Non-finite signed estimate"); file<<','<<value; }
                for (const auto& values : {mixed,picCv,mixedCv})
                    for (double value : values) {if (!std::isfinite(value)) throw std::runtime_error("Non-finite mixed estimate");file<<','<<value;}
                file<<','<<group.size(ParticleKind::Positive)<<','<<group.size(ParticleKind::Negative)<<','<<nf<<','
                    <<vd/(vd+vf)<<','<<initialization<<','<<collisionTime<<','<<diagnosticTime<<'\n';
                if (step==steps) break;
                if (mode=="hdp" && nd>nf)
                    throw std::runtime_error("No-resampling window exhausted: signed count exceeds full count at replica "+std::to_string(replica)+", step "+std::to_string(step));
                const auto collisionStart=Clock::now();
                if (strength>0) {
                    if (mode=="hdp") NegativeParticleCollisions(grid,parameters,random).collideHomogeneous(group);
                    else CollisionOperator(parameters,random).collideHomogeneous(group.list(ParticleKind::Full),nf,1.);
                }
                collisionTime+=std::chrono::duration<double>(Clock::now()-collisionStart).count();
            }
            if ((replica+1)%100==0) std::cout<<"replica "<<replica+1<<'/'<<replicas<<std::endl;
        }
        if (!file) throw std::runtime_error("Output write failed");
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<std::endl;
        return 2;
    }
}
