// Diagnostic only: isolate legacy source sampling and collideWithFull stages.
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "Collisions.h"
#include "Grid.h"
#include "NegativeParticleCollisions.h"
#include "NegativeParticleSampling.h"
#include "ParticleGroup.h"
#include "RandomSampling.h"
#include "SimulationConfig.h"
using namespace coulomb;

std::array<double,5> moments(const NeParticleGroup& group, double weight) {
    std::array<double,5> result{};
    for (auto kind : {ParticleKind::Positive,ParticleKind::Negative})
        for (const auto& particle : group.list(kind)) {
            const double w=kind==ParticleKind::Positive?weight:-weight;
            result[0]+=w;
            for(int j=0;j<3;++j) {
                result[1+j]+=w*particle.velocity(j);
                result[4]+=w*particle.velocity(j)*particle.velocity(j);
            }
        }
    return result;
}

int main(int argc,char**argv) {
    try {
        if(argc<5) throw std::runtime_error("kernel input output dt; or stages output dt replicas");
        ParaClass parameters;
        parameters.dt=std::stod(argv[4]); parameters.coeffBinaryColl=5.;
        NeParticleGroup base; base.rhoM=base.rhoF=1.; base.tprtM=1.; base.setXRange(0.,1.);
        NegativeParticleSampling sampler;
        if(std::string(argv[1])=="kernel") {
            std::ifstream input(argv[2]); std::ofstream output(argv[3]);
            output<<std::setprecision(17)<<"source,maxwellian\n";
            double x,y,z,a,b,c;
            while(input>>x>>y>>z>>a>>b>>c)
                output<<sampler.evaluateSource({x,y,z},{a,b,c},base,parameters)<<','
                      <<sampler.evaluateMaxwellian({x,y,z},base)<<'\n';
            return 0;
        }
        // stages output dt replicas: unlike kernel, argument 4 is replicas.
        parameters.dt=std::stod(argv[3]);
        int replicas=std::stoi(argv[4]);
        NumericGridClass grid(1);grid.dx=1.;grid.neff=.3/128;grid.neffF=1./2048;
        grid.dt=parameters.dt;
        sampler.updateBounds(base,parameters);
        std::ofstream output(argv[2]);output<<std::setprecision(17);
        output<<"replica,fresh_cache,source_mass,source_px,source_py,source_pz,source_v2,transport_mass,transport_px,transport_py,transport_pz,transport_v2,added_positive,added_negative\n";
        for(int r=0;r<replicas;++r) {
            RandomContext initial;initial.reseed(2100000000u+static_cast<unsigned>(r)*7919u);
            RandomSampling draw(initial);
            auto group=base;
            auto sample=[&](bool b) {return Particle1D3D({draw.normal()*std::sqrt(b?2.:1.),draw.normal()*std::sqrt(b?.5:1.),draw.normal()*std::sqrt(b?.5:1.)});};
            for(int j=0;j<2048;++j)group.pushBack(sample(draw.uniform()<.3),ParticleKind::Full);
            for(int j=0;j<128;++j) {group.pushBack(sample(true),ParticleKind::Positive);group.pushBack(sample(false),ParticleKind::Negative);}
            for(int fresh=0;fresh<2;++fresh) {
                auto test=group;
                if(fresh)test.computeMoments();
                RandomContext random;random.reseed(2300000000u+static_cast<unsigned>(r)*7919u);
                NeParticleGroup source;
                sampler.sampleDelta(test,source,parameters,grid.neff,random);
                auto sm=moments(source,grid.neff),before=moments(test,grid.neff);
                NegativeParticleCollisions(grid,parameters,random).collideWithFull(test);
                auto after=moments(test,grid.neff);
                output<<r<<','<<fresh;
                for(double v:sm)output<<','<<v;
                for(int j=0;j<5;++j)output<<','<<after[j]-before[j];
                output<<','<<source.size(ParticleKind::Positive)<<','<<source.size(ParticleKind::Negative)<<'\n';
            }
        }
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 2;}
}
