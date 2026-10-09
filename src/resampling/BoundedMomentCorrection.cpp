#include "BoundedMomentCorrection.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include "RandomSampling.h"

namespace coulomb::resampling {
SignedLowMoments BoundedMomentCorrection::moments(const NeParticleGroup& source,double weight) {
    if (!std::isfinite(weight) || !(weight>0.)) throw std::invalid_argument("Invalid moment weight");
    SignedLowMoments out{};
    out[0]=weight*(static_cast<double>(source.size(ParticleKind::Positive))-source.size(ParticleKind::Negative));
    for(auto kind:{ParticleKind::Positive,ParticleKind::Negative}) {
        const double w=kind==ParticleKind::Positive?weight:-weight;
        for(const auto& p:source.list(kind)) {
            for(std::size_t j=0;j<3;++j) {
                const double v=p.velocity(static_cast<int>(j));
                if(!std::isfinite(v)) throw std::invalid_argument("Nonfinite moment velocity");
                out[1+j]+=w*v; out[4+j]+=w*v*v;
            }
        }
    }
    return out;
}

MomentCorrectionResult BoundedMomentCorrection::apply(NeParticleGroup& core,
    const SignedLowMoments& target,double weight,const std::array<double,3>& center,
    double radius,RandomContext& random,MomentCorrectionConfig config) const {
    if(!std::isfinite(weight) || !(weight>0.) || !std::isfinite(radius) || !(radius>0.)
        || config.maxIterations==0 || !std::isfinite(config.relativeTolerance) || !(config.relativeTolerance>0.)
        || !std::isfinite(config.maxRmsDisplacementFraction) || !(config.maxRmsDisplacementFraction>0.))
        throw std::invalid_argument("Invalid bounded correction configuration");
    for(double x:target) if(!std::isfinite(x)) throw std::invalid_argument("Nonfinite target moment");
    for(double x:center) if(!std::isfinite(x)) throw std::invalid_argument("Nonfinite correction center");
    MomentCorrectionResult result;
    const double massCount=target[0]/weight, desired=std::round(massCount);
    const double massTolerance=std::min(1e-6,128.*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(massCount)));
    if(!std::isfinite(massCount) || std::abs(massCount-desired)>massTolerance
        || std::abs(desired)>std::numeric_limits<int>::max()) {
        result.status=MomentCorrectionStatus::IncompatibleMass; return result;
    }
    auto candidate=core;
    const int np=candidate.size(ParticleKind::Positive),nn=candidate.size(ParticleKind::Negative);
    const double remove=std::abs(static_cast<double>(np)-nn-desired);
    const auto excess=static_cast<double>(np)-nn>desired?ParticleKind::Positive:ParticleKind::Negative;
    if(remove>candidate.size(excess)) {
        result.status=MomentCorrectionStatus::InsufficientParticles; return result;
    }
    result.removed=static_cast<int>(remove);
    for(int k=0;k<result.removed;++k) {
        const int index=static_cast<int>(RandomSampling(random).uniform()*candidate.size(excess));
        candidate.erase(index,excess);
    }
    struct Node { std::array<double,3> x{},original{}; double sign; Particle1D3D* particle; };
    std::vector<Node> nodes;
    for(auto kind:{ParticleKind::Positive,ParticleKind::Negative}) for(auto& p:candidate.list(kind)) {
        Node node; node.sign=kind==ParticleKind::Positive?1.:-1.; node.particle=&p;
        double r2=0.;
        for(std::size_t j=0;j<3;++j) {
            node.x[j]=(p.velocity(static_cast<int>(j))-center[j])/radius;
            if(!std::isfinite(node.x[j])) throw std::invalid_argument("Nonfinite correction velocity");
            r2+=node.x[j]*node.x[j];
        }
        if(r2>1.+1e-12) throw std::invalid_argument("Correction particle outside core sphere");
        node.original=node.x; nodes.push_back(node);
    }
    std::array<double,6> goal{};
    for(std::size_t j=0;j<3;++j) {
        goal[j]=(target[1+j]-center[j]*target[0])/weight/radius;
        goal[3+j]=(target[4+j]-2.*center[j]*target[1+j]+center[j]*center[j]*target[0])/weight/radius/radius;
    }
    for(double x:goal) if(!std::isfinite(x)) throw std::invalid_argument("Normalized target overflow");
    const double scale=std::max(1.,static_cast<double>(nodes.size()));
    const auto residual=[&](const std::vector<Node>& values) {
        auto r=goal;
        for(const auto& n:values) for(std::size_t j=0;j<3;++j) {
            r[j]-=n.sign*n.x[j]; r[3+j]-=n.sign*n.x[j]*n.x[j];
        }
        return r;
    };
    const auto norm=[](const std::array<double,6>& r) {
        double out=0.; for(double x:r) out=std::max(out,std::abs(x)); return out;
    };
    const auto displacement=[&]() {
        double sum=0.,maximum=0.;
        for(const auto& n:nodes) {
            double d2=0.; for(std::size_t j=0;j<3;++j) d2+=(n.x[j]-n.original[j])*(n.x[j]-n.original[j]);
            sum+=d2; maximum=std::max(maximum,d2);
        }
        result.rmsDisplacement=nodes.empty()?0.:radius*std::sqrt(sum/nodes.size());
        result.maxDisplacement=radius*std::sqrt(maximum);
    };
    for(std::size_t iteration=0;iteration<=config.maxIterations;++iteration) {
        const auto r=residual(nodes); const double error=norm(r);
        result.normalizedResidual=error/scale; result.iterations=iteration;
        if(error<=config.relativeTolerance*scale) {
            displacement();
            if(result.rmsDisplacement>radius*config.maxRmsDisplacementFraction) {
                result.status=MomentCorrectionStatus::DisplacementLimit; return result;
            }
            for(auto& n:nodes) {
                std::array<double,3> velocity{};
                for(std::size_t j=0;j<3;++j) velocity[j]=center[j]+radius*n.x[j];
                n.particle->setVelocity(velocity);
            }
            candidate.computeMoments();
            core=std::move(candidate); result.status=MomentCorrectionStatus::Success; return result;
        }
        if(iteration==config.maxIterations) break;
        if(nodes.empty()) { result.status=MomentCorrectionStatus::InsufficientParticles; return result; }
        std::vector<bool> free(nodes.size(),true);
        std::vector<std::array<double,3>> delta(nodes.size());
        double alpha=1.; bool ready=false;
        // Freeze nodes whose sphere constraint blocks the local minimum-norm
        // step, then recompute it on the remaining free nodes.
        for(std::size_t activePass=0;activePass<=nodes.size();++activePass) {
            const auto count=std::count(free.begin(),free.end(),true);
            if(count==0) break;
            bool singular=false;
            for(std::size_t j=0;j<3;++j) {
                double mean=0.;
                for(std::size_t i=0;i<nodes.size();++i) if(free[i]) mean+=nodes[i].x[j];
                mean/=static_cast<double>(count);
                double variance=0.;
                for(std::size_t i=0;i<nodes.size();++i) if(free[i]) variance+=(nodes[i].x[j]-mean)*(nodes[i].x[j]-mean);
                if(variance<1e-14*count) { singular=true; break; }
                const double b=(r[3+j]-2.*mean*r[j])/(4.*variance);
                for(std::size_t i=0;i<nodes.size();++i) delta[i][j]=free[i]
                    ?nodes[i].sign*(r[j]/static_cast<double>(count)+2.*b*(nodes[i].x[j]-mean)):0.;
            }
            if(singular) { result.status=MomentCorrectionStatus::Singular; displacement(); return result; }
            alpha=1.; bool froze=false;
            for(std::size_t i=0;i<nodes.size();++i) if(free[i]) {
                double a=0.,b=0.,c=-1.;
                for(std::size_t j=0;j<3;++j) { a+=delta[i][j]*delta[i][j]; b+=nodes[i].x[j]*delta[i][j]; c+=nodes[i].x[j]*nodes[i].x[j]; }
                if(!std::isfinite(a) || !std::isfinite(b)) { result.status=MomentCorrectionStatus::Singular; return result; }
                if(a==0.) continue;
                c=std::min(0.,c);
                const double root=std::sqrt(std::max(0.,b*b-a*c));
                const double limit=b>0.?-c/(b+root):(-b+root)/a;
                if(limit<.1) { free[i]=false; delta[i]={}; froze=true; }
                else alpha=std::min(alpha,.99*limit);
            }
            if(!froze) { ready=true; break; }
        }
        if(!ready) { result.status=MomentCorrectionStatus::SupportLimited; displacement(); return result; }
        bool accepted=false;
        for(int backtrack=0;backtrack<24;++backtrack) {
            auto trial=nodes; bool supported=true;
            for(std::size_t i=0;i<trial.size();++i) {
                double r2=0.;
                for(std::size_t j=0;j<3;++j) { trial[i].x[j]+=alpha*delta[i][j]; r2+=trial[i].x[j]*trial[i].x[j]; }
                supported=supported && std::isfinite(r2) && r2<=1.+1e-12;
            }
            if(supported && norm(residual(trial))<error) { nodes=std::move(trial); accepted=true; break; }
            alpha*=.5;
        }
        if(!accepted) { result.status=MomentCorrectionStatus::SupportLimited; displacement(); return result; }
        displacement();
        if(result.rmsDisplacement>radius*config.maxRmsDisplacementFraction) {
            result.status=MomentCorrectionStatus::DisplacementLimit; return result;
        }
    }
    displacement(); result.status=MomentCorrectionStatus::IterationLimit; return result;
}
} // namespace coulomb::resampling
