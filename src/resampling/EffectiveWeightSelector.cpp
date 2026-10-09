#include "EffectiveWeightSelector.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace coulomb {
EffectiveWeightSelection EffectiveWeightSelector::select(
    double ws, double wf, double wsMin, double wsMax, double wfMin, double wfMax,
    double gamma0, double gamma1, double dt, double collision,
    const std::vector<EffectiveWeightCell>& cells) const {
    const double inputs[] = {ws,wf,wsMin,wsMax,wfMin,wfMax,gamma0,gamma1,dt,collision};
    for (double value : inputs)
        if (!std::isfinite(value)) throw std::invalid_argument("Non-finite effective-weight input");
    if (ws<=0 || wf<=0 || wsMin<=0 || wfMin<=0 || wsMin>wsMax || wfMin>wfMax)
        throw std::invalid_argument("Invalid effective-weight bounds");
    if (cells.empty()) return {ws,wf};
    // x=ws/newSignedWeight, y=wf/newFullWeight. All constraints are linear.
    struct Constraint { double a,b,c; }; // a*x+b*y >= c
    std::vector<Constraint> constraints{{1,0,ws/wsMax},{-1,0,-ws/wsMin},
                                       {0,1,wf/wfMax},{0,-1,-wf/wfMin}};
    double signedCost=0, fullCost=0;
    const double alpha=1+std::max(0.0,gamma0)+std::max(0.0,gamma1)*std::max(0.0,dt)*std::max(0.0,collision);
    for (const auto& cell : cells) {
        if (!std::isfinite(cell.omega) || !std::isfinite(cell.macroMass) ||
            cell.macroMass<0 || cell.signedParticleCount<0 || cell.fullParticleCount<0)
            throw std::invalid_argument("Invalid effective-weight cell");
        const double nd=cell.signedParticleCount, nf=cell.macroMass/wf;
        signedCost+=alpha*nd;
        fullCost+=0.5*nf;
        constraints.push_back({-nd,nf,0});
        if (nd>0 && cell.fullParticleCount>0) {
            const double omega=std::clamp(cell.omega,0.0,1.0);
            constraints.push_back({1-omega,omega,1});
        }
    }
    double bestX=1, bestY=1, bestCost=std::numeric_limits<double>::infinity();
    const auto consider=[&](double x,double y) {
        if (!std::isfinite(x) || !std::isfinite(y) || x<=0 || y<=0) return;
        for (const auto& c : constraints) {
            const double lhs=c.a*x+c.b*y;
            if (lhs+1e-10*std::max({1.0,std::abs(lhs),std::abs(c.c)})<c.c) return;
        }
        const double cost=signedCost*x+fullCost*y;
        if (cost<bestCost) {bestX=x; bestY=y; bestCost=cost;}
    };
    consider(1,1);
    // A bounded two-variable linear program attains its minimum at a vertex.
    for (std::size_t i=0;i<constraints.size();++i)
        for (std::size_t j=i+1;j<constraints.size();++j) {
            const auto& p=constraints[i]; const auto& q=constraints[j];
            const double determinant=p.a*q.b-q.a*p.b;
            if (std::abs(determinant)<=1e-14*std::max({1.0,std::abs(p.a*q.b),std::abs(q.a*p.b)})) continue;
            consider((p.c*q.b-q.c*p.b)/determinant,(p.a*q.c-q.a*p.c)/determinant);
        }
    if (!std::isfinite(bestCost)) throw std::runtime_error("Effective-weight constraints infeasible within configured bounds");
    return {std::clamp(ws/bestX,wsMin,wsMax),std::clamp(wf/bestY,wfMin,wfMax)};
}
} // namespace coulomb
