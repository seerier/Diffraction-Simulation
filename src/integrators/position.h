
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_POSITION_H
#define PBRT_INTEGRATORS_POSITION_H

// integrators/position.h*
#include "integrator.h"
#include "pbrt.h"
#include "scene.h"

namespace pbrt {

// PositionIntegrator Declarations
class PositionIntegrator : public SamplerIntegrator {
  public:
    // PositionIntegrator Public Methods
    PositionIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
                     std::shared_ptr<Sampler> sampler,
                     const Bounds2i &pixelBounds)
        : SamplerIntegrator(camera, sampler, pixelBounds), maxDepth(maxDepth) {}
    Spectrum Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const;

  private:
    // PositionIntegrator Private Data
    const int maxDepth;
};

PositionIntegrator *CreatePositionIntegrator(const ParamSet &params,
                                         std::shared_ptr<Sampler> sampler,
                                         std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_POSITION_H
