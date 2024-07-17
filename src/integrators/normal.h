
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_NORMAL_H
#define PBRT_INTEGRATORS_NORMAL_H

// integrators/normal.h*
#include "integrator.h"
#include "pbrt.h"
#include "scene.h"

namespace pbrt {

// NormalIntegrator Declarations
class NormalIntegrator : public SamplerIntegrator {
  public:
    // NormalIntegrator Public Methods
    NormalIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
                     std::shared_ptr<Sampler> sampler,
                     const Bounds2i &pixelBounds)
        : SamplerIntegrator(camera, sampler, pixelBounds), maxDepth(maxDepth) {}
    Spectrum Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const;

  private:
    // NormalIntegrator Private Data
    const int maxDepth;
};

NormalIntegrator *CreateNormalIntegrator(const ParamSet &params,
                                         std::shared_ptr<Sampler> sampler,
                                         std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_Normal_H
