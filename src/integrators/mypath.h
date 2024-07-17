#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_MYPATH_H
#define PBRT_INTEGRATORS_MYPATH_H

// integrators/mypath.h*
#include "integrator.h"
#include "lightdistrib.h"
#include "pbrt.h"

namespace pbrt {

// MyPathIntegrator Declarations
class MyPathIntegrator : public SamplerIntegrator {
  public:
    // MyPathIntegrator Public Methods
    MyPathIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
                   std::shared_ptr<Sampler> sampler,
                   const Bounds2i &pixelBounds, Float rrThreshold = 1,
                   const std::string &lightSampleStrategy = "spatial");

    void Preprocess(const Scene &scene, Sampler &sampler);
    Spectrum Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const;

  private:
    // MyPathIntegrator Private Data
    const int maxDepth;
    const Float rrThreshold;
    const std::string lightSampleStrategy;
    std::unique_ptr<LightDistribution> lightDistribution;
};

MyPathIntegrator *CreateMyPathIntegrator(const ParamSet &params,
                                     std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_MyPATH_H