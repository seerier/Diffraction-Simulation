#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_MYVOLPATH_H
#define PBRT_INTEGRATORS_MYVOLPATH_H

// integrators/myvolpath.h*
#include "pbrt.h"
#include "integrator.h"
#include "lightdistrib.h"

namespace pbrt {

    // MyVolPathIntegrator Declarations
    class MyVolPathIntegrator : public SamplerIntegrator {
    public:
        // MyVolPathIntegrator Public Methods
        MyVolPathIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
            std::shared_ptr<Sampler> sampler,
            const Bounds2i &pixelBounds, Float rrThreshold = 1,
            const std::string &lightSampleStrategy = "spatial")
            : SamplerIntegrator(camera, sampler, pixelBounds),
            maxDepth(maxDepth),
            rrThreshold(rrThreshold),
            lightSampleStrategy(lightSampleStrategy) { }
        void Preprocess(const Scene &scene, Sampler &sampler);
        Spectrum Li(const RayDifferential &ray, const Scene &scene,
            Sampler &sampler, MemoryArena &arena, int depth) const;

    private:
        // MyVolPathIntegrator Private Data
        const int maxDepth;
        const Float rrThreshold;
        const std::string lightSampleStrategy;
        std::unique_ptr<LightDistribution> lightDistribution;
    };

    MyVolPathIntegrator *CreateMyVolPathIntegrator(
        const ParamSet &params, std::shared_ptr<Sampler> sampler,
        std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_MYVOLPATH_H