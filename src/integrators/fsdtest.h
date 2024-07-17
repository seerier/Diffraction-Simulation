
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_FSDTEST_H
#define PBRT_INTEGRATORS_FSDTEST_H

// integrators/fsdtest.h*
#include "integrator.h"
#include "pbrt.h"
#include "scene.h"

namespace pbrt {

    // FsdTestIntegrator Declarations
    class FsdTestIntegrator : public SamplerIntegrator {
    public:
        // FsdTestIntegrator Public Methods
        FsdTestIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
            std::shared_ptr<Sampler> sampler,
            const Bounds2i &pixelBounds)
            : SamplerIntegrator(camera, sampler, pixelBounds), maxDepth(maxDepth) {}
        Spectrum Li(const RayDifferential &ray, const Scene &scene,
            Sampler &sampler, MemoryArena &arena, int depth) const;

    private:
        // FsdTestIntegrator Private Data
        const int maxDepth;
    };

    FsdTestIntegrator *CreateFsdTestIntegrator(const ParamSet &params,
        std::shared_ptr<Sampler> sampler,
        std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_FSDTEST_H
