#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_SIMPLE_H
#define PBRT_INTEGRATORS_SIMPLE_H


// integrators/simple.h*
//#include "pbrt.h"
#include "integrator.h"
#include "scene.h"

namespace pbrt {
	//A try to make a simple integrator
class SimpleIntegrator : public SamplerIntegrator {
  public:
    //SimpleIntegrator(std::shared_ptr<const Camera> camera,
    //                         std::shared_ptr<Sampler> sampler,
    //                         const Bounds2i &pixelBounds)
    //    : SamplerIntegrator(camera, sampler, pixelBounds) {}

    using SamplerIntegrator::SamplerIntegrator; //Inherited Constructor.

    Spectrum Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena,
                int depth = 0) const override {
        return Spectrum(0.5f);
    }
};

	SimpleIntegrator *CreateSimpleIntegrator(
            const ParamSet &params, std::shared_ptr<Sampler> sampler,
            std::shared_ptr<const Camera> camera);

        

}//namespace pbrt

#endif  //PBRT_INTEGRATORS_SIMPLE_H
