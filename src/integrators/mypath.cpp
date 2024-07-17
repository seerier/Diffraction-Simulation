// integrators/mypath.cpp*
#include "bssrdf.h"
#include "camera.h"
#include "film.h"
#include "integrators/mypath.h"
#include "interaction.h"
#include "paramset.h"
#include "scene.h"
#include "stats.h"

namespace pbrt {

STAT_PERCENT("Integrator/Zero-radiance paths", zeroRadiancePaths, totalPaths);
STAT_INT_DISTRIBUTION("Integrator/Path length", pathLength);

// MyPathIntegrator Method Definitions
MyPathIntegrator::MyPathIntegrator(int maxDepth,
                               std::shared_ptr<const Camera> camera,
                               std::shared_ptr<Sampler> sampler,
                               const Bounds2i &pixelBounds, Float rrThreshold,
                               const std::string &lightSampleStrategy)
    : SamplerIntegrator(camera, sampler, pixelBounds),
      maxDepth(maxDepth),
      rrThreshold(rrThreshold),
      lightSampleStrategy(lightSampleStrategy) {}

void MyPathIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
    lightDistribution =
        CreateLightSampleDistribution(lightSampleStrategy, scene);
}

Spectrum MyPathIntegrator::Li(const RayDifferential &r, const Scene &scene,
    Sampler &sampler, MemoryArena &arena,
                              int depth) const {
    
    Spectrum L(0.), beta(1.);
    RayDifferential ray(r);
    bool sampledspecular = false;

    for (int bounces = 0;; ++bounces) {

        SurfaceInteraction intr;

        bool foundIntersection = scene.Intersect(ray, &intr);

        if (bounces == 0 || sampledspecular) {

            if (foundIntersection) {
                L += beta * intr.Le(-ray.d);
            } 
            else {
                for (const auto &light : scene.infiniteLights) {
                    L += beta * light->Le(ray);
                }
            }
        }

        if (!foundIntersection) {
            return L;
        }

        //const Normal3f n = intr.shading.n;

        intr.ComputeScatteringFunctions(ray, arena, true);
        if (!intr.bsdf) {
            --bounces;
            ray = intr.SpawnRay(ray.d);
            continue;
        }

        if (bounces >= maxDepth) return L;

        const Distribution1D *distribution = lightDistribution->Lookup(intr.p);

        if (intr.bsdf->NumComponents(BxDFType(BSDF_ALL & ~BSDF_SPECULAR)) > 0) {

            Spectrum Ld = UniformSampleOneLight(intr, scene, arena, sampler, false, distribution);
            L += Ld * beta;

        }



        Vector3f wo = -ray.d, wi;
        Float pdf;
        BxDFType sampledType;

        Spectrum f = intr.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, BSDF_ALL, &sampledType);
        if (f.IsBlack() || pdf == 0.f) {
            break;
        }

        ray = intr.SpawnRay(wi);
        beta *= f * AbsDot(intr.shading.n, wi) / pdf;
        sampledspecular = (sampledType & BSDF_SPECULAR) != 0;

        if (bounces >= 3) {
            Float q = std::min(beta.MaxComponentValue(), 0.95f);
            if (sampler.Get1D() > q) {
                return L;
            }
            Float inv_q = 1.0 / q;
            beta *= inv_q;
        }
        


        //Vector3f wo = -ray.d, wi;

        //intr.bsdf->Sample_f()

    }

    return L;

}

MyPathIntegrator *CreateMyPathIntegrator(const ParamSet &params,
                                     std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    int np;
    const int *pb = params.FindInt("pixelbounds", &np);
    Bounds2i pixelBounds = camera->film->GetSampleBounds();
    if (pb) {
        if (np != 4)
            Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                  np);
        else {
            pixelBounds = Intersect(pixelBounds,
                                    Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
            if (pixelBounds.Area() == 0)
                Error("Degenerate \"pixelbounds\" specified.");
        }
    }
    Float rrThreshold = params.FindOneFloat("rrthreshold", 1.);
    std::string lightStrategy =
        params.FindOneString("lightsamplestrategy", "spatial");
    return new MyPathIntegrator(maxDepth, camera, sampler, pixelBounds,
                              rrThreshold, lightStrategy);
}

}  // namespace pbrt
