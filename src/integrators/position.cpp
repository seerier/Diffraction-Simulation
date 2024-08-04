// integrators/position.cpp*
#include "camera.h"
#include "film.h"
#include "integrators/position.h"
#include "interaction.h"
#include "paramset.h"

namespace pbrt {

// PositionIntegrator Method Definitions
Spectrum PositionIntegrator::Li(const RayDifferential &ray, const Scene &scene,
                              Sampler &sampler, MemoryArena &arena,
                              int depth) const {
    Spectrum L(0.);

    SurfaceInteraction intr;
    if (!scene.Intersect(ray, &intr)) {
        return Spectrum(0.);
    }

    Point3f p = intr.p;
    //n = Normalize(n);
    p = (Point3f(1, 1, 1) + p) / 2.f;
    Float rgb[3] = {p.x, p.y, p.z};
    // return Spectrum(CoefficientSpectrum<3>(n.x, n.y, n.z));
    return RGBSpectrum::FromRGB(rgb);

    /*





    Spectrum L(0.);
    // Find closest ray intersection or return background radiance
    SurfaceInteraction isect;
    if (!scene.Intersect(ray, &isect)) {
        for (const auto &light : scene.lights) L += light->Le(ray);
        return L;
    }

    // Compute emitted and reflected light at ray intersection point

    // Initialize common variables for Position integrator
    const Position3f &n = isect.shading.n;
    Vector3f wo = isect.wo;

    // Compute scattering functions for surface interaction
    isect.ComputeScatteringFunctions(ray, arena);
    if (!isect.bsdf)
        return Li(isect.SpawnRay(ray.d), scene, sampler, arena, depth);

    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Add contribution of each light source
    for (const auto &light : scene.lights) {
        Vector3f wi;
        Float pdf;
        VisibilityTester visibility;
        Spectrum Li =
            light->Sample_Li(isect, sampler.Get2D(), &wi, &pdf, &visibility);
        if (Li.IsBlack() || pdf == 0) continue;
        Spectrum f = isect.bsdf->f(wo, wi);
        if (!f.IsBlack() && visibility.Unoccluded(scene))
            L += f * Li * AbsDot(wi, n) / pdf;
    }
    if (depth + 1 < maxDepth) {
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, isect, scene, sampler, arena, depth);
        L += SpecularTransmit(ray, isect, scene, sampler, arena, depth);
    }
    return L;




    */
}

PositionIntegrator *CreatePositionIntegrator(const ParamSet &params,
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
    return new PositionIntegrator(maxDepth, camera, sampler, pixelBounds);
}

}  // namespace pbrt
