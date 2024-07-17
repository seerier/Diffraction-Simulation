// integrators/fsdtest.cpp*
#include "camera.h"
#include "film.h"
#include "integrators/fsdtest.h"
#include "interaction.h"
#include "paramset.h"
#include "accelerators/kdtreeaccel.h"
#include "shapes/triangle.h"

namespace pbrt {

    // FsdTestIntegrator Method Definitions
    Spectrum FsdTestIntegrator::Li(const RayDifferential &ray, const Scene &scene,
        Sampler &sampler, MemoryArena &arena,
        int depth) const {
        Spectrum L(0.);

        SurfaceInteraction intr;
        if (!scene.Intersect(ray, &intr)) {
            return Spectrum(1.);
        }

        std::shared_ptr<KdTreeAccel> tree = std::dynamic_pointer_cast<KdTreeAccel>(scene.aggregate);
        if (!tree) return Spectrum(1.);

        bool flag = false;
        std::vector<Triangle> triangles = tree->RadiusSearch(intr.p, .2f, flag);
        

        int n = triangles.size();
        if (n == 0) return Spectrum(0.f);
        return Spectrum(1.f / n);

        //return L;
        
    }

    FsdTestIntegrator *CreateFsdTestIntegrator(const ParamSet &params,
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
                    Bounds2i{ {pb[0], pb[2]}, {pb[1], pb[3]} });
                if (pixelBounds.Area() == 0)
                    Error("Degenerate \"pixelbounds\" specified.");
            }
        }
        return new FsdTestIntegrator(maxDepth, camera, sampler, pixelBounds);
    }

}  // namespace pbrt
