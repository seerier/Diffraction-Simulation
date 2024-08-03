// integrators/fsdtest.cpp*
#include "camera.h"
#include "film.h"
#include "integrators/fsdtest.h"
#include "interaction.h"
#include "paramset.h"
#include "accelerators/kdtreeaccel.h"
#include "shapes/triangle.h"

#include <random>

namespace pbrt {

// from zgx
inline float random_float() {
    //Generate a float in [0, 1).
    static std::uniform_real_distribution<float> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

    // FsdTestIntegrator Method Definitions
    Spectrum FsdTestIntegrator::Li(const RayDifferential &ray, const Scene &scene,
        Sampler &sampler, MemoryArena &arena,
        int depth) const {
        Spectrum L(0.);

        SurfaceInteraction intr;
        if (!scene.Intersect(ray, &intr)) {
            return Spectrum(0.);
        }

        std::shared_ptr<KdTreeAccel> tree = std::dynamic_pointer_cast<KdTreeAccel>(scene.aggregate);
        if (!tree) return Spectrum(0.);

        //Float r = 1e-3f;
        Float wavelength = random_float() * 300e-9f + 400e-9f;
        Float beam_sigma = wavelength * 25.f;
        Float search_radius = 3.f * beam_sigma;
        std::vector<Triangle> triangles = tree->RadiusSearch(intr.p, search_radius);

        int n = triangles.size();
        if (n == 0) return Spectrum(0.f);
        return Spectrum(1.0f / n);

        /*
        Vector3f zDir = Vector3f(1.f, 0.f, 0.f);
        Vector3f xDir, yDir;
        CoordinateSystem(zDir, &xDir, &yDir);
        const Point3f &p = intr.p;
        Float area = 0;
        
        for (const Triangle &triangle : triangles) {
            // get vertex positions
            const Point3f &p0 = triangle.mesh->p[triangle.v[0]];
            const Point3f &p1 = triangle.mesh->p[triangle.v[1]];
            const Point3f &p2 = triangle.mesh->p[triangle.v[2]];

            // TODO: continue when triangle is back facing
            Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
            Normal3f n = Normal3f(Normalize(Cross(dp02, dp12)));
            if (triangle.reverseOrientation ^ triangle.transformSwapsHandedness)
                n = -n;
            if (Dot(n, zDir) > 0) continue;

            // project triangle to virtual screen
            Vector3f dp0 = p0 - p;
            Vector3f dp1 = p1 - p;
            Vector3f dp2 = p2 - p;
            Point2f a(Dot(dp0, xDir), Dot(dp0, yDir));
            Point2f b(Dot(dp1, xDir), Dot(dp1, yDir));
            Point2f c(Dot(dp2, xDir), Dot(dp2, yDir));

            // compute the area that triangle covers
            // S_abc = S_aob+S_boc+S_coa
            area += areaCircleTri(r, a, b, c);

        }

        area /= Pi * r * r;
        if (area < 1e-3 || area>1 - 1e-3) return Spectrum(.0f);
        
        return Spectrum(area);
        */
        
        
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
