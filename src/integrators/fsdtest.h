
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





    inline float areaTri(const Point2f &a, const Point2f &b) { return .5f * (a.x * b.y - b.x * a.y); }
    inline int intersectCircleLine(const float r, const Point2f &a, const Point2f &b, Point2f &p1, Point2f &p2, float &t1, float &t2) {
        const auto d = b - a;
        const auto adotd = Dot((Vector2f)a, d);
        const float delta2 = 4 * adotd * adotd - 4 * d.LengthSquared() * (a.x * a.x + a.y * a.y - r * r);
        if (delta2 <= 0) { t1 = 0; t2 = 1; return 0; }

        const auto delta = sqrtf(delta2);
        t1 = (-adotd + .5f * delta) / d.LengthSquared();
        t2 = (-adotd - .5f * delta) / d.LengthSquared();
        if (t2 < t1) std::swap(t1, t2);
        p1 = a + t1 * d;
        p2 = a + t2 * d;

        const bool p1v = t1 > 0 && t1 < 1;
        const bool p2v = t2 > 0 && t2 < 1;
        if (!p1v && !p2v) return 0;
        if (p1v && p2v) return 2;
        if (!p1v) p1 = p2;
        return 1;
    }
    inline int intersectCircleLine(const float r, const Point2f &a, const Point2f &b, Point2f &p1, Point2f &p2) {
        Float t1, t2;
        return intersectCircleLine(r, a, b, p1, p2, t1, t2);
    }
    inline int intersectCircleLine(const float r, const Point2f &a, const Point2f &b, Float &t1, Float &t2) {
        Point2f p1, p2;
        return intersectCircleLine(r, a, b, p1, p2, t1, t2);
    }
    inline float areaCircSector(const float r, float theta) {
        theta = fabsf(theta);
        if (theta > Pi) theta = 2 * Pi - theta;
        return .5f * r * r * theta;
    }
    inline float areaCircSectorLine(const float r, const Point2f &a, const Point2f &b) {
        Point2f p1 = { 0,0 }, p2 = { 0,0 };
        const auto ps = intersectCircleLine(r, a, b, p1, p2);

        const auto thetaa = atan2f(a.y, a.x);
        const auto thetab = atan2f(b.y, b.x);
        const auto thetap1 = ps > 0 ? atan2f(p1.y, p1.x) : .0f;
        const auto thetap2 = ps > 1 ? atan2f(p2.y, p2.x) : .0f;
        auto wnd = thetaa - thetab;
        if (wnd < -Pi) wnd += 2 * Pi;
        if (wnd > Pi)  wnd -= 2 * Pi;
        const auto sgn = wnd < .0f ? -1.f : 1.f;

        if (ps == 2) {
            return sgn * (fabsf(areaTri(p1, p2)) + areaCircSector(r, thetaa - thetap1) + areaCircSector(r, thetap2 - thetab));
        } else if (ps == 1) {
            if (a.x * a.x + a.y * a.y < r * r) return sgn * (fabsf(areaTri(a, p1)) + areaCircSector(r, thetab - thetap1));
            else                 return sgn * (fabsf(areaTri(p1, b)) + areaCircSector(r, thetap1-thetaa)); //gxzhao thinks this is incorrect
        }
        const auto thetaab = thetaa - thetab;
        if (a.x * a.x + a.y * a.y < r * r) return sgn * fabsf(areaTri(a, b));
        else                 return sgn * areaCircSector(r, thetaab);
    }
    inline float areaCircleTri(const float r, const Point2f &a, const Point2f &b, const Point2f &c) {
        return fabsf(areaCircSectorLine(r, a, b) +
            areaCircSectorLine(r, b, c) +
            areaCircSectorLine(r, c, a));
    }

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_FSDTEST_H
