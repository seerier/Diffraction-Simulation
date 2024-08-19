// materials/fsd.cpp*
#include "materials/fsd.h"
#include "paramset.h"
#include "reflection.h"
#include "interaction.h"
#include "texture.h"
#include "interaction.h"

#include "matte.h"
#include "textures/constant.h"

#include "accelerators/kdtreeaccel.h"
#include "shapes/triangle.h"
#include "scene.h"
#include "integrators/fsdtest.h"
#include "sampling.h"
#include "samplers/random.h"

#include <fstream>
#include <array>
#include <random>


namespace pbrt {

// from zgx
inline float random_float() {
    //Generate a float in [0, 1).
    static std::uniform_real_distribution<float> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

// from fsdBSDF source code
// from zgx
inline Float sqr(Float t) {
    return t * t;
}
// From boost
template<typename T>
inline T sinc(const T x) {
    T const taylor_0_bound = std::numeric_limits<T>::epsilon();
    T const taylor_2_bound = static_cast<T>(0.00034526698300124390839884978618400831996329879769945L);
    T const taylor_n_bound = static_cast<T>(0.018581361171917516667460937040007436176452688944747L);

    if (abs(x) >= taylor_n_bound)
        return sin(x) / x;
    // approximation by taylor series in x at 0 up to order 0
    T    result = static_cast<T>(1);
    if (abs(x) >= taylor_0_bound) {
        T    x2 = x * x;
        // approximation by taylor series in x at 0 up to order 2
        result -= x2 / static_cast<T>(6);
        if (abs(x) >= taylor_2_bound)
            result += (x2 * x2) / static_cast<T>(120);
    }

    return result;
}
// Alpha functions and Psi
static inline Float alpha1(Float x, Float y) { return x == 0 ? 0 : y / (x * x + y * y) * Pi * (cos(x / 2) - sinc(x / 2)) / (2 * x); }
static inline Float alpha2(Float x, Float y) { return x == 0 ? 0 : y / (x * x + y * y) * Pi * sinc(x / 2) / 2; }
static inline Float chi(Float r2) { return sqrtf(1 - exp(-.5f * r2 / 3)); }
static inline c_t Psihat(c_t a, c_t b, Vector2f e, Point2f v, Float k, Vector2f xi) {
    const Float ee = e.Length();
    const Float vxi = Dot((Vector2f)v, xi);
    const auto m = Vector2f{ e.y,-e.x };

    xi = k * Vector2f{ Dot(e,xi), Dot(m,xi) };
    const float chi0 = chi(Dot(xi, xi));
    const c_t a1 = (a - b) * alpha1(xi.x, xi.y);
    const c_t a2 = (a + b) / 2.f * alpha2(xi.x, xi.y);

    return k * ee * ee * chi0 * std::polar(1.f, -k * vxi) * (a1 + c_t{ 0,1 }*a2);
}
static inline Float Psihat2(c_t a, c_t b, Vector2f e, Float k, Vector2f xi) {
    const Float ee = e.Length();
    const auto m = Vector2f{ e.y,-e.x };

    xi = k * Vector2f{ Dot(e,xi), Dot(m,xi) };
    const float chi0 = chi(Dot(xi, xi));
    const c_t a1 = (a - b) * alpha1(xi.x, xi.y);
    const c_t a2 = (a + b) / 2.f * alpha2(xi.x, xi.y);

    return sqr(k * ee * ee * chi0) * std::norm(a1 + c_t{ 0,1 }*a2);
}

// Powers
static inline Float Pt(Point2f u1, Point2f u2, Point2f u3, Float ph1, Float ph2, Float ph3) {
    return (
        fabs(-u1.y * u2.x + u1.x * u2.y + u1.y * u3.x - u2.y * u3.x - u1.x * u3.y + u2.x * u3.y) *
        (sqr(ph3) + sqr(ph2) + sqr(ph1) + ph2 * ph3 + ph1 * ph2 + ph1 * ph3)
        ) / Float(12);
}
static inline Float Pj(c_t a, c_t b, Vector2f e) {
    return e.LengthSquared() * (std::norm(a - b) * 0.004973515f + std::norm(a + b) / 4.f * 14.2569085397f);
}
static inline Float Pjhat(c_t a, c_t b, Vector2f e) {
    return e.LengthSquared() * (std::norm(a - b) * 0.0045255085f + std::norm(a + b) / 4.f * 0.114875434f);
}

static inline Float Psi0t(Point2f u1, Point2f u2, Point2f u3, Float ph1, Float ph2, Float ph3) {
    return ((ph1 + ph2 + ph3) * fabs(-(u1.y * u2.x) + u1.x * u2.y + u1.y * u3.x - u2.y * u3.x - u1.x * u3.y + u2.x * u3.y)) / Float(6);
}
// 0-th order peak covariance
static inline Vector3f Sigmat(Point2f u1, Point2f u2, Point2f u3, Float ph1, Float ph2, Float ph3) {
    const auto u1x = u1.x, u1y = u1.y, u2x = u2.x, u2y = u2.y, u3x = u3.x, u3y = u3.y;
    const auto a = (
        ((3 * ph1 + ph2 + ph3) * sqr(u1x) + (ph1 + 3 * ph2 + ph3) * sqr(u2x) + (ph1 + 2 * (ph2 + ph3)) * u2x * u3x + (ph1 + ph2 + 3 * ph3) * sqr(u3x) + u1x * ((2 * (ph1 + ph2) + ph3) * u2x + (2 * ph1 + ph2 + 2 * ph3) * u3x)) * fabs(-(u1y * u2x) + u1x * u2y + u1y * u3x - u2y * u3x - u1x * u3y + u2x * u3y)
        ) / Float(60);
    const auto b = (
        (u1x * (2 * (3 * ph1 + ph2 + ph3) * u1y + (2 * (ph1 + ph2) + ph3) * u2y + (2 * ph1 + ph2 + 2 * ph3) * u3y) + u3x * ((2 * ph1 + ph2 + 2 * ph3) * u1y + (ph1 + 2 * (ph2 + ph3)) * u2y + 2 * (ph1 + ph2 + 3 * ph3) * u3y) + u2x * ((2 * (ph1 + ph2) + ph3) * u1y + 2 * (ph1 + 3 * ph2 + ph3) * u2y + (ph1 + 2 * (ph2 + ph3)) * u3y)) * fabs(-(u1y * u2x) + u1x * u2y + u1y * u3x - u2y * u3x - u1x * u3y + u2x * u3y)
        ) / Float(120);
    const auto c = (
        ((3 * ph1 + ph2 + ph3) * sqr(u1y) + (ph1 + 3 * ph2 + ph3) * sqr(u2y) + (ph1 + 2 * (ph2 + ph3)) * u2y * u3y + (ph1 + ph2 + 3 * ph3) * sqr(u3y) + u1y * ((2 * (ph1 + ph2) + ph3) * u2y + (2 * ph1 + ph2 + 2 * ph3) * u3y)) * fabs(-(u1y * u2x) + u1x * u2y + u1y * u3x - u2y * u3x - u1x * u3y + u2x * u3y)
        ) / Float(60);

    return { Float(a),Float(b),Float(c) };
}



inline Float FsdBxDF::eval(const Vector2f &xi) const {
    if (Ppl_A <= 0)
        return 0;

    const auto cosine = Float(1) / sqrtf(1 + Dot(xi, xi));
    c_t bsdf = {};
    for (int i = 0; i < edgeNum; ++i) {
        const Edge &e = *edges[i];
        bsdf += Psihat(e.a, e.b, e.e, e.v, 2 * Pi / wavelength, xi);
    }
    return 1.f / (cosine * Ppl_A) * std::norm(bsdf);
}

inline Float FsdBxDF::evalPdf(const Vector2f &xi) const {
    if (sumPhat_j == 0)
        return 0;

    Float pdf = 0;
    for (int i = 0; i < edgeNum; ++i) {
        const Edge &e = *edges[i];
        pdf += Psihat2(e.a, e.b, e.e, 2 * Pi / wavelength, xi);
    }
    return pdf / sumPhat_j;
}





/*
inline Float FsdBxDF::eval(const Vector2f &xi) const {
    if (Ppl_A <= 0)
        return 0;

    const auto cosine = Float(1) / sqrtf(1 + Dot(xi, xi));
    c_t bsdf = {};
    for (const auto &eP : edges) {
        const Edge &e = *eP;
        bsdf += Psihat(e.a, e.b, e.e, e.v, 2 * Pi / wavelength, xi);
    }
    return 1.f / (cosine * Ppl_A) * std::norm(bsdf);
}

inline Float FsdBxDF::evalPdf(const Vector2f &xi) const {
    if (sumPhat_j == 0)
        return 0;

    Float pdf = 0;
    for (const auto &eP : edges) {
        const Edge &e = *eP;
        pdf += Psihat2(e.a, e.b, e.e, 2 * Pi / wavelength, xi);
    }
    return pdf / sumPhat_j;
}
*/


struct fsdPrecomputedTables {
    inline Vector2f importanceSampleCDF1(const Point3f &rand3) const {
        return importanceSampleCDF(rand3, iCDFtheta1, iCDF1);
    }
    inline Vector2f importanceSampleCDF2(const Point3f &rand3) const {
        return importanceSampleCDF(rand3, iCDFtheta2, iCDF2);
    }

private:
    static constexpr std::size_t Nsamples = 1024, Msamples = 1024;
    //std::array<Float, Nsamples> iCDFtheta1, iCDFtheta2;
    //std::array<std::array<Float, Msamples>, Msamples> iCDF1, iCDF2;
    std::unique_ptr<std::array<Float, Nsamples>> iCDFtheta1, iCDFtheta2;
    std::unique_ptr<std::array<std::array<Float, Msamples>, Msamples>> iCDF1, iCDF2;

    template <std::size_t S>
    static inline Float lerp(Float x, const std::array<Float, S> &iCDFtheta) {
        x *= S;
        const std::size_t l = std::min((std::size_t)x, S - 1);
        const auto h = std::min(l + 1, S - 1);
        const auto f = std::max<Float>(1, x - (Float)l);
        return f * iCDFtheta[h] + (1 - f) * iCDFtheta[l];
    }
    template <std::size_t S>
    static inline Float lerp(Float theta, Float rx, const std::array<std::array<Float, S>, S> &iCDF) {
        const auto x = theta * 2 / Pi * S;
        const std::size_t l = std::min((std::size_t)x, S - 1);
        const auto h = std::min(l + 1, S - 1);
        const auto f = std::max<Float>(1, x - (Float)l);
        return f * lerp(rx, iCDF[h]) + (1 - f) * lerp(rx, iCDF[l]);
    }

    static inline Vector2f importanceSampleCDF(const Point3f &rand3,
        const std::unique_ptr<std::array<Float, Nsamples>> &iCDFtheta,
        const std::unique_ptr<std::array<std::array<Float, Msamples>, Msamples>> &iCDF) {
        const auto theta = lerp(rand3.x, *iCDFtheta);
        const auto r = std::max((Float)0, lerp(theta, rand3.y, *iCDF));
        const auto q = std::min(3, (int)(rand3.z * 4));

        auto xi = r * Vector2f{ cosf(theta), sinf(theta) };
        xi.x *= Float(((q + 1) / 2) % 2 == 0 ? 1 : -1);
        xi.y *= Float((q / 2) % 2 == 0 ? 1 : -1);
        return xi;
    }


    /*
    static inline std::array<Float, Nsamples> loadiCDFtheta(const std::string &path) {
        std::array<Float, Nsamples> data = {};

        std::ifstream f(path, std::ios::in);
        if (!f) {
            std::cerr << "Could not open \"" << path << "\"" << std::endl;
            return data;
        }

        std::string str;
        for (auto &l : data) {
            if (!std::getline(f, str)) {
                assert(false);
                break;
            }
            std::stringstream conv(str);
            conv >> l;
        }

        return data;
    }
    static inline std::array<std::array<Float, Msamples>, Msamples> loadiCDF(const std::string &path) {
        //std::array<std::array<Float, Msamples>, Msamples> data = {};
        std::unique_ptr<std::array<std::array<Float, Msamples>, Msamples>> data = std::make_unique<std::array<std::array<Float, Msamples>, Msamples>>();

        std::ifstream f(path, std::ios::in);
        if (!f) {
            std::cerr << "Could not open \"" << path << "\"" << std::endl;
            return data;
        }

        std::string line, str;
        for (auto &row : data) {
            if (!std::getline(f, line)) {
                assert(false);
                break;
            }
            std::stringstream ss(line);
            for (auto &cell : row) {
                if (!std::getline(ss, str, ',')) {
                    assert(false);
                    break;
                }
                std::stringstream conv(str);
                conv >> cell;
            }
        }

        return data;
    }
    */


    static inline std::unique_ptr<std::array<Float, Nsamples>> loadiCDFtheta(const std::string &path) {
        std::unique_ptr<std::array<Float, Nsamples>> data = std::make_unique<std::array<Float, Nsamples>>();

        std::ifstream f(path, std::ios::in);
        if (!f) {
            std::cerr << "Could not open \"" << path << "\"" << std::endl;
            return data;
        }

        std::string str;
        for (auto &l : *data) {
            if (!std::getline(f, str)) {
                assert(false);
                break;
            }
            std::stringstream conv(str);
            conv >> l;
        }

        return data;
    }
    static inline std::unique_ptr<std::array<std::array<Float, Msamples>, Msamples>> loadiCDF(const std::string &path) {
        //std::array<std::array<Float, Msamples>, Msamples> data = {};
        std::unique_ptr<std::array<std::array<Float, Msamples>, Msamples>> data = std::make_unique<std::array<std::array<Float, Msamples>, Msamples>>();

        std::ifstream f(path, std::ios::in);
        if (!f) {
            std::cerr << "Could not open \"" << path << "\"" << std::endl;
            return data;
        }

        std::string line, str;
        for (auto &row : *data) {
            if (!std::getline(f, line)) {
                assert(false);
                break;
            }
            std::stringstream ss(line);
            for (auto &cell : row) {
                if (!std::getline(ss, str, ',')) {
                    assert(false);
                    break;
                }
                std::stringstream conv(str);
                conv >> cell;
            }
        }

        return data;
    }


    void loadiCDFs() {
        std::string basePath("E:/Coding/github_repo/fsd/fsd-pbrt/precompiled_fsd_tables");
        iCDFtheta1 = loadiCDFtheta(basePath + "/iCDFa1theta.csv");
        iCDFtheta2 = loadiCDFtheta(basePath + "/iCDFa2theta.csv");
        iCDF1 = loadiCDF(basePath + "/iCDFa1.csv");
        iCDF2 = loadiCDF(basePath + "/iCDFa2.csv");
    }

public:
    fsdPrecomputedTables() { loadiCDFs(); }
};


/*
// Utility Functions
Float EdgeCoverArea(const Point2f &a, const Point2f &b, Float r) {
    bool aIn = a.x * a.x + a.y * a.y < r *r *(1 + 1e-5);
    bool bIn = b.x * b.x + b.y * b.y < r *r *(1 + 1e-5);

    Float thetaA = atan2f(a.y, a.x);
    Float thetaB = atan2f(b.y, b.x);

    
    if (aIn && bIn) {
        return abs(0.5f * (a.x * b.y - a.y * b.x));
    } else if (aIn && !bIn) {

    }
}


inline Float TriangleCoverArea(const Point2f &a, const Point2f &b, const Point2f &c, Float r) {
    return EdgeCoverArea(a, b, r) + EdgeCoverArea(b, c, r) + EdgeCoverArea(c, a, r);
}
*/






/*
// struct Edge definition
struct Edge {
    Point3f v;
    Vector3f e;
    Vector3f m;
    Float a, b;
};
*/


// FsdMaterial Method Definitions
void FsdMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
    MemoryArena &arena,
    TransportMode mode,
    bool allowMultipleLobes, const Scene &scene) const {
    /*
    // Perform bump mapping with _bumpMap_, if present
    if (bumpMap) Bump(bumpMap, si);

    // Evaluate textures for _FsdMaterial_ material and allocate BRDF
    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
    Spectrum r = Kd->Evaluate(*si).Clamp();
    Float sig = Clamp(sigma->Evaluate(*si), 0, 90);
    if (!r.IsBlack()) {
        if (sig == 0)
            si->bsdf->Add(ARENA_ALLOC(arena, LambertianReflection)(r));
        else
            si->bsdf->Add(ARENA_ALLOC(arena, OrenNayar)(r, sig));
    }
    */

    //Spectrum s = scale->Evaluate(*si).Clamp();

    

    m_Material->ComputeScatteringFunctions(si, arena, mode, allowMultipleLobes, scene);

    BSDF *fsdBSDF = ARENA_ALLOC(arena, FsdBSDF)(*si, scene, si->bsdf, arena);

    //BSDF *fsdBSDF = new FsdBSDF(*si, scene, si->bsdf, arena);
    si->bsdf = fsdBSDF;//->scaledBxDF->bxdf
    //delete fsdBSDF;


}

FsdMaterial *CreateFsdMaterial(const TextureParams &mp, const std::shared_ptr<Material> &material) {

    /*
    std::shared_ptr<Texture<Spectrum>> Kd =
        mp.GetSpectrumTexture("Kd", Spectrum(0.5f));
    std::shared_ptr<Texture<Float>> sigma = mp.GetFloatTexture("sigma", 0.f);
    std::shared_ptr<Texture<Float>> bumpMap =
        mp.GetFloatTextureOrNull("bumpmap");
    return new FsdMaterial(Kd, sigma, bumpMap);
    */


    //return new FsdMaterial(std::make_shared<MatteMaterial>(std::make_shared<ConstantTexture<Spectrum>>(0.5f),
    //    std::make_shared<ConstantTexture<Float>>(0.f), std::make_shared<ConstantTexture<Float>>(0.f)));

    std::shared_ptr<Texture<Spectrum>> scale = mp.GetSpectrumTexture("scale", Spectrum(1.0f));
    

    return new FsdMaterial(material, scale);
}


// FsdBxDF Methods
FsdBxDF::FsdBxDF(const SurfaceInteraction &intr, const Scene &scene, MemoryArena &arena)
    : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_DIFFUSE)) {
    //LOG(INFO) << "edges.size() = " << edges.size();
    //edges.reserve(512);
    //edges.shrink_to_fit();
    build(intr, scene, arena);
}


/*
void FsdBxDF::build(const SurfaceInteraction &intr, const Scene &scene, MemoryArena &arena) {
    //wavelength = random_float() * 300e-9f + 400e-9f;
    //wavelength = 5e-5f;
    //wavelength = 1e-4f;
    wavelength = 5.5e-5f;


    // sample wavelength and handle related computations
    //spectrumIndex = static_cast<int>(59.9999 * random_float());
    //wavelength = 400e-7f + 5e-7f * spectrumIndex;
    //wavelength = 400e-7f + 5e-7f * spectrumIndex;

    Float beam_sigma = wavelength * 25.f;
    Float search_radius = 3.f * beam_sigma;
    Float k = 2 * Pi / wavelength;

    // util
    //RandomSampler r(1);

    Point3f p = intr.p;
    zDir = Normalize(-intr.wo);
    CoordinateSystem(zDir, &xDir, &yDir);

    //
    //isect = p;
    //if (abs(zDir.x) < 0.5f) return;

    std::shared_ptr<KdTreeAccel> tree = std::dynamic_pointer_cast<KdTreeAccel>(scene.aggregate);
    bool flag = false;
    std::vector<Triangle> triangles;
    if (tree) {
        triangles = std::move(tree->RadiusSearch(intr.p, search_radius));
    }

    Float area = 0.f;

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
        ////LOG(INFO) << "Back Facing Test: Dot(n, zDir) = " << Dot(n, zDir);
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
        area += areaCircleTri(search_radius, a, b, c);

    }

    area /= Pi * search_radius * search_radius;





    // for some test
    if (area < 1e-6 || area > 2 - 1e-3) {
        return;
    }

    ////LOG(INFO) << "Computed area: " << area;
    //LOG(INFO) << "Computed area: " << area << ", wo = " << intr.wo << ", zDir = " << zDir << ", p = " << p << ", distance = " << Distance(p, Point3f(-0.4f, 1.f, 0.f));




    // from fsdBSDF source code
    Vector3f Sigmat_0 = { 0,0,0 };
    Float psi0 = 0;
    Float eavg = 0;
    const auto winding = [](const Vector2f &e, const Point2f &v, const Point2f &c) -> Float {
        const auto m = Vector2f{ e.y,-e.x };
        return Dot(m, v - c) > 0 ? 1 : -1;
    };
    const auto addedge = [this, &eavg, &winding, k, &arena](const Point2f &u1, const Point2f &u2, Float z1, Float z2, Float a, Float b, const Point2f &tric) {
        const auto v = (u1 + u2) / Float(2);
        auto e = u2 - u1;
        if (winding(e, v, tric) == -1) {
            e = -e;
            std::swap(a, b);
            std::swap(z1, z2);
        }

        const auto ca = std::polar(a, k * z1);
        const auto cb = std::polar(b, k * z2);

        const auto P = Pjhat(ca, cb, e);
        //LOG(INFO) << "I ENTERED EDGE ADDING!!!";
        if (P == 0) return;
        this->sumPhat_j += P;
        this->edges.push_back(ARENA_ALLOC(arena, Edge)(e, v, ca, cb, P, this->sumPhat_j));
        //LOG(INFO) << "ADDED EDGE!!!";
        eavg += P * e.Length();
    };

    const auto phi = [sigma = beam_sigma](const Point2f p, const Float z) -> Float
    { return expf(-.25f * (p.x * p.x + p.y * p.y + sqr(z)) / sqr(sigma)) / (sqrtf(2 * Pi) * sigma); };

    const auto maxlength2 = sqr(beam_sigma);
    std::function<void(bool, bool, bool, const Point2f &, const Point2f &, const Point2f &, Float, Float, Float, int)> addtri;
    addtri = [&, beam_sigma, search_radius, maxlength2, max_depth = 5](bool edge12, bool edge13, bool edge23, const Point2f &u1, const Point2f &u2, const Point2f &u3, Float z1, Float z2, Float z3, int recr_depth) -> void {
        if (areaCircleTri(search_radius, u1, u2, u3) < 1e-10f)
            return;

        const auto c = (u1 + u2 + u3) / 3.f;
        const auto z0 = (z1 + z2 + z3) / 3.f;

        const bool subdivide12 = (u2 - u1).LengthSquared() > maxlength2;
        const bool subdivide13 = (u3 - u1).LengthSquared() > maxlength2;
        const bool subdivide23 = (u3 - u2).LengthSquared() > maxlength2;
        const auto sss = ((int)subdivide12) + ((int)subdivide13) + ((int)subdivide23) == 1;

        // Subdivide triangles if needed. If only one edge is too long, subdivide along that edge only, otherwise split in the centre as well.
        if (recr_depth < max_depth && (subdivide12 || subdivide13 || subdivide23)) {
            if (subdivide12) {
                addtri(edge12, sss && edge13, false, u1, (u1 + u2) / 2, sss ? u3 : c, z1, (z1 + z2) / 2, sss ? z3 : z0, recr_depth + 1);
                addtri(edge12, false, sss && edge23, (u1 + u2) / 2, u2, sss ? u3 : c, (z1 + z2) / 2, z2, sss ? z3 : z0, recr_depth + 1);
            } else if (!sss)
                addtri(edge12, false, false, u1, u2, c, z1, z2, z0, recr_depth + 1);
            if (subdivide13) {
                addtri(sss && edge12, edge13, false, u1, sss ? u2 : c, (u1 + u3) / 2, z1, sss ? z2 : z0, (z1 + z3) / 2, recr_depth + 1);
                addtri(false, edge13, sss && edge23, (u1 + u3) / 2, sss ? u2 : c, u3, (z1 + z3) / 2, sss ? z2 : z0, z3, recr_depth + 1);
            } else if (!sss)
                addtri(false, edge13, false, u1, c, u3, z1, z0, z3, recr_depth + 1);
            if (subdivide23) {
                addtri(false, sss && edge13, edge23, sss ? u1 : c, (u2 + u3) / 2, u3, sss ? z1 : z0, (z2 + z3) / 2, z3, recr_depth + 1);
                addtri(sss && edge12, false, edge23, sss ? u1 : c, u2, (u2 + u3) / 2, sss ? z1 : z0, z2, (z2 + z3) / 2, recr_depth + 1);
            } else if (!sss)
                addtri(false, false, edge23, c, u2, u3, z0, z2, z3, recr_depth + 1);
            return;
        }

        const auto ph1 = phi(u1, z1);
        const auto ph2 = phi(u2, z2);
        const auto ph3 = phi(u3, z3);

        if (edge12) addedge(u2, u1, z2, z1, ph2, ph1, c);
        if (edge13) addedge(u1, u3, z1, z3, ph1, ph3, c);
        if (edge23) addedge(u3, u2, z3, z2, ph3, ph2, c);


        // Bookkeeping
        Ppl_A += Pt(u1, u2, u3, ph1, ph2, ph3);
        Sigmat_0 += Sigmat(u1, u2, u3, ph1, ph2, ph3);
        psi0 += Psi0t(u1, u2, u3, ph1, ph2, ph3);
    };







    // this method might be incorrect
    // lambda funstion which determines if the edge is boundary
    auto edgeBoundary = [](const std::vector<Triangle> &triangles, int v1, int v2, const Vector3f zDir, const int *v, const std::shared_ptr<TriangleMesh> &mesh) ->bool {
        // handle a special case for cornell box scene
        Point3f &mv1 = mesh->p[v1];
        Point3f &mv2 = mesh->p[v2];

        if (fabs(mv1.x) > 0.9f || fabs(mv2.x) > 0.9f) {
            //LOG(INFO) << "hitten the box face. invalid!" << " : fabs(mesh->p[v1].x) = " << fabs(mv1.x) << ", fabs(mesh->p[v2].x) = " << fabs(mv2.x);
            return false;
        }


        // main logic
        for (const Triangle &tri : triangles) {
            int t1 = tri.v[0];
            int t2 = tri.v[1];
            int t3 = tri.v[2];
            Point3f &mt1 = tri.mesh->p[t1];
            Point3f &mt2 = tri.mesh->p[t2];
            Point3f &mt3 = tri.mesh->p[t3];
            if ((mt1 == mv1 && mt2 == mv2) ||
                (mt1 == mv1 && mt3 == mv2) ||
                (mt2 == mv1 && mt1 == mv2) ||
                (mt2 == mv1 && mt3 == mv2) ||
                (mt3 == mv1 && mt1 == mv2) ||
                (mt3 == mv1 && mt2 == mv2)) {
                if (tri.v == v) {
                    //LOG(INFO) << "ZGX SKIPPED!!!";
                    continue;
                }

                //LOG(INFO) << "Not SKIPPED!";

                const Point3f &p0 = tri.mesh->p[t1];
                const Point3f &p1 = tri.mesh->p[t2];
                const Point3f &p2 = tri.mesh->p[t3];

                if (fabs(p0.x) > 0.9f || fabs(p1.x) > 0.9f || fabs(p2.x) > 0.9f) {
                    //LOG(INFO) << "In loop: hitten the box face. invalid!" << "  : p0 = " << p0 << ", p1 = " << p1 << ", p2 = " << p2;
                    continue;
                }

                Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
                Normal3f n = Normal3f(Normalize(Cross(dp02, dp12)));
                if (tri.reverseOrientation ^ tri.transformSwapsHandedness)
                    n = -n;

                //LOG(INFO) << "In edgeBoundary: n = " << n << ", zDir = " << zDir << ", p0 = " << p0 << ", p1 = " << p1 << ", p2 = " << p2 << ", Dot(zDir, n) = " << Dot(zDir, n);

                if (Dot(zDir, n) > 0) {
                    //LOG(INFO) << "Is edge boundary!!!";
                    return true;
                } else return false;
            }
        }
        return false;
    };

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

        // along z-axis
        Float z0 = Dot(dp0, zDir);
        Float z1 = Dot(dp1, zDir);
        Float z2 = Dot(dp2, zDir);

        // determine which edge is boundary
        bool b0 = edgeBoundary(triangles, triangle.v[1], triangle.v[2], zDir, triangle.v, triangle.mesh);
        bool b1 = edgeBoundary(triangles, triangle.v[2], triangle.v[0], zDir, triangle.v, triangle.mesh);
        bool b2 = edgeBoundary(triangles, triangle.v[0], triangle.v[1], zDir, triangle.v, triangle.mesh);

        //LOG(INFO) << "b0 = " << b0 << ", b1 = " << b1 << ", b2 = " << b2 << ", n = " << n << ", zDir = " << zDir << ", p0 = " << p0 << ", p1 = " << p1 << ", p2 = " << p2
        //    << ", dp0 = " << dp0 << ", dp1 = " << dp1 << ", dp2 = " << dp2 << ", a = " << a << ", b = " << b << ", c = " << c;

        // add triangle
        addtri(b2, b1, b0, a, b, c, z0, z1, z2, 0);
    }

    //LOG(INFO) << "Edges size: " << edges.size() << ", Ppl_A: " << Ppl_A;
    if (edges.size() == 0 || Ppl_A < 1e-2f) {
        Ppl_A = 0;
        return;
    }

    enabled = true;

    // from fsdBSDF source code
    // Power in 0-th order lobe
    eavg /= sumPhat_j;
    const Float sigma_xi = sqrt(3) / (k * eavg);
    Sigmat_0 *= 6 * k * k / psi0;
    const auto Sigma0 = Sigmat_0 + Vector3f{ 1 / sqr(sigma_xi),0,1 / sqr(sigma_xi) };
    const auto detSigma0 = std::max<Float>(0, Sigma0.x * Sigma0.z - sqr(Sigma0.y));
    Ppl_0 = k * k / (18 * Pi) / sqrtf(detSigma0) * sqr(psi0);
}
*/




void FsdBxDF::build(const SurfaceInteraction &intr, const Scene &scene, MemoryArena &arena) {
    //wavelength = random_float() * 300e-9f + 400e-9f;
    //wavelength = 5e-5f;
    //wavelength = 1e-4f;
    //wavelength = 5.5e-5f;

    //LOG(INFO) << "FsdBxDF::build()";

#ifdef PBRT_SAMPLED_SPECTRUM
    // sample wavelength and handle related computations
    spectrumIndex = static_cast<int>(59.9999 * random_float());
    wavelength = 400e-7f + 5e-7f * spectrumIndex;
    //wavelength = 400e-7f + 5e-7f * spectrumIndex;
#else
    wavelength = 5.5e-5f;
#endif

    Float beam_sigma = wavelength * 25.f;
    Float search_radius = 3.f * beam_sigma;
    Float k = 2 * Pi / wavelength;

    //edges.reserve(512);

    // util
    //RandomSampler r(1);

    Point3f p = intr.p;
    zDir = Normalize(-intr.wo);
    CoordinateSystem(zDir, &xDir, &yDir);

    //
    //isect = p;
    //if (abs(zDir.x) < 0.5f) return;

    std::shared_ptr<KdTreeAccel> tree = std::dynamic_pointer_cast<KdTreeAccel>(scene.aggregate);
    bool flag = false;
    std::vector<Triangle> triangles;
    if (tree) {
        triangles = std::move(tree->RadiusSearch(intr.p, search_radius));
    }

    Float area = 0.f;

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
        ////LOG(INFO) << "Back Facing Test: Dot(n, zDir) = " << Dot(n, zDir);
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
        area += areaCircleTri(search_radius, a, b, c);

    }

    area /= Pi * search_radius * search_radius;




    // for some test
    if (area < 1e-6 || area > 2 - 1e-3) {
        return;
    }

    ////LOG(INFO) << "Computed area: " << area;
    //LOG(INFO) << "Computed area: " << area << ", wo = " << intr.wo << ", zDir = " << zDir << ", p = " << p << ", distance = " << Distance(p, Point3f(-0.4f, 1.f, 0.f));




    // from fsdBSDF source code
    Vector3f Sigmat_0 = { 0,0,0 };
    Float psi0 = 0;
    Float eavg = 0;
    const auto winding = [](const Vector2f &e, const Point2f &v, const Point2f &c) -> Float {
        const auto m = Vector2f{ e.y,-e.x };
        return Dot(m, v - c) > 0 ? 1 : -1;
    };
    const auto addedge = [this, &eavg, &winding, k, &arena](const Point2f &u1, const Point2f &u2, Float z1, Float z2, Float a, Float b, const Point2f &tric) {
        const auto v = (u1 + u2) / Float(2);
        auto e = u2 - u1;
        if (winding(e, v, tric) == -1) {
            e = -e;
            std::swap(a, b);
            std::swap(z1, z2);
        }

        const auto ca = std::polar(a, k * z1);
        const auto cb = std::polar(b, k * z2);

        const auto P = Pjhat(ca, cb, e);
        //LOG(INFO) << "I ENTERED EDGE ADDING!!!";
        if (P == 0) return;
        this->sumPhat_j += P;
        //this->edges.emplace_back(Edge{ e,v,ca,cb,P,this->sumPhat_j });
        //this->edges.emplace_back(e, v, ca, cb, P, this->sumPhat_j);
        //this->edges.push_back(std::make_unique<Edge>(e, v, ca, cb, P, this->sumPhat_j));
        //this->edges.push_back(std::make_unique<Edge>(e, v, ca, cb, P, this->sumPhat_j));
        //this->edges.push_back(ARENA_ALLOC(arena, Edge)(e, v, ca, cb, P, this->sumPhat_j));
        if (edgeNum < MaxEdges) {
            edges[edgeNum] = ARENA_ALLOC(arena, Edge)(e, v, ca, cb, P, this->sumPhat_j);
            ++edgeNum;
        }
        //alloca()
        //LOG(INFO) << "ADDED EDGE!!!";
        eavg += P * e.Length();
    };

    const auto phi = [sigma = beam_sigma](const Point2f p, const Float z) -> Float
    { return expf(-.25f * (p.x * p.x + p.y * p.y + sqr(z)) / sqr(sigma)) / (sqrtf(2 * Pi) * sigma); };

    const auto maxlength2 = sqr(beam_sigma);
    std::function<void(bool, bool, bool, const Point2f &, const Point2f &, const Point2f &, Float, Float, Float, int)> addtri;
    addtri = [&, beam_sigma, search_radius, maxlength2, max_depth = 5](bool edge12, bool edge13, bool edge23, const Point2f &u1, const Point2f &u2, const Point2f &u3, Float z1, Float z2, Float z3, int recr_depth) -> void {
        if (areaCircleTri(search_radius, u1, u2, u3) < 1e-10f)
            return;

        const auto c = (u1 + u2 + u3) / 3.f;
        const auto z0 = (z1 + z2 + z3) / 3.f;

        const bool subdivide12 = (u2 - u1).LengthSquared() > maxlength2;
        const bool subdivide13 = (u3 - u1).LengthSquared() > maxlength2;
        const bool subdivide23 = (u3 - u2).LengthSquared() > maxlength2;
        const auto sss = ((int)subdivide12) + ((int)subdivide13) + ((int)subdivide23) == 1;

        // Subdivide triangles if needed. If only one edge is too long, subdivide along that edge only, otherwise split in the centre as well.
        if (recr_depth < max_depth && (subdivide12 || subdivide13 || subdivide23)) {
            if (subdivide12) {
                addtri(edge12, sss && edge13, false, u1, (u1 + u2) / 2, sss ? u3 : c, z1, (z1 + z2) / 2, sss ? z3 : z0, recr_depth + 1);
                addtri(edge12, false, sss && edge23, (u1 + u2) / 2, u2, sss ? u3 : c, (z1 + z2) / 2, z2, sss ? z3 : z0, recr_depth + 1);
            } else if (!sss)
                addtri(edge12, false, false, u1, u2, c, z1, z2, z0, recr_depth + 1);
            if (subdivide13) {
                addtri(sss && edge12, edge13, false, u1, sss ? u2 : c, (u1 + u3) / 2, z1, sss ? z2 : z0, (z1 + z3) / 2, recr_depth + 1);
                addtri(false, edge13, sss && edge23, (u1 + u3) / 2, sss ? u2 : c, u3, (z1 + z3) / 2, sss ? z2 : z0, z3, recr_depth + 1);
            } else if (!sss)
                addtri(false, edge13, false, u1, c, u3, z1, z0, z3, recr_depth + 1);
            if (subdivide23) {
                addtri(false, sss && edge13, edge23, sss ? u1 : c, (u2 + u3) / 2, u3, sss ? z1 : z0, (z2 + z3) / 2, z3, recr_depth + 1);
                addtri(sss && edge12, false, edge23, sss ? u1 : c, u2, (u2 + u3) / 2, sss ? z1 : z0, z2, (z2 + z3) / 2, recr_depth + 1);
            } else if (!sss)
                addtri(false, false, edge23, c, u2, u3, z0, z2, z3, recr_depth + 1);
            return;
        }

        const auto ph1 = phi(u1, z1);
        const auto ph2 = phi(u2, z2);
        const auto ph3 = phi(u3, z3);

        if (edge12) addedge(u2, u1, z2, z1, ph2, ph1, c);
        if (edge13) addedge(u1, u3, z1, z3, ph1, ph3, c);
        if (edge23) addedge(u3, u2, z3, z2, ph3, ph2, c);


        // Bookkeeping
        Ppl_A += Pt(u1, u2, u3, ph1, ph2, ph3);
        Sigmat_0 += Sigmat(u1, u2, u3, ph1, ph2, ph3);
        psi0 += Psi0t(u1, u2, u3, ph1, ph2, ph3);
    };







    // this method might be incorrect
    // lambda funstion which determines if the edge is boundary
    auto edgeBoundary = [](const std::vector<Triangle> &triangles, int v1, int v2, const Vector3f zDir, const int *v, const std::shared_ptr<TriangleMesh> &mesh) ->bool {
        // handle a special case for cornell box scene
        Point3f &mv1 = mesh->p[v1];
        Point3f &mv2 = mesh->p[v2];

        if (fabs(mv1.x) > 0.9f || fabs(mv2.x) > 0.9f) {
            //LOG(INFO) << "hitten the box face. invalid!" << " : fabs(mesh->p[v1].x) = " << fabs(mv1.x) << ", fabs(mesh->p[v2].x) = " << fabs(mv2.x);
            return false;
        }


        // main logic
        for (const Triangle &tri : triangles) {
            int t1 = tri.v[0];
            int t2 = tri.v[1];
            int t3 = tri.v[2];
            Point3f &mt1 = tri.mesh->p[t1];
            Point3f &mt2 = tri.mesh->p[t2];
            Point3f &mt3 = tri.mesh->p[t3];
            if ((mt1 == mv1 && mt2 == mv2) ||
                (mt1 == mv1 && mt3 == mv2) ||
                (mt2 == mv1 && mt1 == mv2) ||
                (mt2 == mv1 && mt3 == mv2) ||
                (mt3 == mv1 && mt1 == mv2) ||
                (mt3 == mv1 && mt2 == mv2)) {
                if (tri.v == v) {
                    //LOG(INFO) << "ZGX SKIPPED!!!";
                    continue;
                }

                //LOG(INFO) << "Not SKIPPED!";

                const Point3f &p0 = tri.mesh->p[t1];
                const Point3f &p1 = tri.mesh->p[t2];
                const Point3f &p2 = tri.mesh->p[t3];

                if (fabs(p0.x) > 0.9f || fabs(p1.x) > 0.9f || fabs(p2.x) > 0.9f) {
                    //LOG(INFO) << "In loop: hitten the box face. invalid!" << "  : p0 = " << p0 << ", p1 = " << p1 << ", p2 = " << p2;
                    continue;
                }

                Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
                Normal3f n = Normal3f(Normalize(Cross(dp02, dp12)));
                if (tri.reverseOrientation ^ tri.transformSwapsHandedness)
                    n = -n;

                //LOG(INFO) << "In edgeBoundary: n = " << n << ", zDir = " << zDir << ", p0 = " << p0 << ", p1 = " << p1 << ", p2 = " << p2 << ", Dot(zDir, n) = " << Dot(zDir, n);

                if (Dot(zDir, n) > 0) {
                    //LOG(INFO) << "Is edge boundary!!!";
                    return true;
                } else return false;
            }
        }
        return false;
    };

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

        // along z-axis
        Float z0 = Dot(dp0, zDir);
        Float z1 = Dot(dp1, zDir);
        Float z2 = Dot(dp2, zDir);

        // determine which edge is boundary
        bool b0 = edgeBoundary(triangles, triangle.v[1], triangle.v[2], zDir, triangle.v, triangle.mesh);
        bool b1 = edgeBoundary(triangles, triangle.v[2], triangle.v[0], zDir, triangle.v, triangle.mesh);
        bool b2 = edgeBoundary(triangles, triangle.v[0], triangle.v[1], zDir, triangle.v, triangle.mesh);

        //LOG(INFO) << "b0 = " << b0 << ", b1 = " << b1 << ", b2 = " << b2 << ", n = " << n << ", zDir = " << zDir << ", p0 = " << p0 << ", p1 = " << p1 << ", p2 = " << p2
        //    << ", dp0 = " << dp0 << ", dp1 = " << dp1 << ", dp2 = " << dp2 << ", a = " << a << ", b = " << b << ", c = " << c;

        // add triangle
        addtri(b2, b1, b0, a, b, c, z0, z1, z2, 0);
    }

    //LOG(INFO) << "Edges size: " << edges.size() << ", Ppl_A: " << Ppl_A;
    if (edgeNum == 0 || Ppl_A < 1e-2f) {
        Ppl_A = 0;
        return;
    }

    //LOG(INFO) << "edgeNum = " << edgeNum << ", when enabled.";

    enabled = true;

    // from fsdBSDF source code
    // Power in 0-th order lobe
    eavg /= sumPhat_j;
    const Float sigma_xi = sqrt(3) / (k * eavg);
    Sigmat_0 *= 6 * k * k / psi0;
    const auto Sigma0 = Sigmat_0 + Vector3f{ 1 / sqr(sigma_xi),0,1 / sqr(sigma_xi) };
    const auto detSigma0 = std::max<Float>(0, Sigma0.x * Sigma0.z - sqr(Sigma0.y));
    Ppl_0 = k * k / (18 * Pi) / sqrtf(detSigma0) * sqr(psi0);


}



Spectrum FsdBxDF::f(const Vector3f &wo, const Vector3f &wiWorld) const {
    const Vector3f wi = Normalize(wiWorld);
    Vector3f wiLocal(Dot(wi, xDir), Dot(wi, yDir), Dot(wi, zDir));
    if (wiLocal.z < 0) return 0;
    Vector2f wiScreen(-wiLocal.x, -wiLocal.y);
    //Vector2f wiScreen(wiLocal.x, wiLocal.y);

    if (Ppl_A <= 0)
        return 0;

    const auto cosine = Float(1) / sqrtf(1 + Dot(wiScreen, wiScreen));
    c_t bsdf = {};
    for (int i = 0; i < edgeNum; ++i) {
        const Edge &e = *edges[i];
        bsdf += Psihat(e.a, e.b, e.e, e.v, 2 * Pi / wavelength, wiScreen);
    }

    // handle sampled spectrum from the specific wavelength
    Float result = std::max(.0f, 1.f / (cosine * Ppl_A) * std::norm(bsdf));

    //LOG(INFO) << "FsdBxDF::f: " << "woWorld = " << wo << ", wiWorld = " << wiWorld << ", wiLocal = " << wiLocal << ", f = " << result;
    //SampledSpectrum ss(.0f);
    //ss.c[spectrumIndex] = 60.f * result;
    //return ss;
#ifdef PBRT_SAMPLED_SPECTRUM
    SampledSpectrum ss(.0f);
    ss.c[spectrumIndex] = 60.f * result;
    return ss;
#else
    return Spectrum(result);
#endif
}



Float FsdBxDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    const Vector3f wolocal(Dot(wi, xDir), Dot(wi, yDir), Dot(wi, zDir));
    if (wolocal.z <= 0) return 0;
    const auto xi = -Vector2f{ wolocal.x,wolocal.y };
    //const auto xi = Vector2f{ wolocal.x,wolocal.y };
    //LOG(INFO) << "FsdBxDF::Pdf: " << evalPdf(xi);
    return evalPdf(xi);
}



fsdPrecomputedTables FsdBxDF::tables;






//from source code
Vector2f FsdBxDF::sampleEdge(const Edge &e, Float &pdf) const {
    //const Matrix2x2 Xi = k * Matrix2x2(e.e.x, e.e.y, e.e.y, -e.e.x);
    //Matrix2x2 invXi = { 0,0,0,0 };
    //Xi.invert2x2(invXi);
    Float denominator = (e.e.x * e.e.x + e.e.y * e.e.y) * 2 * Pi / wavelength;
    Float invDenom = 1.0f / denominator;
    if (denominator == 0.) invDenom = 0.;


    const float A = std::norm(e.a - e.b), B = std::norm(e.a + e.b) / 4.f;
    const bool alpha1 = random_float() * (A + B) <= A;
    const auto rand3 = Point3f{ random_float(),random_float(),random_float() };
    Vector2f xi = alpha1 ?
        tables.importanceSampleCDF1(rand3) :
        tables.importanceSampleCDF2(rand3);
    xi = invDenom * Vector2f(e.e.x * xi.x + e.e.y * xi.y, e.e.y * xi.x - e.e.x * xi.y);
    //const Vector2f xi = invXi * (alpha1 ?
    //    tables.importanceSampleCDF1(rand3) :
    //    tables.importanceSampleCDF2(rand3));

    // PDF
    pdf = evalPdf(xi);
    //LOG(INFO) << "FsdBxDF::sampleEdge: xi = " << xi << ", pdf = " << pdf;
    return xi;
}






Spectrum FsdBxDF::Sample_f(const Vector3f &wo, Vector3f *wi,
    const Point2f &sample, Float *pdf,
    BxDFType *sampledType) const {
    Float bsdf;
    Vector2f xi;

    if (edgeNum == 1) {
        xi = sampleEdge(*edges[0], *pdf);
        bsdf = *pdf > 0 ? eval(xi) : .0f;
    } else {
        // SIR when multiple edges are present to improve IS quality
        static constexpr std::size_t N = 8;
        //Float bsdfs[N];
        //Float aggw[N];
        Float bsdfs[N] = { 0.f };
        Float aggw[N] = { 0.f };
        Vector2f xis[N];

        Float sumw = .0f;
        for (std::size_t n = 0; n < N; ++n) {
            const auto rand = sample.x * sumPhat_j;
            const auto &edgeit = std::lower_bound(edges, edges + edgeNum, rand,
                [](const auto &e, auto v) {
                    return (*e).PhatAccum < v;
                });
            const auto &e = (edgeit != edges + edgeNum) ? *edgeit : *(edges + edgeNum - 1);
            //const auto &e = (edgeit != edges + edgeNum) ? *edgeit : (edgeNum!=0?*(edges + edgeNum - 1):*edgeit);

            Float pdf;
            xis[n] = sampleEdge(*e, pdf);
            bsdfs[n] = pdf > 0 ? eval(xis[n]) : .0f;
            const auto w = pdf > 0 ? bsdfs[n] / pdf : .0f;
            sumw += w;
            aggw[n] = sumw;
        }
        const auto n = std::min<std::size_t>(N - 1,
            std::lower_bound(&aggw[0], &aggw[N], sample.y * sumw) - (&aggw[0])
            );
        xi = xis[n];
        bsdf = bsdfs[n];
        *pdf = sumw > 0 ? Float(N) * bsdf / sumw : .0f;
    }


    if (*pdf > 0) {
        // xi in exit frame
        //*wi = Vector3f{ xi.x,xi.y,sqrtf(std::max<Float>(0,1 - Dot(xi,xi))) };
        //*wi = xi.x * xDir + xi.y * yDir + sqrtf(std::max<Float>(0, 1 - Dot(xi, xi))) * zDir;
        *wi = -xi.x * xDir - xi.y * yDir + sqrtf(std::max<Float>(0, 1 - Dot(xi, xi))) * zDir;

        ////LOG(INFO) << "FsdBxDF::Sample_f: wi = " << *wi << ", pdf = " << *pdf;

        // handle sampled spectrum from the specific wavelength
        Float result = std::max(.0f, bsdf);

        //LOG(INFO) << "FsdBxDF::Sample_f: " << "woWorld = " << wo << ", wiWorld = " << *wi << ", Dot(woWorld, wiWorld) = " << Dot(wo, *wi) << ", f = " << result
        //    << ", xDir = " << xDir << ", yDir = " << yDir << ", zDir = " << zDir << ", p = " << isect;
#ifdef PBRT_SAMPLED_SPECTRUM
        SampledSpectrum ss(.0f);
        ss.c[spectrumIndex] = 60.f * result;
        return ss;
#else
        return Spectrum(result);
#endif
    }

    *wi = Vector3f(1.0f, 0.f, 0.f);
    return 0.;
}



} // namespace pbrt