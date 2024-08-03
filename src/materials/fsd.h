#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_FSD_H
#define PBRT_MATERIALS_FSD_H

// materials/fsd.h*
#include "pbrt.h"
#include "material.h"
#include "reflection.h"
#include "scene.h"

#include "rng.h"


namespace pbrt {


struct Edge;
class KdTreeAccel;
struct fsdPrecomputedTables;


// FsdMaterial Declarations
class FsdMaterial : public Material {
    /*
public:
    // FsdMaterial Public Methods
    FsdMaterial(const std::shared_ptr<Texture<Spectrum>> &Kd,
        const std::shared_ptr<Texture<Float>> &sigma,
        const std::shared_ptr<Texture<Float>> &bumpMap)
        : Kd(Kd), sigma(sigma), bumpMap(bumpMap) {}
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
        TransportMode mode,
        bool allowMultipleLobes, const Scene &scene) const;

private:
    // FsdMaterial Private Data
    std::shared_ptr<Texture<Spectrum>> Kd;
    std::shared_ptr<Texture<Float>> sigma, bumpMap;
    */


public:
    // FsdMaterial Public Methods
    FsdMaterial(const std::shared_ptr<Material> &material,
        const std::shared_ptr<Texture<Spectrum>> &scale)
        :m_Material(material), scale(scale) {}

    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
        TransportMode mode,
        bool allowMultipleLobes, const Scene &scene) const;

private:
    // FsdMaterial Private Data
    std::shared_ptr<Material> m_Material;
    std::shared_ptr<Texture<Spectrum>> scale;
    
};

FsdMaterial *CreateFsdMaterial(const TextureParams &mp, const std::shared_ptr<Material> &material);











// FsdBxDF class which handles BTDF calculation and sampling of diffraction effects.
class FsdBxDF : public BxDF {
public:
    FsdBxDF(const SurfaceInteraction &intr, const Scene &scene);

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const {
        return std::string("I'm FsdBxDF");
    }

    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi,
        const Point2f &sample, Float *pdf,
        BxDFType *sampledType = nullptr) const;

    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;


    // member utility functions
    Float eval(const Vector2f &xi) const;

    Float evalPdf(const Vector2f &xi) const;

    Vector2f sampleEdge(const Edge &e, Float &pdf) const;

    static fsdPrecomputedTables tables;


    std::vector<Edge> edges;
    Float wavelength;
    int spectrumIndex;
    Vector3f xDir, yDir, zDir;
    bool enabled = false;

    Float Ppl_A = .0f, sumPhat_j = .0f, Ppl_0 = .0f;
};













// Also create a FsdBSDF to rewrite functions in BSDF class.
// Rewritten functions include Sample_f, pdf, f .etc
class FsdBSDF : public BSDF {
public:
    // FsdBSDF Public Methods
    FsdBSDF(const SurfaceInteraction &si, const Scene &scene, const BSDF *mbsdf, MemoryArena &arena)
        :BSDF(si, eta), fsdBxDF(std::make_unique<FsdBxDF>(si, scene)) {
        int n = mbsdf->NumComponents();

        Float rgb[3] = { 1.0f, 1.0f, 1.0f };
        Spectrum whiteSpectrum = RGBSpectrum::FromRGB(rgb);
        for (int i = 0; i < n; ++i) {
            Add(ARENA_ALLOC(arena, ScaledBxDF)(mbsdf->bxdfs[i], whiteSpectrum));
        }
    }

    Spectrum f(const Vector3f &woW, const Vector3f &wiW,
        BxDFType flags) const override {
        
        // handle fsdBxDF
        if (fsdBxDF->enabled) { 
            LOG(INFO) << "in FsdBSDF::f enabled: Dot(woWorld, wiWorld) = " << Dot(woW, wiW) << ", f = " << fsdBxDF->f(woW, wiW);
            return fsdBxDF->f(woW, wiW); 
        }
        
        
        
        // normal case
        //ProfilePhase pp(Prof::BSDFEvaluation);
        Vector3f wi = WorldToLocal(wiW), wo = WorldToLocal(woW);
        if (wo.z == 0) return 0.;
        bool reflect = Dot(wiW, ng) * Dot(woW, ng) > 0;
        Spectrum f(0.f);
        for (int i = 0; i < nBxDFs; ++i)
            if (bxdfs[i]->MatchesFlags(flags) &&
                ((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
                    (!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))))
                f += bxdfs[i]->f(wo, wi);
        return f;
        
        
    }

    Spectrum Sample_f(const Vector3f &woWorld, Vector3f *wiWorld,
        const Point2f &u, Float *pdf, BxDFType type,
        BxDFType *sampledType) const override {

        
        // handle fsdBxDF
        if (fsdBxDF->enabled) {
            //Vector3f wi, wo = WorldToLocal(woWorld);
            //if (wo.z == 0) return 0.;
            Vector3f wiW;
            *pdf = 0;
            if (sampledType) *sampledType = BxDFType(BSDF_DIFFUSE | BSDF_TRANSMISSION);
            Spectrum f = fsdBxDF->Sample_f(woWorld, &wiW, u, pdf, sampledType);
            if (*pdf == 0) {
                if (sampledType) *sampledType = BxDFType(0);
                return 0.;
            }
            *wiWorld = wiW;
            LOG(INFO) << "in FsdBSDF::Sample_f enabled: Dot(woWorld, wiWorld) = " << Dot(woWorld, *wiWorld) << ", pdf = " << *pdf << ", f = " << f << ", ratio = " << f.y() / *pdf;
            return f;
        }
        
        
        

        // normal case
        //ProfilePhase pp(Prof::BSDFSampling);
        // Choose which _BxDF_ to sample
        int matchingComps = NumComponents(type);
        if (matchingComps == 0) {
            *pdf = 0;
            if (sampledType) *sampledType = BxDFType(0);
            return Spectrum(0);
        }
        int comp =
            std::min((int)std::floor(u[0] * matchingComps), matchingComps - 1);

        // Get _BxDF_ pointer for chosen component
        BxDF *bxdf = nullptr;
        int count = comp;
        for (int i = 0; i < nBxDFs; ++i)
            if (bxdfs[i]->MatchesFlags(type) && count-- == 0) {
                bxdf = bxdfs[i];
                break;
            }
        CHECK(bxdf != nullptr);
        VLOG(2) << "BSDF::Sample_f chose comp = " << comp << " / matching = " <<
            matchingComps << ", bxdf: " << bxdf->ToString();

        // Remap _BxDF_ sample _u_ to $[0,1)^2$
        Point2f uRemapped(std::min(u[0] * matchingComps - comp, OneMinusEpsilon),
            u[1]);

        // Sample chosen _BxDF_
        Vector3f wi, wo = WorldToLocal(woWorld);
        if (wo.z == 0) return 0.;
        *pdf = 0;
        if (sampledType) *sampledType = bxdf->type;
        Spectrum f = bxdf->Sample_f(wo, &wi, uRemapped, pdf, sampledType);
        VLOG(2) << "For wo = " << wo << ", sampled f = " << f << ", pdf = "
            << *pdf << ", ratio = " << ((*pdf > 0) ? (f / *pdf) : Spectrum(0.))
            << ", wi = " << wi;
        if (*pdf == 0) {
            if (sampledType) *sampledType = BxDFType(0);
            return 0;
        }
        *wiWorld = LocalToWorld(wi);

        // Compute overall PDF with all matching _BxDF_s
        if (!(bxdf->type & BSDF_SPECULAR) && matchingComps > 1)
            for (int i = 0; i < nBxDFs; ++i)
                if (bxdfs[i] != bxdf && bxdfs[i]->MatchesFlags(type))
                    *pdf += bxdfs[i]->Pdf(wo, wi);
        if (matchingComps > 1) *pdf /= matchingComps;

        // Compute value of BSDF for sampled direction
        if (!(bxdf->type & BSDF_SPECULAR)) {
            bool reflect = Dot(*wiWorld, ng) * Dot(woWorld, ng) > 0;
            f = 0.;
            for (int i = 0; i < nBxDFs; ++i)
                if (bxdfs[i]->MatchesFlags(type) &&
                    ((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
                        (!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))))
                    f += bxdfs[i]->f(wo, wi);
        }
        VLOG(2) << "Overall f = " << f << ", pdf = " << *pdf << ", ratio = "
            << ((*pdf > 0) ? (f / *pdf) : Spectrum(0.));
        return f;
    }

    Float Pdf(const Vector3f &woWorld, const Vector3f &wiWorld,
        BxDFType flags) const override {
        
        // handle fsdBxDF
        if (fsdBxDF->enabled) {
            //Vector3f wo = WorldToLocal(woWorld), wi = WorldToLocal(wiWorld);
            //if (wo.z == 0) return 0.;
            Float pdf = 0.f;
            pdf = fsdBxDF->Pdf(woWorld, wiWorld);
            LOG(INFO) << "in FsdBSDF::Pdf enabled: Dot(woWorld, wiWorld) = " << Dot(woWorld, wiWorld) << ", pdf = " << pdf;
            return pdf;
        }
        
        

        // normal case
        //ProfilePhase pp(Prof::BSDFPdf);
        if (nBxDFs == 0.f) return 0.f;
        Vector3f wo = WorldToLocal(woWorld), wi = WorldToLocal(wiWorld);
        if (wo.z == 0) return 0.;
        Float pdf = 0.f;
        int matchingComps = 0;
        for (int i = 0; i < nBxDFs; ++i)
            if (bxdfs[i]->MatchesFlags(flags)) {
                ++matchingComps;
                pdf += bxdfs[i]->Pdf(wo, wi);
            }
        Float v = matchingComps > 0 ? pdf / matchingComps : 0.f;
        return v;
    }


    // FsdBSDF Private Data
    std::unique_ptr<FsdBxDF> fsdBxDF;
};
















// Utility Functions






} // namespace pbrt




#endif  // PBRT_MATERIALS_Fsd_H