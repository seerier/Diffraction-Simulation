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

namespace pbrt {

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
            bool allowMultipleLobes) const;

    private:
        // FsdMaterial Private Data
        std::shared_ptr<Texture<Spectrum>> Kd;
        std::shared_ptr<Texture<Float>> sigma, bumpMap;
        */


    public:
        // FsdMaterial Public Methods
        FsdMaterial(const std::shared_ptr<Material> &material,
                    const std::shared_ptr<Texture<Spectrum>> &scale)
            :m_Material(material),scale(scale){}

        void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
            TransportMode mode,
            bool allowMultipleLobes) const;

    private:
        // FsdMaterial Private Data
        std::shared_ptr<Material> m_Material;
        std::shared_ptr<Texture<Spectrum>> scale;
    };

    FsdMaterial *CreateFsdMaterial(const TextureParams &mp, const std::shared_ptr<Material> &material);











    // FsdBxDF class which handles BTDF calculation and sampling of diffraction effects.
    class FsdBxDF : public BxDF {

    };













    // Also create a FsdBSDF to rewrite functions in BSDF class.
    // Rewrited functions include Sample_f, pdf, f .etc
    class FsdBSDF : public BSDF {
        // FsdBSDF Public Methods
        FsdBSDF(const SurfaceInteraction &si, const std::shared_ptr<FsdBxDF> &fsdBxDF, Float eta = 1)
            :BSDF(si, eta), fsdBxDF(fsdBxDF){}


        // FsdBSDF Private Data
        std::shared_ptr<FsdBxDF> fsdBxDF;
    };

}  // namespace pbrt

#endif  // PBRT_MATERIALS_Fsd_H