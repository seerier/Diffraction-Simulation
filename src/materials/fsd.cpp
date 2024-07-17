// materials/fsd.cpp*
#include "materials/fsd.h"
#include "paramset.h"
#include "reflection.h"
#include "interaction.h"
#include "texture.h"
#include "interaction.h"

#include "matte.h"
#include "textures/constant.h"

namespace pbrt {

    // FsdMaterial Method Definitions
    void FsdMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
        MemoryArena &arena,
        TransportMode mode,
        bool allowMultipleLobes) const {
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

        Spectrum s = scale->Evaluate(*si).Clamp();

        m_Material->ComputeScatteringFunctions(si, arena, mode, allowMultipleLobes);

        int n = si->bsdf->NumComponents();
        for (int i = 0; i < n; ++i) {
            si->bsdf->bxdfs[i] = ARENA_ALLOC(arena, ScaledBxDF)(si->bsdf->bxdfs[i], s);
        }

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

        std::shared_ptr<Texture<Spectrum>> scale = mp.GetSpectrumTexture("scale", Spectrum(0.3f));

        return new FsdMaterial(material, scale);
    }

}  // namespace pbrt