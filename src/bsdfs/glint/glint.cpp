//
// Created by lxs on 17-6-12.
//

#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/plugin.h>
#include "../microfacet.h"
#include "../ior.h"
#include "perspective_glint.h"
#include "plane.h"

MTS_NAMESPACE_BEGIN

class Glint : public BSDF
{
public:
    Glint(const Properties &prop): BSDF(prop)
    {
        ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

        m_specularReflectance = new ConstantSpectrumTexture(
                prop.getSpectrum("specularReflectance", Spectrum(1.0f)));

        std::string materialName = prop.getString("material", "Cu");

        Spectrum intEta, intK;
        if (boost::to_lower_copy(materialName) == "none") {
            intEta = Spectrum(0.0f);
            intK = Spectrum(1.0f);
        } else {
            intEta.fromContinuousSpectrum(InterpolatedSpectrum(
                    fResolver->resolve("data/ior/" + materialName + ".eta.spd")));
            intK.fromContinuousSpectrum(InterpolatedSpectrum(
                    fResolver->resolve("data/ior/" + materialName + ".k.spd")));
        }

        Float extEta = lookupIOR(prop, "extEta", "air");

        m_eta = prop.getSpectrum("eta", intEta) / extEta;
        m_k   = prop.getSpectrum("k", intK) / extEta;

        MicrofacetDistribution distribution(prop, MicrofacetDistribution::EBeckmann);
        m_type = distribution.getType();
        if(m_type != MicrofacetDistribution::EBeckmann)
            SLog(EError, "Glint's microfacet must be Beckmann microfacet distribution.");
        m_sampleVisible = distribution.getSampleVisible();
        m_alphaU = new ConstantFloatTexture(distribution.getAlphaU());
        if (distribution.getAlphaU() == distribution.getAlphaV())
            m_alphaV = m_alphaU;
        else
            m_alphaV = new ConstantFloatTexture(distribution.getAlphaV());


//        Properties new_prop = Properties(prop);
//        new_prop.setPluginName("perspective");
//        camera_temp = static_cast<PerspectiveCamera*>(PluginManager::getInstance()->
//                createObject(MTS_CLASS(PerspectiveCamera), new_prop));
//        std::ostringstream oss;
//        oss << "Camera: " << camera_temp->toString() << endl;
//        SLog(EError, oss.str().c_str());

        camera = new PerspectiveCameraGlint(prop);
        camera->configure();

//        ref<const AnimatedTransform> trans =  prop.getAnimatedTransform("toWorld", Transform());
//        std::ostringstream oss;
//        oss << "to world: " << trans.toString() << "sample2camera: "
//            << endl << camera->m_sampleToCamera.toString() << endl;
//        SLog(EInfo, oss.str().c_str());
    }

    Glint(Stream *stream, InstanceManager *manager) : BSDF(stream, manager)
    {
        m_type = (MicrofacetDistribution::EType) stream->readUInt();
        m_sampleVisible = stream->readBool();
        m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
        m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
        m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_eta = Spectrum(stream);
        m_k = Spectrum(stream);

        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const
    {
        BSDF::serialize(stream, manager);

        stream->writeUInt((uint32_t) m_type);
        stream->writeBool(m_sampleVisible);
        manager->serialize(stream, m_alphaU.get());
        manager->serialize(stream, m_alphaV.get());
        manager->serialize(stream, m_specularReflectance.get());
        m_eta.serialize(stream);
        m_k.serialize(stream);
    }

    void configure() {
        unsigned int extraFlags = 0;
        if (m_alphaU != m_alphaV)
            extraFlags |= EAnisotropic;

        if (!m_alphaU->isConstant() || !m_alphaV->isConstant() ||
            !m_specularReflectance->isConstant())
            extraFlags |= ESpatiallyVarying;

        m_components.clear();
        m_components.push_back(EGlossyReflection | EFrontSide | extraFlags);

        /* Verify the input parameters and fix them if necessary */
        m_specularReflectance = ensureEnergyConservation(
                m_specularReflectance, "specularReflectance", 1.0f);

        m_usesRayDifferentials =
                m_alphaU->usesRayDifferentials() ||
                m_alphaV->usesRayDifferentials() ||
                m_specularReflectance->usesRayDifferentials();

        BSDF::configure();
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const
    {
        return Spectrum(0.5f);
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const
    {
        return 1.0f;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const
    {
        SLog(EInfo, "Sample Intersection: ", bRec.its.toString());
        MicrofacetDistribution distr(
                m_type,
                m_alphaU->eval(bRec.its).average(),
                m_alphaV->eval(bRec.its).average(),
                m_sampleVisible
        );
        Float pdf;
        Normal m = distr.sample(bRec.wi, sample, pdf);
        bRec.wo = (2 * dot(bRec.wi, m) * Vector(m) - bRec.wi);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EGlossyReflection;
        return Spectrum(0.5f);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const
    {
//        mitsuba::Thread *thread = mitsuba::Thread::getThread();
//        int thread_id = thread->getID() - 3;
        std::ostringstream oss;


        const Intersection &its = bRec.its;
        const TriMesh *mesh = static_cast<const TriMesh*>(its.shape);
        Triangle tri = mesh->getTriangles()[its.primIndex];

        Vector ray_d = its.toWorld(-its.wi), d_diff[4];
        Point ray_o = its.p - its.t * ray_d, p_film = ray_o + camera->getNearClip() * ray_d,
                sample_diff[4], intersect_diff[4];
        Point2 tex_diff[4];
        Ray *ray_diff[4];
        camera->get_sample_differential(p_film, sample_diff);
        bool all_intersect = true;
        TexPlane plane = TexPlane(its.p, Normal(its.toWorld(Vector(0, 0, 1.f))));
        // Sample texture coord for intersect diff
        plane.add_tri(mesh->getVertexPositions(), mesh->getVertexTexcoords(), tri.idx);
        for(int i = 0; i < 4; i++)
        {
            float t;
            d_diff[i] = (sample_diff[i] - ray_o) / camera->getNearClip();
            ray_diff[i] = new Ray(ray_o, d_diff[i], 0.f);
            if(plane.intersect(*ray_diff[i], t))
                intersect_diff[i] = (*ray_diff[i])(t);
            else {
                SLog(EWarn, "Shape is parallel with ray");
                return Spectrum(0.5f);
            }
            plane.sample_tex_coord(intersect_diff[i], tex_diff[i]);
        }
        AABB2 tex_box = AABB2(tex_diff[0]);
        tex_box.expandBy(tex_diff[1]);
        tex_box.expandBy(tex_diff[2]);
        tex_box.expandBy(tex_diff[3]);

        oss << "Hit: " << its.uv.toString() << endl
            << "Box: " << endl
            << tex_box.toString() << endl;
        SLog(EInfo, oss.str().c_str());


        MicrofacetDistribution distr(
                m_type,
                m_alphaU->eval(bRec.its).average(),
                m_alphaV->eval(bRec.its).average(),
                m_sampleVisible
        );
        Normal m = distr.sample(bRec.wi, sample, pdf);
        bRec.wo = (2 * dot(bRec.wi, m) * Vector(m) - bRec.wi);

        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EGlossyReflection;

        // Clear variable
        for(int i = 0; i < 4; i++)
            delete ray_diff[i];

        return Spectrum(0.5f);
    }

    void addChild(const std::string &name, ConfigurableObject *child)
    {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
            if (name == "alpha")
                m_alphaU = m_alphaV = static_cast<Texture *>(child);
            else if (name == "alphaU")
                m_alphaU = static_cast<Texture *>(child);
            else if (name == "alphaV")
                m_alphaV = static_cast<Texture *>(child);
            else if (name == "specularReflectance")
                m_specularReflectance = static_cast<Texture *>(child);
            else
                BSDF::addChild(name, child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    Float getRoughness(const Intersection &its, int component) const
    {
        return 0.5f * (m_alphaU->eval(its).average()
                       + m_alphaV->eval(its).average());
    }

    ~Glint()
    {
        delete camera;
    }

    std::string toString() const
    {
        std::ostringstream oss;
        oss << "Glint Microfacet[" << endl
            << "id = \"" << getID() << "\"," << endl
            << "distributon = \"EBeckmann\"]" << endl;
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;
    MTS_DECLARE_CLASS()
private:
    MicrofacetDistribution::EType m_type;
    ref<Texture> m_specularReflectance;
    ref<Texture> m_alphaU, m_alphaV;
    bool m_sampleVisible;
    Spectrum m_eta, m_k;

    PerspectiveCameraGlint *camera;

    PerspectiveCamera *camera_temp;
};


class GlintShader : public Shader
{
public:
    GlintShader(Renderer *renderer, const Texture *specularReflectance,
                const Texture *alphaU, const Texture *alphaV, const Spectrum &eta,
                const Spectrum &k) : Shader(renderer, EBSDFShader),
                                     m_specularReflectance(specularReflectance),
                                     m_alphaU(alphaU), m_alphaV(alphaV) {
        m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
        m_alphaUShader = renderer->registerShaderForResource(m_alphaU.get());
        m_alphaVShader = renderer->registerShaderForResource(m_alphaV.get());

        /* Compute the reflectance at perpendicular incidence */
        m_R0 = fresnelConductorExact(1.0f, eta, k);
    }
    bool isComplete() const {
        return m_specularReflectanceShader.get() != NULL &&
               m_alphaUShader.get() != NULL &&
               m_alphaVShader.get() != NULL;
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_specularReflectanceShader.get());
        deps.push_back(m_alphaUShader.get());
        deps.push_back(m_alphaVShader.get());
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_specularReflectance.get());
        renderer->unregisterShaderForResource(m_alphaU.get());
        renderer->unregisterShaderForResource(m_alphaV.get());
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_R0", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
        program->setParameter(parameterIDs[0], m_R0);
    }

    void generateCode(std::ostringstream &oss,
                      const std::string &evalName,
                      const std::vector<std::string> &depNames) const {
        oss << "uniform vec3 " << evalName << "_R0;" << endl
            << endl
            << "float " << evalName << "_D(vec3 m, float alphaU, float alphaV) {" << endl
            << "    float ct = cosTheta(m), ds = 1-ct*ct;" << endl
            << "    if (ds <= 0.0)" << endl
            << "        return 0.0f;" << endl
            << "    alphaU = 2 / (alphaU * alphaU) - 2;" << endl
            << "    alphaV = 2 / (alphaV * alphaV) - 2;" << endl
            << "    float exponent = (alphaU*m.x*m.x + alphaV*m.y*m.y)/ds;" << endl
            << "    return sqrt((alphaU+2) * (alphaV+2)) * 0.15915 * pow(ct, exponent);" << endl
            << "}" << endl
            << endl
            << "float " << evalName << "_G(vec3 m, vec3 wi, vec3 wo) {" << endl
            << "    if ((dot(wi, m) * cosTheta(wi)) <= 0 || " << endl
            << "        (dot(wo, m) * cosTheta(wo)) <= 0)" << endl
            << "        return 0.0;" << endl
            << "    float nDotM = cosTheta(m);" << endl
            << "    return min(1.0, min(" << endl
            << "        abs(2 * nDotM * cosTheta(wo) / dot(wo, m))," << endl
            << "        abs(2 * nDotM * cosTheta(wi) / dot(wi, m))));" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_schlick(float ct) {" << endl
            << "    float ctSqr = ct*ct, ct5 = ctSqr*ctSqr*ct;" << endl
            << "    return " << evalName << "_R0 + (vec3(1.0) - " << evalName << "_R0) * ct5;" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "   if (cosTheta(wi) <= 0 || cosTheta(wo) <= 0)" << endl
            << "    	return vec3(0.0);" << endl
            << "   vec3 H = normalize(wi + wo);" << endl
            << "   vec3 reflectance = " << depNames[0] << "(uv);" << endl
            << "   float alphaU = max(0.2, " << depNames[1] << "(uv).r);" << endl
            << "   float alphaV = max(0.2, " << depNames[2] << "(uv).r);" << endl
            << "   float D = " << evalName << "_D(H, alphaU, alphaV)" << ";" << endl
            << "   float G = " << evalName << "_G(H, wi, wo);" << endl
            << "   vec3 F = " << evalName << "_schlick(1-dot(wi, H));" << endl
            << "   return reflectance * F * (D * G / (4*cosTheta(wi)));" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
            << "    	return vec3(0.0);" << endl
            << "    return " << evalName << "_R0 * inv_pi * inv_pi * cosTheta(wo);"<< endl
            << "}" << endl;
    }
    MTS_DECLARE_CLASS()
private:
    ref<const Texture> m_specularReflectance;
    ref<const Texture> m_alphaU;
    ref<const Texture> m_alphaV;
    ref<Shader> m_specularReflectanceShader;
    ref<Shader> m_alphaUShader;
    ref<Shader> m_alphaVShader;
    Spectrum m_R0;
};

Shader *Glint::createShader(Renderer *renderer) const {
    return new GlintShader(renderer,m_specularReflectance.get(), m_alphaU.get(),
                           m_alphaV.get(), m_eta, m_k);
}

MTS_IMPLEMENT_CLASS(GlintShader, false, Shader);
MTS_IMPLEMENT_CLASS_S(Glint, false, BSDF);
MTS_EXPORT_PLUGIN(Glint, "Glint BRDF");
MTS_NAMESPACE_END
