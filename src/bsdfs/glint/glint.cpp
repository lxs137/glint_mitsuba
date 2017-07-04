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

#define SET_P_SIZE 10000
#define SPATIAL_SAMPLE_NUM 100000000
#define DIRECTION_SAMPLE_NUM 100000000
#define EPSILON_COUNT 0.1f
#define GAMMA_DIRECTION 0.087f
#define PI 3.141593f
#define PI_2 1.570795f
#define _2_PI 6.28318f
#define Spherical2Cartesian(dir) Point3(cosf(dir.x)*sinf(dir.y),sinf(dir.x)*sinf(dir.y),cosf(dir.y))


MTS_NAMESPACE_BEGIN

struct DirectionTri
{
    Point2 vertices[3];
    DirectionTri(const Point2 &p1, const Point2 &p2, const Point2 &p3) {
        vertices[0] = p1, vertices[1] = p2, vertices[2] = p3;
    }
    Point2 edge_center(int index) const
    {
        int i = index, j = (index + 1) % 3;
        if(theta_is_PI_2(i))
            return Point2(vertices[j].x, 0.5f * (vertices[i].y + vertices[j].y));
        else if(theta_is_PI_2(j))
            return Point2(vertices[i].x, 0.5f * (vertices[i].y + vertices[j].y));
        else
            return Point2(0.5f * (vertices[i].x + vertices[j].x), 0.5f * (vertices[i].y + vertices[j].y));
    }
    bool theta_is_PI_2(int index) const {
        return fabsf(vertices[index].y - PI_2) < 1e-4;
    }
    Point2 operator[](int i) const
    {
        if(i < 0 || i > 3)
            return vertices[0];
        else
            return vertices[i];
    }
};

typedef std::pair<AABB2, uint32_t> SpatialNode;
typedef std::pair<DirectionTri, uint32_t> DirectionNode;
struct ConicQuery{
    Vector wi, wo;
    Float Gamma;
    Matrix3x3 C;

    ConicQuery(const Vector &_wi, const Vector &_wo, Float _Gamma): wi(_wi), wo(_wo), Gamma(_Gamma)
    {
        Vector x = normalize(cross(wi, wo)), y = normalize(wi - wo), z = normalize(wi + wo);
        Float lambda1 = (dot(wi, wo) + cosf(Gamma)) / (1 - cosf(Gamma)), lambda2 = 1.f / (tanf(Gamma/2)*tanf(Gamma/2));
        Matrix3x3 Q(x, y, z), A(lambda1, 0.f, 0.f, 0.f, lambda2, 0.f, 0.f, 0.f, -1.f), Q_transpose;
        Q.transpose(Q_transpose);
        C = Q;
        C *= A;
        C *= Q_transpose;
    }

    bool is_in(const Point2 &direction) const
    {
        Vector _m = Vector(Spherical2Cartesian(direction)), temp = C.preMult(_m);
        return dot(temp, _m) <= 0.f;
    }

    bool is_intersect(const DirectionTri &_tri) const
    {
        Point3 tri[3] = {Spherical2Cartesian(_tri[0]), Spherical2Cartesian(_tri[1]), Spherical2Cartesian(_tri[2])};
        Vector c, d;
        Float param_a, param_b, param_c, axis;
        for(int i = 0; i < 3; i++)
        {
            c = tri[i] + tri[(i + 1) % 3] - Point3(0.f, 0.f, 0.f);
            d = tri[i] - tri[(i + 1) % 3];
            param_a = dot(C.preMult(d), d);
            param_b = 2 * dot(C.preMult(d), c);
            param_c = dot(C.preMult(c), c);
            if((param_a - param_b + param_c) * (param_a + param_b + param_c) < 0.f)
                return true;
            else
            {
                axis = - param_b / (2 * param_a);
                if(axis < -1.f || axis > 2.f || param_a * param_c - param_b * param_b * 0.25f >= 0.f)
                    continue;
                else
                    return true;
            }
        }
        return false;
    }

    bool is_overlap(const DirectionTri &_tri) const
    {
        if(is_intersect(_tri))
            return true;
        else
            return is_in(_tri[0]) && is_in(_tri[1]) && is_in(_tri[2]);
    }

    bool is_contain(const DirectionTri &_tri) const
    {
        return !is_intersect(_tri) && (is_in(_tri[0]) && is_in(_tri[1]) && is_in(_tri[2]));
    }


};

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

        ref<Random> random = new Random((uint64_t)time(NULL));
        Float x, y, z, w, u, variance;
        std::ostringstream oss;
        for(int i = 0; i < SET_P_SIZE;)
        {
            x = random->nextFloat();
            y = (1 - x) * random->nextFloat();
            z = (1 - x - y) * random->nextFloat();
            w = 1 - x - y - z;
            u = (x + y + z + w) / 4;
            variance = ((x - u)*(x - u) + (y - u)*(y - u) + (z - u)*(z - u) + (w - u)*(w - u)) / 4;
            Float max = 1.f - 0.96f * powf(1.1f, -i);
            if(variance < max) {
                set_p.push_back(Point4(x, y, z, w));
                i++;
            }
        }
        SLog(EInfo, "Generate set points");

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
        /* Stop if this component was not requested */
        if (measure != ESolidAngle ||
            Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 ||
            ((bRec.component != -1 && bRec.component != 0) ||
             !(bRec.typeMask & EGlossyReflection)))
            return Spectrum(0.0f);


        std::ostringstream oss;
        const Intersection &its = bRec.its;
        const TriMesh *mesh = static_cast<const TriMesh*>(its.shape);
        Triangle tri = mesh->getTriangles()[its.primIndex];

        Vector ray_d = its.toWorld(-its.wi), d_diff[4];
        Point ray_o = its.p - its.t * ray_d, p_film = ray_o + camera->getNearClip() * ray_d,
                sample_diff[4], intersect_diff[4];
        Point2 tex_diff[4];
        Ray ray_diff[4];
        camera->get_sample_differential(p_film, sample_diff);
        TexPlane plane = TexPlane(its.p, Normal(its.toWorld(Vector(0, 0, 1.f))));
        // Sample texture coord for intersect diff
        plane.add_tri(mesh->getVertexPositions(), mesh->getVertexTexcoords(), tri.idx);
        for(int i = 0; i < 4; i++)
        {
            float t;
            d_diff[i] = (sample_diff[i] - ray_o) / camera->getNearClip();
            ray_diff[i] = Ray(ray_o, d_diff[i], 0.f);
            if(plane.intersect(ray_diff[i], t))
                intersect_diff[i] = (ray_diff[i])(t);
            else {
                SLog(EWarn, "Triangle is parallel with ray.");
                return Spectrum(0.5f);
            }
            plane.sample_tex_coord(intersect_diff[i], tex_diff[i]);
        }
        AABB2 tex_box = AABB2(tex_diff[0]);
        tex_box.expandBy(tex_diff[1]);
        tex_box.expandBy(tex_diff[2]);
        tex_box.expandBy(tex_diff[3]);

        /* Calculate the reflection half-vector */
        Vector H = normalize(bRec.wo+bRec.wi);

        /* Construct the microfacet distribution matching the
           roughness values at the current surface position. */
        MicrofacetDistribution distr(
                m_type,
                m_alphaU->eval(bRec.its).average(),
                m_alphaV->eval(bRec.its).average(),
                m_sampleVisible
        );
//        Float D = count_spatial(tex_box);
//        D *= count_direction(bRec.wi, bRec.wo);
        Float D = count_direction(bRec.wi, bRec.wo);
        oss << D <<std::endl;
        SLog(EInfo, oss.str().c_str());

        /* Evaluate the microfacet normal distribution */
        if (D == 0)
            return Spectrum(0.0f);

        /* Fresnel factor */
        const Spectrum F = fresnelConductorExact(dot(bRec.wi, H), m_eta, m_k) *
                           m_specularReflectance->eval(bRec.its);

        /* Smith's shadow-masking function */
        const Float G = distr.G(bRec.wi, bRec.wo, H);

        /* Calculate the total amount of reflection */
        Float model = dot(bRec.wi, H) * D * G /
                (tex_box.getSurfaceArea() * PI * (1 - cosf(GAMMA_DIRECTION)) * Frame::cosTheta(bRec.wi));

        return F * model;
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
//        std::ostringstream oss;
//
//        const Intersection &its = bRec.its;
//        const TriMesh *mesh = static_cast<const TriMesh*>(its.shape);
//        Triangle tri = mesh->getTriangles()[its.primIndex];
//
//        Vector ray_d = its.toWorld(-its.wi), d_diff[4];
//        Point ray_o = its.p - its.t * ray_d, p_film = ray_o + camera->getNearClip() * ray_d,
//                sample_diff[4], intersect_diff[4];
//        Point2 tex_diff[4];
//        Ray ray_diff[4];
//        camera->get_sample_differential(p_film, sample_diff);
//        TexPlane plane = TexPlane(its.p, Normal(its.toWorld(Vector(0, 0, 1.f))));
//        // Sample texture coord for intersect diff
//        plane.add_tri(mesh->getVertexPositions(), mesh->getVertexTexcoords(), tri.idx);
//        for(int i = 0; i < 4; i++)
//        {
//            float t;
//            d_diff[i] = (sample_diff[i] - ray_o) / camera->getNearClip();
//            ray_diff[i] = Ray(ray_o, d_diff[i], 0.f);
//            if(plane.intersect(ray_diff[i], t))
//                intersect_diff[i] = (ray_diff[i])(t);
//            else {
//                SLog(EWarn, "Triangle is parallel with ray.");
//                return Spectrum(0.5f);
//            }
//            plane.sample_tex_coord(intersect_diff[i], tex_diff[i]);
//        }
//        AABB2 tex_box = AABB2(tex_diff[0]);
//        tex_box.expandBy(tex_diff[1]);
//        tex_box.expandBy(tex_diff[2]);
//        tex_box.expandBy(tex_diff[3]);

        MicrofacetDistribution distr(
                m_type,
                m_alphaU->eval(bRec.its).average(),
                m_alphaV->eval(bRec.its).average(),
                m_sampleVisible
        );
        Normal m = distr.sample(bRec.wi, sample, pdf);

        bRec.wo = 2 * dot(bRec.wi, m) * Vector(m) - bRec.wi;
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EGlossyReflection;

//        Float D = count_spatial(tex_box);
//        D *= count_direction(bRec.wi, bRec.wo);

//        return Spectrum(0.2f * (tex_box.getSurfaceArea() / 0.005f));
        return Spectrum(0.2f);
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

    Float count_spatial(const AABB2 &query) const
    {
        uint32_t count = 0, split_times = 0;
        std::vector<SpatialNode> queue;
        queue.push_back(std::make_pair(AABB2(Point2(0.f, 0.f), Point2(1.f, 1.f)), SPATIAL_SAMPLE_NUM));
        SpatialNode node;
        AABB2 clip_box;
        Float p;

        while(!queue.empty())
        {
            node = queue.back();
            queue.pop_back();
            if(!node.first.overlaps(query) || node.second <= 0) {
                // Pass
            }
            else if(query.contains(node.first)) {
                count += node.second;
            }
            else {
                p = count / SPATIAL_SAMPLE_NUM;
                if(sqrt(node.second * p * (1 - p)) < EPSILON_COUNT * count) {
                    clip_box = node.first;
                    clip_box.clip(query);
                    count += (uint32_t)(node.second * clip_box.getSurfaceArea() / node.first.getSurfaceArea());
                }
                else {
                    assert(split_times < SET_P_SIZE);
                    Point2 center = node.first.getCenter(), min = node.first.min, max = node.first.max;
                    queue.push_back(std::make_pair(
                            AABB2(min, center),
                            (uint32_t)(node.second * set_p[split_times].x)));
                    queue.push_back(std::make_pair(
                            AABB2(Point2(center.x, min.y), Point2(max.x, center.y)),
                            (uint32_t)(node.second * set_p[split_times].y)));
                    queue.push_back(std::make_pair(
                            AABB2(center, max),
                            (uint32_t)(node.second * set_p[split_times].z)));
                    queue.push_back(std::make_pair(
                            AABB2(Point2(min.x, center.y), Point2(center.x, max.y)),
                            (uint32_t)(node.second * set_p[split_times].w)));
                    split_times++;
                }
            }
        }

        return (Float)count / SPATIAL_SAMPLE_NUM;
    }

    Float count_direction(const Vector &wi, const Vector &wo) const
    {
        uint32_t count = 0, split_times = 1;
        std::vector<DirectionNode> queue;
        ConicQuery query(wi, wo, GAMMA_DIRECTION);
        queue.push_back(std::make_pair(DirectionTri(Point2(0.f, 0.f), Point2(PI_2, 0.f), Point2(0.f, PI_2)),
                                       (uint32_t)(DIRECTION_SAMPLE_NUM * set_p[0].x)));
        queue.push_back(std::make_pair(DirectionTri(Point2(PI_2, 0.f), Point2(PI, 0.f), Point2(0.f, PI_2)),
                                       (uint32_t)(DIRECTION_SAMPLE_NUM * set_p[0].y)));
        queue.push_back(std::make_pair(DirectionTri(Point2(PI, 0.f), Point2(PI_2 + PI, 0.f), Point2(0.f, PI_2)),
                                       (uint32_t)(DIRECTION_SAMPLE_NUM * set_p[0].z)));
        queue.push_back(std::make_pair(DirectionTri(Point2(PI + PI_2, 0.f), Point2(_2_PI, 0.f), Point2(0.f, PI_2)),
                                       (uint32_t)(DIRECTION_SAMPLE_NUM * set_p[0].w)));
        bool is_intersect;
        Float p;

        while(!queue.empty())
        {
            DirectionNode node = queue.back();
            queue.pop_back();
            is_intersect = query.is_intersect(node.first);
            if((!is_intersect && !(query.is_in(node.first[0]) && query.is_in(node.first[1]) && query.is_in(node.first[2])))
               || node.second <= 0) {
                // pass
            }
            else if(!is_intersect && (query.is_in(node.first[0]) && query.is_in(node.first[1]) && query.is_in(node.first[2]))) {
                count += node.second;
            }
            else
            {
                p = count / SPATIAL_SAMPLE_NUM;
                if(sqrt(node.second * p * (1 - p)) < EPSILON_COUNT * count) {
                    count += node.second;
                }
                else {
                    assert(split_times < SET_P_SIZE);
                    queue.push_back(std::make_pair(
                            DirectionTri(node.first.edge_center(0), node.first.edge_center(1), node.first.edge_center(2)),
                            (uint32_t)(node.second * set_p[split_times].x)));
                    queue.push_back(std::make_pair(
                            DirectionTri(node.first.edge_center(0), node.first.edge_center(2), node.first[0]),
                            (uint32_t)(node.second * set_p[split_times].y)));
                    queue.push_back(std::make_pair(
                            DirectionTri(node.first.edge_center(0), node.first.edge_center(1), node.first[1]),
                            (uint32_t)(node.second * set_p[split_times].z)));
                    queue.push_back(std::make_pair(
                            DirectionTri(node.first.edge_center(1), node.first.edge_center(2), node.first[2]),
                            (uint32_t)(node.second * set_p[split_times].w)));
                    split_times++;
                }
            }
        }
        return (Float)count / DIRECTION_SAMPLE_NUM;
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
    std::vector<Point4> set_p;
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
