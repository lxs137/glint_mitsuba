//
// Created by lxs on 17-6-14.
//
#pragma once
#if !defined(__MITSUBA_PERSPECTIVE_H_)
#define __MITSUBA_PERSPECTIVE_H_

#include <mitsuba/render/common.h>
#include <mitsuba/render/film.h>
#include <mitsuba/render/emitter.h>

MTS_NAMESPACE_BEGIN


class PerspectiveCameraImpl : public PerspectiveCamera {
public:
    PerspectiveCameraImpl(const Properties &props);

    PerspectiveCameraImpl(Stream *stream, InstanceManager *manager);

    void configure();

    /**
     * \brief Compute the directional sensor response function
     * of the camera multiplied with the cosine foreshortening
     * factor associated with the image plane
     *
     * \param d
     *     A normalized direction vector from the aperture position to the
     *     reference point in question (all in local camera space)
     */
    inline Float importance(const Vector &d);

    Spectrum sampleRay(Ray &ray, const Point2 &pixelSample,
                       const Point2 &otherSample, Float timeSample);

    Spectrum sampleRayDifferential(RayDifferential &ray, const Point2 &pixelSample,
                                   const Point2 &otherSample, Float timeSample);

    Spectrum samplePosition(PositionSamplingRecord &pRec,
                            const Point2 &sample, const Point2 *extra);

    Spectrum evalPosition(const PositionSamplingRecord &pRec);

    Float pdfPosition(const PositionSamplingRecord &pRec);

    Spectrum sampleDirection(DirectionSamplingRecord &dRec,
                             PositionSamplingRecord &pRec,
                             const Point2 &sample, const Point2 *extra);

    Float pdfDirection(const DirectionSamplingRecord &dRec,
                       const PositionSamplingRecord &pRec) const {
        if (dRec.measure != ESolidAngle)
            return 0.0f;

        const Transform &trafo = m_worldTransform->eval(pRec.time);

        return importance(trafo.inverse()(dRec.d));
    }

    Spectrum evalDirection(const DirectionSamplingRecord &dRec,
                           const PositionSamplingRecord &pRec) const {
        if (dRec.measure != ESolidAngle)
            return Spectrum(0.0f);

        const Transform &trafo = m_worldTransform->eval(pRec.time);

        return Spectrum(importance(trafo.inverse()(dRec.d)));
    }

    bool getSamplePosition(const PositionSamplingRecord &pRec,
                           const DirectionSamplingRecord &dRec, Point2 &samplePosition) const {
        Transform invTrafo = m_worldTransform->eval(pRec.time).inverse();
        Point local(Point(invTrafo(dRec.d)));

        if (local.z <= 0)
            return false;

        Point screenSample = m_cameraToSample(local);
        if (screenSample.x < 0 || screenSample.x > 1 ||
            screenSample.y < 0 || screenSample.y > 1)
            return false;

        samplePosition = Point2(
                screenSample.x * m_resolution.x,
                screenSample.y * m_resolution.y);

        return true;
    }

    Spectrum sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample) const {
        const Transform &trafo = m_worldTransform->eval(dRec.time);

        /* Transform the reference point into the local coordinate system */
        Point refP = trafo.inverse().transformAffine(dRec.ref);

        /* Check if it is outside of the clip range */
        if (refP.z < m_nearClip || refP.z > m_farClip) {
            dRec.pdf = 0.0f;
            return Spectrum(0.0f);
        }

        Point screenSample = m_cameraToSample(refP);
        dRec.uv = Point2(screenSample.x, screenSample.y);
        if (dRec.uv.x < 0 || dRec.uv.x > 1 ||
            dRec.uv.y < 0 || dRec.uv.y > 1) {
            dRec.pdf = 0.0f;
            return Spectrum(0.0f);
        }

        dRec.uv.x *= m_resolution.x;
        dRec.uv.y *= m_resolution.y;

        Vector localD(refP);
        Float dist = localD.length(),
                invDist = 1.0f / dist;
        localD *= invDist;

        dRec.p = trafo.transformAffine(Point(0.0f));
        dRec.d = (dRec.p - dRec.ref) * invDist;
        dRec.dist = dist;
        dRec.n = trafo(Vector(0.0f, 0.0f, 1.0f));
        dRec.pdf = 1;
        dRec.measure = EDiscrete;

        return Spectrum(
                importance(localD) * invDist * invDist);
    }

    Float pdfDirect(const DirectSamplingRecord &dRec) const {
        return (dRec.measure == EDiscrete) ? 1.0f : 0.0f;
    }

    Transform getProjectionTransform(const Point2 &apertureSample,
                                     const Point2 &aaSample) const {
        Float right = std::tan(m_xfov * M_PI/360) * m_nearClip, left = -right;
        Float top = right / m_aspect, bottom = -top;

        Vector2 offset(
                (right-left)/m_film->getSize().x * (aaSample.x-0.5f),
                (top-bottom)/m_film->getSize().y * (aaSample.y-0.5f));

        return m_clipTransform *
               Transform::glFrustum(left+offset.x, right+offset.x,
                                    bottom+offset.y, top+offset.y, m_nearClip, m_farClip);
    }

    AABB getAABB() const {
        return m_worldTransform->getTranslationBounds();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "PerspectiveCamera[" << endl
            << "  fov = [" << getXFov() << ", " << getYFov() << "]," << endl
            << "  nearClip = " << m_nearClip << "," << endl
            << "  farClip = " << m_farClip << "," << endl
            << "  worldTransform = " << indent(m_worldTransform.toString()) << "," << endl
            << "  sampler = " << indent(m_sampler->toString()) << "," << endl
            << "  film = " << indent(m_film->toString()) << "," << endl
            << "  medium = " << indent(m_medium.toString()) << "," << endl
            << "  shutterOpen = " << m_shutterOpen << "," << endl
            << "  shutterOpenTime = " << m_shutterOpenTime << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    Transform m_cameraToSample;
    Transform m_sampleToCamera;
    Transform m_clipTransform;
    AABB2 m_imageRect;
    Float m_normalization;
    Vector m_dx, m_dy;
};

MTS_NAMESPACE_END

#endif //MITSUBA_PERSPECTIVE_H
