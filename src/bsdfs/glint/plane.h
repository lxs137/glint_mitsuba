//
// Created by lxs on 17-6-15.
//

#pragma once
#ifndef __MITSUBA_GLINT_PLANE_H
#define __MITSUBA_GLINT_PLANE_H

#define PLANE_EPSILON 1e-3

MTS_NAMESPACE_BEGIN

// texture所在的平面
class TexPlane {
public:
    TexPlane(): epsilon(PLANE_EPSILON) {
        p = Point();
        n = Normal();
    }
    TexPlane(const Point &pp, const Normal &nn): epsilon(PLANE_EPSILON) {
        p = pp, n = nn;
    }
    TexPlane(const Point &p1, const Point &p2, const Point &p3): epsilon(PLANE_EPSILON) {
        p = p1;
        Vector v1 = p1 - p2, v2 = p1 - p3;
        n = Normal(cross(v1, v2));
    }
    TexPlane(const TexPlane *plane): epsilon(PLANE_EPSILON) {
        p = plane->p, n = plane->n;
    }
    void add_tri(const Point *vertexs, const Point2 *texs, uint32_t *index)
    {
        tri_on_plane[0] = vertexs[index[0]];
        tri_on_plane[1] = vertexs[index[1]];
        tri_on_plane[2] = vertexs[index[2]];
        tri_tex[0] = texs[index[0]];
        tri_tex[1] = texs[index[1]];
        tri_tex[2] = texs[index[2]];
    }
    bool intersectP(const Ray &ray) const
    {
        float parallel_test = dot(ray.d, this->n);
        if(std::fabs(parallel_test) < 1e-5)
            return false;
        else
        {
            float t = dot((this->p - ray.o), this->n) / dot(ray.d, this->n);
            return (t > this->epsilon && t > ray.mint && t < ray.maxt);
        }
    }
    bool intersect(const Ray &ray, float &hit_t) const
    {
        float parallel_test = absDot(ray.d, this->n);
        // 通过平面法向量和光线的方向向量的点乘结果,判断两者是否平行
        if(std::fabs(parallel_test) < 1e-5)
            return false;
        else
        {
            float t = dot((this->p - ray.o), this->n) / dot(ray.d, this->n);
            if(t > this->epsilon && t > ray.mint && t < ray.maxt)
            {
                hit_t = t;
                return true;
            }
            else
                return false;
        }
    }
    void sample_tex_coord(Point &p, Point2 &tex) const
    {
        Vector v1 = tri_on_plane[1] - tri_on_plane[0], v2 = tri_on_plane[2] - tri_on_plane[0];
        float D = v1.x * v2.y - v1.y * v2.x,
              u = ((p.x - tri_on_plane[0].x) * v2.y - (p.y - tri_on_plane[0].y) * v2.x)/D,
              v = ((p.y - tri_on_plane[0].y) * v1.x - (p.x - tri_on_plane[0].x) * v1.y)/D;
        tex = tri_tex[0] + u * (tri_tex[1] - tri_tex[0]) + v * (tri_tex[2] - tri_tex[0]);
    }
private:
    Point p, tri_on_plane[3];
    Point2 tri_tex[3];
    Normal n;
    const float epsilon;
};

MTS_NAMESPACE_END

#endif //__MITSUBA_GLINT_PLANE_H
