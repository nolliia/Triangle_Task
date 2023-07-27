#include "kacem_triangle.h"
#include <iostream>
#include <cmath>
#include <limits>

namespace kacem {

    bool arePointsEqual(const double p1[3], const double p2[3]) {
        for (int i = 0; i < 3; ++i) {
            if (p1[i] != p2[i])
                return false;
        }
        return true;
    }

    void crossProduct(const double a[3], const double b[3], double result[3]) {
        result[0] = a[1] * b[2] - a[2] * b[1];
        result[1] = a[2] * b[0] - a[0] * b[2];
        result[2] = a[0] * b[1] - a[1] * b[0];
    }

    void subtractVectors(const double a[3], const double b[3], double result[3]) {
        for (int i = 0; i < 3; ++i) {
            result[i] = a[i] - b[i];
        }
    }

    bool isPointInsideTriangle2D(const double p[2], const double t[6]) {
        double x1 = t[0];
        double y1 = t[1];
        double x2 = t[2];
        double y2 = t[3];
        double x3 = t[4];
        double y3 = t[5];

        double denominator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));

        // Check if the denominator is close to zero (indicating a degenerate triangle)
        if (std::abs(denominator) < std::numeric_limits<double>::epsilon()) {
            return false;
        }

        double lambda1 = ((y2 - y3) * (p[0] - x3) + (x3 - x2) * (p[1] - y3)) / denominator;
        double lambda2 = ((y3 - y1) * (p[0] - x3) + (x1 - x3) * (p[1] - y3)) / denominator;
        double lambda3 = 1.0 - lambda1 - lambda2;

        // Check if the barycentric coordinates are within the valid range
        if (lambda1 >= 0.0 && lambda1 <= 1.0 &&
            lambda2 >= 0.0 && lambda2 <= 1.0 &&
            lambda3 >= 0.0 && lambda3 <= 1.0) {
            return true;
        }

        return false;
    }

    bool isPointInsideTriangle3D(const double p[3], const double t[9]) {
        double v0[3], v1[3], v2[3], n[3], u[3], v[3];

        // Extract the vertices of the triangle
        for (int i = 0; i < 3; ++i) {
            v0[i] = t[i];
            v1[i] = t[3 + i];
            v2[i] = t[6 + i];
        }

        // Compute the triangle normal
        subtractVectors(v1, v0, u);
        subtractVectors(v2, v0, v);
        crossProduct(u, v, n);

        // Check if the triangle is degenerate (zero area)
        double area = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if (area < std::numeric_limits<double>::epsilon()) {
            return false;
        }

        // Compute the barycentric coordinates of the point
        double uDotU = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        double uDotV = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        double vDotV = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        double uDotP = u[0] * (p[0] - v0[0]) + u[1] * (p[1] - v0[1]) + u[2] * (p[2] - v0[2]);
        double vDotP = v[0] * (p[0] - v0[0]) + v[1] * (p[1] - v0[1]) + v[2] * (p[2] - v0[2]);

        double denom = uDotU * vDotV - uDotV * uDotV;
        double lambda1 = (vDotV * uDotP - uDotV * vDotP) / denom;
        double lambda2 = (uDotU * vDotP - uDotV * uDotP) / denom;
        double lambda3 = 1.0 - lambda1 - lambda2;

        // Check if the barycentric coordinates are within the valid range
        if (lambda1 >= 0.0 && lambda1 <= 1.0 &&
            lambda2 >= 0.0 && lambda2 <= 1.0 &&
            lambda3 >= 0.0 && lambda3 <= 1.0) {
            return true;
        }

        return false;
    }

    bool segmentsIntersect(const double p1[3], const double p2[3], const double q1[3], const double q2[3]) {
        double p[3], q[3];
        for (int i = 0; i < 3; ++i) {
            p[i] = p2[i] - p1[i];
            q[i] = q2[i] - q1[i];
        }

        double epsilon = std::numeric_limits<double>::epsilon();
        double r[3];
        kacem::crossProduct(p, q, r);

        // Check if the line segments are parallel
        if (std::abs(r[0]) < epsilon && std::abs(r[1]) < epsilon && std::abs(r[2]) < epsilon) {
            // Check if the line segments are collinear
            double t = ((q1[0] - p1[0]) * p[0] + (q1[1] - p1[1]) * p[1] + (q1[2] - p1[2]) * p[2]) /
                (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            double u = ((q1[0] - p1[0]) * q[0] + (q1[1] - p1[1]) * q[1] + (q1[2] - p1[2]) * q[2]) /
                (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);

            if ((t >= 0.0 && t <= 1.0) || (u >= 0.0 && u <= 1.0)) {
                return true;
            }
        }
        else {
            double pq[3];
            kacem::subtractVectors(p1, q1, pq);
            double pqCrossP[3];
            kacem::crossProduct(pq, p, pqCrossP);
            double pqCrossQ[3];
            kacem::crossProduct(pq, q, pqCrossQ);
            double t = (pqCrossQ[0] * q[0] + pqCrossQ[1] * q[1] + pqCrossQ[2] * q[2]) /
                (r[0] * q[0] + r[1] * q[1] + r[2] * q[2]);
            double u = (pqCrossP[0] * p[0] + pqCrossP[1] * p[1] + pqCrossP[2] * p[2]) /
                (r[0] * p[0] + r[1] * p[1] + r[2] * p[2]);

            if (t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0) {
                return true;
            }
        }

        return false;
    }




    bool isTriangleInsideTriangle(const double t1[9], const double t2[9]) {
        // Check if all vertices of t1 are inside t2
        for (int i = 0; i < 3; ++i) {
            double p[3] = { t1[i * 3], t1[i * 3 + 1], t1[i * 3 + 2] };
            if (!isPointInsideTriangle3D(p, t2))
                return false;
        }

        return true;
    }

    bool intersection(double t1[9], double t2[9]) {
        // Check if the triangles share any vertices
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (arePointsEqual(&t1[i * 3], &t2[j * 3]))
                    return true;
            }
        }

        // Check if the triangles' edges intersect
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double p1[3] = { t1[i * 3], t1[i * 3 + 1], t1[i * 3 + 2] };
                double p2[3] = { t1[(i * 3 + 3) % 9], t1[(i * 3 + 4) % 9], t1[(i * 3 + 5) % 9] };
                double q1[3] = { t2[j * 3], t2[j * 3 + 1], t2[j * 3 + 2] };
                double q2[3] = { t2[(j * 3 + 3) % 9], t2[(j * 3 + 4) % 9], t2[(j * 3 + 5) % 9] };

                // Check if the line segments (p1, p2) and (q1, q2) intersect
                // Line Segment Intersection using Oriented Segments (LSIOS) test
                if (segmentsIntersect(p1, p2, q1, q2))
                    return true;
            }
        }

        // Check if a triangle lies completely inside the other triangle
        if (isTriangleInsideTriangle(t1, t2) || isTriangleInsideTriangle(t2, t1))
            return true;

        return false;
    }


} // namespace kacem

int main() {
    // Test the intersection function with sample triangles
    double t1[9] = { 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0, 0.0 };
    double t2[9] = { 3.0, 3.0, 0.0, 5.0, 3.0, 0.0, 3.0, 5.0, 1.0 };

    bool isIntersecting = kacem::intersection(t1, t2);
    std::cout << "Triangles intersect: " << (isIntersecting ? "true" : "false") << std::endl;

    return 0;
}
