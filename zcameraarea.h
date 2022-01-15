/*-
    Porting of CameraCalculator.java

    This code is a C++ porting of the java code:
        https://github.com/zelenmi6/thesis/blob/master/src/geometry/CameraCalculator.java
        @author Milan Zelenka

    using Eigen
*/
#ifndef ZCAMERAAREA_H
#define ZCAMERAAREA_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>

#include <Eigen/Dense>

namespace zCameraAreaTester
{
    class zCameraArea
    {
    public:
        /**
         * Get corners of the polygon captured by the camera on the ground. 
         * The calculations are performed in the axes origin (0, 0, altitude)
         * and the points are not yet translated to camera's X-Y coordinates.
         * @param FOVh Horizontal field of view in degree
         * @param FOVv Vertical field of view in degree
         * @param altitude Altitude of the camera in meters
         * @param heading Heading of the camera (z axis) in degree
         * @param roll Roll of the camera (x axis) in degree
         * @param pitch Pitch of the camera (y axis) in degree
         */
        zCameraArea(
            double hFOV, double vFOV, double anglelimit,
            double roll, double pitch, double yaw,
            double altitude )
        {

            vFOV = vFOV*M_PI/180.0;
            hFOV = hFOV*M_PI/180.0;

            Eigen::Vector3d ray1(tan(hFOV/2), tan(vFOV/2), -1);
            Eigen::Vector3d ray2(tan(hFOV/2), -tan(vFOV/2), -1);
            Eigen::Vector3d ray3(-tan(hFOV/2), -tan(vFOV/2), -1);
            Eigen::Vector3d ray4(-tan(hFOV/2), tan(vFOV/2), -1);
            ray1.normalize();
            ray2.normalize();
            ray3.normalize();
            ray4.normalize();

            Eigen::Vector3d rotatedVectors[4];
            Eigen::Vector3d rays[4];
            rays[0] = ray1;
            rays[1] = ray2;
            rays[2] = ray3;
            rays[3] = ray4;
            roll = limit(roll, -anglelimit, anglelimit);
            pitch = limit(pitch, -anglelimit, anglelimit);
            rotateRays(rays,roll,pitch,yaw,rotatedVectors);

            Eigen::Vector3d origin(0,0,altitude);
            Eigen::Vector3d interx[4];
            getRayGroundIntersections(rotatedVectors,origin,interx,4);

        }
        
        void printVector(Eigen::Vector3d v, const char *name)
        {
            printf("%s: v.x %lf v.y %lf v.z %lf\n",
                name, v.x(),v.y(),v.z());

        }

        void rotateRays(
            Eigen::Vector3d *rays,
            double roll, double pitch, double yaw,
            Eigen::Vector3d *rayArray) {

            roll = roll*M_PI/180.0;
            pitch = pitch*M_PI/180.0;
            yaw = yaw*M_PI/180.0;

            double sinAlpha = sin(yaw);
            double sinBeta = sin(pitch);
            double sinGamma = sin(roll);
            double cosAlpha = cos(yaw);
            double cosBeta = cos(pitch);
            double cosGamma = cos(roll);
            double m00 = cosAlpha * cosBeta;
            double m01 = cosAlpha * sinBeta * sinGamma - sinAlpha * cosGamma;
            double m02 = cosAlpha * sinBeta * cosGamma + sinAlpha * sinGamma;
            double m10 = sinAlpha * cosBeta;
            double m11 = sinAlpha * sinBeta * sinGamma + cosAlpha * cosGamma;
            double m12 = sinAlpha * sinBeta * cosGamma - cosAlpha * sinGamma;
            double m20 = -sinBeta;
            double m21 = cosBeta * sinGamma;
            double m22 = cosBeta * cosGamma;
            
            Eigen::MatrixXd mRot(3,3);
            Rot << m00,m01,m02,m10,m11,m12,m20,m21,m22;

            rayArray[0] = (mRot*rays[0]);
            rayArray[1] = (mRot*rays[1]);
            rayArray[2] = (mRot*rays[2]);
            rayArray[3] = (mRot*rays[3]);

        }

        /**
         * Finds a ray-vector's intersection with the ground approximated by a plane
         * @param ray Ray-vector
         * @param origin Camera's position
         * @return
         */
        Eigen::Vector3d findRayGroundIntersection(Eigen::Vector3d ray, Eigen::Vector3d origin) {
            
            // Parametric form of an equation
            // P = origin + vector * t
            Eigen::Vector2d x(origin.x(),ray.x());
            Eigen::Vector2d y(origin.y(),ray.y());
            Eigen::Vector2d z(origin.z(),ray.z());
            
            // Equation of the horizontal plane (ground)
            // -z = 0
            
            // Calculate t by substituting z
            double t = - (z.x() / z.y());
            
            // Substitute t in the original parametric equations to get points of intersection
            Eigen::Vector3d ret(x.x() + x.y() * t, y.x() + y.y() * t, z.x() + z.y() * t);
            
            return ret;
        }

        /**
         * Finds the intersections of the camera's ray-vectors 
         * and the ground approximated by a horizontal plane
         * @param rays Array of 4 ray-vectors
         * @param origin Position of the camera. The computation were developed 
         * assuming the camera was at the axes origin (0, 0, altitude) and the 
         * results translated by the camera's real position afterwards.
         * @return
         */
        void getRayGroundIntersections(
            Eigen::Vector3d *rays, Eigen::Vector3d origin,
            Eigen::Vector3d *intersections,
            int length) {
            for (int i = 0; i < length; i ++) {
                intersections[i] = findRayGroundIntersection(rays[i], origin);
            }
        }

        double limit(double value, double min, double max)
        {
            if (value < min) { return min; }
            if (value > max) { return max; }
            return value;
        }

    };
}


#endif

