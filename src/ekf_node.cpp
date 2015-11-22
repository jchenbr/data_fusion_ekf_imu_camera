#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/Range.h>
#include <nav_msgs/Odometry.h>
#include <Eigen/Eigen>

#include <cmath>

using namespace std;
using namespace Eigen;

ros::Publisher odom_pub;
ros::Publisher msm_pub;
const double _PI = acos(-1);

#ifndef _EKF_STATE_DIM_NAME
#define _EKF_STATE_DIM_NAME

#define _STT_LEN 15

#define _POS_BEG 0
#define _POS_LEN 3
#define _POS_END 3

#define _ORT_BEG 3
#define _ORT_LEN 3
#define _ORT_END 6

#define _VEL_BEG 6
#define _VEL_LEN 3
#define _VEL_END 9

#define _B_G_BEG 9
#define _B_G_LEN 3
#define _B_G_END 12

#define _B_A_BEG 12
#define _B_A_LEN 3
#define _B_A_END 15

#define _DIM_X 0
#define _DIM_Y 1
#define _DIM_Z 2
#define _DIM_W 3
#define _DIM_LEN 3


#define _P_M_BEG 0
#define _P_M_LEN 3
#define _P_M_END 3

#define _O_M_BEG 3
#define _O_M_LEN 3
#define _O_M_END 6

#define _N_G_BEG 0
#define _N_G_LEN 3
#define _N_A_BEG 3
#define _N_A_LEN 3
#define _NBG_BEG 6
#define _NBG_LEN 3
#define _NBA_BEG 9
#define _NBA_LEN 3

#endif 

VectorXd x = VectorXd::Zero(_STT_LEN);
MatrixXd cov_x = MatrixXd::Identity(_STT_LEN, _STT_LEN) * 10;
bool _is_first_imu = true;
ros::Time _last_imu_stamp;

MatrixXd RPYtoR(double roll, double pitch, double yaw)
{
    Matrix3d R(_DIM_LEN, _DIM_LEN);
    R << 
         cos(pitch)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw), -cos(roll)*sin(yaw), cos(yaw)*sin(pitch) + cos(pitch)*sin(roll)*sin(yaw),
         cos(pitch)*sin(yaw) + cos(yaw)*sin(pitch)*sin(roll),  cos(roll)*cos(yaw), sin(pitch)*sin(yaw) - cos(pitch)*cos(yaw)*sin(roll),
                                       -cos(roll)*sin(pitch),           sin(roll),                                cos(pitch)*cos(roll);
    return R;
}

MatrixXd RPYtoR(const VectorXd &x)
{
    return RPYtoR(x(0), x(1), x(2));
}

Matrix3d RPYtoG(double roll, double pitch, double yaw)
{
    MatrixXd G(_DIM_LEN, _DIM_LEN);
    G <<  
         cos(pitch), 0, -cos(roll)*sin(pitch),
                  0, 1,             sin(roll),
         sin(pitch), 0,  cos(pitch)*cos(roll);
    return G;
}

Matrix3d RPYtoG(const VectorXd &x)
{
    return RPYtoG(x(0), x(1), x(2));
}

Matrix3d RPYtoInvG(const VectorXd &x)
{
    double roll = x(0), pitch = x(1);
    MatrixXd D(_DIM_LEN, _DIM_LEN);
    D << 
                               cos(pitch), 0,                        sin(pitch),
         (sin(pitch)*sin(roll))/cos(roll), 1, -(cos(pitch)*sin(roll))/cos(roll),
                    -sin(pitch)/cos(roll), 0,              cos(pitch)/cos(roll);
    return D;
}

double sqr(double x)
{
    return x * x;
}

MatrixXd RPYtoInvGdR(const VectorXd &x)
{   
    double roll = x(0), pitch = x(1);
    MatrixXd D(_DIM_LEN, _DIM_LEN);
    D << 
                                                   0, 0,                                     0,
                    -sin(pitch)/(sqr(sin(roll)) - 1), 0,            -cos(pitch)/sqr(cos(roll)),
         (sin(pitch)*sin(roll))/(sqr(sin(roll)) - 1), 0, (cos(pitch)*sin(roll))/sqr(cos(roll));
    return D;
}

MatrixXd RPYtoInvGdP(const VectorXd &x)
{
    double roll = x(0), pitch = x(1);
    MatrixXd D(_DIM_LEN, _DIM_LEN);
    D << 
                              -sin(pitch), 0,                       cos(pitch),
         (cos(pitch)*sin(roll))/cos(roll), 0, (sin(pitch)*sin(roll))/cos(roll),
                    -cos(pitch)/cos(roll), 0,            -sin(pitch)/cos(roll);
    return D;
}

MatrixXd RPYtoInvGdY(const VectorXd &x)
{
    return Matrix3d::Zero(_DIM_LEN, _DIM_LEN);
}

MatrixXd RPYtoRdR(const VectorXd &x)
{
    double roll = x(0), pitch = x(1), yaw = x(2);
    MatrixXd D(_DIM_LEN, _DIM_LEN);
    D << 
         -cos(roll)*sin(pitch)*sin(yaw),  sin(roll)*sin(yaw),  cos(pitch)*cos(roll)*sin(yaw),
          cos(roll)*cos(yaw)*sin(pitch), -cos(yaw)*sin(roll), -cos(pitch)*cos(roll)*cos(yaw),
                   sin(pitch)*sin(roll),           cos(roll),          -cos(pitch)*sin(roll);
    return D;
}

MatrixXd RPYtoRdP(const VectorXd &x)
{
    double roll = x(0), pitch = x(1), yaw = x(2);
    MatrixXd D(_DIM_LEN, _DIM_LEN);
    D << 
         - cos(yaw)*sin(pitch) - cos(pitch)*sin(roll)*sin(yaw), 0, cos(pitch)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw),
           cos(pitch)*cos(yaw)*sin(roll) - sin(pitch)*sin(yaw), 0, cos(pitch)*sin(yaw) + cos(yaw)*sin(pitch)*sin(roll),
                                         -cos(pitch)*cos(roll), 0,                               -cos(roll)*sin(pitch);
    return D;
}

MatrixXd RPYtoRdY(const VectorXd &x)
{
    double roll = x(0), pitch = x(1), yaw = x(2);
    MatrixXd D(_DIM_LEN, _DIM_LEN);
    D << 
         - cos(pitch)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll), -cos(roll)*cos(yaw), cos(pitch)*cos(yaw)*sin(roll) - sin(pitch)*sin(yaw),
           cos(pitch)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw), -cos(roll)*sin(yaw), cos(yaw)*sin(pitch) + cos(pitch)*sin(roll)*sin(yaw),
                                                             0,                   0,                                                   0;
    return D;
}

VectorXd RotMtoRPY(const MatrixXd &m)
{
    VectorXd rpy(3);
    double yaw = atan2(m(1,0), m(0,0));
    double pitch = atan2(-m(2,0), m(0,0)*cos(yaw) * m(1,0)*sin(yaw));
    double roll = atan2(m(0,2)*sin(yaw)-m(1,2)*cos(yaw),-m(0,1)*sin(yaw)+m(1,1)*cos(yaw));

    roll = asin(m(1, 2));
    yaw = atan2(-m(1, 0)/cos(roll), m(1, 1)/cos(roll));
    pitch = atan2(-m(0, 2)/cos(roll), m(2, 2)/cos(roll));

    rpy << roll, pitch, yaw;
    return rpy;
}

void pubOdometry()
{
    nav_msgs::Odometry odom;
    {
        odom.header.stamp = _last_imu_stamp;
        odom.header.frame_id = "map";

        odom.pose.pose.position.x = x(0);
        odom.pose.pose.position.y = x(1);
        odom.pose.pose.position.z = x(2);

        Quaterniond q;
        q = RPYtoR(x(3), x(4), x(5)).block<3,3>(0, 0);
        odom.pose.pose.orientation.x = q.x();
        odom.pose.pose.orientation.y = q.y();
        odom.pose.pose.orientation.z = q.z();
        odom.pose.pose.orientation.w = q.w();
    }

    //ROS_WARN("[update] publication done");
    ROS_WARN_STREAM("[final] b_g = " << x.segment(_B_G_BEG, _B_G_LEN).transpose());
    ROS_WARN_STREAM("[final] b_a = " << x.segment(_B_A_BEG, _B_A_LEN).transpose() << endl);
    ///ROS_WARN_STREAM("[final] cov_x = " << endl << cov_x << endl);
    odom_pub.publish(odom);
}

void imu_callback(const sensor_msgs::Imu::ConstPtr &msg)
{

    // remember the first stamp;
    if (_is_first_imu)
    {
        _is_first_imu = false;
        _last_imu_stamp = msg->header.stamp;
        return ;
    }
    double dt = (msg->header.stamp - _last_imu_stamp).toSec();
    _last_imu_stamp = msg->header.stamp;
    //ROS_WARN("[propagate] time stamp done");

    VectorXd w_m(_DIM_LEN), a_m(_DIM_LEN);
    MatrixXd Q_t = MatrixXd::Zero(_DIM_LEN << 2, _DIM_LEN << 2);
    { // get imu measurement
        w_m << 
            msg->angular_velocity.x, 
            msg->angular_velocity.y, 
            msg->angular_velocity.z;

        a_m << 
            msg->linear_acceleration.x, 
            msg->linear_acceleration.y, 
            msg->linear_acceleration.z;

        //ROS_WARN("[propagate] w_m, a_m done");
        
        Q_t(0, 0) = (5.0/180.0)*acos(-1)*(5.0/180.0)*acos(-1);
        Q_t(1, 1) = (5.0/180.0)*acos(-1)*(5.0/180.0)*acos(-1);
        Q_t(2, 2) = (5.0/180.0)*acos(-1)*(5.0/180.0)*acos(-1);
        Q_t(3, 3) = 1 * 1;
        Q_t(4, 4) = 1 * 1;
        Q_t(5, 5) = 1 * 1;
        Q_t(6, 6) = 0.5*0.5;
        Q_t(7, 7) = 0.5*0.5;
        Q_t(8, 8) = 0.5*0.5;
        Q_t(9, 9) = 0.5*0.5;
        Q_t(10, 10) = 0.5*0.5;
        Q_t(11, 11) = 0.5*0.5;
        //ROS_WARN("[propagate] Q_t done");
    }

    //ROS_WARN("[propagate] preparation done");

    { // update the mean value
        VectorXd x_dot = VectorXd::Zero(_STT_LEN);

        // update position
        x_dot.segment(_POS_BEG, _POS_LEN) = x.segment(_VEL_BEG, _VEL_LEN);

        // update orientation
        x_dot.segment(_ORT_BEG, _ORT_LEN) =  
            RPYtoInvG(x.segment(_ORT_BEG, _ORT_LEN)) * (w_m - x.segment(_B_G_BEG, _B_G_LEN));

        // update velocity
        VectorXd g = VectorXd::Zero(_DIM_LEN);
        g(2) = 9.8;
        x_dot.segment(_VEL_BEG, _VEL_LEN) = g + 
           RPYtoR(x.segment(_ORT_BEG, _ORT_LEN)) * (a_m - x.segment(_B_A_BEG, _B_A_LEN)); 

        x += x_dot * dt;
    }

    //ROS_WARN("[propagate] mean value propagation done");
    
    // update the covariance 
    {
        //> 1. compute the A_t, F_t matrix
        MatrixXd A_t = MatrixXd::Zero(_STT_LEN, _STT_LEN);

        // first block row 
        A_t.block(_POS_BEG, _VEL_BEG, _POS_LEN, _VEL_LEN) << Matrix3d::Identity();

        //ROS_WARN("[propagate] first block row done");
        
        // second block row
        auto rpy = x.segment(_ORT_BEG, _ORT_LEN);
        auto b_g = x.segment(_B_G_BEG, _B_G_LEN);
        A_t.block(_ORT_BEG, _ORT_BEG, _ORT_LEN, _ORT_LEN) << 
            RPYtoInvGdR(rpy) * (w_m - b_g), 
            RPYtoInvGdP(rpy) * (w_m - b_g),
            RPYtoInvGdY(rpy) * (w_m - b_g);  
        A_t.block(_ORT_BEG, _B_G_BEG, _ORT_LEN, _B_G_LEN) <<
            -RPYtoInvG(rpy);
        
        //ROS_WARN("[propagate] second block row done");

        // third block row
        auto b_a = x.segment(_B_A_BEG, _B_A_LEN);
        A_t.block(_VEL_BEG, _ORT_BEG, _VEL_LEN, _ORT_LEN) << 
            RPYtoRdR(rpy) * (a_m - b_a),
            RPYtoRdP(rpy) * (a_m - b_a),
            RPYtoRdY(rpy) * (a_m - b_a);
        A_t.block(_VEL_BEG, _B_A_BEG, _VEL_LEN, _B_A_LEN) << 
            -RPYtoR(rpy);
        
        //ROS_WARN("[propagate] A_t done");

        //> 2. compute the U_t, V_t matrix 
        MatrixXd U_t = MatrixXd::Zero(_STT_LEN, _DIM_LEN << 2);

        U_t.block(_ORT_BEG, _N_G_BEG, _ORT_LEN, _N_G_LEN) <<
            -RPYtoInvG(rpy);
        U_t.block(_VEL_BEG, _N_A_BEG, _VEL_LEN, _N_A_LEN) <<
            -RPYtoR(rpy);
        U_t.block(_B_G_BEG, _NBG_BEG, _B_G_LEN, _NBG_LEN) <<
            Matrix3d::Identity();
        U_t.block(_B_A_BEG, _NBA_BEG, _B_A_LEN, _NBA_LEN) <<
            Matrix3d::Identity();


        //ROS_WARN("[propagate] U_t done");
        MatrixXd F_t = A_t * dt, V_t = U_t * dt;
        for (int i = 0; i < _STT_LEN; ++i) F_t(i, i) += 1.0;

        //ROS_WARN_STREAM("A_t = endl" << endl << A_t);
        //ROS_WARN_STREAM("F_t = endl" << endl << F_t);
        //ROS_WARN_STREAM("U_t = endl" << endl << U_t);
        //ROS_WARN_STREAM("V_t = endl" << endl << V_t);

        //ROS_WARN_STREAM("cov_x = " << endl << cov_x);
        //ROS_WARN_STREAM("cov_x_1 = " << (F_t * cov_x * F_t.transpose()));
        //ROS_WARN_STREAM("cov_x_2 = " << (V_t * Q_t * V_t.transpose()));
        cov_x = F_t * cov_x * F_t.transpose() + V_t * Q_t * V_t.transpose();
    }
    pubOdometry();
    //ROS_WARN("[propagate] covariance propagation done");
}

//Rotation from the camera frame to the IMU frame

void odom_callback(const nav_msgs::Odometry::ConstPtr &msg)
{
    //your code for update
    //camera position in the IMU frame = (0, -0.05, +0.02)
    //camera orientaion in the IMU frame = Quaternion(0, 0, -1, 0); w x y z, respectively
    VectorXd z = VectorXd::Zero(_DIM_LEN << 1);
    MatrixXd R = MatrixXd::Zero(_DIM_LEN << 1, _DIM_LEN << 1);
    { // get the measurement from the msg
        auto q = msg->pose.pose.orientation;
        auto p = msg->pose.pose.position;
        Quaterniond q_ic(0, 0, -1, 0);
        Quaterniond q_cw(q.w, q.x, q.y, q.z);

        MatrixXd g_cw = MatrixXd::Zero(4, 4), g_ic = MatrixXd::Zero(4, 4);

        g_cw.block(0, 0, 3, 3) << q_cw.toRotationMatrix();
        g_cw.block(0, 3, 3, 1) << p.x, p.y, p.z;
        g_cw(3, 3) = 1.0;

        g_ic.block(0, 0, 3, 3) << q_ic.toRotationMatrix();
        g_ic.block(0, 3, 3, 1) << 0, -0.05, +0.02;
        g_ic(3, 3) = 1.0;

        MatrixXd g_wi = (g_ic * g_cw).inverse();
        //ROS_WARN_STREAM("g_ic = " << endl << g_ic);
        //ROS_WARN_STREAM("g_cw = " << endl << g_cw);
        //ROS_WARN_STREAM("g_wi = " << endl << g_wi);


        z << 
            g_wi.block(0, 3, 3, 1),
            RotMtoRPY(g_wi.block(0, 0, 3, 3));

        nav_msgs::Odometry odom = *msg;
        Quaterniond qq ;
        qq = RPYtoR(z.segment(3, 3)).block<3, 3>(0, 0);
        odom.pose.pose.position.x = z(0);
        odom.pose.pose.position.y = z(1);
        odom.pose.pose.position.z = z(2);
        odom.pose.pose.orientation.w = qq.w();
        odom.pose.pose.orientation.x = qq.x();
        odom.pose.pose.orientation.y = qq.y();
        odom.pose.pose.orientation.z = qq.z();
        msm_pub.publish(odom);

#if 1
        R(0, 0) = 0.02 * 0.02;
        R(1, 1) = 0.02 * 0.02;
        R(2, 2) = 0.02 * 0.02;
        R(3, 3) = (2.0/180.0)*acos(-1) * (2.0/180.0)*acos(-1);
        R(4, 4) = (2.0/180.0)*acos(-1) * (2.0/180.0)*acos(-1);
        R(5, 5) = (2.0/180.0)*acos(-1) * (2.0/180.0)*acos(-1);
#endif
    }

    //ROS_WARN("[update] preparation done");

    MatrixXd C_t = MatrixXd::Zero(_DIM_LEN << 1, _STT_LEN);
    MatrixXd K_t = MatrixXd::Zero(_STT_LEN, _DIM_LEN << 1);
    { // compute the C_t and K_t
        C_t.block(_P_M_BEG, _POS_BEG, _P_M_LEN, _POS_LEN) << 
            Matrix3d::Identity();
        C_t.block(_O_M_BEG, _ORT_BEG, _O_M_LEN, _ORT_LEN) <<
            Matrix3d::Identity();
        K_t = cov_x * C_t.transpose() * (C_t * cov_x * C_t.transpose() + R).inverse();
    }

    //ROS_WARN("[update] C_t K_t done");
#if 1
    {
        VectorXd error = z - C_t * x;
        for (int i = _O_M_BEG; i < _O_M_END; ++i)
        {
            if (error(i) < -_PI) 
            {
                //ROS_WARN_STREAM("error = " << error.transpose() << endl);
                error(i) += 2.0 * _PI;
                //ROS_WARN_STREAM("ERROR = " << error.transpose() << endl);
            }
            if (error(i) > _PI) 
            {
                //ROS_WARN_STREAM("error = " << error.transpose() << endl);
                error(i) -= 2.0 * _PI;
                //ROS_WARN_STREAM("ERROR = " << error.transpose() << endl);
            }

        }
        x = x + K_t * error;
        cov_x = cov_x - K_t * C_t * cov_x;
    }
#endif
    //ROS_WARN("[update] x cov_x done");

    //pubOdometry();
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "ekf");
    ros::NodeHandle n("~");

    ros::Subscriber s1 = n.subscribe("imu", 100, imu_callback);
    ros::Subscriber s2 = n.subscribe("tag_odom", 100, odom_callback);

    odom_pub = n.advertise<nav_msgs::Odometry>("ekf_odom", 100);
    msm_pub = n.advertise<nav_msgs::Odometry>("rcv_odom", 100);

    ros::spin();
}
