/*
 * @Author: dongtaihong 2396203400@qq.com
 * @Date: 2024-06-16 23:12:51
 * @LastEditTime: 2024-06-24 20:48:57
 * @FilePath: /AlgExper/common/cartesian_frenet_conversion.cc
 * @Description:
 */
#include "cartesian_frenet_conversion.h"

#include <iostream>

void FrenetToCartesian(std::vector<TrajectoryPoint>& frenet_path,
                       std::vector<TrajectoryPoint>& ref_cartesian_path) {
  TrajectoryPoint proj_site;
  for (unsigned int iter_frenet = 0; iter_frenet < frenet_path.size();++iter_frenet) 
  {
    //根据s匹配最近的参考点
    proj_site = FindFrenetProjPoint(ref_cartesian_path,
                                    frenet_path[iter_frenet].frenet_info.s);

    FrenetToCartesian(proj_site.length, proj_site.xg, proj_site.yg,
                      proj_site.global_angle,
                      frenet_path[iter_frenet].frenet_info.s,
                      frenet_path[iter_frenet].frenet_info.l,
                      frenet_path[iter_frenet].xg, frenet_path[iter_frenet].yg);
  }
  std::cout << "proj_site xg: " << proj_site.xg << " yg: " << proj_site.yg
            << " s: " << frenet_path.back().frenet_info.s << " s_ref: "
            << ref_cartesian_path[ref_cartesian_path.size() - 2].frenet_info.s
            << " s_ref2: " << ref_cartesian_path.back().frenet_info.s
            << std::endl;
}

//不完整，仅sl，配合ComputePathProfile可以得到完整的笛卡尔坐标
void FrenetToCartesian(const double& rs, const double& rx, const double& ry,
                       const double& rtheta, const double& s_condition,
                       const double& d_condition, double& x, double& y) {
  const double cos_theta_r = std::cos(rtheta);
  const double sin_theta_r = std::sin(rtheta);

  x = rx - sin_theta_r * d_condition;
  y = ry + cos_theta_r * d_condition;
}
//输入：一条带有 Frenet s 值的参考轨迹 + 某个目标 s
//输出：该 s 在轨迹上的对应点的全局坐标和朝向（插值获得）
TrajectoryPoint FindFrenetProjPoint(
    const std::vector<TrajectoryPoint>& ref_cartesian_path, const double& s) {
  double min_dis = 9999.0;
  TrajectoryPoint proj_site;//初始化一个投影点，用来保存最终结果
  //逐对轨迹点进行检查，寻找目标 s 落在哪两点之间
  for (unsigned int index = 0; index < ref_cartesian_path.size() - 1; ++index) {
    //如果 s 已经超过轨迹末尾，直接返回末尾点当作投影点
    if (s > ref_cartesian_path.back().frenet_info.s) {
      proj_site = std::move(ref_cartesian_path.back());
    }
    //确定 s 介于 index 与 index+1 间，接下来开始算精确位置
    if (s >= ref_cartesian_path[index].frenet_info.s &&
        s < ref_cartesian_path[index + 1].frenet_info.s) {
      proj_site.frenet_info.s = s;
      //proportion 是在区间中的相对长度权重，范围 0 到 1： 0 表示完全靠前一个点 1 表示完全靠后一个点
      double proportion = (s - ref_cartesian_path[index].frenet_info.s) /
                          (ref_cartesian_path[index + 1].frenet_info.s -
                           ref_cartesian_path[index].frenet_info.s);
      //线性插值计算全局坐标
      proj_site.xg = (1 - proportion) * ref_cartesian_path[index].xg +
                     proportion * ref_cartesian_path[index + 1].xg;
      proj_site.yg = (1 - proportion) * ref_cartesian_path[index].yg +
                     proportion * ref_cartesian_path[index + 1].yg;
      //线性插值计算朝向角，注意处理角度跳变问题，如果差超过 2π，直接取后一个角，当成未跳变处理【可以进一步优化】
      if (std::fabs(ref_cartesian_path[index].global_angle -
                    ref_cartesian_path[index + 1].global_angle) >= 6.28) {
        proj_site.global_angle = ref_cartesian_path[index + 1].global_angle;
      } else {
        proj_site.global_angle =
            (1 - proportion) * ref_cartesian_path[index].global_angle +
            proportion * ref_cartesian_path[index + 1].global_angle;
      }
      break;
    }
  }
  return proj_site;//投影点已经获得 s + xg + yg + global_angle 等信息
}

void CartesianToFrenet(std::vector<TrajectoryPoint>& cartesian_path,
                       std::vector<TrajectoryPoint>& ref_cartesian_path) {
  for (auto& point : cartesian_path) {
    CartesianToFrenet(point, ref_cartesian_path);
  }
}

void CartesianToFrenet(TrajectoryPoint& cardesian_point,
                       std::vector<TrajectoryPoint>& ref_cartesian_path) {
  int match_index = getMatchPoint(cardesian_point, ref_cartesian_path);
  TrajectoryPoint proj_point;
  getProjectPoint(cardesian_point, ref_cartesian_path[match_index], proj_point);
  CartesianToFrenet(cardesian_point, proj_point);
}

void CartesianToFrenet(TrajectoryPoint& cardesian_point,
                       const TrajectoryPoint& project_point) {
  Eigen::Vector2d r_h = {cardesian_point.xg, cardesian_point.yg};
  Eigen::Vector2d r_r = {project_point.xg, project_point.yg};
  Eigen::Vector2d tao_r = {cos(project_point.global_angle),
                           sin(project_point.global_angle)};
  Eigen::Vector2d n_r = {-sin(project_point.global_angle),
                         cos(project_point.global_angle)};
  Eigen::Vector2d n_h = {-sin(cardesian_point.global_angle),
                         cos(cardesian_point.global_angle)};

  cardesian_point.frenet_info.s = project_point.frenet_info.s;
  cardesian_point.frenet_info.l = (r_h - r_r).dot(n_r);

  // TODO: 这里的vx,vy是全局坐标系下的，还是车辆坐标系下的？
  Eigen::Vector2d v = {cardesian_point.v * cos(cardesian_point.global_angle),
                       cardesian_point.v * sin(cardesian_point.global_angle)};
  cardesian_point.frenet_info.l_d = v.dot(n_r);
  std::cout << "cardesian2Frenet: v" << cardesian_point.v << "----"
            << cardesian_point.frenet_info.l_d << std::endl;
  cardesian_point.frenet_info.s_d =
      (1 / (1 - cardesian_point.frenet_info.l * project_point.curvature)) *
      v.dot(tao_r);
  if (fabs(cardesian_point.frenet_info.s_d) < 1e-8) {
    cardesian_point.frenet_info.l_ds = 0;
  } else {
    cardesian_point.frenet_info.l_ds =
        cardesian_point.frenet_info.l_d / cardesian_point.frenet_info.s_d;
  }

  Eigen::Vector2d a = {cardesian_point.a * cos(cardesian_point.global_angle),
                       cardesian_point.a * sin(cardesian_point.global_angle)};

  cardesian_point.frenet_info.l_dd =
      a.dot(n_r) -
      project_point.curvature *
          (1 - project_point.curvature * cardesian_point.frenet_info.l) *
          pow(cardesian_point.frenet_info.s_d, 2);
  // because d_k/ds is small, use 0 replace
  cardesian_point.frenet_info.s_dd =
      (1 / (1 - cardesian_point.frenet_info.l * project_point.curvature)) *
      (a.dot(tao_r) + 2 * project_point.curvature *
                          cardesian_point.frenet_info.l_ds *
                          pow(cardesian_point.frenet_info.s_d, 2));
  if (fabs(cardesian_point.frenet_info.s_d) < 1e-8) {
    cardesian_point.frenet_info.l_dds = 0;
  } else {
    cardesian_point.frenet_info.l_dds =
        (1 / pow(cardesian_point.frenet_info.s_d, 2)) *
        (cardesian_point.frenet_info.l_dd -
         cardesian_point.frenet_info.l_ds * cardesian_point.frenet_info.s_dd);
  }
}

int getMatchPoint(TrajectoryPoint hostPoint,
                  const std::vector<TrajectoryPoint>& vecTraj) {
  std::vector<double> distVec;
  for (int i = 0; i < vecTraj.size() - 1; ++i) {
    double tempDist = pow(hostPoint.xg - vecTraj.at(i).xg, 2) +
                      pow(hostPoint.yg - vecTraj.at(i).yg, 2);
    distVec.push_back(tempDist);
  }
  auto itr = std::min_element(distVec.begin(), distVec.end());
  return static_cast<int>(std::distance(distVec.begin(), itr));
}

void getProjectPoint(TrajectoryPoint hostPoint, TrajectoryPoint& matchPoint,
                     TrajectoryPoint& projectPoint) {
  Eigen::Vector2d hostPos = {hostPoint.xg, hostPoint.yg},
                  matchPointPos = {matchPoint.xg, matchPoint.yg};
  Eigen::Vector2d d = hostPos - matchPointPos;
  Eigen::Vector2d tao = {cos(matchPoint.global_angle),
                         sin(matchPoint.global_angle)};
  Eigen::Vector2d projectPose = matchPointPos + tao * (d.transpose() * tao);
  projectPoint.curvature = matchPoint.curvature;
  projectPoint.global_angle =
      matchPoint.global_angle +
      matchPoint.curvature * static_cast<double>((d.transpose() * tao));
  projectPoint.xg = projectPose[0];
  projectPoint.yg = projectPose[1];

  // s有方向的概念
  Eigen::Vector2d dir = {projectPoint.xg - matchPoint.xg,
                         projectPoint.yg - matchPoint.yg};
  double s_temp = std::sqrt(pow(projectPoint.xg - matchPoint.xg, 2) +
                            pow(projectPoint.yg - matchPoint.yg, 2));
  projectPoint.frenet_info.s = dir.dot(tao) > 0
                                   ? matchPoint.frenet_info.s + s_temp
                                   : matchPoint.frenet_info.s - s_temp;
}
//将参考轨迹 Frenet s 坐标系以车辆当前位置为原点重新校准，让路径规划可以基于相对 s 而不是全局 s
void updateRefLineS(TrajectoryPoint& start_point,
                    std::vector<TrajectoryPoint>& ref_path) {
  int match_index = getMatchPoint(start_point, ref_path);
  TrajectoryPoint proj_point;
  getProjectPoint(start_point, ref_path[match_index], proj_point);//把车辆在轨迹上的位置精确定位到参考路径上的“最近点”，得到对应的 Frenet s
  for (auto& p : ref_path) {
    p.frenet_info.s -= proj_point.frenet_info.s;//把车辆当前位置作为 Frenet s 的 0 点
  }
}