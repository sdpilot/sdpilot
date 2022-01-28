#include "car.h"

namespace {
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 5.991464547107981;

/******************************************************************************
 *                       Code generated with sympy 1.8                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_2216711326364192709) {
   out_2216711326364192709[0] = delta_x[0] + nom_x[0];
   out_2216711326364192709[1] = delta_x[1] + nom_x[1];
   out_2216711326364192709[2] = delta_x[2] + nom_x[2];
   out_2216711326364192709[3] = delta_x[3] + nom_x[3];
   out_2216711326364192709[4] = delta_x[4] + nom_x[4];
   out_2216711326364192709[5] = delta_x[5] + nom_x[5];
   out_2216711326364192709[6] = delta_x[6] + nom_x[6];
   out_2216711326364192709[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_7483620347635243954) {
   out_7483620347635243954[0] = -nom_x[0] + true_x[0];
   out_7483620347635243954[1] = -nom_x[1] + true_x[1];
   out_7483620347635243954[2] = -nom_x[2] + true_x[2];
   out_7483620347635243954[3] = -nom_x[3] + true_x[3];
   out_7483620347635243954[4] = -nom_x[4] + true_x[4];
   out_7483620347635243954[5] = -nom_x[5] + true_x[5];
   out_7483620347635243954[6] = -nom_x[6] + true_x[6];
   out_7483620347635243954[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_2162824484004701577) {
   out_2162824484004701577[0] = 1.0;
   out_2162824484004701577[1] = 0.0;
   out_2162824484004701577[2] = 0.0;
   out_2162824484004701577[3] = 0.0;
   out_2162824484004701577[4] = 0.0;
   out_2162824484004701577[5] = 0.0;
   out_2162824484004701577[6] = 0.0;
   out_2162824484004701577[7] = 0.0;
   out_2162824484004701577[8] = 0.0;
   out_2162824484004701577[9] = 1.0;
   out_2162824484004701577[10] = 0.0;
   out_2162824484004701577[11] = 0.0;
   out_2162824484004701577[12] = 0.0;
   out_2162824484004701577[13] = 0.0;
   out_2162824484004701577[14] = 0.0;
   out_2162824484004701577[15] = 0.0;
   out_2162824484004701577[16] = 0.0;
   out_2162824484004701577[17] = 0.0;
   out_2162824484004701577[18] = 1.0;
   out_2162824484004701577[19] = 0.0;
   out_2162824484004701577[20] = 0.0;
   out_2162824484004701577[21] = 0.0;
   out_2162824484004701577[22] = 0.0;
   out_2162824484004701577[23] = 0.0;
   out_2162824484004701577[24] = 0.0;
   out_2162824484004701577[25] = 0.0;
   out_2162824484004701577[26] = 0.0;
   out_2162824484004701577[27] = 1.0;
   out_2162824484004701577[28] = 0.0;
   out_2162824484004701577[29] = 0.0;
   out_2162824484004701577[30] = 0.0;
   out_2162824484004701577[31] = 0.0;
   out_2162824484004701577[32] = 0.0;
   out_2162824484004701577[33] = 0.0;
   out_2162824484004701577[34] = 0.0;
   out_2162824484004701577[35] = 0.0;
   out_2162824484004701577[36] = 1.0;
   out_2162824484004701577[37] = 0.0;
   out_2162824484004701577[38] = 0.0;
   out_2162824484004701577[39] = 0.0;
   out_2162824484004701577[40] = 0.0;
   out_2162824484004701577[41] = 0.0;
   out_2162824484004701577[42] = 0.0;
   out_2162824484004701577[43] = 0.0;
   out_2162824484004701577[44] = 0.0;
   out_2162824484004701577[45] = 1.0;
   out_2162824484004701577[46] = 0.0;
   out_2162824484004701577[47] = 0.0;
   out_2162824484004701577[48] = 0.0;
   out_2162824484004701577[49] = 0.0;
   out_2162824484004701577[50] = 0.0;
   out_2162824484004701577[51] = 0.0;
   out_2162824484004701577[52] = 0.0;
   out_2162824484004701577[53] = 0.0;
   out_2162824484004701577[54] = 1.0;
   out_2162824484004701577[55] = 0.0;
   out_2162824484004701577[56] = 0.0;
   out_2162824484004701577[57] = 0.0;
   out_2162824484004701577[58] = 0.0;
   out_2162824484004701577[59] = 0.0;
   out_2162824484004701577[60] = 0.0;
   out_2162824484004701577[61] = 0.0;
   out_2162824484004701577[62] = 0.0;
   out_2162824484004701577[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_1607279644594056566) {
   out_1607279644594056566[0] = state[0];
   out_1607279644594056566[1] = state[1];
   out_1607279644594056566[2] = state[2];
   out_1607279644594056566[3] = state[3];
   out_1607279644594056566[4] = state[4];
   out_1607279644594056566[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_1607279644594056566[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_1607279644594056566[7] = state[7];
}
void F_fun(double *state, double dt, double *out_7079384484779116511) {
   out_7079384484779116511[0] = 1;
   out_7079384484779116511[1] = 0;
   out_7079384484779116511[2] = 0;
   out_7079384484779116511[3] = 0;
   out_7079384484779116511[4] = 0;
   out_7079384484779116511[5] = 0;
   out_7079384484779116511[6] = 0;
   out_7079384484779116511[7] = 0;
   out_7079384484779116511[8] = 0;
   out_7079384484779116511[9] = 1;
   out_7079384484779116511[10] = 0;
   out_7079384484779116511[11] = 0;
   out_7079384484779116511[12] = 0;
   out_7079384484779116511[13] = 0;
   out_7079384484779116511[14] = 0;
   out_7079384484779116511[15] = 0;
   out_7079384484779116511[16] = 0;
   out_7079384484779116511[17] = 0;
   out_7079384484779116511[18] = 1;
   out_7079384484779116511[19] = 0;
   out_7079384484779116511[20] = 0;
   out_7079384484779116511[21] = 0;
   out_7079384484779116511[22] = 0;
   out_7079384484779116511[23] = 0;
   out_7079384484779116511[24] = 0;
   out_7079384484779116511[25] = 0;
   out_7079384484779116511[26] = 0;
   out_7079384484779116511[27] = 1;
   out_7079384484779116511[28] = 0;
   out_7079384484779116511[29] = 0;
   out_7079384484779116511[30] = 0;
   out_7079384484779116511[31] = 0;
   out_7079384484779116511[32] = 0;
   out_7079384484779116511[33] = 0;
   out_7079384484779116511[34] = 0;
   out_7079384484779116511[35] = 0;
   out_7079384484779116511[36] = 1;
   out_7079384484779116511[37] = 0;
   out_7079384484779116511[38] = 0;
   out_7079384484779116511[39] = 0;
   out_7079384484779116511[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_7079384484779116511[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_7079384484779116511[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7079384484779116511[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7079384484779116511[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_7079384484779116511[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_7079384484779116511[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_7079384484779116511[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_7079384484779116511[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_7079384484779116511[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_7079384484779116511[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7079384484779116511[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7079384484779116511[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_7079384484779116511[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_7079384484779116511[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_7079384484779116511[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7079384484779116511[56] = 0;
   out_7079384484779116511[57] = 0;
   out_7079384484779116511[58] = 0;
   out_7079384484779116511[59] = 0;
   out_7079384484779116511[60] = 0;
   out_7079384484779116511[61] = 0;
   out_7079384484779116511[62] = 0;
   out_7079384484779116511[63] = 1;
}
void h_25(double *state, double *unused, double *out_1447412912867946261) {
   out_1447412912867946261[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3370942617829710619) {
   out_3370942617829710619[0] = 0;
   out_3370942617829710619[1] = 0;
   out_3370942617829710619[2] = 0;
   out_3370942617829710619[3] = 0;
   out_3370942617829710619[4] = 0;
   out_3370942617829710619[5] = 0;
   out_3370942617829710619[6] = 1;
   out_3370942617829710619[7] = 0;
}
void h_24(double *state, double *unused, double *out_796968104925971332) {
   out_796968104925971332[0] = state[4];
   out_796968104925971332[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3252620506830277588) {
   out_3252620506830277588[0] = 0;
   out_3252620506830277588[1] = 0;
   out_3252620506830277588[2] = 0;
   out_3252620506830277588[3] = 0;
   out_3252620506830277588[4] = 1;
   out_3252620506830277588[5] = 0;
   out_3252620506830277588[6] = 0;
   out_3252620506830277588[7] = 0;
   out_3252620506830277588[8] = 0;
   out_3252620506830277588[9] = 0;
   out_3252620506830277588[10] = 0;
   out_3252620506830277588[11] = 0;
   out_3252620506830277588[12] = 0;
   out_3252620506830277588[13] = 1;
   out_3252620506830277588[14] = 0;
   out_3252620506830277588[15] = 0;
}
void h_30(double *state, double *unused, double *out_5084862441121968611) {
   out_5084862441121968611[0] = state[4];
}
void H_30(double *state, double *unused, double *out_6242956293119472661) {
   out_6242956293119472661[0] = 0;
   out_6242956293119472661[1] = 0;
   out_6242956293119472661[2] = 0;
   out_6242956293119472661[3] = 0;
   out_6242956293119472661[4] = 1;
   out_6242956293119472661[5] = 0;
   out_6242956293119472661[6] = 0;
   out_6242956293119472661[7] = 0;
}
void h_26(double *state, double *unused, double *out_8020014396701552442) {
   out_8020014396701552442[0] = state[7];
}
void H_26(double *state, double *unused, double *out_4677178389784551764) {
   out_4677178389784551764[0] = 0;
   out_4677178389784551764[1] = 0;
   out_4677178389784551764[2] = 0;
   out_4677178389784551764[3] = 0;
   out_4677178389784551764[4] = 0;
   out_4677178389784551764[5] = 0;
   out_4677178389784551764[6] = 0;
   out_4677178389784551764[7] = 1;
}
void h_27(double *state, double *unused, double *out_1589373138576386424) {
   out_1589373138576386424[0] = state[3];
}
void H_27(double *state, double *unused, double *out_7530538280956097973) {
   out_7530538280956097973[0] = 0;
   out_7530538280956097973[1] = 0;
   out_7530538280956097973[2] = 0;
   out_7530538280956097973[3] = 1;
   out_7530538280956097973[4] = 0;
   out_7530538280956097973[5] = 0;
   out_7530538280956097973[6] = 0;
   out_7530538280956097973[7] = 0;
}
void h_29(double *state, double *unused, double *out_8360208364926737360) {
   out_8360208364926737360[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3020725773613447789) {
   out_3020725773613447789[0] = 0;
   out_3020725773613447789[1] = 1;
   out_3020725773613447789[2] = 0;
   out_3020725773613447789[3] = 0;
   out_3020725773613447789[4] = 0;
   out_3020725773613447789[5] = 0;
   out_3020725773613447789[6] = 0;
   out_3020725773613447789[7] = 0;
}
void h_28(double *state, double *unused, double *out_1950553220345849760) {
   out_1950553220345849760[0] = state[5];
   out_1950553220345849760[1] = state[6];
}
void H_28(double *state, double *unused, double *out_8176291088642779162) {
   out_8176291088642779162[0] = 0;
   out_8176291088642779162[1] = 0;
   out_8176291088642779162[2] = 0;
   out_8176291088642779162[3] = 0;
   out_8176291088642779162[4] = 0;
   out_8176291088642779162[5] = 1;
   out_8176291088642779162[6] = 0;
   out_8176291088642779162[7] = 0;
   out_8176291088642779162[8] = 0;
   out_8176291088642779162[9] = 0;
   out_8176291088642779162[10] = 0;
   out_8176291088642779162[11] = 0;
   out_8176291088642779162[12] = 0;
   out_8176291088642779162[13] = 0;
   out_8176291088642779162[14] = 1;
   out_8176291088642779162[15] = 0;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_2216711326364192709) {
  err_fun(nom_x, delta_x, out_2216711326364192709);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7483620347635243954) {
  inv_err_fun(nom_x, true_x, out_7483620347635243954);
}
void car_H_mod_fun(double *state, double *out_2162824484004701577) {
  H_mod_fun(state, out_2162824484004701577);
}
void car_f_fun(double *state, double dt, double *out_1607279644594056566) {
  f_fun(state,  dt, out_1607279644594056566);
}
void car_F_fun(double *state, double dt, double *out_7079384484779116511) {
  F_fun(state,  dt, out_7079384484779116511);
}
void car_h_25(double *state, double *unused, double *out_1447412912867946261) {
  h_25(state, unused, out_1447412912867946261);
}
void car_H_25(double *state, double *unused, double *out_3370942617829710619) {
  H_25(state, unused, out_3370942617829710619);
}
void car_h_24(double *state, double *unused, double *out_796968104925971332) {
  h_24(state, unused, out_796968104925971332);
}
void car_H_24(double *state, double *unused, double *out_3252620506830277588) {
  H_24(state, unused, out_3252620506830277588);
}
void car_h_30(double *state, double *unused, double *out_5084862441121968611) {
  h_30(state, unused, out_5084862441121968611);
}
void car_H_30(double *state, double *unused, double *out_6242956293119472661) {
  H_30(state, unused, out_6242956293119472661);
}
void car_h_26(double *state, double *unused, double *out_8020014396701552442) {
  h_26(state, unused, out_8020014396701552442);
}
void car_H_26(double *state, double *unused, double *out_4677178389784551764) {
  H_26(state, unused, out_4677178389784551764);
}
void car_h_27(double *state, double *unused, double *out_1589373138576386424) {
  h_27(state, unused, out_1589373138576386424);
}
void car_H_27(double *state, double *unused, double *out_7530538280956097973) {
  H_27(state, unused, out_7530538280956097973);
}
void car_h_29(double *state, double *unused, double *out_8360208364926737360) {
  h_29(state, unused, out_8360208364926737360);
}
void car_H_29(double *state, double *unused, double *out_3020725773613447789) {
  H_29(state, unused, out_3020725773613447789);
}
void car_h_28(double *state, double *unused, double *out_1950553220345849760) {
  h_28(state, unused, out_1950553220345849760);
}
void car_H_28(double *state, double *unused, double *out_8176291088642779162) {
  H_28(state, unused, out_8176291088642779162);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
