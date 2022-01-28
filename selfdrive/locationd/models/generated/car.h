#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_2216711326364192709);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7483620347635243954);
void car_H_mod_fun(double *state, double *out_2162824484004701577);
void car_f_fun(double *state, double dt, double *out_1607279644594056566);
void car_F_fun(double *state, double dt, double *out_7079384484779116511);
void car_h_25(double *state, double *unused, double *out_1447412912867946261);
void car_H_25(double *state, double *unused, double *out_3370942617829710619);
void car_h_24(double *state, double *unused, double *out_796968104925971332);
void car_H_24(double *state, double *unused, double *out_3252620506830277588);
void car_h_30(double *state, double *unused, double *out_5084862441121968611);
void car_H_30(double *state, double *unused, double *out_6242956293119472661);
void car_h_26(double *state, double *unused, double *out_8020014396701552442);
void car_H_26(double *state, double *unused, double *out_4677178389784551764);
void car_h_27(double *state, double *unused, double *out_1589373138576386424);
void car_H_27(double *state, double *unused, double *out_7530538280956097973);
void car_h_29(double *state, double *unused, double *out_8360208364926737360);
void car_H_29(double *state, double *unused, double *out_3020725773613447789);
void car_h_28(double *state, double *unused, double *out_1950553220345849760);
void car_H_28(double *state, double *unused, double *out_8176291088642779162);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}