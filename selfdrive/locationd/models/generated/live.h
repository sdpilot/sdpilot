#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_411448996248703735);
void live_err_fun(double *nom_x, double *delta_x, double *out_1085822966685210764);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_7262348582700839748);
void live_H_mod_fun(double *state, double *out_3888029891761610350);
void live_f_fun(double *state, double dt, double *out_6542603415890197301);
void live_F_fun(double *state, double dt, double *out_4350196572466825462);
void live_h_4(double *state, double *unused, double *out_7074488734674611985);
void live_H_4(double *state, double *unused, double *out_2717965093872703048);
void live_h_9(double *state, double *unused, double *out_8822594196648213134);
void live_H_9(double *state, double *unused, double *out_4569253841391744422);
void live_h_10(double *state, double *unused, double *out_2903305059900435345);
void live_H_10(double *state, double *unused, double *out_5349451036793780592);
void live_h_12(double *state, double *unused, double *out_2881046531679803007);
void live_H_12(double *state, double *unused, double *out_9099223470915436044);
void live_h_31(double *state, double *unused, double *out_2501729476602464410);
void live_H_31(double *state, double *unused, double *out_7694726252134761153);
void live_h_32(double *state, double *unused, double *out_8542714594002408647);
void live_H_32(double *state, double *unused, double *out_6666282133477867363);
void live_h_13(double *state, double *unused, double *out_4969815048066148841);
void live_H_13(double *state, double *unused, double *out_1152157836340279365);
void live_h_14(double *state, double *unused, double *out_8822594196648213134);
void live_H_14(double *state, double *unused, double *out_4569253841391744422);
void live_h_33(double *state, double *unused, double *out_4563887877865566276);
void live_H_33(double *state, double *unused, double *out_7601460816935932859);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}