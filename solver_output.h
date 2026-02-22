#ifndef SOLVER_OUTPUT_H
#define SOLVER_OUTPUT_H

#ifdef __cplusplus
extern "C"{
#endif

/*
============================================================
初期化

・出力ディレクトリが存在するか確認
・無ければ mkdir
・solverのmainで1回だけ呼ぶ
============================================================
*/
int solver_output_init(const char* dir);

/*
============================================================
BIN書き出し（軽量版）

前提：
solver_output_init が呼ばれていること
============================================================
*/
void write_frame_bin(
    const char* dir,
    int step,
    int N,
    const double* x,
    const double* r,
    double length_factor
);

#ifdef __cplusplus
}
#endif

#endif
