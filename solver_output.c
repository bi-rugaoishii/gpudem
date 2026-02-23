#include "solver_output.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

/*
============================================================
ディレクトリ作成（初期化時のみ）
============================================================
*/
int solver_output_init(const char* dir)
{
    struct stat st;

    /* 既に存在するか確認 */
    if(stat(dir,&st)==0)
    {
        if(S_ISDIR(st.st_mode))
        {
            printf("Output dir exists: %s\n",dir);
            return 0;
        }
    }

    /* 無ければ作成 */
    if(mkdir(dir,0755)==0)
    {
        printf("Created output dir: %s\n",dir);
        return 0;
    }

    if(errno==EEXIST)
        return 0;

    printf("ERROR: cannot create dir %s\n",dir);
    return -1;
}

/*
============================================================
ファイル名生成
============================================================
*/
static void build_filename(
    char* out,
    const char* dir,
    int step)
{
    sprintf(out,"%s/result_%04d.bin",dir,step);
}

/*
============================================================
BIN出力

フォーマット:
[int N]
[float pos[N*3]]
[float radius[N]]
============================================================
*/
void write_frame_bin(
    const char* dir,
    int step,
    int N,
    const double* x,
    const double* r,const double length_factor)
{
    char filename[256];

    build_filename(filename,dir,step);

    FILE* fp=fopen(filename,"wb");

    if(fp==NULL)
    {
        printf("ERROR: cannot open %s\n",filename);
        return;
    }

    /* 粒子数 */
    fwrite(&N,sizeof(int),1,fp);

    /* float変換バッファ */
    float* pos_f=(float*)malloc(sizeof(float)*N*3);
    float* r_f  =(float*)malloc(sizeof(float)*N);

    int i;
    for(i=0;i<N;i++)
    {
        pos_f[i*3+0]=(float)(x[i*3+0]*length_factor);
        pos_f[i*3+1]=(float)(x[i*3+1]*length_factor);
        pos_f[i*3+2]=(float)(x[i*3+2]*length_factor);

        r_f[i]=(float)(r[i]*length_factor);
    }

    fwrite(pos_f,sizeof(float),N*3,fp);
    fwrite(r_f,sizeof(float),N,fp);

    free(pos_f);
    free(r_f);

    fclose(fp);
}
