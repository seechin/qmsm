#include    <stdio.h>
#include    <stdlib.h>
#include    <math.h>
#include    <stdint.h>
#include    <string.h>
#include    <stdlib.h>
#include    <sys/types.h>
#include    <sys/wait.h>
#include    <sys/stat.h>
#include    <sys/mman.h>
#include    <sys/resource.h>
#include    <signal.h>
#include    <fcntl.h>
#include    <ctype.h>
#include    <sys/time.h>
#include    <pthread.h>

#include    "StringX.cpp"
#include    "HeapSort.cpp"
#include    "matrics_dyson_kernel.cpp"

const char * software_name = "dyson-kernel";
const char * software_version = "(c) Cao Siqin, 2018.6.15";
const char * szHelp = "\
    [-f]                TPM file to handle\n\
    -fs, -list          TPM file list to handle\n\
    -dt 1               time step\n\
    -silent             don't display iteration details\n\
    -show-...           K (default), iK, miK, dT. Will turn off -predict\n\
        the parameter list following -show-miK is the stationary populations\n\
    -%9.3g              output format\n\
    -predict[-to]       time of prediction. No less than training data\n\
    -skip               time of skip at the beginning.\n\
    -predict-*[-to]     -predict with: eigen_value(*=ev), ITS(*=its)\n\
    -debug[-level 0]    0/1/2/3, show debug information\n\
";
const char * default_real_output_format = "%9.3g";
//-[no-]precise-dT0   use (T(2)^0.0001-I)/0.0001 as dT(0), default on

#define MAX_FILE 1024
#define MAX_WORD 20000
#define MAX_BUFFER (MAX_WORD*24)
#define MACHINE_REASONABLE_ERROR 1e-12
#define OVERFLOWSIZE 0

extern "C" void dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);
bool eigen_decompose(int n, int nv, double * matrix, double * evl_real, double * evl_image, double * ev_left, double * ev_right, double * work, int nwork){
    int ldvl = n; int ldvr = n; int ret = 0;
    dgeev_("V","V", &n, matrix, &nv, evl_real, evl_image, ev_left, &ldvl, ev_right, &ldvr, work, &nwork, &ret);
    return ret==0;
}
bool eigen_value(int n, int nv, double * matrix, double * evl_real, double * evl_image, double * work, int nwork){
    int ldvl = n; int ldvr = n; int ret = 0;
    dgeev_("N","N", &n, matrix, &nv, evl_real, evl_image, NULL, &ldvl, NULL, &ldvr, work, &nwork, &ret);
    return ret==0;
}

typedef bool (FunctionOfMatrix) (int,double*,double param);
bool matrix_pow(int dim, double * data, double param){
    bool success = true;
    for (int i=0; i<dim; i++){
        if (data[i*dim+i]<0) success = false;
        data[i*dim+i] = pow(data[i*dim+i], param);
    }
    return success;
}
void print_matrix_data(double * m, int dim, const char * sep, const char * format, bool multiline_display){
    if (!format||!format[0]) format = default_real_output_format;
    for (int a=0; a<dim; a++){
        for (int b=0; b<dim; b++){
            printf(format, fabs(m[a*dim+b])<1e-99?0:m[a*dim+b]);
            printf(" ");
        }; printf("%s", sep); if (multiline_display) printf("\n");
    }
    printf("\n");
}
void print_matrix(MatrixNS::matrix * m, const char * sep, const char * format, bool multiline_display){
    if (!format||!format[0]) format = default_real_output_format;
    for (int a=0; a<m->n; a++){
        for (int b=0; b<m->n; b++){
            printf(format, fabs(m->m[a*m->n+b])<1e-99?0:m->m[a*m->n+b]);
            printf(" ");
        }; printf("%s", sep); if (multiline_display) printf("\n");
    }
    printf("\n");
}
void print_memory_allocation(long int __global_memory_allocated, char * out, int size){
    if (__global_memory_allocated>1000000000) snprintf(out, size, "%.2f GB", ((double)__global_memory_allocated)/1024/1024/1024);
    else if (__global_memory_allocated>1000000) snprintf(out, size, "%.2f MB", ((double)__global_memory_allocated)/1024/1024);
    else if (__global_memory_allocated>1000) snprintf(out, size, "%.2f KB", ((double)__global_memory_allocated)/1024);
    else snprintf(out, size, "%ld B", __global_memory_allocated);
}
void print_matrix_array(double ** out_Matrix, int begin, int end, int dimension, MatrixNS::matrix mm[4], MatrixNS::MatrixDataContainer * mc, double * eigen_work, double dt, double t_skip, int output_type, const char * output_format, bool multiline_display){
    for (int i=begin; i<end; i++){
        if (output_type==1 || output_type==2){ // ev or ITS
            MatrixNS::matrix M; M.harness(out_Matrix[i], dimension, mc);
            mm[0] = M.m;
            //eigen_decompose(M.n, M.n, mm[0].m, mm[1].m, mm[2].m, mm[3].m, mm[4].m, eigen_work, 4*M.n);
            eigen_value(M.n, M.n, mm[0].m, mm[1].m, mm[2].m, eigen_work, 4*M.n);
            mm[3] = mm[1].m; mm[3] *= -1;
            HeapSortNS::HeapSortHeap <double, double> (M.n, mm[3].m, mm[1].m).HeapSort();
            for (int j=output_type==1?0:1; j<M.n; j++){
                double value = output_type==1? mm[1].m[j] : -(i*dt+t_skip+dt) / log(mm[1].m[j]);
                printf(output_format, fabs(value)<1e-99?0:value);
                printf(j+1<M.n? " " : "\n");
            };
        } else {
            print_matrix_data(&out_Matrix[i][0], dimension, "", output_format, multiline_display);
        }
    }
}
/*bool perform_matrix_function(MatrixNS::matrix & A, MatrixNS::matrix & ret, MatrixNS::matrix * temp[6], double * eigen_work, FunctionOfMatrix * matrix_function, double param=0){ // perform A * f(A)
    bool success = true;
    memset(eigen_work, 0, sizeof(double) * 4 * A.n);
    MatrixNS::matrix * ev = temp[0];
    MatrixNS::matrix * evi = temp[1];
    MatrixNS::matrix * vlt = temp[2];
    MatrixNS::matrix * vr = temp[3];
    ret = A;
  // eigen decomposition first
    eigen_decompose(A.n, A.n, ret.m, ev->m, evi->m, vr->m, vlt->m, eigen_work, 4*A.n);
    vr->transpose();
  // convert eigenvalue array into eigenvalue matrix
    for (int i=1; i<A.n; i++){ int k = i*A.n+i;
        ev->m[k] = ev->m[i]; ev->m[i] = 0;
        evi->m[k] = ev->m[i]; evi->m[i] = 0;
    }

//    printf("A:\n"); print_matrix_data(A.m, A.n, "\n");
//    printf("ev:\n"); print_matrix_data(ev->m, A.n, "\n");
//    printf("vlt:\n"); print_matrix_data(vlt->m, A.n, "\n");
//    printf("vr:\n"); print_matrix_data(vr->m, A.n, "\n");

  // get normalized eigenvectors
    *temp[4] = 0.0; for (int i=0; i<A.n; i++){
        double rescale = 0; for (int j=0; j<A.n; j++) rescale += vlt->m[i*A.n+j] * vr->m[j*A.n+i];
        temp[4]->m[i*A.n+i] = sqrt(fabs(1/rescale));
    }
    * vlt = *temp[4] * *vlt; * vr = *vr * *temp[4];

  // do function
    success &= matrix_function(A.n, ev->m, param);

    //ret = *vlt * A; ret = ret * *vr; // test L^T*A*V = lambda
    //ret = *vr * *ev; ret = ret * *vlt; // test V*lambda*L^T = A

  // reconstructing the result matrix
    ret = *vr * *ev; ret = ret * *vlt;

    //printf("vlt:\n"); print_matrix_data(vlt->m, A.n, "\n");
    //printf("vr:\n"); print_matrix_data(vr->m, A.n, "\n");
    //printf("ret:\n"); print_matrix_data(ret.m, A.n, "\n");

    return success;
}*/

long int __global_memory_allocated = 0; int debug_level = 0;
void print_current_memory_allocation(long int len, bool already_allocated = true){
    char memory_out_buffer[2][128];
    print_memory_allocation(len, memory_out_buffer[0], sizeof(memory_out_buffer[0]));
    print_memory_allocation(__global_memory_allocated, memory_out_buffer[1], sizeof(memory_out_buffer[1]));
    if (already_allocated){
        fprintf(stderr, "%s : debug : %s allocated, totally %s\n", software_name, memory_out_buffer[0], memory_out_buffer[1]);
    } else {
        fprintf(stderr, "%s : debug : allocating %s, totally %s\n", software_name, memory_out_buffer[0], memory_out_buffer[1]);
    }
}

double ** init_matrix(int m, int n, int overflow_chars=0){
    long int lenh = sizeof(double*) * m; long int len = lenh + sizeof(double) * m * n + overflow_chars;
    if (debug_level>=3) print_current_memory_allocation(len, false);
    char * p = (char*) malloc(len); if (!p) return NULL; __global_memory_allocated += len;
    if (debug_level>=2) print_current_memory_allocation(len);
    memset(p, 0, len); double * d = (double*)(p + lenh);
    double ** a = (double**) p; for (int i=0; i<m; i++) a[i] = &d[i*n];
    return a;
}
void cp_double_array(double * dst, double * src, int n){
    for (int i=0; i<n; i++) dst[i] = src[i];
}
void set_matrix(MatrixNS::matrix * m, double * data, int n){
    cp_double_array(&m->m[0], data, n);
}
void get_matrix(MatrixNS::matrix * m, double * out, int n){
    cp_double_array(out, &m->m[0], n);
}
void set_temp_matrix(MatrixNS::matrix * m, double * array, int dimension){
    m->n = dimension; m->m = array;
}

int main(int argc, char * argv[]){
    bool success = true; bool silent = false; double dt = 1; double dt2 = 1; double t_skip = 0; int n_skip = 0;
    double t_predict = -1; int n_predict = 0; int predict_output = 0; // output type: 0-TPM, 1-ITS
    int show_K = 0;
    char * filename = NULL; bool file_is_list = false;
    const char * output_format = default_real_output_format;
    double * stationary_population = NULL; int n_stationary_population = 0;

    if (argc<2){ printf("%s %s\n%s", software_name, software_version, szHelp); return 0; }
    for (int i=1; i<argc; i++){
        if (argv[i][0] == '-'){
            StringNS::string key = argv[i];
            if (key=="-h" || key=="--h" || key=="-help" || key=="--help"){
                printf("%s %s\n%s", software_name, software_version, szHelp); success = false;
            } else if (key=="-f" || key=="--f"){
                if (i+1<argc && argv[i+1][0]!='-'){ i++; filename = argv[i]; file_is_list = false; }
            } else if (key=="-fs" || key=="--fs" || key=="-list" || key=="--list"){
                if (i+1<argc && argv[i+1][0]!='-'){ i++; filename = argv[i]; file_is_list = true; }
            } else if (key=="-dt" || key=="--dt"){
                if (i+1<argc && StringNS::is_string_number(argv[i+1])){
                    i++; dt = atof(argv[i]); if (dt<MACHINE_REASONABLE_ERROR) dt = MACHINE_REASONABLE_ERROR;
                    dt2 = dt * dt;
                }
            } else if (key=="-silent" || key=="--silent"){
                silent = true;
            } else if (key=="-debug" || key=="--debug" || key=="-debug-level" || key=="--debug-level" || key=="-debug_level" || key=="--debug_level"){
                debug_level = 1;
                if (i+1<argc && StringNS::is_string_number(argv[i+1])){
                    i++; debug_level = atoi(argv[i]);
                }
            } else if (key=="-skip" || key=="--skip" || key=="-skip-time" || key=="--skip-time" || key=="-skip_time" || key=="--skip_time"){
                if (i+1<argc && StringNS::is_string_number(argv[i+1])){
                    i++; t_skip = atof(argv[i]);
                }
            } else if (key=="-predict" || key=="--predict" || key=="-predict-to" || key=="--predict-to" || key=="-predict_to" || key=="--predict_to"){
                if (i+1<argc && StringNS::is_string_number(argv[i+1])){
                    i++; t_predict = atof(argv[i]);
                } else t_predict = MACHINE_REASONABLE_ERROR;
                predict_output = 0;
            } else if (key=="-predict-ev" || key=="--predict-ev" || key=="-predict_ev" || key=="--predict_ev" || key=="-predict-ev-to" || key=="--predict-ev-to" || key=="-predict_ev_to" || key=="--predict_ev_to"){
                if (i+1<argc && StringNS::is_string_number(argv[i+1])){
                    i++; t_predict = atof(argv[i]);
                } else t_predict = MACHINE_REASONABLE_ERROR;
                predict_output = 1;
            } else if (key=="-predict-its" || key=="--predict-its" || key=="-predict_its" || key=="--predict_its" || key=="-predict-its-to" || key=="--predict-its-to" || key=="-predict_its_to" || key=="--predict_its_to"){
                if (i+1<argc && StringNS::is_string_number(argv[i+1])){
                    i++; t_predict = atof(argv[i]);
                } else t_predict = MACHINE_REASONABLE_ERROR;
                predict_output = 2;
            } else if (key=="-show-dT" || key=="--show-dT" || key=="-show_dT" || key=="--show_dT"){
                show_K = 1;
            } else if (key=="-show-K" || key=="--show-K" || key=="-show_K" || key=="--show_K"){
                show_K = 0;
            } else if (key=="-show-iK" || key=="--show-iK" || key=="-show_iK" || key=="--show_iK"){
                show_K = 2;
            } else if (key=="-show-miK" || key=="--show-miK" || key=="-show_miK" || key=="--show_miK" || key=="-show-m0iK" || key=="--show-m0iK" || key=="-show_m0iK" || key=="--show_m0iK"){
                show_K = 3;
            } else if (key=="-show-m1iK" || key=="--show-m1iK" || key=="-show_m1iK" || key=="--show_m1iK" || key=="-show-m2iK" || key=="--show-m2iK" || key=="-show_m2iK" || key=="--show_m2iK"){
                if (key=="-show-m1iK" || key=="--show-m1iK" || key=="-show_m1iK" || key=="--show_m1iK"){
                    show_K = 4;
                } else {
                    show_K = 5;
                }
                int isp_start = i+1; while (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; }
                int isp_end = i+1;
                if (isp_end > isp_start){
                    n_stationary_population = isp_end - isp_start;
                    //printf("sps: %d, %d ~ %d\n", n_stationary_population, isp_start, isp_end);
                    stationary_population = (double*)malloc(sizeof(double)*n_stationary_population);
                    for (int k=0; k<n_stationary_population; k++) stationary_population[k] = atof(argv[isp_start+k]);
                }
                //for (int k=0; k<n_stationary_population; k++) printf("sp[%d] = %g\n", k+1, stationary_population[k]);
            } else if (key.text[0]=='-' && key.text[1]=='%'){
                output_format = &key.text[1];
            } else if (key.text[0]=='-' && key.text[1]=='-' && key.text[2]=='%'){
                output_format = &key.text[2];
            } else {
                fprintf(stderr, "%s : error : unrecognizable parameter %s\n", software_name, argv[i]); success = false;
            }
        } else {
            filename = argv[i]; file_is_list = false;
        }
    }
    if (!success) return 0;
    if (!filename){ fprintf(stderr, "%s : please specify a file to handle.\n", software_name); success = false; }
    FILE * file = NULL;
    if (success && filename){
       file = fopen(filename, "r"); if (!file){ fprintf(stderr, "%s : error : cannot open %s\n", software_name, filename); success = false; }
    }

    n_skip = (int)(t_skip/dt);

  // get dimensions: the matrix dimension (square root of column number), and the number of matrixes (line number)
    int line_num = 0; int col_num = 0; int dimension = 0;
    if (success && file){ int inner_line = 0;
if (debug_level>=1) fprintf(stderr, "%s : debug : preprocessing of %s\n", software_name, filename);
        if (file_is_list){
            char fninput[MAX_FILE];
            fseek(file, 0, SEEK_SET); while (fgets(fninput, sizeof(fninput), file)) {
                StringNS::string sl[MAX_WORD];
                int nw = StringNS::analysis_line(fninput, sl, MAX_WORD, true);
                int ni = 0; if (sl[0].text[0]=='#') {
                    if (nw>=2 && sl[0]=="#dimension"){
                        dimension = atof(sl[1].text); col_num = dimension*dimension;
                    } else if (nw>=3 && sl[0]=="#" && sl[1]=="dimension"){
                        dimension = atof(sl[2].text); col_num = dimension*dimension;
                    }
                    continue;
                } else for (int i=0; i<nw; i++){
                    if (sl[i].text[0]=='#' || sl[i].text[0]==';' || sl[i].text[0]=='@' || (sl[i].text[0]=='/'&&sl[i].text[1]=='/')) break;
                    ni ++;
                }
                if (ni<=0) continue;
                inner_line ++; if (inner_line<=n_skip) continue;
                line_num ++; //printf("file at %d : %s\n", inner_line, sl[0].text);
            }
            if (dimension<=0){ fprintf(stderr, "%s : error : please specify #dimension in \"%s\"\n", software_name, filename); success = false; }
        } else {
            char input[MAX_BUFFER];
            fseek(file, 0, SEEK_SET); while (fgets(input, sizeof(input), file)) {
                StringNS::string sl[MAX_WORD];
                int nw = StringNS::analysis_line(input, sl, MAX_WORD);
                int ni = 0; for (int i=0; i<nw; i++){
                    if (sl[i].text[0]=='#' || sl[i].text[0]==';' || sl[i].text[0]=='@' || (sl[i].text[0]=='/'&&sl[i].text[1]=='/')) break;
                    ni ++;
                }
                if (ni<=0) continue;
                inner_line ++; if (inner_line<=n_skip) continue;
                line_num ++; if (col_num<ni) col_num = ni;
            }
            dimension = (int) sqrt(fabs(col_num));
            if (dimension*dimension != col_num){
                int new_col_num = dimension*dimension;
                fprintf(stderr, "%s : warning : file %s contains %d column%s, only first %d column%s are considered\n", software_name, filename, col_num, col_num>1?"s":"", new_col_num, new_col_num>1?"s":"");
                col_num = new_col_num;
            }
            //printf("file %s contains %d lines, %dx%d, requires %d columns\n", filename, line_num, dimension, dimension, col_num);
        }
    }
    if (n_stationary_population>0){
        if (n_stationary_population<dimension){
            fprintf(stderr, "%s : error : %d stationary population%s specified, but %d required\n", software_name, n_stationary_population, dimension>1?"s":"", dimension); success = false;
        } else {
            double sum = 0; for (int i=0; i<dimension; i++) sum += stationary_population[i];
            if (fabs(sum-1)>MACHINE_REASONABLE_ERROR){
                for (int i=0; i<dimension; i++) stationary_population[i] /= sum;
                fprintf(stderr, "%s : warning : stationary population rescaled to ", software_name);
                for (int j=0; j<dimension; j++) fprintf(stderr, j+1<dimension?"%g ":"%g\n", stationary_population[j]);
            }
        }
    }
  // initialize memory: array and matrix
    MatrixNS::MatrixDataContainer mc; mc.init(sizeof(double)*dimension*dimension*7); MatrixNS::matrix mm[5]; MatrixNS::matrix m_dT0, m_T0inv;
    double ** M_in = NULL; double ** M_dt = NULL; double ** M_predict = NULL; int len_M_predict = line_num + 1;
    double ** K = NULL;
    double * eigen_work = NULL;
    double ** iK = NULL;
    if (success){
if (debug_level>=1) fprintf(stderr, "%s : debug : allocate memory\n", software_name);
        m_dT0.init(dimension, &mc);
        m_T0inv.init(dimension, &mc);
        M_in = init_matrix(line_num, col_num, OVERFLOWSIZE);
        M_dt = init_matrix(line_num, col_num, OVERFLOWSIZE);
        K = init_matrix(line_num, col_num, OVERFLOWSIZE);
        iK = init_matrix(2, col_num, OVERFLOWSIZE);
        n_predict = (size_t)((t_predict-t_skip)/dt);
        if (debug_level>=3) print_current_memory_allocation(sizeof(double) * 4 * dimension, false);
        eigen_work = (double*)malloc(sizeof(double) * 4 * dimension); __global_memory_allocated += sizeof(double) * 4 * dimension;
        if (debug_level>=2) print_current_memory_allocation(sizeof(double) * 4 * dimension);
        if (t_predict>0){
            if (n_predict<line_num){ n_predict = line_num;
                if (t_predict>MACHINE_REASONABLE_ERROR) fprintf(stderr, "%s : prediction extend to -predict %g\n", software_name, n_predict*dt+t_skip);
            }
        }
        if (n_predict>0){
            M_predict = init_matrix(len_M_predict, col_num, OVERFLOWSIZE);
        }
        if (!M_in||!M_dt||(n_predict>0&&!M_predict)||!K||!iK || !mc.root || !eigen_work){ fprintf(stderr, "%s : error : malloc failure\n", software_name); success = false; }
        if (success){ for (int i=0; i<sizeof(mm)/sizeof(mm[0]); i++) mm[i].init(dimension, &mc); }
    }
  // read input data: TPM in function of time, store into M_in[time][col*row]
    if (success){
if (debug_level>=1) fprintf(stderr, "%s : debug : read TPM inputs\n", software_name);
        int iline = 0; char input[MAX_BUFFER]; int inner_line = 0;
        if (file_is_list){
            char fninput[MAX_FILE];
            fseek(file, 0, SEEK_SET); while (fgets(fninput, sizeof(fninput), file)) {
                StringNS::string sl[MAX_WORD];
                int nw = StringNS::analysis_line(fninput, sl, MAX_WORD, true);
                int ni = 0; if (sl[0].text[0]=='#') {
                    continue;
                } else for (int i=0; i<nw; i++){
                    if (sl[i].text[0]=='#' || sl[i].text[0]==';' || sl[i].text[0]=='@' || (sl[i].text[0]=='/'&&sl[i].text[1]=='/')) break;
                    ni ++;
                }
                if (ni<=0) continue;
                inner_line ++; if (inner_line<=n_skip) continue;
                FILE * filei = fopen(fninput, "r"); int icolnow = 0; if (filei){
                    while (fgets(input, sizeof(input), filei)) {
                        StringNS::string sl[MAX_WORD];
                        int nw = StringNS::analysis_line(input, sl, MAX_WORD, true);
                        int ni = 0; for (int i=0; i<nw; i++){
                            if (sl[i].text[0]=='#' || sl[i].text[0]==';' || sl[i].text[0]=='@' || (sl[i].text[0]=='/'&&sl[i].text[1]=='/')) break;
                            ni ++;
                        }
                        if (ni<=0) continue;
                        for (int i=0; icolnow<col_num && i<nw; i++) M_in[iline][icolnow++] = atof(sl[i].text);
                    }
                    fclose(filei);
                } else {
                    fprintf(stderr, "%s[%d] : cannot open \"%s\"\n", filename, inner_line+1, fninput); success = false;
                }
                iline ++;
            }
            if (dimension<=0){ fprintf(stderr, "%s : error : please specify #dimension in \"%s\"\n", software_name, filename); success = false; }
        } else {
            fseek(file, 0, SEEK_SET); while (fgets(input, sizeof(input), file)) {
                StringNS::string sl[MAX_WORD];
                int nw = StringNS::analysis_line(input, sl, MAX_WORD, true);
                int ni = 0; for (int i=0; i<nw; i++){
                    if (sl[i].text[0]=='#' || sl[i].text[0]==';' || sl[i].text[0]=='@' || (sl[i].text[0]=='/'&&sl[i].text[1]=='/')) break;
                    ni ++;
                }
                if (ni<=0) continue;
                inner_line ++; if (inner_line<=n_skip) continue;
                for (int i=0; i<col_num && i<nw; i++) M_in[iline][i] = atof(sl[i].text);
                iline ++;
            };
        }
        //for (int j=0; j<iline; j++){ for (int i=0; i<col_num; i++) printf("%8.3f", M_in[j][i]); printf("\n"); }
    }


  // calculate time derivatives of input TPMs
    bool continue_work = true;
    if (success && continue_work){
if (debug_level>=1) fprintf(stderr, "%s : debug : calculate derivatives\n", software_name);
        mc.push_ip();
      // 1. dTPM/dt and d^2TPM/dt^2
        for (int row=0; row<line_num; row++){
            for (int col=0; col<col_num; col++){
                if (row==0){
                    M_dt[row][col] = (M_in[row+1][col]-M_in[row][col])/dt;
                } else if (row+1==line_num){
                    M_dt[row][col] = (M_in[row][col]-M_in[row-1][col])/dt;
                } else {
                    M_dt[row][col] = (M_in[row+1][col]-M_in[row][col])/dt;
                }
            }
        };
      // m_T0inv (T(0)^-1) and dT0 (dTo/dt, not dT(0)/dt) from definition or interia requirement
        mm[0] = M_in[0]; mm[1] = M_dt[0];
        m_T0inv = mm[0].inverse();
        m_dT0 = m_T0inv * mm[1];
      // output dT if required
        if (n_predict<=0){
            if (show_K==1) for (int i=0; i<line_num; i++) print_matrix_data(&M_dt[i][0], dimension, "", output_format, file_is_list);
        }
        //printf("Y[0]=\n"); print_matrix_data(y[0], dimension, "\n"); continue_work = false;
      // test
        /*
        print_matrix_data(M_in[0], dimension, "\n"); print_matrix_data(m_T0inv.m, dimension, "\n"); mm[0] = M_in[0]; print_matrix_data((m_T0inv * mm[0]).m, dimension, "\n");
        print_matrix_data(M_dt[0], dimension, "\n"); print_matrix_data(m_dT0.m, dimension, "\n"); mm[0] = M_dt[0]; print_matrix_data((m_dT0 * mm[0]).m, dimension, "\n");
        success = false;//*/
        mc.pop_ip();
    }

  // direct calculate K with greedy algorithm
    if (success && continue_work){
if (debug_level>=1) fprintf(stderr, "%s : debug : calculate memory\n", software_name);
        mc.push_ip();
        MatrixNS::matrix mK, mT, mdT;

        mK.harness(K[0], dimension, &mc); mK = (double)0;

        for (int krow=1; krow<line_num; krow++){
if (debug_level>=2) fprintf(stderr, "%s : debug : calculate K[%d] (out of %d)%s", software_name, krow+1, line_num, krow+1>=line_num?"\n":"\r");
            mc.push_ip();
            mK.harness(K[krow], dimension, &mc); mT.harness(M_in[krow], dimension, &mc); mdT.harness(M_dt[krow], dimension, &mc);

          // mm[0] = T dT0 - T: let mm[0] be the sum of the whole history
            mm[1] = M_in[krow]; mm[1].product(m_dT0, mm[0]); mm[0] -= mdT;
          // mm[0] = T0 K0: take away the previous memory from mm[0]
            for (int i=0; i<krow && i<line_num; i++){
                mm[1] = M_in[krow-i]; mm[2] = K[i];
                mm[1].product(mm[2], mm[3]);
                mm[3] *= dt;
                mm[0] -= mm[3];
            }
          // extract current memory term K
            mm[1] = m_T0inv; mm[1].product(mm[0], mK); mK /= dt;

            mc.pop_ip();
        }

      // print K if required
        if (n_predict<=0){
            if (show_K==0){
                for (int i=0; i<line_num; i++) print_matrix_data(&K[i][0], dimension, "", output_format, file_is_list);
            } else if (show_K==2){  //iK
                for (int i=0; i<2; i++) for (int j=0; j<col_num; j++) iK[i][j] = 0;
                for (int i=0; i<line_num; i++){
                    for (int j=0; j<col_num;j++) iK[0][j] += K[i][j] * dt;
                    //double value = 0; for (int j=0; j<col_num;j++) value += iK[i][j]*iK[i][j]/col_num;
                    print_matrix_data(&iK[0][0], dimension, "", output_format, file_is_list);
                    //printf(output_format, iK); printf("\n");
                }
            } else if (show_K==3){  //miK
                for (int i=0; i<2; i++) for (int j=0; j<col_num; j++) iK[i][j] = 0;
                for (int i=0; i<line_num; i++){
                    for (int j=0; j<col_num;j++) iK[0][j] += K[i][j] * dt;
                    double value = 0; for (int j=0; j<col_num;j++) value += iK[0][j]*iK[0][j]/col_num;
                    printf(output_format, sqrt(value)); printf("\n");
                }
            } else if (show_K==4){  //m1iK
                for (int i=0; i<2; i++) for (int j=0; j<col_num; j++) iK[i][j] = 0;
                for (int i=0; i<line_num; i++){
                    for (int j=0; j<col_num;j++) iK[0][j] += K[i][j] * dt;
                    double value = 0; for (int j=0; j<dimension;j++) for (int k=0; k<dimension;k++) value += iK[0][j*dimension+k]*iK[0][j*dimension+k]*(stationary_population? stationary_population[k] : 1.0/col_num);
                    printf(output_format, sqrt(value)); printf("\n");
                }
            } else if (show_K==5){  //m2iK
                for (int i=0; i<2; i++) for (int j=0; j<col_num; j++) iK[i][j] = 0;
                for (int i=0; i<line_num; i++){
                    for (int j=0; j<col_num;j++) iK[0][j] += K[i][j] * dt;
                    double value = 0; for (int j=0; j<dimension;j++) for (int k=0; k<dimension;k++) value += iK[0][j*dimension+k]*iK[0][j*dimension+k]*(stationary_population? stationary_population[k]*stationary_population[k] : 1.0/col_num);
                    printf(output_format, sqrt(value)); printf("\n");
                }
            }
        }

      // debug: check reconstruction of y
        /*
        for (int row=0; row<line_num; row++){
            mm[1] = 0.0;
            for (int i=0; i<=row && i<line_num; i++){
                MatrixNS::matrix mTnow, mKnow;
                mTnow.harness(M_in[row-i], dimension, &mc); mKnow.harness(K[i], dimension, &mc);
                mm[1] += mTnow * mKnow * dt;
            }
            printf("Debug: reconstruct y at %d\n", row);
            printf("  y:    "); print_matrix_data(y[row], dimension, ",", "%9.4f");
            printf("  p:    "); print_matrix_data(mm[1].m, dimension, ",", "%9.4f");
        }
        continue_work = false;
        //*/
        mc.pop_ip();
    }

  // direct extend the prediction of T to -predict
    if (success && continue_work){
if (debug_level>=1) fprintf(stderr, "%s : debug : predict kinetics\n", software_name);
        mc.push_ip();
      // predict forward, mm[0]=T_now, mm{1]=Tp_now
        mm[0] = M_in[0]; mm[1] = M_dt[0];
        for (int row=0; row<n_predict; row++){
if (debug_level>=2) fprintf(stderr, "%s : debug : %s T[%d] (out of %d)%s", software_name, row+1<=line_num?"recover":"predict", row+1, n_predict, row+1>=n_predict?"\n":"\r");
            mc.push_ip();
            mm[0].product(m_dT0, mm[1]);
            for (int i=0; i<=row && i<line_num; i++){
                mc.push_ip();
                MatrixNS::matrix mTnow, mKnow; mTnow.harness(M_predict[(row-i)%len_M_predict], dimension, &mc); mKnow.harness(K[i], dimension, &mc);
                /*
                if (i<5){
                    printf("Debug at %d out of %4d:\n", i, row);
                    printf("  T:    "); print_matrix(&mTnow, ",", "%9.4f");
                    printf("  K:    "); print_matrix(&mKnow, ",", "%9.4f");
                    mm[5] = mTnow * mKnow;
                    printf("  T*K:  "); print_matrix(&mm[5], ",", "%9.4f");
                    printf("  dT_:  "); print_matrix(&mm[1], ",", "%9.4f");
                    mm[2] = mm[1] - mTnow * mKnow * dt;
                    printf("  dT:   "); print_matrix(&mm[2], ",", "%9.4f");
                }; mm[3] = mm[0];//*/
                mTnow.product(mKnow, mm[2]); mm[2] *= dt;
                mm[1] -= mm[2];
                mc.pop_ip();
            }
            MatrixNS::matrix mMnow; mMnow.harness(M_predict[row%len_M_predict], dimension, &mc);
            mMnow = mm[0].m; mm[1] *= dt; mm[0] += mm[1];
            /*
            printf("dT:     "); print_matrix(&mm[1], ",", "%9.4f");
            printf("T_p:    "); print_matrix(&mm[3], ",", "%9.4f");
            printf("T_in:   "); print_matrix_data(M_in[row], dimension, ",", "%9.4f");
            printf("T_pn:   "); print_matrix(&mm[0], ",", "%9.4f");
            printf("---------------------\n");//*/

            print_matrix_array(&M_predict[row%len_M_predict], 0, 1, dimension, &mm[1], &mc, eigen_work, dt, t_skip+row*dt, predict_output, output_format, file_is_list);

            mc.pop_ip();
        }
        mc.pop_ip();
    }

    if (debug_level>0){
        char memory_out_buffer[64]; print_memory_allocation(__global_memory_allocated, memory_out_buffer, sizeof(memory_out_buffer));
        fprintf(stderr, "%s : debug : memory allocated : %s \n", software_name, memory_out_buffer);
    }

  // done of all
    if (stationary_population) free(stationary_population);
    if (file) fclose(file); if (mc.root) mc.dispose();
    if (M_in) free(M_in); if (M_dt) free(M_dt);
    if (K) free(K); if (eigen_work) free(eigen_work); if (M_predict) free(M_predict);
    return 0;
}
