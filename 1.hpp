void synchronize(int);
int Solve(int, double*, double*, double*, int*,int,int,int*);
double f(int, int,int);
double newyaz(int, double*, double*, double*,int/*,int*/,double);
double err(int, double *);
double norm(int, double*);
double norm_paral(int, double*,int,int,double);
int  matrix_input_file(int,double *, double*, char*);
void solution_print(int, double*, double *, double *, int);                  //можно ещё 1 int
void matrix_input_formula(int,double*,double*);
double newyaz_nepar(int, double*, double *, double*);
struct ARGS_gauss{
    int n;
    double *a;
    double *b;
    double *x;
    int *ind;
    int thread_num;
    int p;
    int* er;
    double time;
};
struct ARGSS_newyaz{
    int n;
    double *a0;
    double *b0;
    double *x;
    int *ind;
    int thread_num;
    int p;
    int* er;
    double time;
double sum_4tread;
};


struct ARGSS_norm{
    int n;
    double *a0;
    double *b0;
    double *x;
    int *ind;
    int thread_num;
    int p;
    int* er;
    double time;
        double nia;
double sum_4tread;
};

