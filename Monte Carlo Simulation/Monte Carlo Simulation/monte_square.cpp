using namespace std;

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define L 100
#define J 1.0
#define T_i 0.5
#define dT 0.5f
#define T_f 10.0
#define AVG_SIZE ceil((T_f - T_i)/dT) + 1
#define SAMPLE_SIZE (L * 1000)

double **monte();
int flip_sign(int[],double);
inline int rand_int(double, double);
inline double rand_double(double, double);
int calc_dE(int[],int);
int sum_int(int[],int);
double sum_double(double[],int);
void data(double **,double **,int);
void data_gnuplot(int,double **,int);
void auto_corr(double *);
void auto_corr_gnuplot();
void linfit_gnuplot();
void std_dev(double *);

int main()
{
    int n = 10;
    double **O;
    
    O = monte();
    auto_corr(O[n]);
    std_dev(O[n]);
    
    return 0;
} // main()

double **monte()
{
    srand((unsigned int)time(0));
    
    int i = 0, j = 0, k, count, N_m = 0;
    double T, **E, **M_sq;
    
    E = new double *[(int)AVG_SIZE];
    M_sq = new double *[(int)AVG_SIZE];
    for(k = 0; k < AVG_SIZE; k++)
    {
        E[k] = new double[SAMPLE_SIZE + 1];
        M_sq[k] = new double[SAMPLE_SIZE + 1];
    } // for
    
    T = T_i;
    
    while(T <= T_f)
    {
        int S[L], E_sum, E_f, dE, r;
        
        for(k = 0; k < L; k++){
            S[k] = 1;
        } // for i
    
        // Initialize E and dE.
        E_sum = 0;
        for(k = 0; k < (L - 1); k++){
            E_sum = E_sum + S[k] * S[k + 1];
        } // for n
        E_sum = E_sum + S[L - 1] * S[0];
        E_f = -J * E_sum;
        
        for(count = 0; count < SAMPLE_SIZE; count++)
        {
            r = flip_sign(S,T);
            
            if(r >= 0){
                dE = calc_dE(S,r);
                E_f = E_f + dE;
                S[r] = -1 * S[r];
            } // if r > 0
            
            E[i][j] = E_f;
            M_sq[i][j] = pow(sum_int(S,L),2);
            
            j++;
            N_m++;
        } // for count
        
        j = 0;
        
        E[i][SAMPLE_SIZE] = sum_double(E[i],SAMPLE_SIZE)/SAMPLE_SIZE;
        M_sq[i][SAMPLE_SIZE] = sum_double(M_sq[i],SAMPLE_SIZE)/SAMPLE_SIZE;
        
        i++;
        
        T = T + dT;
    } // while
    
    data(E,M_sq,AVG_SIZE);
    
    return E;
} // monte()


int flip_sign(int S[],double T)
{
    int r, r1, dE;
    double r2, beta;
    
    r1 = rand_int(0,L);
    dE = calc_dE(S,r1);
    
    if(dE < 0){
        r = r1;
    } // if
    else
    {
        r2 = rand_double(0,1);
        
        beta = 1/T;
    
        if(r2 < exp(-dE * beta)){
            r = r1;
        } // if
        else{
            r = -1;
        } // else
    } //  else
    
    return r;
} // flip_sign


inline int rand_int(double x0, double x1)
{
    double r;
    r = x0 + (x1 - x0) * rand() / RAND_MAX;
    
    return r;
} // rand_int


inline double rand_double(double x0, double x1)
{
    double r;
    r = x0 + (x1 - x0) * (double)rand() / (double)RAND_MAX;
    
    return r;
} // rand_double


int calc_dE(int S[],int r)
{
    int dE = 0;
    
    if(r == 0){
        dE = 2 * J * S[r] * (S[L - 1] + S[r + 1]);
    } // if r == 0
    if(r == L - 1){
        dE = 2 * J * S[r] * (S[r - 1] + S[0]);
    } // if r == L
    if(r != 0 && r != (L - 1)){
        dE = 2 * J * S[r] * (S[r - 1] + S[r + 1]);
    } // if r != 0 && r != L
    
    return dE;
} // calc_dE()


int sum_int(int A[],int size)
{
    int i, s = 0;
    
    for(i = 0; i < size; i++){
        s += A[i];
    } // for i
    
    return s;
} // sum_int()

double sum_double(double A[],int size)
{
    int i, s = 0;
    
    for(i = 0; i < size; i++){
        s += A[i];
    } // for i
    
    return s;
} // sum_double()


void data(double **E,double **M_sq,int size)
{
    int i, signal;
    double T = T_i, E_exp, M_sq_exp, **data;
    
    data = new double *[5];
    for(i = 0; i < 5; i++)
        data[i] = new double[size];
    
    for(i = 0; i < size; i++)
    {
        data[0][i] = T; // Temperature
        data[1][i] = E[i][SAMPLE_SIZE]; // Computed E
        
        if(L == 2)
            E_exp = -L * J * tanh(2 * 1/T * J);
        else
            E_exp = -L * J * (tanh(1/T * J) + pow(tanh(1/T * J),(L - 1))) / (1 + pow(tanh(1/T * J),L));
        
        data[2][i] = E_exp; // Expected E
        data[3][i] = M_sq[i][SAMPLE_SIZE]; // Computed M^2
        
        M_sq_exp = L * (1 + tanh(1/T * J)) / (1 - tanh(1/T * J));
        data[4][i] = M_sq_exp; // Expected M^2
        
        T = T + dT;
    } // for i
    
    /* for(i = 0; i < size; i++)
        printf("T: %f   E_avg: %f   E_exp: %f   Diff: %f   M^2_avg: %f   M^2_exp: %f\n",data[0][i],data[1][i],data[2][i],(double)abs(data[1][i] - data[2][i]),data[3][i],data[4][i]);
    */
    
    signal = 1; // E
    data_gnuplot(signal,data,size);
    signal = 2; // M
    data_gnuplot(signal,data,size);
} // data()

void data_gnuplot(int signal,double **data,int size)
{
    FILE *pipe;
    char *title;
    int i, j;
    
    if (signal == 1) // E
    {
        title = (char * )"E Average vs. T";
        j = 1;
    } // if signal == 1
    else // M
    {
        title = (char *)"M^2 Average vs. T";
        j = 3;
    } // else
    
    pipe = popen("/opt/local/bin/gnuplot","w");
    
    if (pipe) {
        fprintf(pipe,"set title \"%s vs. Temperature\"\n",title);
        fprintf(pipe,"set xrange [%f:%f]\n",T_i,T_f);
        fprintf(pipe,"set xlabel \"Temperature\"\nset ylabel \"%s\"\n",title);
        fprintf(pipe,"plot '-' with lines title 'Computed Values','-' with lines title 'Expected Values'\n");
        
        fflush(pipe);
        
        for (i = 0; i < size; i++)
            fprintf(pipe,"%lf %lf\n",data[0][i],data[j][i]);
        fprintf(pipe,"exit \n");
        for (i = 0; i < size; i++)
            fprintf(pipe,"%lf %lf\n",data[0][i],data[j + 1][i]);
        fprintf(pipe,"exit \n");
    } // if
    else
        printf("Gnuplot not found...");
} // data_gnuplot()


void auto_corr(double *O)
{
    int i, t, tau, N_m = SAMPLE_SIZE, N_p, lim;
    double sum = 0, A, x[N_m], y[N_m], ln_y[N_m], O_avg;
    ofstream file;
    
    O_avg = O[SAMPLE_SIZE];

    for (t = 0; t < N_m; t++)
    {
        x[t] = t;
        y[t] = O[t];
    } // for t
    
    lim = 200;
    for (tau = 1; tau <= lim; tau++)
    {
        // A(tau)
        N_p = N_m - tau;
        sum = 0;
        for(t = 0; t <= N_p; t++)
            sum = sum + (O[t] - O_avg) * (O[t + tau] - O_avg);
    
        A = 1.0 / (double)N_p * sum;
        ln_y[tau] = log(abs(A));
    } // for tau
    
    file.open ("O.txt");
    for (i = 0; i < N_m; i++)
        file << x[i] << " " << y[i] << endl;
    file.close();
    
    file.open ("ln.txt");
    for (i = 0; i <= lim; i++)
        file << x[i] << " " << ln_y[i] << endl;
    file.close();
    
    // auto_corr_gnuplot();
    linfit_gnuplot();
} // auto_corr()


void auto_corr_gnuplot()
{
    FILE *pipe;
    pipe = popen("/opt/local/bin/gnuplot","w");
    
    if (pipe) {
        fprintf(pipe,"set title \"O(t) vs. t\"\n");
        fprintf(pipe,"plot 'O.txt' with lines title 'O(t)''\n");
        
        fflush(pipe);
        fprintf(pipe,"exit \n");
    } // if
    else
        printf("Gnuplot not found...");
} // auto_corr_gnuplot()


void linfit_gnuplot()
{
    FILE *pipe;
    
    pipe = popen("/opt/local/bin/gnuplot","w");
    
    if (pipe) {
        fprintf(pipe,"set title \"ln(A(t)) vs. t\"\n");
        
        fprintf(pipe,"f(x) = m*x + b\n");
        fprintf(pipe,"fit f(x) 'ln.txt' via m,b\n");
        fprintf(pipe,"plot 'ln.txt' with lines title 'ln(A(t))', f(x) title 'Linear Fit'\n");
        
        fflush(pipe);
        fprintf(pipe,"exit \n");
    } // if
    else
        printf("Gnuplot not found...");
} // linfit_gnuplot()


void std_dev(double *O)
{
    int i, N_m = SAMPLE_SIZE, K;
    double T = dT, O_avg, Osq[N_m], Osq_avg, O_std, pct_error, slope, tau;
    
    cin >> slope;
    
    tau = -1.0 / slope;
    K = ceil(double(N_m)/tau);
    
    O_avg = O[SAMPLE_SIZE];
    
    for (i = 0; i < N_m; i++)
        Osq[i] = pow(O[i],2);
    Osq_avg = 1.0 / (double)N_m * sum_double(Osq,N_m);
    
    O_std = sqrt(abs(Osq_avg - pow(O_avg,2))/K);
    pct_error = abs(O_std / O_avg * 100);
    
    cout << "T: " << T << "   O_avg: " << O_avg << "   Osq_avg: " << Osq_avg << "   O_std: " << O_std << "   Pct_error: " << pct_error << endl;
    
    T = T + dT;
} // std_dev()