#include <string.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

int signet(int a, int b);
int testr(int a, int b, int c);
void signe_Rec(int n,vector <double> &pw,vector <double> &wrange);
void wilc_Rec(int m, int n,vector <double> &pw,vector <double> &wrange);
long int wilcoxon_exact(int *m, vector <double> &pw,vector <double> &wrange);

/////////////////////////////////////////////////////////////////////

//Вспомогательная функция двухвыборочного критерия Уилкоксона

int testr(int a, int b, int c) {
    if (c < 0) return 0;
    if (a == 0 || b == 0) {
        if (c == 0) return 1;
        return 0;
    }
    return -1;
}

////////////////////////////////////////////////////////////////////

//Вспомогательная функция критерия знаковых рангов Уилкоксона

int signet(int a, int b) {
    if (b < 0) return 0;
    if (a == 0) {
        if (b == 0) return 1;
        return 0;
    }
    return -1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Точное распределение критерия знаковых рангов Уилкоксона по рекурентной формуле(в книге Хетманспергера ошибка, смотри signet)

void signe_Rec(int n,vector <double> &pw,vector <double> &wrange) {

    //pw - вероятность
    //wrange -статистика критерия

    double z;
    int i,j,k,s,h1,h2,mn;
    mn=n*(n+1)/2+1;
    int**w;

    w = new int* [n+1];
    for (i = 0;i <=n;i++) {
        w[i] = new int[mn+1];
    }
    for (i = 0; i <=n; i++) {
        for (j = 0; j<=mn; j++) w[i][j] = 0;
    }
    s=0;
    z = pow(2, n);
    for (k=0;k<mn;k++) {
        for (i=1;i<= n;i++) {
            h1=signet(i-1,k);
            if (h1== -1) h1=w[i-1][k];
            h2=signet(i-1,k-i);
            if (h2==-1) h2=w[i-1][k-i];
            w[i][k]=h1+h2;
        }
        s+= w[n][k];
        pw.push_back(s/z);
        wrange.push_back(k);
    }

    delete [] w;
}

////////////////////////////////////////////////////////////////////////////////

//Точное распределение критерия Уилкоксона по рекурентной формуле

void wilc_Rec(int m, int n,vector <double> &pw,vector <double> &wrange) {

    //pw - вероятность
    //wrange -статистика критерия

    int mn,i,j,k;
    double s,h1,h2;

    mn=m*n+1;
    double***w;

    w = new double**[m+1];
    for(i=0;i<=m;i++) {
        w[i] = new double*[n+1];
        for (j=0;j<=n;j++) w[i][j]=new double[mn+1];
    }

    for(i=0;i<=m;i++) {
        for(j=0;j<=n;j++) {
            for(k=0;k<=mn;k++) w[i][j][k]= 0;
        }
    }
    
    s = 0;
    for (k=0;k<mn;k++) {
        for (i=1;i<=m;i++) {
            for (j=1;j<= n;j++) {
                h1=testr(i,j-1,k-i);
                if (h1== -1) h1=(w[i][j-1][k-i]);
                h2=testr(i-1,j,k);
                if (h2== -1) h2=(w[i-1][j][k]);
                w[i][j][k]= (h1*j+h2*i)/(1.00*(i+j));
            }
        }
        s+=w[m][n][k];
        pw.push_back(s);
        wrange.push_back(k+m*(m+1)/2);
    }
    delete [] w;
}

/////////////////////////////////////////////////////////////////////////////////////

// Критерий Уилкоксона (точное распределение)*(аналитика)
/*
  AS 62 generates the frequencies for the Mann-Whitney U-statistic.
  Users are much more likely to need the distribution function.
  Code to return the distribution function has been added at the end
  of AS 62 by Alan Miller
*/
long int wilcoxon_exact(int *m, vector <double> &pw,vector <double> &wrange) {

    int minmn, maxmn, inx;
    int n1, kk, k, j, i;
    long int mn;
    double z;
    double *work,*w,sum;

    mn = m[0]*m[1]+1; 
    maxmn = fmax(m[0],m[1]);
    minmn = fmin(m[0],m[1]);
    n1 = maxmn + 1;

    work = new double[mn+10];
    w = new double[mn+10];
   
    for (i = 1; i <= n1; i++)  w[i] = 1;
    n1++;
    for (i = n1; i <= mn; i++) w[i] = 0;
    work[1] = 0; inx = maxmn;

    for (i = 2; i <= minmn; i++) {
        work[i] = 0; inx = inx + maxmn; n1 = inx + 2; kk = 1 + inx / 2; k = i;
        for (j = 1; j <= kk; j++) {
            k++;  n1--; sum = w[j] + work[j];
            w[j] = sum; work[k] = sum - w[n1]; w[n1] = sum;
        }
    }

    for (i = 0; i < mn; i++)  pw.push_back(w[i + 1]);
    sum = 0;
    for (i = 0; i < mn; i++) {
        sum += pw[i];
        z=double(i + minmn * (minmn + 1.) / 2.); //Wilcoxon
        wrange.push_back(z);
        //wrange.push_back((double(i)); //Mann-Whitney
        pw[i] = sum;
    }
    for (i = 0; i < mn; i++) pw[i] = pw[i] / sum;

    delete[] work,w;
    return mn;
}


/////////////////////////////////////////////////////

int main() {

  int i,k,*m;
  long int j,nn;
  string st,ff;
  vector <double> pw;
  vector <double> wrange;

  ifstream inp1("wilcoxon.inp");
  inp1>>ff;
  inp1.close();
 
  ifstream inp("Inp/" + ff + ".inp");
  ofstream out("Out/" + ff + ".out");

  inp>>st;
  inp>>k;
  m=new int[k];
  inp>>st;
  for(i=0;i<k;i++) inp>>m[i];
  inp.close();

  out<<"Criterion:"<<ff<<endl;
  for(i=0;i<k;i++) out<<m[i]<<";";
  out <<endl;


  if(ff=="Wilcoxon_exact") wilc_Rec(m[0],m[1],pw,wrange);
  if(ff=="Wilcoxon_AS62") nn=wilcoxon_exact(m,pw,wrange);
  if(ff=="WilcSigne_exact") signe_Rec(m[0],pw,wrange);


  nn=pw.size();
  out << "Size=" << nn << endl;
  for (j=0;j<nn;j++) out<<(j+1)<<". "<<wrange[j]<<"  "<<pw[j]<<endl;
  out <<endl;
  out.close();

  pw.clear();
  wrange.clear(); 
  delete [] m;
  return 0;
}
