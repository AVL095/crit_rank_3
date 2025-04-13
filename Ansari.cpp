#include <string.h>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

int start1(int n, double* f);
int start2(int n, double* f);
double *frqadd(double *f1,double *f2, int l1in, int l1out, int l2, int nstart);
double *imply(double *f1,double *f2, int l1in, int l1out, int l2, int noff);
int gscale(int test, int other, double *pw);

//***********Точное распределение критерия Ансари-Брэдли*(аналитика)*************

int start1(int n, double *f) {
    int i, lout;

    lout = floor(1 + int(n / 2));
    for (i = 1; i <= lout; i++) f[i] = 2;
    if ((n % 2) == 0) f[lout] = 1;
    return(lout);
}

//////////////////////////////////////////////

int start2(int n, double* f) {
    int one, two, three, four, i, j, a, b, lt1, ndo, nu, lout;

    one = 1; two = 2; three = 3; four = 4;
    nu = n - n % 2;
    j = nu + 1; lout = j; lt1 = lout + 1;
    ndo = floor(lt1 / 2);
    a = one; b = three;
    for (i = 1; i <= ndo; i++) {
        f[i] = a; f[j] = a; j = j - 1; a = a + b; b = four - b;
    }
    if (nu == n) return(lout);
    nu = ndo + 1;
    for (i = nu; i <= lout; i++) f[i] = f[i] + two;
    f[lt1] = two; lout = lt1;
    return(lout);
}

/////////////////////////////////////////////////////////

double* frqadd(double *f1,double *f2, int l1in, int l1out, int l2, int nstart) {

    int i1, i2, nxt;
    double *fadd;

    fadd=new double[1000];

    i2 = 1;
    for (i1 = nstart; i1 <= l1in; i1++) {
        f1[i1] = f1[i1] + 2 * f2[i2];
        i2++;
    }
    nxt = l1in + 1;
    l1out = l2 + nstart - 1;
    for (i1 = nxt; i1 <= l1out; i1++) {
        f1[i1] = 2 * f2[i2];
        i2++;
    }
    nstart++;
    fadd[1] = l1out; fadd[2] = nstart;
    return fadd;
}

/////////////////////////////////////////////////

double *imply(double *f1, double *f2, int l1in, int l1out, int l2, int noff) {
    int  i2, i1, j2, j1, j2min, ndo;
   double sum,diff,*fimply;

    fimply=new double[1000];

    i2 = 1 - noff; j1 = l1out; j2 = l1out - noff; l2 = j2;
    j2min = floor((j2 + 1) / 2);
    ndo = floor((l1out + 1) / 2);

    for (i1 = 1; i1 <= ndo; i1++) {
        if (i2 > 0) {
            sum = f1[i1] + f2[i2];
            f1[i1] = sum;
        }
        else {
            sum = f1[i1];
        }
        i2 = i2 + 1;
        if (j2 >= j2min) {
            if (j1 <= l1in) {
                diff = sum - f1[j1];
            }
            else {
                diff = sum;
            }

            f2[i1] = diff; f2[j2] = diff; j2 = j2 - 1;
        }
        f1[j1] = sum; j1 = j1 - 1;
    }
    fimply[1] = l1out; fimply[2] = l2; fimply[3] = noff;
    return fimply;
}

////////////////////////////////////////////////////////////////////////

int gscale(int test, int other, double* pw) {

    int  i, m, lres, mm1, nm1, nm2, ier, mnow, ks, j, ndo;
    int n, ln1, ln2, nc, l1out, l2out, n2b1, n2b2, ln3, kk,z;
    double *a2, *a3,*fadd, *fimply,ai;
    bool symm;

    ln1 = 0; ln2 = 0; l1out = 0; l2out = 0; ln1 = 0; ln2 = 0; ln3 = 0; n2b1 = 0; n2b2 = 0; kk = 0; nc = 0;
    ndo = 0; ier = 0; ks = 0; j = 0; i = 0;

    m = test;
    if (m < 0) return(0);
    n = other;
    lres = 1 + floor(m*n/2);
    fadd = new double[1500];
    fimply = new double[1500];
    a2 = new double[1500];
    a3 = new double[1500];

    symm = false;
    z = (m + n) % 2;
    if (z == 0) symm = true;
    mm1 = m - 1;
    //*****************************************************         
    if (m <= 2) {
        if (mm1 < 0) {
            pw[1] = 1; return(lres);
        }

        if (mm1 == 0) ln1 = start1(n, pw);
        if (mm1 > 0)  ln1 = start2(n, pw);
        if (symm || (other > test)) return(lres);
        j = lres;
        ndo =floor(lres / 2);
        for (i = 1; i <= ndo; i++) {
            ai = pw[i];
            pw[i] = pw[j];
            pw[j] = ai;
            j = j - 1;
        }
        return(lres);
    }
    //***********************************************************
    nm1 = n - 1; nm2 = n - 2; mnow = 3; nc = 3; ier = 0;
    //**************************************************************
    while (true) {
        if (ier == 0) {
            if ((n % 2) != 1) {
                n2b1 = 3;
                n2b2 = 2;
                ln1 = start2(n, pw);
                ln3 = start2(nm2, a3);
                ln2 = start1(nm1, a2);
                //***********************************************************************          
                fadd = frqadd(a2, a3, ln2, l2out, ln3, n2b2);
                l2out = fadd[1]; n2b2 = fadd[2];
                ln2 = ln2 + nm1;
                fimply = imply(a2, a3, l2out, ln2, j, nc);
                ln2 = fimply[1]; j = fimply[2]; nc = fimply[3];
                //***************************************************************************
                nc = nc + 1;
                if (mnow == m) break;
                mnow = mnow + 1;
            }
            else {
                n2b1 = 2;
                n2b2 = 3;
                ln1 = start1(n, pw);
                ln2 = start2(nm1, a2);
            }
        }
        //*******************************************************************
        fadd = frqadd(pw, a2, ln1, l1out, ln2, n2b1);
        l1out = fadd[1]; n2b1 = fadd[2];
        ln1 = ln1 + n;
        fimply = imply(pw, a3, l1out, ln1, ln3, nc);
        ln1 = fimply[1]; ln3 = fimply[2]; nc = fimply[3];
        nc = nc + 1;
        if (mnow == m) break;
        mnow = mnow + 1;
        fadd = frqadd(a2, a3, ln2, l2out, ln3, n2b2);
        l2out = fadd[1]; n2b2 = fadd[2];
        ln2 = ln2 + nm1;
        fimply = imply(a2, a3, l2out, ln2, j, nc);
        ln2 = fimply[1]; j = fimply[2]; nc = fimply[3];
        //***************************************************************************
        nc = nc + 1;
        if (mnow == m) break;
        mnow = mnow + 1;
        ier = 1;
    }
    //********************************************************************
    if (symm) return(lres);
    ks = floor((m + 3) / 2);
    j = 1;
    for (i = ks; i <= lres; i++) {
        if (i > ln1) {
            pw[i] = a2[j];
        }
        else {
            pw[i] = pw[i] + a2[j];
        }
        j = j + 1;
    }
    if (other < test) return(lres);
    j = lres;
    ndo = floor(lres / 2);
    for (i = 1; i <= ndo; i++) {
        ai = pw[i];
        pw[i] = pw[j];
        pw[j] = ai;
        j = j - 1;
    }
   delete [] fadd, fimply,a2,a3;
    return(lres);
}


/////////////////////////////////////////////////////

int main() {

  int i,k,*m,*wrange;
  string st,ff;
  double *w,sum;
  int min_val, max_val,a0,nrows;

  ff="Ansari_Exact";

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

    min_val = m[0];
    max_val = m[1];
    nrows=1+m[0]*m[1]/2;
    w = new double[nrows+10];
    wrange=new int[nrows+10];

    nrows = gscale(min_val, max_val, w);

    a0 = floor((min_val + 1) / 2) * (1 + floor(min_val / 2));
    sum = 0;
    for (i=1;i<=nrows;i++) {
       wrange[i]=a0+i-1;
       sum+=w[i];
       w[i]=sum;
    }
    
    for (i=1;i<=nrows;i++) {
        w[i]=w[i]/sum;
    }

  out << "Size=" << nrows << endl;
  for (i=1;i<=nrows;i++) out<<(i+1)<<". "<<setprecision(0)<<fixed<<wrange[i]<<"  "<<setprecision(20)<<fixed<<w[i]<<endl;
  out.close();

  delete [] m,w,wrange;
  return 0;
}
