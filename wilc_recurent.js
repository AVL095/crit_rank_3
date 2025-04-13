//***********Точное распределение критерия Уилкоксона******************************************
function wilc_Rec(m,n,range,w,pw) {

/*
  AS 62 generates the frequencies for the Mann-Whitney U-statistic.
  Users are much more likely to need the distribution function.
  Code to return the distribution function has been added at the end
  of AS 62 by Alan Miller
*/


 let work=[],minmn,maxmn,inx,i;
 let n1,l,k,j,sum,mn,s;
  
  mn=m*n+1;
  maxmn=Math.max(m,n);
  minmn=Math.min(m,n);
  n1=maxmn+1;
 
     for(i=1;i<=n1;i++) w[i]=1;
      n1=n1+1;
      for(i=n1;i<=mn;i++) w[i]=0;
  
     work[1]=0;inx=maxmn;

    for(i=2;i<=minmn;i++) {
        work[i]=0;inx=inx+maxmn;n1=inx+2;l=1+inx/2;k=i;

       for(j=1;j<=l;j++) {
          k=k+1;n1=n1-1;
          sum=w[j]+work[j];
          w[j]=sum;work[k]=sum-w[n1];w[n1]=sum;
       }
    }
  
  
  //w[i] - Frequencies
   //for(i=1;i<=mn;i++) pw[i]=w[i];
 
sum=0;
for(i=1;i<=mn;i++) sum=sum+w[i];

    s=0;
     for(i=1;i<=mn;i++) {
        s=s+w[i];
        range[i]=i-1+minmn*(minmn+1)/2; //Wilcoxon
        //range[i]=i-1; //Mann-Whitney
        pw[i]=s/sum;
     }
 
    return(mn);
 }
//*****************Точное распределение критерия Уилкоксона по рекурентной формуле*********************
function wilc_Rec_1(m,n,P) { 
       var W=[];
        w=m*n;s=0;

       for(i=1;i<=m;i++) {
           W[i]=[];
           for(j=1;j<=n;j++) {
             W[i][j]=[];
             for(k=0;k<=w;k++)  W[i][j][k]=0;
              }}

       for(k=0;k<=w;k++) {
         for(i=1;i<=m;i++) {
             for(j=1;j<=n;j++) {
                 h1=test_Rec(i,j-1,k-i);
                 if(h1==-1)  h1=W[i][j-1][k-i];
                 h2=test_Rec(i-1,j,k);
                 if(h2==-1)  h2=W[i-1][j][k];
                 W[i][j][k]=(h1*j+h2*i)/(i+j);
           }}
               s+=W[m][n][k];
              P[k]=s;
        }
 }
//*****************************************************
function test_Rec(a,b,c) {
 if(c<0) return 0;
 if(a==0 || b==0) {
   if(c==0) return 1;
   return 0;
 }
 return -1;
}

//*****Точное распределение критерия знаковых рангов Уилкоксона************************************************

 function wsigne_Rec(n,aa,a,pw) { 
   let k1,k2,k3,i,j,s1,s;

   k1=3;k3=4;
  for(i=1;i<=4;i++) {aa[i]=1;a[i]=1;}
  for(j=3;j<=n;j++) {
      k1+=1;k2=k3; k3+=j;
      for(i=k1;i<=k3;i++) {
        s1=0;
        if(i<=k2) s1=a[i];
        a[i]=aa[i-k1+1]+s1;
      }
     for(i=k1;i<=k3;i++) aa[i]=a[i];
  }
  s=0;s1=Math.pow(2,n);
 for(i=1;i<=k3;i++) {
    s+=a[i]/s1;aa[i]=i-1;pw[i]=s;
   }
}

//**Точное распределение критерия знаковых рангов Уилкоксона по рекурентной формуле(в книгах ошибка, смотри test_signe)***

function signe_Rec(n,P) { 
      let i,k,s;
      let W=[];
      let w=n*(n+1)/2+1;
      let  z=Math.pow(2,n);
      let h1,h2;

     for(i=1;i<=n;i++) {
           W[i]=[];
             for(k=0;k<=w;k++)  W[i][k]=0;
             }
       s=0;
     
        for(k=0;k<=w;k++) {
             for(i=1;i<=n;i++) {
                h1=test_signe(i-1,k);
                if(h1==-1)  h1=W[i-1][k];
                h2=test_signe(i-1,k-i);
                if(h2==-1)  h2=W[i-1][k-i];
                W[i][k]=h1+h2; 
              }
               s+=W[n][k];
               P[k]=s/z;
        }
 }
//*****************************************************
function test_signe(a,b) {
   if(b<0) return 0;
   if(a==0) {
     if(b==0) return 1;
     return 0;
   }
 return -1;
}

//////////////////////////////////////////////////////////

function Test() {

    let m=2;
    let n=2;
    let nsigne=10;
    let wrange=[],w=[],pw=[];

    //crit="wilc_Rec_1"; // Уилкоксон рекурент
    //crit="wilc_Rec";  // Уилкоксон AS62
    crit="wsigne_Rec";  // знаковых рангов Уилкоксона аналитика
    //crit="signe_Rec";  //знаковых рангов Уилкоксона рекурент

     if(crit==="wilc_Rec_1") {
        wilc_Rec_1(m,n,pw);
        let k=pw.length;
        document.write("Size=",k,"<br>");
        document.write("m=",m,"     n=",n,"<br>");
        for (i=0;i<k;i++){
          wrange[i]=parseFloat(m*(m+1)*0.5+i);
          document.write(wrange[i],"  ",pw[i],"<br>");
        }
     }
     if(crit==="wilc_Rec") {
        wilc_Rec(m,n,wrange,w,pw);
        let k=pw.length;
        document.write("Size=",k,"<br>");
        document.write("m=",m,"     n=",n,"<br>");
        for (i=1;i<k;i++) document.write(w[i],"    ",wrange[i],"   ",pw[i],"<br>");
     }

     if(crit==="wsigne_Rec") {
        wsigne_Rec(nsigne,wrange,w,pw);
        let k=pw.length;
        document.write("Size=",k,"<br>");
        document.write("n=",nsigne,"<br>");
        for (i=1;i<k;i++) document.write(w[i],"    ",wrange[i],"   ",pw[i],"<br>");
     }

     if(crit==="signe_Rec") {
        document.write("n=",nsigne,"<br>");
        signe_Rec(nsigne,pw);
        let k=pw.length;
        document.write("Size=",k,"<br>");
        for (i=0;i<k;i++) document.write(pw[i],"<br>");
     }
    
}
