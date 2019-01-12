#include<bits/stdc++.h>
using namespace std;
const int N=300005;
const double eps=1e-10;
inline int dcmp(const double&x){return x>eps?1:(x<-eps?-1:0);}
struct P{
	double x,y;
	inline void in(){scanf("%lf%lf",&x,&y);}
	inline P(double _x=0,double _y=0):x(_x),y(_y){}
	inline double len(){return sqrt(x*x+y*y);}
	inline P operator*(const double&k)const{return P(x*k,y*k);}
	inline P operator-(const P&rhs)const{return P(x-rhs.x,y-rhs.y);}
	inline P operator+(const P&rhs)const{return P(x+rhs.x,y+rhs.y);}
	inline double operator*(const P&rhs)const{return x*rhs.x+y*rhs.y;}
	inline double operator%(const P&rhs)const{return x*rhs.y-rhs.x*y;}
}a[N];
vector<P>v[N];
struct seg{
	P x,y;
	double l,d,ll,dd,be,en,k;
	int key;
	inline seg(){}
	inline seg(P a,P b):x(a),y(b){
		l=(a-b).len();dd=dcmp(l)>0?fabs(a%b/l):a.len();ll=dcmp(l)==0?0:a*(a-b)/l;
		if(dcmp(ll)<=0)d=a.len(),key=1;
			else if(dcmp(ll-l)>=0)d=b.len(),key=-1;
					else d=dd,key=0;
	}
	inline pair<double,double>jiao(const double&r){
		if(dcmp(r-d)==-1)return make_pair(-1,-1);
		if(k==0)return make_pair(be,en);
		double x=sqrt(r*r-dd*dd);
		if(key==-1)return make_pair(be+max(0.0,ll-x)*k,be+l*k);
		if(key==1)return make_pair(be,be+min(l,x+ll)*k);
		return make_pair(be+max(0.0,ll-x)*k,be+(ll+min(l-ll,x))*k);
	}
};
vector<seg>u[N];
int n,m,i,j,k,xb;
double l,r,mid,c;
pair<double,int>b[N];
int buc[N];
int main(){
	scanf("%d%d",&n,&m);
	for(i=1;i<=n;++i)a[i].in();a[n+1]=a[1];
	for(i=1;i<=n;++i)c+=(a[i+1]-a[i]).len();c/=m;
	if(m>n){printf("%.8f\n",c);return 0;}
	l=0;v[xb=1].push_back(a[1]);
	for(i=1;i<=n;++i){
		double d=(a[i+1]-a[i]).len(),td=d;
		while(dcmp(l+td-c)>=0){
			td-=c-l;l=0;
			v[xb].push_back((a[i+1]-a[i])*(1-td/d)+a[i]);
			v[++xb].push_back((a[i+1]-a[i])*(1-td/d)+a[i]);
		}
		l+=td;
		if((a[i+1]-v[xb].back()).len()>1e-7)v[xb].push_back(a[i+1]);
	}
	v[m+1]=v[1];
	for(i=1;i<=m;++i){
		double be1=0,be2=0,l1,l2,dt,totl=0;
		j=k=0;
		for(;j+1<v[i].size() && k+1<v[i+1].size();){
			l1=(v[i][j+1]-v[i][j]).len();
			l2=(v[i+1][k+1]-v[i+1][k]).len();
			dt=min(l1-be1,l2-be2);
			P a=v[i][j],b=v[i][j+1],c=v[i+1][k],d=v[i+1][k+1],e,f,g,h;
			e=a+(b-a)*(be1/l1);f=a+(b-a)*((be1+dt)/l1);g=c+(d-c)*(be2/l2);h=c+(d-c)*((be2+dt)/l2);
			u[i].push_back(seg(g-e,h-f));u[i].back().be=be1+totl;u[i].back().en=be1+totl+dt;
			u[i].back().k=u[i].back().l>eps?dt/u[i].back().l:0;
			be1+=dt;be2+=dt;
			if(dcmp(be1-l1)==0)be1=0,++j,totl+=l1;
			if(dcmp(be2-l2)==0)be2=0,++k;
		}
	}
	for(l=0,r=c;r-l>1e-8;){
		mid=(l+r)/2;xb=0;
		for(i=1;i<=m;++i)
			for(j=0;j<u[i].size();++j){
				pair<double,double>z=u[i][j].jiao(mid);
				if(z.first!=-1){
					b[++xb]=make_pair(z.first-eps,i);
					b[++xb]=make_pair(z.second+eps,-i);
				}
			}
		sort(b+1,b+xb+1);memset(buc+1,0,m<<2);
		int su=0;
		for(i=1;i<=xb && su<m;)j=abs(b[i].second),su-=buc[j]>0,buc[j]+=j==b[i].second?1:-1,su+=buc[j]>0,++i;
		if(su==m)r=mid;else l=mid;
	}
	printf("%.8f\n",l);
	return 0;
}
