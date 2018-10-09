#include<bits/stdc++.h>
using namespace std;
const double PI=acos(-1);
const int INF=1<<20;
struct Point{//定义点 
	double x,y;
	Point (double x=0,double y=0):x(x),y(y){}
}p[10000];
typedef Point Vector;//定义向量 
Vector operator +(Vector A,Vector B){//点+向量=点  向量+向量=向量 
	return Vector(A.x+B.x,A.y+B.y);
}
Vector operator -(Point A,Point B){//点-点=向量 
	return Vector(A.x-B.x,A.y-B.y);
}
Vector operator *(Vector A,double p){//向量*数=向量 
	return Vector(A.x*p,A.y*p);
}
Vector operator /(Vector A,double p){//向量/数=向量 
	return Vector(A.x/p,A.y/p);
}
bool operator <(const Point& a,const Point& b){//判断位置 
	if(a.x==b.x)return a.y<b.y;
	return a.x<b.x;
}
const double eps=1e-10;
int dcmp(double x){//判断x正负性 
	if(fabs(x)<eps)return 0;
	else return x<0?-1:1;
}
bool operator ==(const Point& a,const Point& b){//点位置是否相等 
	return dcmp(a.x-b.x)==0&&dcmp(a.y-b.y)==0;
}
double Dot(Vector A,Vector B){//向量点乘 
	return A.x*B.x+A.y*B.y;
}
double Length(Vector A){//向量的模 
	return sqrt(Dot(A,A));
}
double Angel(Vector A,Vector B){//向量夹角 
	return acos(Dot(A,B)/Length(A)/Length(B));
}
double Cross(Vector A,Vector B){//向量叉乘 
	return A.x*B.y-A.y*B.x;
}
double Area2(Point A,Point B,Point C){//向量叉乘2 
	return Cross(B-A,C-A);
}
Vector Rotate(Vector A,double rad){//向量逆时针旋转rad弧度 
	return Vector(A.x*cos(rad)-A.y*sin(rad),A.x*sin(rad)+A.y*cos(rad));
}
Vector Normal(Vector A){//非零向量A的法向量 
	double L=Length(A);
	return Vector(-A.y/L,A.x/L);
}
Point GetLineIntersection(Point A,Vector v,Point B,Vector w){//直线交点 向量点版 
	Vector u=A-B;
	double t=Cross(w,u)/Cross(v,w);
	return A+v*t;
}
double DistanceToLine(Point P,Point A,Point B){//P点到经过点A、B的直线的距离 
	Vector v1=P-A,v2=B-A;
	return fabs(Cross(v1,v2)/Length(v2));
}
double DistanceToSegment(Point P,Point A,Point B){//P点到线段AB的距离 
	if(A==B)return Length(P-A);
	Vector v1=B-A,v2=P-A,v3=P-B;
	if(dcmp(Dot(v1,v2))<0)return Length(v2);
	else if(dcmp(Dot(v1,v3))>0)return Length(v3);
	else return fabs(Cross(v1,v2)/Length(v1));
}
Point GetLineProjection(Point P,Point A,Point B){//点P在直线AB上的投影 
	Vector v=B-A;
	return A+v*(Dot(v,P-A)/Dot(v,v));
}
bool SegmentProperIntersection(Point a1,Point a2,Point b1,Point b2){//判断线段不在端点处是否相交 
	double c1=Cross(a2-a1,b1-a1),c2=Cross(a2-a1,b2-a1),c3=Cross(b2-b1,a1-b1),c4=Cross(b2-b1,a2-b1);
	return dcmp(c1)*dcmp(c2)<0&&dcmp(c3)*dcmp(c4)<0;
}
bool OnSegment(Point p,Point a1,Point a2){//判断点是否在线段上 
	return dcmp(Cross(a1-p,a2-p))==0&&dcmp(Dot(a1-p,a2-p))<0;
}
double ConvexPolygonArea(Point *p,int n){//n突边形面积 
	double area=0;
	for(int i=1;i<n-1;i++){
		area+=Cross(p[i]-p[0],p[i+1]-p[0]);
	}
	return area/2;
}
double PolgonArea(Point* p,int n){//多边形的有向面积 
	double area=0;
	for(int i = 1;i < n-1;i++){
		area+=Cross(p[i]-p[0],p[i+1]-p[0]);
	}
	return area/2;
} 
Point read_point(){//获得点 
	Point a;
	scanf("%lf%lf",&a.x,&a.y);
	return a;
}
struct Circle{//定义圆 
	Point c;
	double r;
	Circle(Point c,double r):c(c),r(r){}
	Point point(double a){//求圆上的点
		return Point(c.x+r*cos(a),c.y+r*sin(a));
	}
};
struct Line{//定义直线 
	Point p;
	Vector v;
	double ang;
	Line (){}
	Line(Point p,Vector v):p(p),v(v){ang=atan2(v.y,v.x);}
	bool operator <(const Line& L)const{
		return ang<L.ang;
	}
	Point point(double a){
		return p+v*a;
	}
};
int getLineCircleIntersection(Line L,Circle C,double& t1,double &t2,vector<Point>& sol){//直线与圆的交点 
	double a=L.v.x,b=L.p.x-C.c.x,c=L.v.y,d=L.p.y-C.c.y;
	double e=a*a+c*c,f=2*(a*b+c*d),g=b*b+d*d-C.r*C.r;
	double delta=f*f-4*e*g;//判别式
	if(dcmp(delta)<0)return 0;//相离 
	if(dcmp(delta)==0){
		t1=t2=-f/(2*e);
		sol.push_back(L.point(t1));
		return 1;//相切 
	}
	t1=(-f-sqrt(delta))/(2*e);
	sol.push_back(L.point(t1));
	t2=(-f+sqrt(delta))/(2*e);
	sol.push_back(L.point(t2));
	return 2;//相交 
}
double angle(Vector v){//计算向量极角 
	return atan2(v.y,v.x);
}
int getCIrcleCirclrIntersection(Circle C1,Circle C2,vector<Point>& sol){//两圆相交 
	double d=Length(C1.c-C2.c);
	if(dcmp(d)==0){//同心圆 
		if(dcmp(C1.r-C2.r)==0)return -1;//两圆重合
		return 0; //相离 
	}
	if(dcmp(C1.r+C2.r-d)<0)return 0;//半径之和小于圆心距  相离 
	if(dcmp(fabs(C1.r-C2.r)-d)>0)return 0; //半径之差大于圆心距 相离 
	double a=angle(C2.c-C1.c);//向量C1C2极角 
	double da=acos((C1.r*C1.r+d*d-C2.r*C2.r))/(2*C1.r*d);//C1C2到C1P1的角
	Point p1=C1.point(a-da),p2=C1.point(a+da);
	sol.push_back(p1);
	if(p1==p2)return 1;//相切 
	sol.push_back(p2);
	return 2;//相交 
} 
int getTangents(Point p,Circle C,Vector* v){//过点p做圆C的切线
	Vector u=C.c-p;
	double dist=Length(u);
	if(dist<C.r)return 0;//点在圆内无法做切线 
	else if(dcmp(dist-C.r)==0){//点在圆上,只能有一条切线
		v[0]=Rotate(u,PI/2); 
		return 1;
	}
	else{//点在圆外,两条切线 
		double ang=asin(C.r/dist);
		v[0]=Rotate(u,ang);
		v[1]=Rotate(u,-ang);
		return 2;
	}
}
int getTangents(Circle A,Circle B,Point *a,Point *b){//两圆公切线 
	int cnt=0;
	if(A.r<B.r){
		swap(A,B);
		swap(a,b);
	}
	int d2=(A.c.x-B.c.x)*(A.c.x-B.c.x)+(A.c.y-B.c.y)*(A.c.y-B.c.y);
	int rdiff=A.r-B.r;
	int rsum=A.r+B.r;
	if(d2<rdiff*rdiff)return 0;//内含
	double base=atan2(B.c.y-A.c.y,B.c.x-A.c.x);
	if(d2==0&&A.r==B.r)return -1;
	if(d2==rdiff*rdiff){//内切 
		a[cnt]=A.point(base);
		b[cnt]=B.point(base);
		cnt++;
		return 1;
	}
	double ang=acos((A.r-B.r)/sqrt(d2));
	
	if(d2=rsum*rsum){//外切 一条内公切线 
		a[cnt]=A.point(base);
		b[cnt]=B.point(base+PI);
		cnt++;
	}
	else if(d2>rsum*rsum){//两条内公切线 
		double ang=acos((A.r+B.r)/sqrt(d2));
		a[cnt]=A.point(base+ang);
		b[cnt]=B.point(base+ang);
		cnt++;
		a[cnt]=A.point(base-ang);
		b[cnt]=B.point(base-ang);
		cnt++;
	}
	return cnt;
}
Circle CircumscribedCircle(Point A,Point B,Point C){//三角形外接圆 
	double bx=B.x-A.x,by=B.y-A.y;
	double cx=C.x-A.x,cy=C.y-A.y;
	double d=2*(bx*cy-by*cx);
	double px=(cy*(bx*bx+by*by)-by*(cx*cx+cy*cy))/d+A.x; 
	double py=(bx*(cx*cx+cy*cy)-cx*(bx*bx+by*by))/d+A.y; 
	Point p=Point(px,py);
	return Circle(p,Length(A-p));
}
Circle InscribedCircle(Point A,Point B,Point C){//三角形内切圆 
	double a=Length(B-C);
	double b=Length(C-A);
	double c=Length(A-B);
	Point p=(A*a+B*b+C*c)/(a+b+c);
	return Circle(p,DistanceToLine(p,A,B));
}
typedef vector<Point>polygon;
int isPointInPolygon(Point p,polygon poly){//点在多边形内判定 
	int wn=0;
	int n=poly.size();
	for(int i=0;i<n;i++){
		if(OnSegment(p,poly[i],poly[(i+1)%n]))return -1;//在边界上 
		int k=dcmp(Cross(poly[(i+1)%n]-poly[i],p-poly[i]));
		int d1=dcmp(poly[i].y-p.y);
		int d2=dcmp(poly[(i+1)%n].y-p.y);
		if(k>0&&d1<=0&&d2>0)wn++;
		if(k<0&&d2<=0&&d1>0)wn--;
	}
	if(wn!=0)return 1;//内部 
	return 0;//外部 
}
int ConvexHull(Point* p,int n,Point* ch){//计算并返回凸包顶点个数 
	sort(p,p+n);
	int m=0;
	for(int i=0;i<n;i++){
		while(m>1&&Cross(ch[m-1]-ch[m-2],p[i]-ch[m-2])<=0)m--;
		ch[m++]=p[i];
	}
	int k=m;
	for(int i=n-2;i>=0;i--){
		while(m>k&&Cross(ch[m-1]-ch[m-2],p[i]-ch[m-2])<=0)m--;
		ch[m++]=p[i];
	}
	if(n>1)m--;
	return m;
}
vector<Point> ConvexHull(vector<Point>& p){//动态数组凸包 
	sort(p.begin(),p.end());
	p.erase(unique(p.begin(),p.end()),p.end());
	int n=p.size();
	int m=0;
	vector<Point>ch(n+1);
	for(int i=0;i<n;i++){
		while(m>1&&Cross(ch[m-1]-ch[m-2],p[i]-ch[m-2])<=0)m--;
		ch[m++]=p[i];
	}
	int k=m;
	for(int i=n-2;i>=0;i--){
		while(m>k&&Cross(ch[m-1]-ch[m-2],p[i]-ch[m-2])<=0)m--;
		ch[m++]=p[i];
	}
	if(n>1)m--;
	ch.resize(m);
	return ch;
}
int getMaxDirmater(vector<Point>&points){//凸包最长点距(平方) 
	vector<Point>p=ConvexHull(points);
	int n=p.size();
	if(n==1)return 0;
	if(n==2)return (int)Dot(p[0]-p[1],p[0]-p[1]);
	p.push_back(p[0]);
	int ans=0;
	for(int i=0,j=1;i<n;i++){
		while(1){
			int diff=(int)Cross(p[i+1]-p[i],p[j+1]-p[j]);
			if(diff<=0){
				ans=max(ans,(int)Dot(p[i]-p[j],p[i]-p[j]));
				if(diff==0)ans=max(ans,(int)Dot(p[i]-p[j+1],p[i]-p[j+1]));
				break;
			}
			j=(j+1)%n;
		}
	}
	return ans;
}
polygon CutPolygon(polygon poly,Point A,Point B){//有向直线A->B切割多边形poly,返回左侧 
	polygon newpoly;
	int n=poly.size();
	for(int i=0;i<n;i++){
		Point C=poly[i];
		Point D=poly[(i+1)%n];
		if(dcmp(Cross(B-A,C-A))>=0)newpoly.push_back(C);
		if(dcmp(Cross(B-A,C-A))!=0){
			Point ip=GetLineIntersection(A,B-A,C,D-C);
			if(OnSegment(ip,C,D))newpoly.push_back(ip);
		}
	}
	return newpoly;
}
bool OnLeft(Line L,Point p){//点P在直线L左边 
	return Cross(L.v,p-L.p)>0;
}
Point GetIntersection(Line a,Line b){//直线交点直线版 
	Vector u=a.p-b.p;
	double t=Cross(b.v,u)/Cross(a.v,b.v);
	return a.p+a.v*t;
}
int HalfplaneIntersection(Line* L,int n,Point* poly){//半面相交 
	sort(L,L+n);//按极角排序 
	int frist,last;//双端队列的第一个元素和最后一个元素的下标 
	Point *p=new Point[n];//p[i]为q[i]与q[i+1]的交点 
	Line *q=new Line[n];//双端队列 
	q[frist=last=0]=L[0];//双端队列初始化为只有一个半平面L[0] 
	for(int i=1;i<n;i++){
		while(frist<last&&!OnLeft(L[i],p[last-1]))last--;
		while(frist<last&&!OnLeft(L[i],p[frist]))frist++;
		q[++last]=L[i];
		if(fabs(Cross(q[last].v,q[last-1].v))<eps){//两向量平行,选内侧 
			last--;
			if(OnLeft(q[last],L[i].p))q[last]=L[i];
		}
		if(frist<last)p[last-1]=GetIntersection(q[last-1],q[last]);
	}
	while(frist<last&&!OnLeft(q[frist],p[last-1]))last--;//删除无用平面
	if(last-frist<=1)return 0;
	p[last]=GetIntersection(q[last],q[frist]); //计算首尾平面交点
	//从deque复制到输出中
	int m=0;
	for(int i=frist;i<=last;i++)poly[m++]=p[i];
	return m;
}
bool cmp(int& a,int& b){
	return p[a].y<p[b].y;
}
int temp[10000];
double merge(int l,int r){//平面最近点 s排序 
	double d=INF;
	if(l==r)return d;
	if(l+1==r)return Length(p[l]-p[r]);
	int m=(l+r)/2;
	double d1=merge(l,m);
	double d2=merge(m+1,r);
	d=min(d1,d2);
	int i,j,k=0;
	for(i=l;i<=r;i++)if(fabs(p[m].x-p[i].x)<d)temp[k++]=i;
	sort(temp,temp+k,cmp);
	for(i=0;i<k;i++)
		for(j=i+1;j<k&&p[temp[j]].y-p[temp[i]].y<d;j++){
			double d3=Length(p[temp[i]]-p[temp[j]]);
			d=min(d,d3);
		}
	return d;
}
bool cmp2(Point a,Point b){
	return a.y<b.y;
}
vector<Point>q;
double merge2(int l,int r){//平面最近点 p排序 
	double d=INF;
	if(l==r)return d;
	if(l==r-1)return Length(p[r]-p[l]);
	int m=(l+r)/2;
	double d1=merge(l,m);
	double d2=merge(m+1,r);
	d=min(d1,d2);
	q.clear();
	for(int i=l;i<=r;i++)if(fabs(p[i].x-p[m].x)<d)q.push_back(p[i]);
	int n=q.size();
	sort(q.begin(),q.end(),cmp2);
	for(int i=0;i<n-1;i++)
		for(int j=i+1;j<n&&q[j].y-q[i].y<d;j++)d=min(d,Length(q[i]-q[j]));
	return d;
}
int main(){
	return 0;
}
