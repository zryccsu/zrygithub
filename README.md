#include<bits/stdc++.h>  
using namespace std;  
const double PI=acos(-1);  
struct Point3{//定义三维坐标点   
	double x,y,z;  
	Point3(double x=0,double y=0,double z=0):x(x),y(y),z(z){}  
};  
typedef Point3 Vector3;   
Vector3 operator +(Vector3 A,Vector3 B){//向量+向量=向量 点+向量=点   
	return Vector3(A.x+B.x,A.y+B.y,A.z+B.z);  
}  
Vector3 operator -(Point3 A,Point3 B){//点-点=向量   
	return Vector3(A.x-B.x,A.y-B.y,A.z-B.z);  
}  
Vector3 operator *(Vector3 A,double p){//向量*数=向量   
	return Vector3(A.x*p,A.y*p,A.z*p);  
}  
Vector3 operator /(Vector3 A,double p){//向量/数=向量   
	return Vector3(A.x/p,A.y/p,A.z/p);  
}  
const double eps=1e-10;  
int dcmp(double x){  
	if(fabs(x)<eps)return 0;  
	return x<0?-1:1;	  
}  
bool operator ==(const Point3& a,const Point3& b){  
	return dcmp(a.x-b.x)==0&&dcmp(a.y-b.y)==0&&dcmp(a.z-b.z)==0;  
}  
double Dot(Vector3 A,Vector3 B){//三维点乘   
	return A.x*B.x+A.y*B.y+A.z*B.z;  
}  
double Length(Vector3 A){//向量的模   
	return sqrt(Dot(A,A));  
}  
double Angle(Vector3 A,Vector3 B){//向量夹角   
	return acos(Dot(A,B)/Length(A)/Length(B));  
}  
double DistanceToPlane(const Point3& p,const Point3& p0,const Vector3& n){//点p到平面距离p0为平面一点n为单位法向量   
	return fabs(Dot(p-p0,n));  
}  
Point3 GetPlaneProjection(const Point3& p,const Point3& p0,const Vector3& n){//点p在p0-n上的投影   
	return p-n*Dot(p-p0,n);  
}  
Point3 LinePlaneIntersection(Point3 p1,Point3 p2,Point3 p0,Vector3 n){//直线p1-p2与面p0交点   
	Vector3 v=p2-p1;  
	double t=(Dot(n,p0-p1)/Dot(n,p2-p1));  
	return p1+v*t;  
}  
Vector3 Cross(Vector3 A,Vector3 B){//三维向量叉乘   
	return Vector3(A.y*B.z-A.z*B.y,A.z*B.x-A.x*B.z,A.x*B.y-A.y*B.x);  
}  
double Area2(Point3 A,Point3 B,Point3 C){//三角形ABC面积   
	return Length(Cross(B-A,C-A));  
}  
bool PointInTri(Point3 P,Point3 P0,Point3 P1,Point3 P2){//点P是否在三角形P0P1P2中   
	double area1=Area2(P,P0,P1);  
	double area2=Area2(P,P1,P2);  
	double area3=Area2(P,P2,P0);  
	return dcmp(area1+area2+area3-Area2(P0,P1,P2))==0;  
}  
bool TriSegIntersection(Point3 P0,Point3 P1,Point3 P2,Point3 A,Point3 B,Point3& P){//线段AB与三角形P0P1P2是否相交   
	Vector3 n=Cross(P1-P0,P2-P0);  
	if(dcmp(Dot(n,B-A))==0)return false;//线段AB与平面平行或在平面内  
	else{  
		double t=Dot(n,P0-A)/Dot(n,B-A);  
		if(dcmp(t)<0||dcmp(t-1)>0)return false;//交点不在线段AB上   
		P=A+(B-A)*t;  
		return PointInTri(P,P0,P1,P2);  
	}  
}  
double DistTanceceToLine(Point3 P,Point3 A,Point3 B){//点P到直线AB的距离   
	Vector3 v1=B-A,v2=P-A;  
	return Length(Cross(v1,v2))/Length(v1);  
}  
double DistanceToSegment(Point3 P,Point3 A,Point3 B){//点P到线段AB的距离   
	if(A==B)return Length(B-A);  
	Vector3 v1=B-A,v2=P-A,v3=P-B;  
	if(dcmp(Dot(v1,v2))<0)return Length(v2);  
	else if(dcmp(Dot(v1,v3)>0))return Length(v3);  
	else return Length(Cross(v1,v2))/Length(v1);  
}  
double Volume6(Point3 A,Point3 B,Point3 C,Point3 D){//6*V=S*H=(AB×AC)·AD   
	return Dot(D-A,Cross(B-A,C-A));  
}  
struct Face{  
	int v[3];  
	Vector3 normal(Point3 *P)const{  
		return Cross(P[v[1]]-P[v[0]],P[v[2]]-P[v[0]]);  
	}  
	int cansee(Point3 *P,int i)const{  
		return Dot(P[i]-P[v[0]],normal(P))>0?1:0;  
	}  
};  
int vis[505][505];  
vector<Face> CH3D(Point3 *P,int n){//增量法求三维凸包   
	vector<Face>cur;  
	//前三点不共线   
	cur.push_back((Face){{0,1,2}});  
	cur.push_back((Face){{2,1,0}});  
	for(int i=3;i<n;i++){  
		vector<Face>next;  
		//计算每条边"左边"的可能性   
		for(int j=0;j<cur.size();j++){  
			Face& f=cur[j];  
			int res=f.cansee(P,i);  
			if(!res)next.push_back(f);  
			for(int k=0;k<3;k++)vis[f.v[k]][f.v[(k+1)%3]]=res;  
		}  
		for(int j=0;j<cur.size();j++)  
			for(int k=0;k<3;k++){  
				int a=cur[j].v[k],b=cur[j].v[(k+1)%3];  
				if(vis[a][b]!=vis[b][a]&&vis[a][b])//(a,b)是分界线,左边对P[i]可见  
					next.push_back((Face){{a,b,i}});   
			}  
		cur=next;  
	}  
	return cur;  
}  
int main(){  
	memset(vis,0,sizeof(vis));  
	return 0;  
}   
