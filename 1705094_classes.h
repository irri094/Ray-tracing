#ifndef CLASSSES_H
#define CLASSSES_H

struct PointLight;
struct SpotLight;
struct Object;
struct Point;
struct Triangle;

extern vector<Object*> objects;
extern vector<SpotLight> spotLights;
extern vector<PointLight> pointLights;
extern Point bgcolor;


const double PI = acos(-1), EPS = 1e-9;
int dcmp(double x){return abs(x)<EPS?0:(x<0?-1:1);}
struct Point {
  double x, y, z;
  Point() : x(0), y(0), z(0) {}
  Point(double X, double Y, double Z) :
      x(X), y(Y), z(Z) {}
  Point operator + (const Point& u) const {
    return Point(x + u.x, y + u.y, z + u.z); }
  Point operator - (const Point& u) const {
    return Point(x - u.x, y - u.y, z - u.z); }
  Point operator * (const double u) const {
    return Point(x * u, y * u, z * u); }
  Point operator * (const Point &u) const {
    return Point(x * u.x, y * u.y, z * u.z); }
  Point operator / (const double u) const {
    return Point(x / u, y / u, z / u); }
  Point operator- () {
    return Point(-x, -y, -z);
  }
  friend std::ostream &operator << (
            std::ostream &os, const Point &p) {
    return os << p.x << " " << p.y <<" "<<p.z; }
  friend std::istream &operator >>
            (std::istream &is, Point &p) {
    return is >> p.x >> p.y >> p.z; }
};
double dot(Point a, Point b) {
  return a.x * b.x + a.y * b.y + a.z * b.z; }
Point cross(Point a, Point b) {
  return Point(a.y*b.z - a.z*b.y,
     a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }
double lenn(Point a) { return sqrt(dot(a, a));}
void colnormalize(Point &a) {
    if(a.x>1) a.x = 1;
    if(a.y>1) a.y = 1;
    if(a.z>1) a.z = 1;
    if(a.x<0) a.x = 0;
    if(a.y<0) a.y = 0;
    if(a.z<0) a.z = 0;
}
double distance(Point a, Point b) {
  return lenn(a-b); }
Point unit(const Point &p) { return p/lenn(p); }
// Rotate p around axis x, with angle radians.
Point rotate(Point p, Point axis, double angle) {
  axis = unit(axis);Point comp1 = p * cos(angle);
  Point comp2 = axis*(1-cos(angle))*dot(axis,p);
  Point comp3 = cross(axis, p) * sin(angle);
  return comp1 + comp2 + comp3;
}
Point rotatepi(Point p, Point axis) {
  axis = unit(axis);
  Point comp2 = axis*2*dot(axis,p);
  return comp2 - p;
}

struct Ray{
	Point start, dir;
    Ray (Point s, Point d){
        start = s, dir = d;
    }
};

void smallcube(Point cen, Point col){
    glColor3f(col.x, col.y, col.z);
    double d = 0.5;
    glBegin(GL_QUADS);{
		glVertex3f(cen.x-d, cen.y+d, cen.z+d); glVertex3f(cen.x-d, cen.y-d, cen.z+d); glVertex3f(cen.x+d, cen.y-d, cen.z+d); glVertex3f(cen.x+d, cen.y+d, cen.z+d);
		glVertex3f(cen.x-d, cen.y+d, cen.z-d); glVertex3f(cen.x-d, cen.y-d, cen.z-d); glVertex3f(cen.x+d, cen.y-d, cen.z-d); glVertex3f(cen.x+d, cen.y+d, cen.z-d);
		glVertex3f(cen.x-d, cen.y+d, cen.z-d); glVertex3f(cen.x-d, cen.y+d, cen.z+d); glVertex3f(cen.x-d, cen.y-d, cen.z+d); glVertex3f(cen.x-d, cen.y-d, cen.z-d);
		glVertex3f(cen.x+d, cen.y+d, cen.z-d); glVertex3f(cen.x+d, cen.y+d, cen.z+d); glVertex3f(cen.x+d, cen.y-d, cen.z+d); glVertex3f(cen.x+d, cen.y-d, cen.z-d);
		glVertex3f(cen.x-d, cen.y+d, cen.z-d); glVertex3f(cen.x-d, cen.y+d, cen.z+d); glVertex3f(cen.x+d, cen.y+d, cen.z+d); glVertex3f(cen.x+d, cen.y+d, cen.z-d);
		glVertex3f(cen.x-d, cen.y-d, cen.z-d); glVertex3f(cen.x-d, cen.y-d, cen.z+d); glVertex3f(cen.x+d, cen.y-d, cen.z+d); glVertex3f(cen.x+d, cen.y-d, cen.z-d);
	}
	glEnd();
}

struct PointLight{
	Point light_pos;
	Point color;
    friend istream &operator >> (istream &is, PointLight &pl) {
        return is >> pl.light_pos >> pl.color;
    }
    void draw(){
        smallcube(light_pos, color);
    }
};

struct SpotLight{
    PointLight point_light;
	Point light_direction;
	double cutoff_angle;
    friend istream &operator >> (istream &is, SpotLight &sl) {
        is >> sl.point_light >> sl.light_direction >> sl.cutoff_angle;
        sl.cutoff_angle *= PI/180.0;
        return is;
    }
    void draw(){
        smallcube(point_light.light_pos, point_light.color);
    }
};

struct Object{
    Point reference_point;	// should have x, y, z
    double height, width, length;
    Point color ;
    virtual Point getColorAt(Point intersecpoint){
        return color;
    }
    array <double, 4> coEfficients;	// ambient, diffuse, specular, reflection coefficients
    int shine;	// exponent term of specular component
    Object() {}
    virtual void draw(){}
    void setColor() {}
    void setShine() {}
    void setCoEfficients() {}
    virtual double intersect(Ray &r, Point &color, int level){
        return -1;
    }
    virtual bool oneye(Point pt){
        return false;
    }

    Point lightworks(Point &objcolor, Point &intersecpoint, Point &normalvec, Ray &r){

        Point color = objcolor;
        color = color * this->coEfficients[0];

        for(PointLight &pl : pointLights){
            Ray rey(pl.light_pos, unit(intersecpoint-pl.light_pos));
            if( dcmp(dot(rey.dir, normalvec)) > 0 ) continue;
            double kasert = inf;
            for (Object *o : objects){
                Point dummyColor;
                double t = o->intersect(rey, dummyColor, -1);
//                    cout<<"ekta intersection done "<<endl;
                if(t < 0 or t > kasert) continue;
                kasert = t;
            }
            if(dcmp(inf-kasert)==0) continue;
            if( dcmp( lenn(pl.light_pos-intersecpoint ) - kasert ) !=0 ) continue;
            double costheta = -dot(rey.dir, normalvec);
            if(costheta > 0)
                color = color + pl.color * this->coEfficients[1] * costheta * objcolor;
//                cout<<"specular shit done "<<endl;
            Point reflecc = rotatepi(-rey.dir , normalvec);
            double cosphi = dot(reflecc, -r.dir);
            if(cosphi > 0)
                color = color + pl.color * this->coEfficients[2] * pow(cosphi, this->shine) * objcolor;
        }

        for(SpotLight &sl : spotLights){
            Ray rey(sl.point_light.light_pos, unit(intersecpoint-sl.point_light.light_pos));
            if( dcmp(dot(rey.dir, normalvec)) > 0 ) continue;
            double cosbeta =   dot(rey.dir, sl.light_direction);
            if( cos( sl.cutoff_angle ) >= cosbeta ) continue;

            double kasert = inf;
            for (Object *o : objects){
                Point dummyColor;
                double t = o->intersect(rey, dummyColor, -1);
                if(t < 0 or t > kasert) continue;
                kasert = t;
            }
            if(dcmp(inf-kasert)==0) continue;
            if( dcmp( lenn(sl.point_light.light_pos-intersecpoint ) - kasert ) !=0 ) continue;
            double costheta = -dot(rey.dir, normalvec);
            if(costheta > 0)
                color = color + sl.point_light.color * this->coEfficients[1] * costheta * objcolor;
            Point reflecc = rotatepi(-rey.dir , normalvec);
            double cosphi = dot(reflecc, -r.dir);
            if(cosphi > 0)
                color = color + sl.point_light.color * this->coEfficients[2] * pow(cosphi, this->shine) * objcolor;
        }

        return color;
    }

    void reflection(int level, Ray &r, Point &normalvec, Point &intersecpoint, Point &color){
        if(level<=0) return;
        Point reflecc = rotatepi(-r.dir , normalvec);
        Ray bounced(intersecpoint+reflecc*EPS*1000, reflecc);
        double kasert = inf;
        Point rong(bgcolor.x,bgcolor.y,bgcolor.z);
        for (Object *o : objects){
            Point dummyColor;
            double t = o->intersect(bounced, dummyColor, level-1);
            if(t < 0 or t > kasert) continue;
            kasert = t, rong = dummyColor;
        }
        color = color + rong * this->coEfficients[3];
    }

};

struct Sphere : Object{
    Sphere(){}
	Sphere(Point center, double radius){
		reference_point = center;
		length = radius;
	}
    friend istream &operator >> (istream &is, Sphere &S) {
        is >> S.reference_point >> S.length ;
        is >> S.color;
        for(int j=0; j<4; j++) is >> S.coEfficients[j];
        is >> S.shine;
        return is;
    }

    double intersect(Ray &r, Point &color, int level){
        double t = -1;
        Point Ro = r.start - this->reference_point;
        double Ro2 = dot(Ro, Ro);
        double tp = - dot(Ro, r.dir);
        double d2 = Ro2 - tp*tp;
        double r2 = this->length * this->length;
        if(d2 > r2 or (Ro2 >= r2 and tp<0) ) return t;

        double tprime = sqrtl(r2 - d2);
        if(Ro2 > r2) t = tp - tprime;
        else t = tp + tprime;

        if(t<0) return t;

        Point intersecpoint = r.start + r.dir * t;
        Point normalvec = unit(intersecpoint - this->reference_point);
        Point objcolor = this->color;

        if(level!=-1)
            color = lightworks(objcolor, intersecpoint, normalvec, r);

        reflection(level, r, normalvec, intersecpoint, color);

        colnormalize(color);

        return t;
    }
    void draw(){
        vector< vector<Point> >points(stacknum+1, vector<Point>(slicenum+1));
        for(int i=0;i<=stacknum;i++){
            double h=length*sin(((double)i/(double)stacknum)*(PI/2));
            double r=length*cos(((double)i/(double)stacknum)*(PI/2));
            for(int j=0;j<=slicenum;j++){
                points[i][j].x=r*cos(((double)j/(double)slicenum)*2*PI);
                points[i][j].y=r*sin(((double)j/(double)slicenum)*2*PI);
                points[i][j].z=h;
            }
        }
        glPushMatrix();
        glTranslatef(reference_point.x, reference_point.y, reference_point.z);
        Point c = this->color;
        glColor3f(c.x, c.y, c.z);
        //draw quads using generated points
        for(int i=0;i<stacknum;i++){
            for(int j=0;j<slicenum;j++){
                glBegin(GL_QUADS);{
                    //upper hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }glEnd();
            }
        }
        glPopMatrix();
    }
    bool oneye(Point pt){
        return bool( dcmp( lenn(pt-reference_point) - length ) == 0 );
    }
};

struct Floor: Object {
    Floor(int floorWidth, int tileWidth) {
        reference_point = Point(-floorWidth/2,-floorWidth/2,0);
        length=tileWidth;
    }
    Point getColorAt(Point intersecpoint){
        int i = (intersecpoint.x + reference_point.x + EPS) / length;
        int j = (intersecpoint.y + reference_point.y + EPS) / length;
        if((i+j)%2) return Point(1,1,1);
        else return Point(0,0,0);
    }

    double intersect(Ray &r, Point &color, int level){
        double t = -r.start.z/r.dir.z;
        Point intersecpoint = r.start + r.dir * t;
        if(t<0) return t;
        if( intersecpoint.x > -reference_point.x or intersecpoint.x < reference_point.x or intersecpoint.y > -reference_point.y or intersecpoint.y < reference_point.y)
            return -1;

        Point normalvec(0,0,1);
        if(r.start.z<0) normalvec = -normalvec;
        Point objcolor = getColorAt(r.start + r.dir * t);
        if(level!=-1)
            color = lightworks(objcolor, intersecpoint, normalvec, r);

        reflection(level, r, normalvec, intersecpoint, color);

        colnormalize(color);

        return t;
    }
    void draw(){
        int tilecount = int((-reference_point.x+EPS)/length)*2;
        for(int i=0; i<tilecount; i++){
            for(int j=0; j<tilecount; j++){
                if( (i+j)%2 )
                    glColor3f(1,1,1);
                else
                    glColor3f(0,0,0);

                glBegin(GL_QUADS);
                glVertex3f(reference_point.x + i*length, reference_point.y + j*length, 0);
                glVertex3f(reference_point.x + (i+1)*length, reference_point.y + j*length, 0);
                glVertex3f(reference_point.x + (i+1)*length, reference_point.y + (j+1)*length, 0);
                glVertex3f(reference_point.x + i*length, reference_point.y + (j+1)*length, 0);
                glEnd();
            }
        }
    }
    bool oneye(Point pt){
        if(dcmp(pt.z)==0){
            if( pt.x > -reference_point.x or pt.x < reference_point.x or pt.y > -reference_point.y or pt.y < reference_point.y)
                return false;
            return true;
        }
        return false;
    }
};

/// is p1 and p2 in same side of ab
bool sameSide(Point &p1, Point &p2, Point a, Point &b){
    Point c1 = cross(b-a, p1-a), c2 = cross(b-a, p2-a);
    if( dcmp(dot(c1, c2)) != -1) return true;
    return false;
}

bool inTriangle(Point &a, Point &b, Point &c, Point &p){
    if(sameSide(p, a, b, c) and sameSide(p, b, a, c) and sameSide(p, c, b, a) )
        return true;
    return false;
}

struct Triangle: Object {
    Point P[3];
    Triangle(){}
    Triangle(Point p1, Point p2, Point p3){
        P[0] = p1, P[1] = p2, P[2] = p3;
    }
    friend istream &operator >> (istream &is, Triangle &T) {
        is >> T.P[0] >> T.P[1] >> T.P[2];
        is >> T.color;
        for(int j=0; j<4; j++) is >> T.coEfficients[j];
        is >> T.shine;
        return is;
    }
    double intersect(Ray &r, Point &color, int level){
        double t = -1;
        Point normalvec = unit( cross(P[2]-P[0], P[1]-P[0]) );
        int tmp = dcmp( dot(normalvec, r.dir) ) ;
        if(!tmp) return t;
        if(tmp>0) normalvec = -normalvec;

        t = dot(normalvec, P[0]-r.start) / dot(normalvec, r.dir);
        Point intersecpoint = r.start + r.dir * t;

        if(!inTriangle(P[0], P[1], P[2], intersecpoint)) return -1;

        Point objcolor = this->color;
        if(level!=-1)
            color = lightworks(objcolor, intersecpoint, normalvec, r);

        reflection(level, r, normalvec, intersecpoint, color);

        colnormalize(color);

        return t;
    }
    void draw(){
        Point c = this->color;
        glColor3f(c.x, c.y, c.z);

        glBegin(GL_TRIANGLES);
        glVertex3f(P[0].x, P[0].y, P[0].z);
        glVertex3f(P[1].x, P[1].y, P[1].z);
        glVertex3f(P[2].x, P[2].y, P[2].z);
        glEnd();
    }
    bool onplane(Point pt){
        return bool(dcmp( dot( cross(P[2]-P[0], P[1]-P[0]), pt-P[0] ) )==0);
    }
    bool oneye(Point pt){
        if(onplane(pt) and inTriangle(P[0], P[1], P[2], pt))
            return true;
        return false;
    }
};

bool validate(Point &p, Point &rf, double dx, double dy, double dz){
    if(dcmp(dx)>0){
        if(dcmp(p.x-rf.x)<0 or dcmp(p.x-rf.x-dx)>0) return false;
    }
    if(dcmp(dy)>0){
        if(dcmp(p.y-rf.y)<0 or dcmp(p.y-rf.y-dy)>0) return false;
    }
    if(dcmp(dz)>0){
        if(dcmp(p.z-rf.z)<0 or dcmp(p.z-rf.z-dz)>0) return false;
    }
    return true;
}

struct Surf: Object {
    double A,B,C,D,E,F,G,H,I,J;
    Surf(){}
    Surf(double A1, double B1, double C1, double D1, double E1, double F1, double G1, double H1, double I1, double J1) {
        A=A1, B=B1, C=C1, D=D1, E=E1, F=F1, G=G1, H=H1, I=I1, J=J1;
    }
    friend istream &operator >> (istream &is, Surf &S) {
        is>>S.A>>S.B>>S.C>>S.D>>S.E>>S.F>>S.G>>S.H>>S.I>>S.J;
        is >> S.reference_point >> S.length >> S.width >> S.height;
        is >> S.color;
        for(int j=0; j<4; j++) is >> S.coEfficients[j];
        is >> S.shine;
        return is;
    }
    double intersect(Ray &r, Point &color, int level){
        double t = -1;

        double sx = r.start.x, sy = r.start.y, sz = r.start.z;
        double dx = r.dir.x, dy = r.dir.y, dz = r.dir.z;

        double a = A*dx*dx + B*dy*dy + C*dz*dz + D*dx*dy + E*dx*dz + F*dy*dz;
        double b = 2*A*sx*dx + 2*B*sy*dy + 2*C*sz*dz + D*(sx*dy + sy*dx) + E*(sz*dx + sx*dz) + F*(sy*dz + sz*dy) + G*dx + H*dy + I*dz;
        double c = A*sx*sx + B*sy*sy + C*sz*sz + D*sx*sy + E*sx*sz + F*sy*sz + G*sx + H*sy + I*sz + J;

        double Det = b*b-4*a*c;
        if(Det<0) return t;

        double tbig = (-b+sqrt(Det))/(2*a), tsmol = (-b-sqrt(Det))/(2*a);

        Point farpoint = r.start + r.dir*tbig;
        Point clospoint = r.start + r.dir*tsmol;

//        cout<<farpoint<<"\n"<<clospoint<<endl;

        bool farvisible = (dcmp(tbig)>0 and validate(farpoint, reference_point, length, width, height) );
        bool closvisible = (dcmp(tsmol)>0 and validate(clospoint, reference_point, length, width, height) );

//        cout<<farvisible<<" "<<closvisible<<endl;

        Point intersecpoint;
        if(!farvisible and !closvisible) return t;
        if(closvisible) intersecpoint = clospoint, t = tsmol;
        else intersecpoint = farpoint, t = tbig;

//        cout<<intersecpoint<<endl;

        Point normalvec = Point(
            2*A*intersecpoint.x + D*intersecpoint.y + E*intersecpoint.z + G,
            2*B*intersecpoint.y + D*intersecpoint.x + F*intersecpoint.z + H,
            2*C*intersecpoint.z + E*intersecpoint.x + F*intersecpoint.y + I
        );
        normalvec = unit(normalvec);
        if(dcmp(dot(normalvec, r.dir)) > 0) normalvec = -normalvec;

        Point objcolor = this->color;
        if(level!=-1)
            color = lightworks(objcolor, intersecpoint, normalvec, r);

        reflection(level, r, normalvec, intersecpoint, color);

        colnormalize(color);

        return t;
    }
    void draw(){}
    bool oneye(Point p){
        double val = A*p.x*p.x + B*p.y*p.y + C*p.z*p.z + D*p.x*p.y + E*p.x*p.z + F*p.z*p.y + G*p.x + H*p.y + I*p.z + J;
        return bool(dcmp(val)==0);
    }
};

#endif // CLASSSES_H
