#include <GL/glut.h>
#include<bits/stdc++.h>
using namespace std;

#include "bitmap_image.hpp"

const double inf = 1e9;

int width = 600, height = 600;
int imageWidth = 600, imageHeight = 600;
int imgcount = 11;
int reclevel = 0;
const double viewAngle = 80;
const int slicenum = 90, stacknum = 90;

#include "1705094_classes.h"

Point bgcolor;
Point campos, camu, caml, camr;
vector <Object*> objects;
vector <PointLight> pointLights;
vector <SpotLight> spotLights;

void capture(){
    bitmap_image image(imageWidth,imageHeight);
    for(int i=0; i<imageHeight; i++) for(int j=0; j<imageWidth; j++)
        image.set_pixel(j, i, 0, 0, 0);

    double planeDistance  = (height/2.0) / tan(PI*viewAngle/180.0/2.0);
    Point topleft = campos + caml*planeDistance - camr*width/2 + camu*height/2;
    double du = 1.0*width/imageWidth;
    double dv = 1.0*height/imageHeight;
    topleft = topleft + camr*(0.5*du) - camu*(0.5*dv);

    cout<<"capture begin"<<endl;

    /// check for direct intersection
    bool ondho = false;

    for (Object *o : objects){
        if(  o->oneye(campos)  ) {
            ondho = true;
            Point rong = o->getColorAt(campos);
            rong = rong * 255;
            for(int j=0; j<imageWidth; j++) for(int i=0; i<imageHeight; i++)
                image.set_pixel(j, i, int(rong.x), int(rong.y), int(rong.z));
            break;
        }
    }

    if(!ondho){
        for (int j=0; j<imageWidth; j++) {
            for(int i=0; i<imageHeight; i++){
                Point curPixel = topleft + camr*j*du + camu*-i*dv;
                Ray rey(campos, unit(curPixel - campos) );
                double kasert = inf;
                Point rong(bgcolor.x,bgcolor.y,bgcolor.z);
                for (Object *o : objects){
                    Point dummyColor;
                    double t = o->intersect(rey, dummyColor, reclevel);
                    if(t < 0 or t > kasert) continue;
                    kasert = t, rong = dummyColor;
                }
                rong = rong * 255;
                image.set_pixel(j, i, int(rong.x), int(rong.y), int(rong.z));
            }
        }
    }
    cout<<"capture ended"<<endl;

    string cntstr = to_string(imgcount++);
    image.save_image("output_"+cntstr+".bmp");
}

void quit(){
    glutDestroyWindow(glutGetWindow());
    for(Object *ob : objects)
        delete ob;
    objects.clear();
    pointLights.clear();
    spotLights.clear();
    cout<<"memory cleared"<<endl;
}

/// this is in radian
const double rot_angle = 0.05;

/// rotate a and b about c
void lookaround(Point &a, Point &b, Point &c, bool klokwise){
    double theta = rot_angle;
    if(klokwise) theta *= -1;
    a = rotate(a, c, theta);
    b = rotate(b, c, theta);
    a = unit(a), b = unit(b);
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
		case '1':
		    lookaround(caml, camr, camu, false);
			break;
		case '2':
		    lookaround(caml, camr, camu, true);
			break;
		case '3':
		    lookaround(camu, caml, camr, false);
			break;
		case '4':
		    lookaround(camu, caml, camr, true);
			break;
		case '5':
		    lookaround(camu, camr, caml, false);
			break;
		case '6':
		    lookaround(camu, camr, caml, true);
			break;
		case '0':
		    capture();
			break;
        case 'q':
            quit();
            break;
        case 'k':
            cout<<"level komano "<<endl;
            reclevel = max(reclevel-1, 0);
            break;
        case 'l':
            cout<<"level barano "<<endl;
            reclevel++;
            break;
		default:
			break;
	}
}

const double muvlen = 3;
void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:
            campos = campos - caml*muvlen;
			break;
		case GLUT_KEY_UP:
		    campos = campos + caml*muvlen;
			break;
		case GLUT_KEY_RIGHT:
			campos = campos + camr*muvlen;
			break;
		case GLUT_KEY_LEFT:
			campos = campos - camr*muvlen;
			break;
		case GLUT_KEY_PAGE_UP:
			campos = campos + camu*muvlen;
			break;
		case GLUT_KEY_PAGE_DOWN:
			campos = campos - camu*muvlen;
			break;
		default:
			break;
	}
}

void loaddata(){
    if(!freopen("scene.txt", "r", stdin)){
        cout<<"no file";
        exit(0);
    }
    cin>>reclevel;
    cin>>imageWidth;
    imageHeight = imageWidth;
    int objcount;
    cin>>objcount;
    for(int i=0; i<objcount; i++){
        string typ;
        cin>>typ;
        if(typ=="sphere"){
            Sphere *sp = new Sphere();
            cin >> *sp;
            objects.push_back(sp);
        }
        else if(typ=="triangle"){
            Triangle *tr = new Triangle();
            cin>> *tr;
            objects.push_back(tr);
        }
        else {
            Surf *sf = new Surf();
            cin >> *sf;
            objects.push_back(sf);
        }
    }

    int pLcount;
    cin>>pLcount;
    for(int i=0; i<pLcount; i++){
        PointLight pl;
        cin>>pl.light_pos>>pl.color;
        pointLights.push_back(pl);
    }
    int sLcount;
    cin>>sLcount;
    for(int i=0; i<sLcount; i++){
        SpotLight sl;
        cin>>sl.point_light.light_pos>>sl.point_light.color;
        cin>>sl.light_direction>>sl.cutoff_angle;
        sl.cutoff_angle *= PI / 180.0;
        spotLights.push_back(sl);
    }

    Floor *phlor = new Floor(1000, 20);
    phlor->coEfficients = {0.2, 0.3, 0.3, 0.6};
    phlor->shine = 5;
    objects.push_back(phlor);

    cout<<"loading done "<<endl;
}

void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	/********************
	/ set-up camera here
	********************/
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(campos.x, campos.y, campos.z, caml.x+campos.x, caml.y+campos.y, caml.z+campos.z, camu.x, camu.y, camu.z);
	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);

	/****************************
	/ Add your objects from here
	****************************/

    for(Object *ob : objects)
        ob->draw();

    for(PointLight pl : pointLights)
        pl.draw();

    for(SpotLight sl : spotLights)
        sl.draw();

	glutSwapBuffers();
}

void animate(){
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	campos = Point(100, 100, 0);
	camu = Point(0, 0, 1);
	caml = Point(-1/sqrt(2), -1/sqrt(2), 0);
    camr = Point(-1/sqrt(2), 1/sqrt(2), 0);
	//clear the screen
	glClearColor(0,0,0,0);
	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);
	//initialize the matrix
	glLoadIdentity();
	gluPerspective(viewAngle, 1, 1,	1000.0);
}

int main(int argc, char **argv){

    loaddata();

	glutInit(&argc,argv);
	glutInitWindowSize(width, height);
	glutInitWindowPosition(500, 200);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("raytracing");
	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing
	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
