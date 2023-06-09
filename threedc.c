/* threed.c   */

/* Rotating cube demo*/
/* demonstrates many OpenGL features including */
/* smooth shading, antialiasing, texture mapping, */
/* color interpolation */

/* uses a trackball interface */

/* E. Angel, Interactive Computer Graphics */
/* A Top-Down Approach with OpenGL, Third Edition */
/* Addison-Wesley Longman, 2003 */


/* Prints out instructions */

/*Both normals and colors are assigned to the vertices */
/*Cube is centered at origin so (unnormalized) normals
are the same as the vertex values */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
#include <GL/glu.h>

#define bool int
#define false 0
#define true 1
/* #define M_PI 3.14159 */

/*
** Global settings controlled by keyboard input.
*/
/*bool texEnabled 		= false;
bool mipmapEnabled 		= false;
bool fastTexture 		= true;
bool fogEnabled 		= false;
bool depthEnabled 		= true;
bool lineAAEnabled 		= false;
bool lightingEnabled 		= true;
bool smoothEnabled 		= false;
bool drawLines 			= false;
bool idleSpin 			= false;
bool perspectiveXform 		= true;
bool TopOn			= true;*/

bool texEnabled;
bool mipmapEnabled;
bool fastTexture;
bool fogEnabled;
bool depthEnabled;
bool lineAAEnabled;
bool lightingEnabled;
bool smoothEnabled;
bool drawLines;
bool idleSpin;
bool perspectiveXform;
bool TopOn;




/*
** Global settings.
*/

float 	nnear = 3.0;	/* near clipping plane in eye coords */
float 	nfar = 7.0;	/* far clipping plane in eye coords */
float 	viewxform_z = -5.0;

int 	winWidth, winHeight;

float 	angle = 0.0, axis[3], trans[3];
bool 	trackballEnabled = true;
bool 	trackballMove = false;
bool 	trackingMouse = false;
bool 	redrawContinue = false;

/*GLfloat lightXform[4][4] = {
    {1.0, 0.0, 0.0, 0.0},
    {0.0, 1.0, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0},
    {0.0, 0.0, 0.0, 1.0}
};
GLfloat objectXform[4][4] = {
    {1.0, 0.0, 0.0, 0.0},
    {0.0, 1.0, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0},
    {0.0, 0.0, 0.0, 1.0}
};*/
GLfloat lightXform[4][4];
GLfloat objectXform[4][4];
GLfloat *trackballXform = (GLfloat*)objectXform;

/* global arrays containing object */
/* ------------------------------------ */
int numnp,numel;
float xxx[79999],yyy[79999],zzz[79999][2];
int kx[79999][4];
GLfloat vertice1[79999][3];
GLfloat vertice2[79999][3];
GLfloat normal1[79999][3];
GLfloat normal2[79999][3];
GLfloat fnormal1[79999][3];
GLfloat fnormal2[79999][3];
GLfloat vcolor1[79999][3];
GLfloat vcolor2[79999][3];
float colors[256][3];
int ncolors;
/* ------------------------------------ */
void display(void);
void spinCube(void);
void formnormals(void);
void setMenuEntries(bool);
/*----------------------------------------------------------------------*/
/*
** Draw the cube.
*/

/*----------------------------------------------------------------------*/
/* 
** Materials setup.
*/
typedef struct materialStruct {
   GLfloat ambient[4];
   GLfloat diffuse[4];
   GLfloat specular[4];
   GLfloat shininess;
} materialStruct;

materialStruct brassMaterials = {
    {0.33F, 0.22F, 0.03F, 1.0F},
    {0.78F, 0.57F, 0.11F, 1.0F},
    {0.99F, 0.91F, 0.81F, 1.0F},
    27.8F
};
materialStruct redPlasticMaterials = {
    {0.3F, 0.0F, 0.0F, 1.0F},
    {0.6F, 0.0F, 0.0F, 1.0F},
    {0.0F, 0.0F, 0.0F, 1.0F},
    128.0F
};
materialStruct colorCubeMaterials = {
    { 1.0F, 1.0F, 1.0F, 1.0F },
    { 1.0F, 1.0F, 1.0F, 1.0F },
    { 1.0F, 1.0F, 1.0F, 1.0F },
    100.0F
};

materialStruct *currentMaterials = &brassMaterials;

void materials( materialStruct *materials)
{
    /* define material proerties for front face of all polygons */
    glMaterialfv(GL_FRONT, GL_AMBIENT, materials->ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, materials->diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, materials->specular);
    glMaterialf(GL_FRONT, GL_SHININESS, materials->shininess);
}

/*----------------------------------------------------------------------*/
/*
** Lighting setup.
*/
typedef struct lightingStruct {
    GLfloat ambient[4];
    GLfloat diffuse[4];
    GLfloat specular[4];
} lightingStruct;

lightingStruct whiteLighting = {
    {0.5F, 0.5F, 0.5F, 1.0F},
    {1.0F, 1.0F, 1.0F, 1.0F},
    {1.0F, 1.0F, 1.0F, 1.0F}
};
lightingStruct colorCubeLighting = {
    {0.2F, 0.0F, 0.0F, 1.0F},
    {0.0F, 1.0F, 0.0F, 1.0F},
    {0.0F, 0.0F, 1.0F, 1.0F}
};

lightingStruct *currentLighting = &whiteLighting;

void lighting(lightingStruct *lightSettings)
{
    /* set up light 0 ambient, diffuse, specular, and spotlight */

    glLightfv(GL_LIGHT0, GL_AMBIENT, lightSettings->ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightSettings->diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightSettings->specular);

    glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, 0.0F);
    glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 180.0F);
    glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1.0F);
    glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.0F);
    glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.0F);

    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
}

/*----------------------------------------------------------------------*/
/*
** Light position setup.
*/
void setLightPos()
{
    /*GLfloat light0_pos[4] = { 0.90F, 0.90F, 2.25F, 0.00F };
    GLfloat light0_spotDir[3] = {0.0F, 0.0F, -1.0F};*/
    GLfloat light0_pos[4] = {-5.00F,-5.00F,-5.00F, 0.00F };
    GLfloat light0_spotDir[3] = {-1.0F, -1.0F,-1.0F};
    glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light0_spotDir);
}

/*----------------------------------------------------------------------*/
/*
** Initial texture settings.
*/
void texture()
{
   GLubyte image[64][64][3];
   int i, j, r, c;
   for(i=0;i<64;i++)
   {
      for(j=0;j<64;j++)
      {
	 c = ((((i&0x8)==0)^((j&0x8))==0))*255;
	 image[i][j][0]= (GLubyte) c;
	 image[i][j][1]= (GLubyte) c;
	 image[i][j][2]= (GLubyte) c;
      }
   }
   glPixelStorei(GL_UNPACK_ALIGNMENT,1);
   glTexImage2D(GL_TEXTURE_2D,0,3,64,64,0,GL_RGB,GL_UNSIGNED_BYTE, image);
   glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
   glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
   glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
   glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
   gluBuild2DMipmaps(GL_TEXTURE_2D,3,64,64,GL_RGB,GL_UNSIGNED_BYTE, image);
}

/*
** Set initial state.
*/
void initSettings(void)
{
    texture();
    glLineWidth(3.0);
    setMenuEntries(true);
}

/*----------------------------------------------------------------------*/
/*
** Set state according to user interaction.
*/
void userSettings(void) 
{
    lighting(currentLighting);
    materials(currentMaterials);

    if (lightingEnabled) {
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
    } else {
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
    }

    if (smoothEnabled) {
	glShadeModel(GL_SMOOTH);
    } else {
	glShadeModel(GL_FLAT);
    }

    if(idleSpin) {
        glutIdleFunc(spinCube);
    } else {
	glutIdleFunc(NULL);
    }
    if (texEnabled) {
	glEnable(GL_TEXTURE_2D);
    } else {
	glDisable(GL_TEXTURE_2D);
    }
    if (fastTexture) {
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);
    } else {
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    }
    if (mipmapEnabled) {
	glTexParameterf(GL_TEXTURE_2D,
	    GL_TEXTURE_MIN_FILTER,GL_NEAREST_MIPMAP_NEAREST);
    } else {
	glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    }

    if (fogEnabled) {

	float fogColor[] = {0.7, 0.6, 0.6, 1.0};

	glClearColor(fogColor[0], fogColor[1], fogColor[2], fogColor[3]);
	glEnable(GL_FOG);
	glFogi(GL_FOG_MODE, GL_LINEAR);
	glFogf(GL_FOG_DENSITY, 1.0);
	glFogf(GL_FOG_START, nnear);
	glFogf(GL_FOG_END, nfar);
	glFogfv(GL_FOG_COLOR, fogColor);

    } else {
	glDisable(GL_FOG);
	glClearColor(0.0,0.0,0.0,1.0);
    }

    if (lineAAEnabled) {
	glEnable(GL_BLEND);
    } else {
	glDisable(GL_BLEND);
    }

    if (lineAAEnabled) {
	glEnable(GL_LINE_SMOOTH);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    } else {
	glDisable(GL_LINE_SMOOTH);
    }

    if (depthEnabled) {
	glEnable(GL_DEPTH_TEST); /* Enable hidden--surface--removal */
    } else {
	glDisable(GL_DEPTH_TEST);
    }

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if(perspectiveXform) {
	glFrustum(-1.25, 1.25, -1.25, 1.25, nnear, nfar);
	viewxform_z = -5.0;
    } else {
	glOrtho(-2.0, 2.0, -2.0, 2.0, nnear, nfar);
	viewxform_z = -5.0;
    }
    glMatrixMode(GL_MODELVIEW);
}


void polygon(int a, int b, int c , int d, int face)
{

    /* draw a polygon via list of vertices */
    if(drawLines) {
        glColor3f(1.0,1.0,1.0);
        if(TopOn) {
          glBegin(GL_LINE_LOOP);
  	    glVertex3fv(vertice1[a]);
  	    glVertex3fv(vertice1[b]);
  	    glVertex3fv(vertice1[c]);
  	    glVertex3fv(vertice1[d]);
          glEnd();
        }
        glBegin(GL_LINE_LOOP);
  	  glVertex3fv(vertice2[a]);
  	  glVertex3fv(vertice2[b]);
  	  glVertex3fv(vertice2[c]);
  	  glVertex3fv(vertice2[d]);
        glEnd();
    }
    else {
        if(TopOn) {
  	  glNormal3fv(fnormal1[face]);
          glBegin(GL_POLYGON);
  	    glColor3fv(vcolor1[a]);
  	    /* glNormal3fv(normals[a]); */
  	    glTexCoord2f(0.0,0.0);
  	    glVertex3fv(vertice1[a]);
  	    glColor3fv(vcolor1[b]);
  	    /* glNormal3fv(normals[b]); */
  	    glTexCoord2f(0.0,1.0);
  	    glVertex3fv(vertice1[b]);
  	    glColor3fv(vcolor1[c]);
  	    /* glNormal3fv(normals[c]); */
  	    glTexCoord2f(1.0,1.0);
  	    glVertex3fv(vertice1[c]);
  	    glColor3fv(vcolor1[d]);
  	    /* glNormal3fv(normals[d]); */
  	    glTexCoord2f(1.0,0.0);
  	    glVertex3fv(vertice1[d]);
          glEnd();
        }
	glNormal3fv(fnormal2[face]);
        glBegin(GL_POLYGON);
  	  glColor3fv(vcolor2[a]);
  	  /* glNormal3fv(normals[a]); */
  	  glTexCoord2f(0.0,0.0);
  	  glVertex3fv(vertice2[a]);
  	  glColor3fv(vcolor2[b]);
  	  /* glNormal3fv(normals[b]); */
  	  glTexCoord2f(0.0,1.0);
  	  glVertex3fv(vertice2[b]);
  	  glColor3fv(vcolor2[c]);
  	  /* glNormal3fv(normals[c]); */
  	  glTexCoord2f(1.0,1.0);
  	  glVertex3fv(vertice2[c]);
  	  glColor3fv(vcolor2[d]);
  	  /* glNormal3fv(normals[d]); */
  	  glTexCoord2f(1.0,0.0);
  	  glVertex3fv(vertice2[d]);
        glEnd();
    }
}

void colorcube(void)
{
    int i;
    /* map vertices to faces */
    for(i=1;i<=numel;i++)
    {  polygon(kx[i][1],kx[i][2],kx[i][3],kx[i][4],i-1);
    /* printf("%d %d %d %d %d\n",kx[i][1],kx[i][2],kx[i][3],kx[i][4],i-1); */

    }
}

/*----------------------------------------------------------------------*/
/* 
** These functions implement a simple trackball-like motion control.
*/

float lastPos[3] = {0.0F, 0.0F, 0.0F};
int curx, cury;
int startX, startY;

void
trackball_ptov(int x, int y, int width, int height, float v[3])
{
    float d, a;

    /* project x,y onto a hemi-sphere centered within width, height */
    v[0] = (2.0F*x - width) / width;
    v[1] = (height - 2.0F*y) / height;
    d = (float) sqrt(v[0]*v[0] + v[1]*v[1]);
    v[2] = (float) cos((M_PI/2.0F) * ((d < 1.0F) ? d : 1.0F));
    a = 1.0F / (float) sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] *= a;
    v[1] *= a;
    v[2] *= a;
}


void
mouseMotion(int x, int y)
{
    float curPos[3], dx, dy, dz;

    trackball_ptov(x, y, winWidth, winHeight, curPos);

    dx = curPos[0] - lastPos[0];
    dy = curPos[1] - lastPos[1];
    dz = curPos[2] - lastPos[2];

    if (dx || dy || dz) {
	angle = 90.0F * sqrt(dx*dx + dy*dy + dz*dz);

	axis[0] = lastPos[1]*curPos[2] - lastPos[2]*curPos[1];
	axis[1] = lastPos[2]*curPos[0] - lastPos[0]*curPos[2];
	axis[2] = lastPos[0]*curPos[1] - lastPos[1]*curPos[0];

	lastPos[0] = curPos[0];
	lastPos[1] = curPos[1];
	lastPos[2] = curPos[2];
    }

    glutPostRedisplay();
}

void
startMotion(long time, int button, int x, int y)
{
    if (!trackballEnabled) return;

    trackingMouse = true;
    redrawContinue = false;
    startX = x; startY = y;
    curx = x; cury = y;
    trackball_ptov(x, y, winWidth, winHeight, lastPos);
    trackballMove = true;
}

void
stopMotion(long time, int button, int x, int y)
{
    if (!trackballEnabled) return;

    trackingMouse = false;

    if (startX != x || startY != y) {
	redrawContinue = true;
    } else {
	angle = 0.0F;
	redrawContinue = false;
	trackballMove = false;
    }
}

/*----------------------------------------------------------------------*/

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    /* view transform */
    glLoadIdentity();
    glTranslatef(0.0,0.0,viewxform_z);

    if (trackballMove) {
	glPushMatrix();
	    glLoadIdentity();
	    glRotatef(angle, axis[0], axis[1], axis[2]);
	    glMultMatrixf((GLfloat *) trackballXform);
	    glGetFloatv(GL_MODELVIEW_MATRIX, trackballXform);
	glPopMatrix();

    }
    glPushMatrix();
	glMultMatrixf((GLfloat*)lightXform);
	setLightPos();
    glPopMatrix();
    glPushMatrix();
	glMultMatrixf((GLfloat *) objectXform);
	colorcube();
    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}

/*----------------------------------------------------------------------*/

void mouseButton(int button, int state, int x, int y)
{
    switch (button) {
      case GLUT_LEFT_BUTTON:
	 trackballXform = (GLfloat*)objectXform;
	 break;
      case GLUT_RIGHT_BUTTON:
	 trackballXform = (GLfloat*)lightXform;
	 break;
    }
    switch(state) {
    case GLUT_DOWN:
	startMotion(0,1, x,y);
	break;
    case GLUT_UP:
	stopMotion(0,1, x,y);
	break;
    } 
}

void myReshape(int w, int h)
{
    glViewport(0, 0, w, h);
    winWidth = w;
    winHeight = h;
}

void spinCube()
{
    if (redrawContinue) glutPostRedisplay();
}
void backup()
{   FILE *fileid; int i,j; float tmp;
    fileid=fopen("backup.data","w");
    /*printf("backup.data\n");*/
    /*printf(" %d %d %d %d %d %d %d %d %d %d %d %d\n",texEnabled,
            mipmapEnabled, fastTexture, fogEnabled, depthEnabled,
            lineAAEnabled, lightingEnabled, smoothEnabled, drawLines,
            idleSpin, perspectiveXform,TopOn);*/
    fprintf(fileid," %d %d %d %d %d %d %d %d %d %d %d %d\n",texEnabled,
            mipmapEnabled, fastTexture, fogEnabled, depthEnabled,
            lineAAEnabled, lightingEnabled, smoothEnabled, drawLines,
            idleSpin, perspectiveXform,TopOn);
    /*printf(" %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n",
            objectXform[0][0],
            objectXform[0][1],
            objectXform[0][2],
            objectXform[0][3],
            objectXform[1][0],
            objectXform[1][1],
            objectXform[1][2],
            objectXform[1][3],
            objectXform[2][0],
            objectXform[2][1],
            objectXform[2][2],
            objectXform[2][3],
            objectXform[3][0],
            objectXform[3][1],
            objectXform[3][2],
            objectXform[3][3]);*/
    fprintf(fileid," %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n",
            objectXform[0][0],
            objectXform[0][1],
            objectXform[0][2],
            objectXform[0][3],
            objectXform[1][0],
            objectXform[1][1],
            objectXform[1][2],
            objectXform[1][3],
            objectXform[2][0],
            objectXform[2][1],
            objectXform[2][2],
            objectXform[2][3],
            objectXform[3][0],
            objectXform[3][1],
            objectXform[3][2],
            objectXform[3][3]);
    /*printf(" %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n",
            lightXform[0][0],
            lightXform[0][1],
            lightXform[0][2],
            lightXform[0][3],
            lightXform[1][0],
            lightXform[1][1],
            lightXform[1][2],
            lightXform[1][3],
            lightXform[2][0],
            lightXform[2][1],
            lightXform[2][2],
            lightXform[2][3],
            lightXform[3][0],
            lightXform[3][1],
            lightXform[3][2],
            lightXform[3][3]);*/
    fprintf(fileid," %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n",
            lightXform[0][0],
            lightXform[0][1],
            lightXform[0][2],
            lightXform[0][3],
            lightXform[1][0],
            lightXform[1][1],
            lightXform[1][2],
            lightXform[1][3],
            lightXform[2][0],
            lightXform[2][1],
            lightXform[2][2],
            lightXform[2][3],
            lightXform[3][0],
            lightXform[3][1],
            lightXform[3][2],
            lightXform[3][3]);
    fclose(fileid);
}
void readbackup()
{   FILE *fileid; int i,j; float tmp;
    fileid=fopen("backedup.data","r");
    /*printf("in readbackup %d\n",fileid);*/
    j=fscanf(fileid," %d %d %d %d %d %d %d %d %d %d %d %d\n",&texEnabled,
            &mipmapEnabled, &fastTexture, &fogEnabled, &depthEnabled,
            &lineAAEnabled, &lightingEnabled, &smoothEnabled, &drawLines,
            &idleSpin, &perspectiveXform,&TopOn);
    /*printf("EOF %d\n",j);*/
    /*printf(" %d %d %d %d %d %d %d %d %d %d %d %d\n",texEnabled,
            mipmapEnabled, fastTexture, fogEnabled, depthEnabled,
            lineAAEnabled, lightingEnabled, smoothEnabled, drawLines,
            idleSpin, perspectiveXform,TopOn);*/
    fscanf(fileid," %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n",
            &objectXform[0][0],
            &objectXform[0][1],
            &objectXform[0][2],
            &objectXform[0][3],
            &objectXform[1][0],
            &objectXform[1][1],
            &objectXform[1][2],
            &objectXform[1][3],
            &objectXform[2][0],
            &objectXform[2][1],
            &objectXform[2][2],
            &objectXform[2][3],
            &objectXform[3][0],
            &objectXform[3][1],
            &objectXform[3][2],
            &objectXform[3][3]);
    /*printf(" %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n",
            objectXform[0][0],
            objectXform[0][1],
            objectXform[0][2],
            objectXform[0][3],
            objectXform[1][0],
            objectXform[1][1],
            objectXform[1][2],
            objectXform[1][3],
            objectXform[2][0],
            objectXform[2][1],
            objectXform[2][2],
            objectXform[2][3],
            objectXform[3][0],
            objectXform[3][1],
            objectXform[3][2],
            objectXform[3][3]);*/
    fscanf(fileid," %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n",
            &lightXform[0][0],
            &lightXform[0][1],
            &lightXform[0][2],
            &lightXform[0][3],
            &lightXform[1][0],
            &lightXform[1][1],
            &lightXform[1][2],
            &lightXform[1][3],
            &lightXform[2][0],
            &lightXform[2][1],
            &lightXform[2][2],
            &lightXform[2][3],
            &lightXform[3][0],
            &lightXform[3][1],
            &lightXform[3][2],
            &lightXform[3][3]);
    /*printf(" %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n",
            lightXform[0][0],
            lightXform[0][1],
            lightXform[0][2],
            lightXform[0][3],
            lightXform[1][0],
            lightXform[1][1],
            lightXform[1][2],
            lightXform[1][3],
            lightXform[2][0],
            lightXform[2][1],
            lightXform[2][2],
            lightXform[2][3],
            lightXform[3][0],
            lightXform[3][1],
            lightXform[3][2],
            lightXform[3][3]);*/
    fclose(fileid);
}

void userEventAction(int key) {
    int i;
    switch(key) {
    case '0':				/* wire frame/polygon  */
      drawLines = !drawLines;
      break;
    case '1':
      smoothEnabled = !smoothEnabled;
      break;
    case '2':				/* lighting */
      lightingEnabled = !lightingEnabled;
      break;
    case '3':				/* texture */
      texEnabled = !texEnabled;
      break;
    case '4':				/* fog */
      fogEnabled = !fogEnabled;
      break;
    case '5':				/* HSR */
      depthEnabled = !depthEnabled;
      break;
    case '6':				/* line aa */
      lineAAEnabled = !lineAAEnabled;
      break;
    case '7':				/* perpective texture */
      fastTexture = !fastTexture;
      break;
    case 'b':
      currentMaterials = &brassMaterials;
      break;
    case 'c':
      currentMaterials = &colorCubeMaterials;
      break;
    case 'C':
      currentLighting = &colorCubeLighting;
      break;
    case 't':
      TopOn = !TopOn;
      break;
    case 'i':
      idleSpin = !idleSpin;
      break;
    case 'm':				/* mipmapped texture */
      mipmapEnabled = !mipmapEnabled;
      break;
    case 'p':				/* perspective/ortho */
      perspectiveXform = !perspectiveXform;
      break;
    case 'r':
      currentMaterials = &redPlasticMaterials;
      break;
    case 'w':
      currentLighting = &whiteLighting;
      break;
    case 's':
      angle=0;
      break;
    case 'x':
      for(i=1;i<=numnp;i++)
      {  vertice1[i][0]=vertice1[i][0]*3./4.;  
         vertice2[i][0]=vertice2[i][0]*3./4.; }
      formnormals();
      break;
    case 'X':
      for(i=1;i<=numnp;i++)
      {  vertice1[i][0]=vertice1[i][0]*4./3.;  
         vertice2[i][0]=vertice2[i][0]*4./3.; }
      formnormals();
      break;
    case 'y':
      for(i=1;i<=numnp;i++)
      {  vertice1[i][1]=vertice1[i][1]*3./4.;  
         vertice2[i][1]=vertice2[i][1]*3./4.; }
      formnormals();
      break;
    case 'Y':
      for(i=1;i<=numnp;i++)
      {  vertice1[i][1]=vertice1[i][1]*4./3.;  
         vertice2[i][1]=vertice2[i][1]*4./3.; }
      formnormals();
      break;
    case 'z':
      for(i=1;i<=numnp;i++)
      {  vertice1[i][2]=vertice1[i][2]*3./4.;  
         vertice2[i][2]=vertice2[i][2]*3./4.; }
      formnormals();
      break;
    case 'Z':
      for(i=1;i<=numnp;i++)
      {  vertice1[i][2]=vertice1[i][2]*4./3.;  
         vertice2[i][2]=vertice2[i][2]*4./3.; }
      formnormals();
      break;
    case 27:
      backup();
      exit(0);
    default:
      break;
    }
    userSettings();
    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
    userEventAction(key);
}

/*----------------------------------------------------------------------*/

typedef struct menuEntryStruct {
    char *label;
    char key;
} menuEntryStruct;

static menuEntryStruct mainMenu[] = {
    "lines/polygons", 		'0',
    "flat/smooth", 		'1',
    "lighting", 		'2',
    "texture", 			'3',
    "fog", 			'4',
    "HSR", 			'5',
    "line smooth", 		'6',
    "motion", 			'i',
    "ortho/perspective", 	'p',
    "quit", 			27,
};
int mainMenuEntries = sizeof(mainMenu)/sizeof(menuEntryStruct);

void selectMain(int choice)
{
    userEventAction(mainMenu[choice].key);
}

static menuEntryStruct materialsMenu[] = {
    "brass", 		'b',
    "white",		'c',
    "red plastic", 	'r',
};
int materialsMenuEntries = sizeof(materialsMenu)/sizeof(menuEntryStruct);

void selectMaterials(int choice)
{
    userEventAction(materialsMenu[choice].key);
}


static menuEntryStruct lightingMenu[] = {
    "white", 		'w',
    "color",		'C',
};
int lightingMenuEntries = sizeof(lightingMenu)/sizeof(menuEntryStruct);

void selectLighting(int choice)
{
    userEventAction(lightingMenu[choice].key);
}

void setMenuEntries(bool init)
{
    int i, sub1, sub2;

    if (init) {
	sub1 = glutCreateMenu(selectMaterials);
	for (i=0; i < materialsMenuEntries; i++) {
	    glutAddMenuEntry(materialsMenu[i].label, i);
	}
	sub2 = glutCreateMenu(selectLighting);
	for (i=0; i < lightingMenuEntries; i++) {
	    glutAddMenuEntry(lightingMenu[i].label, i);
	}
	glutCreateMenu(selectMain);
	for (i=0; i < mainMenuEntries; i++) {
	    glutAddMenuEntry(mainMenu[i].label, i);
	}
	glutAddSubMenu("materials", sub1);
	glutAddSubMenu("lighting", sub2);
	glutAttachMenu(GLUT_MIDDLE_BUTTON);
    } else {
    }
}
float mymax(float aaa,float bbb)
{  if(aaa>bbb) return aaa; else return bbb;
}
float mymin(float aaa,float bbb)
{  if(aaa<bbb) return aaa; else return bbb;
}
float myscale(float aaa,float amin,float amax)
{  float tmp;
   if(amax==amin) 
     tmp=1; 
   else
     tmp=2*(aaa-amin)/(amax-amin)-1;
   return tmp;
}
void formnormals()
{   int i,lm[4];
    float dx13,dx24,dy13,dy24,dz13,dz24,fnx,fny,fnz,fmag;
    for(i=1;i<=numel;i++)
    { lm[1]=kx[i][1];lm[2]=kx[i][2];lm[3]=kx[i][3];lm[4]=kx[i][4];
      dx13=vertice1[lm[3]][0]-vertice1[lm[1]][0];
      dy13=vertice1[lm[3]][1]-vertice1[lm[1]][1];
      dz13=vertice1[lm[3]][2]-vertice1[lm[1]][2];
      dx24=vertice1[lm[4]][0]-vertice1[lm[2]][0];
      dy24=vertice1[lm[4]][1]-vertice1[lm[2]][1];
      dz24=vertice1[lm[4]][2]-vertice1[lm[2]][2];
      fnx=dz13*dy24-dz24*dy13;
      fny=dx24*dy13-dx13*dy24;
      fnz=dx13*dz24-dx24*dz13;
      fmag=sqrt(fnx*fnx+fny*fny+fnz*fnz);
      fnx=fnx/fmag; fny=fny/fmag; fnz=fnz/fmag;
      fnormal1[i][0]=fnx; fnormal1[i][1]=fny; fnormal1[i][2]=fnz;
      dx13=vertice2[lm[3]][0]-vertice2[lm[1]][0];
      dy13=vertice2[lm[3]][1]-vertice2[lm[1]][1];
      dz13=vertice2[lm[3]][2]-vertice2[lm[1]][2];
      dx24=vertice2[lm[4]][0]-vertice2[lm[2]][0];
      dy24=vertice2[lm[4]][1]-vertice2[lm[2]][1];
      dz24=vertice2[lm[4]][2]-vertice2[lm[2]][2];
      fnx=dz13*dy24-dz24*dy13;
      fny=dx24*dy13-dx13*dy24;
      fnz=dx13*dz24-dx24*dz13;
      fmag=sqrt(fnx*fnx+fny*fny+fnz*fnz);
      fnx=fnx/fmag; fny=fny/fmag; fnz=fnz/fmag;
      fnormal2[i][0]=fnx; fnormal2[i][1]=fny; fnormal2[i][2]=fnz;
      /*printf("%d %f %f %f %f\n",i,fnx,fny,fnz,fmag);*/
    }
}

void loadvertices() 
{   int i,j,icolor;
    float xmin,xmax,ymin,ymax,zmin,zmax,xdel,ydel,xmid,ymid;
    float ztmin,ztmax,zbmin,zbmax;
    xmin=1e30;xmax=-xmin;
    ymin=1e30;ymax=-ymin;
    zmin=1e30;zmax=-zmin;
    for(i=1;i<=numnp;i++)
    {  xmin=mymin(xmin,xxx[i]); xmax=mymax(xmax,xxx[i]);
       ymin=mymin(ymin,yyy[i]); ymax=mymax(ymax,yyy[i]);
       zmin=mymin(zmin,zzz[i][0]); zmax=mymax(zmax,zzz[i][0]);
       ztmin=mymin(ztmin,zzz[i][0]); ztmax=mymax(ztmax,zzz[i][0]);
       if(zzz[i][1]!=-9999.) {
         zmin=mymin(zmin,zzz[i][1]); zmax=mymax(zmax,zzz[i][1]);
         zbmin=mymin(zbmin,zzz[i][1]); zbmax=mymax(zbmax,zzz[i][1]);
       }
    }
    xdel=xmax-xmin;ydel=ymax-ymin;
    xmid=0.5*(xmax+xmin);ymid=0.5*(ymax+ymin);
    if(xdel > ydel)  
    { ymin=ymid-xdel/2.; ymax=ymid+xdel/2.;
    }
    else
    { xmin=xmid-ydel/2.; xmax=xmid+ydel/2.;
    }
    printf("%f %f %f %f\n",xmin,xmax,xmid,xdel);
    printf("%f %f %f %f\n",ymin,ymax,ymid,ydel);
    printf("%f %f \n",zmin,zmax);
    for(i=1;i<=numnp;i++)
    {  vertice1[i][0]= myscale(xxx[i],xmin,xmax);
       vertice1[i][1]= myscale(yyy[i],ymin,ymax);
       vertice1[i][2]= myscale(zzz[i][0],zmin,zmax);
       normal1[i][0]=vertice1[i][0];
       normal1[i][1]=vertice1[i][1];
       normal1[i][2]=vertice1[i][2];
       icolor=ncolors-0.5*(myscale(zzz[i][0],ztmin,ztmax)+1)*(ncolors-1);
      /*printf("%f %d \n",vertice1[i][2],icolor);*/
       vcolor1[i][0]=colors[icolor][0];
       vcolor1[i][1]=colors[icolor][1];
       vcolor1[i][2]=colors[icolor][2];
       vertice1[i][2]= 0.2*vertice1[i][2];
      /* printf("%f %f %f \n",xxx[i],yyy[i],zzz[i][0]); */
       vertice2[i][0]= myscale(xxx[i],xmin,xmax);
       vertice2[i][1]= myscale(yyy[i],ymin,ymax);
       vertice2[i][2]= myscale(zzz[i][1],zmin,zmax);
       normal2[i][0]=vertice2[i][0];
       normal2[i][1]=vertice2[i][1];
       normal2[i][2]=vertice2[i][2];
       icolor=ncolors-0.5*(myscale(zzz[i][1],zbmin,zbmax)+1)*(ncolors-1);
      /*printf("%f %d \n",vertice2[i][2],icolor);*/
       vcolor2[i][0]=colors[icolor][0];
       vcolor2[i][1]=colors[icolor][1];
       vcolor2[i][2]=colors[icolor][2];
       vertice2[i][2]= 0.2*vertice2[i][2];
      /* printf("%f %f %f \n",xxx[i],yyy[i],zzz[i][1]); */
    }
    formnormals();
}
void mycolors() /* reads from file written by caller */
{   FILE *fileid; int i,j,itmp; float tmp;
    fileid=fopen("tmp.color","r");
    fscanf(fileid," %d",&ncolors);
    for(i=1;i<=ncolors;i++)
    {  fscanf(fileid," %f %f %f\n",&colors[i][0],
                                   &colors[i][1],
                                   &colors[i][2]);
    }
    fclose(fileid);
}
void myreader() /* reads from file written by caller */
{   FILE *fileid; int i,j,itmp; float tmp;
    fileid=fopen("tmp.data","r");
    fscanf(fileid,"%d %d",&numnp,&numel);
    for(i=1;i<=numnp;i++)
    {  fscanf(fileid," %f",&tmp); xxx[i]=tmp; }
    for(i=1;i<=numnp;i++)
    {  fscanf(fileid," %f",&tmp); yyy[i]=tmp;
       /*printf("%d %f %f \n",i,xxx[i],yyy[i]); */
    }
    for(i=1;i<=numnp;i++)
    {  fscanf(fileid," %f %f",&zzz[i][0],&zzz[i][1]);
       /* printf("%d %f %f %f\n",i,xxx[i],yyy[i],zzz[i][0]); */
    }
    for(i=1;i<=numel;i++)
    {  for(j=1;j<=4;j++)
       {  fscanf(fileid," %d",&itmp); kx[i][j]=itmp; }
       /*printf("%d %d %d %d %d\n",i,kx[i]);*/
    }
    fclose(fileid);
}

/*----------------------------------------------------------------------*/
int
main(int argc, char **argv)
{   
    mycolors();
    myreader();
    loadvertices();
    glutInit(&argc, argv); 
    readbackup();
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1500, 1500);
    glutCreateWindow("colorcube");
    glutReshapeFunc(myReshape);
    glutDisplayFunc(display);
    glutIdleFunc(spinCube);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMotion);
    glutKeyboardFunc(keyboard);
    initSettings();
    userSettings();
    glutMainLoop();
}
