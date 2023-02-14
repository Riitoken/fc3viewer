// fc3viewer.cpp

/*
class        license
type         GLORYWARE
ipname       FC3VIEWER
trademark    FARCRAFT®
author       Ray Edward Bornert II
date         2020-SEP-22 TUE
royalty      Free
doc          https://docs.google.com/document/d/1xAZ-WAHxwBuu1H-LszPGuMQXLryebBYzBdrHzV2Gom0
*/

#include <windows.h>		/* must include this before GL/gl.h */
#include <GL/gl.h>			/* OpenGL header file */
#include <GL/glu.h>			/* OpenGL utilities header file */
#include <stdio.h>
#include <sys/timeb.h>
#include <string>
using namespace std;

#include "fc3.h"
fc3_s	the_fc3;

string		the_fn			= ""	;
string		the_cmd_line	= ""	;
unsigned	the_dlist		= 1		; // we choose 1
unsigned	the_htex		= 0		; // opengl will choose
double		the_extent		= 1		; // default on load
int			tilt_mode		= 0		;

HDC			the_hDC			= NULL	; // device context
HINSTANCE	the_hInstance	= 0		; // instance id
HINSTANCE	the_prev_inst	= 0		;
HINSTANCE	the_curr_inst	= 0		;
HWND		the_hwnd		= NULL	;

void glerr( const char* tag=NULL )
{
	GLenum E=GL_NO_ERROR;
	GLenum e=GL_NO_ERROR;
	while(1)
	{
		E = e;
		e = glGetError();
		if(E==e || GL_NO_ERROR==e)
			break;
		char msg[128];msg[0]=0;
		sprintf_s(msg,sizeof(msg),"OpenGL error %d %s",e,tag?tag:"");
		MessageBox(NULL,msg, "OpenGL Error", MB_OK);
	}
}

void fc3__gl( const fc3_vertf_s& vert )
{
	glTexCoord2f( vert.t.x,vert.t.y );
	glNormal3f( vert.n.x,vert.n.y,vert.n.z );
	glVertex3f( vert.v.x,vert.v.y,vert.v.z );
}

int fc3__compile( const fc3_s& fc3, const unsigned int dlist )
{
	if(0==dlist)return -1;
	if(NULL==fc3.pv)return -2;
	if(NULL==fc3.pt)return -3;

	glNewList( dlist, GL_COMPILE );
	glBegin(GL_TRIANGLES);
	for(unsigned i=0;i<fc3.h.ntris;i++)
	{
		fc3_tri_s& t = (fc3.pt[i]);
		fc3__gl( fc3.pv[ t.a ] );
		fc3__gl( fc3.pv[ t.b ] );
		fc3__gl( fc3.pv[ t.c ] );
	}
	glEnd(); // triangles
	glEndList(); // list
	glerr();
	return 0;
}

unsigned consume_texture()
{
	unsigned htex=0;
	glGenTextures(1,&htex);
	glerr();
	return htex;
}

void fc3__create_texture( const fc3_s& fc3, unsigned& htex )
{
	if(0==htex)
		htex=consume_texture();

	glBindTexture(GL_TEXTURE_2D, htex);
	glerr();

	// now create the 2d tex object
	glTexImage2D
	(
		 GL_TEXTURE_2D		// target
		,0					// level
		,4					// number of components
		,fc3.h.cwidth		// width in pixels
		,fc3.h.cheight		// height in pixels
		,0					// border width
		,GL_BGRA_EXT		// format
		,GL_UNSIGNED_BYTE	// type
		,fc3.pc
	);
	glerr();
}

struct material_s
{
	float ma[4];
	float md[4];
	float ms[4];
	float me[4];
	float sh;

	void zero()
	{
		memset(this,0,sizeof(material_s));
	}
	void reset()
	{
		zero();
	}

	void redbook() // page 180
	{
		float fa[4]={ 1.0f,0.5f,0.8f,1.0f };
		float fd[4]={ 1.0f,0.5f,0.8f,1.0f };
		float fs[4]={ 1.0f,1.0f,1.0f,1.0f };
		float fe[4]={ 0.3f,0.2f,0.2f,0.0f };
		sh=5.0;
		memcpy(ma,fa,sizeof(ma));
		memcpy(md,fd,sizeof(md));
		memcpy(ms,fs,sizeof(ms));
		memcpy(me,fe,sizeof(me));
	}

	void mytest()
	{
		float Ns = 128;
		float fs[4]={ 0.25,0.50,0.25,1.0 };
		float fe[4]={ 0.0,0.0,0.0,1.0 };
		float fa[4]={ 0.0,0.0,0.0,0.0 };
		float fd[4]={ 0.25,0.50,0.25,1.0 };

		sh=Ns;
		memcpy(ma,fa,sizeof(ma));
		memcpy(md,fd,sizeof(md));
		memcpy(ms,fs,sizeof(ms));
		memcpy(me,fe,sizeof(me));
	}

	void mymat()
	{
		float Ns = 96.078431f;
		float fs[4]={ 0.5,0.5,0.5,1.0 };
		float fe[4]={ 0.0,0.0,0.0,1.0 };
		float fa[4]={ 1.0,1.0,1.0,1.0 };
		float fd[4]={ 0.5,0.5,0.5,1.0 };

		sh=Ns;
		memcpy(ma,fa,sizeof(ma));
		memcpy(md,fd,sizeof(md));
		memcpy(ms,fs,sizeof(ms));
		memcpy(me,fe,sizeof(me));
	}

	void select()
	{
		glDisable(GL_COLOR_MATERIAL);
		glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT		, ma );
		glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE		, md );
		glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR	, ms );
		glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION	, me );
		glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS	, sh );
	}
};
material_s mat;

void select_material(const bool v)
{
	if(false==v)
	{
		mat.zero();
	}
	else
	{
		mat.mymat();
	}
	mat.select();
}

void select_texture(const bool v, const bool wrap )
{
	if(false==v)
	{
		glDisable(GL_TEXTURE_2D);
	}
	else
	{
		glEnable(GL_TEXTURE_2D);

		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER	, GL_LINEAR		);
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER	, GL_LINEAR		);
		if(wrap)
		{
			glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S		, GL_REPEAT		);
			glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T		, GL_REPEAT		);
		}
		else
		{
			glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S		, GL_CLAMP		);
			glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T		, GL_CLAMP		);
		}
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
	}
}

void extra()
{
	// CULLING
	//glEnable		(GL_CULL_FACE		);
	//glCullFace	(GL_BACK			);

	// SHADER
	//glShadeModel	(GL_SMOOTH			); //set the shader to smooth shader
	//glShadeModel	(GL_FLAT			); //set the shader to smooth shader

	// COLOR MATERIAL
	//glEnable		(GL_COLOR_MATERIAL	);
	//glDisable		(GL_COLOR_MATERIAL	);

	// COLOR
	//glColor4f(1,1,1,1);

	// BLEND
	//glDisable		(GL_BLEND			);
}

void enable()
{
	// DEPTH
	glEnable		(GL_DEPTH_TEST		); //enable the depth testing

	// COLOR
	glColor4f(0.5,0.5,0.5,1.0);

	// MATERIAL
	//glEnable(GL_COLOR_MATERIAL);
	glDisable(GL_COLOR_MATERIAL);
	select_material(true);

	// SHADER
	glShadeModel	(GL_SMOOTH			); //set the shader to smooth shader

	// LIGHTING
	glEnable		(GL_LIGHTING);
	glLightModeli	(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE );
	glLightModeli	(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );

	const float amb=1.0f/4;
	float a[4]={amb,amb,amb,1};
	glLightModelfv  ( GL_LIGHT_MODEL_AMBIENT, a );

	// LAMP0
	glEnable		(GL_LIGHT0);

	// TEXTURE
	glBindTexture	(GL_TEXTURE_2D, the_htex );
	select_texture(true,true);
}

void on_resize( const int dx, const int dy )
{
	glViewport(0, 0, dx,dy );
	float aspect = (float)dx/dy;
	//aspect = max(aspect,1.0);

	// PERSPECTIVE
	glMatrixMode	(GL_PROJECTION);
	glLoadIdentity	();
	gluPerspective	( 45.0, aspect, 1.0, 8*(the_extent*sqrt(3.0)) );
	glMatrixMode	(GL_MODELVIEW);
	glLoadIdentity	();
}

void select_vertical_sync()
{
	const char* fname = "wglSwapIntervalEXT";
	PROC p = wglGetProcAddress(fname);
	void* pvoid = (void*)p;
	if(NULL!=pvoid)
	{
		typedef BOOL (WINAPI *pp)(int);
		pp _p = (pp)pvoid;
		BOOL bresult = (*_p)(1);
		if(!bresult)
			MessageBox(NULL,"Vertical Sync request failed","INFO",MB_OK);
	}
}

string getfc3()
{
	OPENFILENAME ofn;       // common dialog box structure
	TCHAR szFile[260]	= { 0 };       // if using TCHAR macros
	memset(&ofn,0,sizeof(ofn));
	ofn.lStructSize		= sizeof(ofn);
	ofn.hwndOwner		= NULL;
	ofn.lpstrFile		= szFile;
	ofn.nMaxFile		= sizeof(szFile);
	ofn.lpstrFilter		= "FC3 Files (*.fc3)\0*.fc3\0";
	ofn.nFilterIndex	= 1;
	ofn.lpstrFileTitle	= NULL;
	ofn.nMaxFileTitle	= 0;
	ofn.lpstrInitialDir	= "./";
	ofn.Flags			= OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
	string fn			= "";

	if(the_cmd_line > " ")
		fn=the_cmd_line;
	else
	if (GetOpenFileName(&ofn) == TRUE)
		fn=ofn.lpstrFile;

	the_cmd_line="";

	return fn;
}

void reload(const string& fn)
{
	if(fn <= " ")
		return;

	the_fc3.reset();
	int e = the_fc3.load(fn.c_str());
	if(fc3_ok != e)
	{
		char msg[128]; msg[0]=0;
		sprintf_s(msg,sizeof(msg),"FC3 loader error: %d",e);
		MessageBox(NULL,msg,"FC3VIWER LOADER ERROR",MB_OK);
		return;
	}

	the_fc3.convert_axis_to_opengl();
	the_fc3.denormalize_t(); // textures only
	if(0>=the_fc3.h.get_npix())
	{
		fc3_s skin;
		if(0<=skin.load("skin.fc3"))
		{
			the_fc3.set_image( (unsigned long int*) skin.pc, skin.h.cwidth, skin.h.cheight );
		}
		else
		{
			unsigned long int pix = 0xff7f7f7f;
			the_fc3.set_image(&pix,1,1);
		}
	}
	fc3__compile( the_fc3, the_dlist );
	fc3__create_texture( the_fc3, the_htex );
}

void reset_window_title( const char units )
{
	string fn = the_fn;
	const double mpi = 0.0254			; // meters per inch
	const double mpf = mpi * 12			; // meters per foot
	const double fpm = 1.0 / mpf		; // feet per meters
	const double ipm = 1.0 / mpi		; // inches per meters
	double rad = pow(2.0,the_fc3.h.vscale) * the_fc3.h.unitlen; // recover original bounding radius in meters
	double w; // width
	double h; // height
	double l; // length
	double ex;
	the_fc3.get_dimensions(w,h,l); // will be normalized [-1,+1]
	w*=rad; h*=rad; l*=rad; // recover original metric dimensions
	the_extent = ex = the_fc3.get_radius(); // in meters
	fc3_vecf_t vhi,vlo;
	the_fc3.find_bounds(vhi,vlo);

	ex*=rad; // want actual dimensions for title text
	vhi *= (float)rad;
	vlo *= (float)rad;

	string sunit = "meters";
	double cu=1;
	switch(units)
	{
	case 'i': cu = ipm; sunit = "inches";break;
	case 'f': cu = fpm; sunit = "feet"  ;break;
	}
	w*=cu; h*=cu; l*=cu; ex*=cu;
	vlo*=(float)cu;
	vhi*=(float)cu;

	char metrics[256]; metrics[0]=0;
	double mdiam = 2*rad;
//	sprintf_s(metrics,sizeof(metrics),"Box %gm - Size (%s) %3.1f(width) %3.1f(height) %3.1f(length) %3.1f(extent)",mdiam,sunit.c_str(), w,h,l,ex);
//	sprintf_s(metrics,sizeof(metrics),"Box %gm - Size (%s)  width:%g  height:%g  length:%g",mdiam,sunit.c_str(), w,h,l);
	sprintf_s(metrics,sizeof(metrics),"Box %gm - Size (%s)  R:%+g U:%+g B:%+g  L:%+g D:%+g F:%+g",mdiam,sunit.c_str()
		,vhi.x,vhi.y,vhi.z
		,vlo.x,vlo.y,vlo.z
	);
	char poly[256];
	sprintf_s(poly,sizeof(poly),"Poly %df/%dv",the_fc3.h.ntris,the_fc3.h.nverts);

	string title = "FC3 Viewer - ";
	title += fn;
	title += " - ";
	title += metrics;
	title += " - ";
	title += poly;
	SetWindowText( the_hwnd, title.c_str() );
}

void reload()
{
	the_fn = getfc3();
	reload(the_fn);
	reset_window_title('f');
}

void init()
{
	if(false==fc3t::is_valid_types())
	{
		// This is a fatal error
		MessageBox(NULL
			,"Failed sanity check for internal data types.\r\n"
			 "FC3Viewer was compiled incorrectly.\r\n"
			 "Press OK to gracefully exit.\r\n"
			,"FATAL ERROR"
			,MB_OK
		);
		exit(0);
	}

	select_vertical_sync();

	reload();

	enable();
}

bool is_cooldown()
{
	static double t0 = 0;
	timeb t;
	ftime(&t);
	double t1 = ((double)t.time)*1000 + t.millitm;
	double dt = t1-t0;
	const double FPS = 64.0; // assume monitor is 60hz even
	const double cusp = 1000.0 / FPS; // 15.625 nice and smooth
	if(cusp > dt)
		return false;
	t0 = t1;
	return true;
}

void display()
{
	if(false==is_cooldown())
		return;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	///////////////////////////////////////////////////

	glPushMatrix();
	glTranslated( 0,0,0.0-the_extent*3);
	static int deg=0;deg++;
	glRotatef((float)(deg%360),0,1,0);
	float tilt_degrees = 30.0f * ((tilt_mode%3)-1); 
	glRotatef(tilt_degrees,1,0,0);
	glCallList(the_dlist);
	glPopMatrix();

	///////////////////////////////////////////////////
	SwapBuffers(the_hDC);
	glFlush();

}

void process_keystroke( const char ks )
{
	char c = (char)tolower(ks);
	switch(c)
	{
	case 'f':
	case 'i':
	case 'm':
		reset_window_title(c);
		break;
	default:
		the_fn="";
		reload();
		break;
	}
}

LONG WINAPI WindowProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	switch(uMsg)
	{
	case WM_PAINT:
		static PAINTSTRUCT ps;
		display();
		BeginPaint(hWnd, &ps);
		EndPaint(hWnd, &ps);
		break;
	case WM_SIZE:
		on_resize( LOWORD(lParam), HIWORD(lParam) );
		break;
	case WM_LBUTTONDOWN:
		tilt_mode++;
		break;
	case WM_CHAR:
		process_keystroke((char)wParam);
		break;
	default:
		return (LONG)DefWindowProc(hWnd, uMsg, wParam, lParam);
		break;
	}
	return 0;
} 

HWND CreateOpenGLWindow()
{
	int						pf;
	HWND					hWnd;
	WNDCLASS				wc;
	PIXELFORMATDESCRIPTOR	pfd;

	if (!the_hInstance)
	{
		the_hInstance		= GetModuleHandle(NULL);
		wc.style			= CS_OWNDC;
		wc.lpfnWndProc		= (WNDPROC)WindowProc;
		wc.cbClsExtra		= 0;
		wc.cbWndExtra		= 0;
		wc.hInstance		= the_hInstance;
		wc.hIcon			= LoadIcon(NULL, IDI_WINLOGO);
		wc.hCursor			= LoadCursor(NULL, IDC_ARROW);
		wc.hbrBackground	= NULL;
		wc.lpszMenuName		= NULL;
		wc.lpszClassName	= "OpenGL";

		if (!RegisterClass(&wc))
		{
			MessageBox(NULL, "RegisterClass() failed", "Error", MB_OK);
			return NULL;
		}
	}

	int dx = GetSystemMetrics( SM_CXFULLSCREEN );
	int dy = GetSystemMetrics( SM_CYFULLSCREEN );
//	int dx = GetSystemMetrics( SM_CXSCREEN );
//	int dy = GetSystemMetrics( SM_CYSCREEN );


	hWnd = CreateWindow
	(
			"OpenGL"
		, "FC3 Viewer"
		, 0
			| WS_OVERLAPPEDWINDOW
			| WS_MAXIMIZE
		, 0,0,dx,dy
		, NULL, NULL, the_hInstance, NULL
	);

	if (hWnd == NULL)
	{
		MessageBox(NULL, "CreateWindow() failed", "Error", MB_OK);
		return NULL;
	}

	the_hDC = GetDC(hWnd);
 
	memset(&pfd, 0, sizeof(pfd));
	pfd.nSize		= sizeof(pfd);
	pfd.nVersion	= 1;
	pfd.dwFlags		= PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER | PFD_GENERIC_ACCELERATED;
	pfd.iPixelType	= PFD_TYPE_RGBA;
	pfd.cColorBits	= 32;

	pf = ChoosePixelFormat(the_hDC, &pfd);
	if (pf == 0)
	{
		MessageBox(NULL, "ChoosePixelFormat() failed", "Error", MB_OK); 
		return 0;
	} 
 
	if (SetPixelFormat(the_hDC, pf, &pfd) == FALSE)
	{
		MessageBox(NULL, "SetPixelFormat() failed", "Error", MB_OK);
		return 0;
	} 

	DescribePixelFormat(the_hDC, pf, sizeof(PIXELFORMATDESCRIPTOR), &pfd);

	ReleaseDC(hWnd,the_hDC);

	return hWnd;
}

// windows will put the file name argument in double quotes
string extract_file_name( LPSTR lpszCmdLine )
{
	if(NULL!=lpszCmdLine)
	{
		char* pa = lpszCmdLine;
		char* pb = pa + strlen(pa)-1;
		if(pa<pb)
		{
			if('"'==*pa)pa++;
			if('"'==*pb)*pb=0;
			return (string)pa;
		}
	}
	return "";
}

int WINAPI WinMain
(
	 _In_ HINSTANCE hCurrentInst
	,_In_opt_ HINSTANCE hPreviousInst
	,_In_ LPSTR lpszCmdLine
	,_In_ int nCmdShow
)
{
	HDC		hDC	; // device context
	HGLRC	hRC ; // opengl context
	HWND	hWnd; // window
	MSG		msg ; // message

	memset(&msg,0,sizeof(msg));
	the_prev_inst = hPreviousInst;
	the_curr_inst = hCurrentInst;
	the_cmd_line = extract_file_name(lpszCmdLine);

	hWnd = CreateOpenGLWindow();
	if (hWnd == NULL)
		return (1);

	the_hwnd = hWnd;
	hDC = GetDC(hWnd);
	hRC = wglCreateContext(hDC);
	wglMakeCurrent(hDC, hRC);
	ShowWindow(hWnd, nCmdShow=SW_MAXIMIZE);

	init(); // once

	while(IsWindow(hWnd))
	{
		while(PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
		display();
	}

	wglMakeCurrent(NULL, NULL);
	ReleaseDC(hWnd,hDC);
	wglDeleteContext(hRC);
	DestroyWindow(hWnd);

	return (int)(msg.wParam);
}
