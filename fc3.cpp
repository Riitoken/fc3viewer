// fc3.cpp

/*
class        license
type         GLORYWARE
ipname       FC3
trademark    FARCRAFT®
author       Ray Edward Bornert II
date         2020-SEP-22 TUE
royalty      Free
doc          https://docs.google.com/document/d/1xAZ-WAHxwBuu1H-LszPGuMQXLryebBYzBdrHzV2Gom0
endclass
*/

#include "fc3.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

#pragma warning( disable : 4018 ) // '<': signed/unsigned mismatch
#pragma warning( disable : 4244 ) // conversion from <big> to <small>, possible loss of data
#pragma warning( disable : 4305 ) // truncation from 'fc3t::intd' to 'double'
#pragma warning( disable : 4389 ) // '==' and '!=': signed/unsigned mismatch

namespace fc3t {

bool is_valid_types()
{
  return(true
  && ( 1 == sizeof(inta) )
  && ( 2 == sizeof(intb) )
  && ( 4 == sizeof(intc) )
  && ( 8 == sizeof(intd) )
//&& (16 == sizeof(inte) )
  && ( 1 == sizeof(unta) )
  && ( 2 == sizeof(untb) )
  && ( 4 == sizeof(untc) )
  && ( 8 == sizeof(untd) )
//&& (16 == sizeof(unte) )
//&& ( 1 == sizeof(floa) )
//&& ( 2 == sizeof(flob) )
  && ( 4 == sizeof(floc) )
  && ( 8 == sizeof(flod) )
//&& (16 == sizeof(floe) )
  );
}

}; // end namespace fc3t

void fc3_header_s::reset()
{
	memset(this,0,sizeof(fc3_header_s));
}

fc3_header_s::fc3_header_s()
{
	reset();
}

bool fc3_header_s::is_valid_file_signature()const
{
	return ( true
		&& 'F' == sigF
		&& 'C' == sigC
		&& '3' == sig3
		&& isalpha(sigv)
		);
}
void fc3_header_s::init()
{
	reset();

	sigF = 'F';
	sigC = 'C';
	sig3 = '3';
	sigv = 'a';

	xaxis = 'R';
	yaxis = 'U';
	zaxis = 'B';
	format = 'b'; // default is 16bit

	endian = (((short)('e')) << 8) + 'E';
	vscale = 0;
	tscale = 0;

	unitlen = 1; // default is one meter

	nverts = 0;
	ntris = 0;
	cwidth = 0;
	cheight=0;
}

bool fc3_header_s::is_endian_valid() const
{
	unsigned char lsb = (endian >> 0) & 0xff;
	unsigned char msb = (endian >> 8) & 0xff;
	bool valid = (msb != lsb);
	return valid;
}

bool fc3_header_s::is_endian_correct() const
{
	unsigned char lsb = (endian >> 0) & 0xff;
	unsigned char msb = (endian >> 8) & 0xff;
	bool correct = (msb > lsb);
	return correct;
}

bool fc3_header_s::is_do_endian() const
{
	return is_endian_valid() && !is_endian_correct();
}

void fc3_header_s::set_scaling_exponents( const char v, const char t )
{
	vscale=v;
	tscale=t;
}

size_t fc3_header_s::get_npix() const
{
	return (size_t)cwidth * (size_t)cheight;
}

bool fc3_header_s::is_image_hash_one() const
{
	return(true
		&& 0==cheight
		&& 1==cwidth
	);
}

bool fc3_header_s::is_image_hash() const
{
	return is_image_hash_one();
}

size_t fc3_header_s::calc_file_size() const
{
	int dz = 1 << (tolower(format) - 'a'); // data size
	size_t z = 0
		+ sizeof(*this)
		+ dz * 8 * nverts
		+ sizeof(fc3_tri_s::a) * 3 * ntris
		+ sizeof(fc3t::untc) * cwidth * cheight
		;
	return z;
}

struct unit_mcoeff_s
{
	char name[16]; // unit name
	double meters; // meters per
};

unit_mcoeff_s unit_mcoeff[32] =
{
	 {"ANGSTROM"	,(    0.0000000001     )}
	,{"MICRON"		,(    0.000001         )}
	,{"POINT"		,(    0.00035277777778 )}
	,{"INCH"		,(    0.0254           )}
	,{"HAND"		,(    0.1016           )}
	,{"SPAN"		,(    0.2286           )}
	,{"FOOT"		,(    0.3048           )}
	,{"CUBIT"		,(    0.4572           )}
	,{"YARD"		,(    0.9144           )}
	,{"METER"		,(    1                )}
	,{"FATHOM"		,(    1.8288           )}
	,{"ROD"			,(    5.0292           )}
	,{"CHAIN"		,(   20.1168           )}
	,{"FURLONG"		,(  201.168            )}
	,{"CABLE"		,(  182.88             )} //600 Feet (100 Fathoms)
	,{"CABLE_UK"	,(  185.3184           )} //608 Feet (1/10 UK nautical mile)
	,{"CABLE_US"	,(  219.456            )} //720 Feet (120 Fathoms)
	,{"MILE"		,( 1609.344            )}
	,{"NAUT"		,( 1852                )} //International unit 
	,{"NAUT_UK"		,( 1853.184            )} //British unit == 6080 FEET
	,{"NAUT_US"		,( 1853.248            )} //American unit== 6080.2 FEET (expired 1959-JUL-1)
	,{"LEAGUE"		,( 4828.032            )} //1 League = about 3 miles
	,{"AU"			,( 1459.0e+8           )} //Astronomical Unit
	,{"LY"			,(    9.2324061e+15    )} //Light Year
	,{"PARSEC"		,(    3.0857e+16       )}
};

double fc3_header_s::get_unit_meters( const string& name )
{
	for(int i=0;i<32;i++)
		if(name == unit_mcoeff[i].name)
			return unit_mcoeff[i].meters;
	return 0;
}

const char* fc3_header_s::get_unit_name( const double meters )
{
	for(int i=0;i<32;i++)
		if(meters == unit_mcoeff[i].meters)
			return unit_mcoeff[i].name;
	return "";
}

void fc3_reverse_bytes( void* v, const int nbytes )
{
	if(NULL==v)return;
	char* a = (char*)v;
	for(int i=0;i<nbytes/2;i++)
	{
		swap(a[i],a[nbytes-1-i]);
	}
}

void fc3_reverse_bytes( fc3_header_s& h )
{
	fc3_reverse_bytes( &h.endian , sizeof(h.endian ));
	fc3_reverse_bytes( &h.cwidth , sizeof(h.cwidth ));
	fc3_reverse_bytes( &h.cheight, sizeof(h.cheight));
	fc3_reverse_bytes( &h.nverts , sizeof(h.nverts ));
	fc3_reverse_bytes( &h.ntris  , sizeof(h.ntris  ));
	fc3_reverse_bytes( &h.unitlen, sizeof(h.unitlen));
}

template<typename T,typename U>
struct fc3_file_vertex_template
{
	T vx;
	T vy;
	T vz;
	T nx;
	T ny;
	T nz;
	T tx;
	T ty;

	T get_unit() const
	{
		return (T)    ( ((U)(~0)) >> 1 )      ;
	}

	void getv( fc3_vecd_t& v )
	{
		const T unit = get_unit();
		v.x = vx; v.x /= unit;
		v.y = vy; v.y /= unit;
		v.z = vz; v.z /= unit;
	}
	void setv( const fc3_vecd_t& v )
	{
		const T unit = get_unit();
		vx = (T)(v.x * unit);
		vy = (T)(v.y * unit);
		vz = (T)(v.z * unit);
	}

	void get( fc3_vec_t& v, fc3_vec_t& n, fc3_vec_t& t )
	{
		const T unit = get_unit();
		v.x = (float)vx; v.x /= unit;
		v.y = (float)vy; v.y /= unit;
		v.z = (float)vz; v.z /= unit;

		n.x = (float)nx; n.x /= unit;
		n.y = (float)ny; n.y /= unit;
		n.z = (float)nz; n.z /= unit;

		t.x = (float)tx; t.x /= unit;
		t.y = (float)ty; t.y /= unit;
		t.z = 0;
	}

	void set( const fc3_vec_t& v, const fc3_vec_t& n, const fc3_vec_t& t )
	{
		const T unit = get_unit();
		vx = (T)(v.x * unit);
		vy = (T)(v.y * unit);
		vz = (T)(v.z * unit);
		nx = (T)(n.x * unit);
		ny = (T)(n.y * unit);
		nz = (T)(n.z * unit);
		tx = (T)(t.x * unit);
		ty = (T)(t.y * unit);
	}

	void reverse_bytes()
	{
		const size_t n = sizeof(T);
		fc3_reverse_bytes(&vx,n);
		fc3_reverse_bytes(&vy,n);
		fc3_reverse_bytes(&vz,n);
		fc3_reverse_bytes(&nx,n);
		fc3_reverse_bytes(&ny,n);
		fc3_reverse_bytes(&nz,n);
		fc3_reverse_bytes(&tx,n);
		fc3_reverse_bytes(&ty,n);
	}

};

typedef fc3_file_vertex_template<signed          char, unsigned          char> fc3_vert_a;
typedef fc3_file_vertex_template<signed short     int, unsigned short     int> fc3_vert_b;
typedef fc3_file_vertex_template<signed long      int, unsigned long      int> fc3_vert_c;
typedef fc3_file_vertex_template<signed long long int, unsigned long long int> fc3_vert_d;

fc3_s::fc3_s() : pt(NULL), pc(NULL), pv(NULL)
{
}

// copy constructor
fc3_s::fc3_s( const fc3_s& rval ) :  pt(NULL), pc(NULL), pv(NULL)
{
	h = rval.h;
	heap();
	for(int i=0;i<h.nverts;i++)
		pv[i] = rval.pv[i];
	for(int i=0;i<h.ntris;i++)
		pt[i] = rval.pt[i];
	set_image( (unsigned long*)rval.pc, rval.h.cwidth, rval.h.cheight );
}

void fc3_s::reset()
{
	unheap();
	h.reset();
}

//DESTRUCTOR
fc3_s::~fc3_s()
{
	reset();
}

// a data is 1 bytes
int fc3_s::loadva( FILE* fp )
{
	fc3_vert_a v;
	for(unsigned i=0;i<h.nverts;i++)
	{
		if(!fread(&v,sizeof(v),1,fp))
			return fc3_fail_read;
		// 1 bytes so no endian processing needed
		v.get( pv[i].v, pv[i].n, pv[i].t );
	}
	return fc3_ok;
}

// b data is 2 bytes
int fc3_s::loadvb( FILE* fp, const bool doend )
{
	fc3_vert_b v;
	for(unsigned i=0;i<h.nverts;i++)
	{
		if(!fread(&v,sizeof(v),1,fp))
			return fc3_fail_read;
		// 2 bytes so endian check needed
		if(doend)
			v.reverse_bytes();
		v.get( pv[i].v, pv[i].n, pv[i].t );
	}
	return fc3_ok;
}

// c data is 4 bytes
int fc3_s::loadvc( FILE* fp, const bool doend )
{
	fc3_vert_c v;
	for(unsigned i=0;i<h.nverts;i++)
	{
		if(!fread(&v,sizeof(v),1,fp))
			return fc3_fail_read;
		if(doend)
			v.reverse_bytes();
		v.get( pv[i].v, pv[i].n, pv[i].t );
	}
	return fc3_ok;
}

// d data is 8 bytes
int fc3_s::loadvd( FILE* fp, const bool doend )
{
	fc3_vert_d v;
	for(unsigned i=0;i<h.nverts;i++)
	{
		if(!fread(&v,sizeof(v),1,fp))
			return fc3_fail_read;
		if(doend)
			v.reverse_bytes();
		v.get( pv[i].v, pv[i].n, pv[i].t );
	}
	return fc3_ok;
}

#if 0
// e data is 16 bytes
int fc3_s::loadvd( FILE* fp ){}
#endif

// a data is 1 bytes
int fc3_s::saveva( FILE* fp )
{
	fc3_vert_a v;
	for(unsigned i=0;i<h.nverts;i++)
	{
		v.set( pv[i].v, pv[i].n, pv[i].t );
		if(!fwrite(&v,sizeof(v),1,fp))
			return fc3_fail_write;
	}
	return fc3_ok;
}

// b data is 2 bytes
int fc3_s::savevb( FILE* fp )
{
	fc3_vert_b v;
	for(unsigned i=0;i<h.nverts;i++)
	{
		v.set( pv[i].v, pv[i].n, pv[i].t );
		if(!fwrite(&v,sizeof(v),1,fp))
			return fc3_fail_write;
	}
	return fc3_ok;
}

// c data is 4 bytes
int fc3_s::savevc( FILE* fp )
{
	fc3_vert_c v;
	for(unsigned i=0;i<h.nverts;i++)
	{
		v.set( pv[i].v, pv[i].n, pv[i].t );
		if(!fwrite(&v,sizeof(v),1,fp))
			return fc3_fail_write;
	}
	return fc3_ok;
}

// d data is 8 bytes
int fc3_s::savevd( FILE* fp )
{
	fc3_vert_d v;
	for(unsigned i=0;i<h.nverts;i++)
	{
		v.set( pv[i].v, pv[i].n, pv[i].t );
		if(!fwrite(&v,sizeof(v),1,fp))
			return fc3_fail_write;
	}
	return fc3_ok;
}

#if 0
// e data is 16 bytes
int fc3_s::saveve( FILE* fp ){}
#endif

int fc3_s::loadh( FILE* fp )
{
	h.reset();
	size_t z = sizeof(h);
	size_t n = fread( &h, z,1, fp );
	if(1!=n)
	{
		return fc3_fail_read;
	}
	if(false==h.is_valid_file_signature())
	{
		return fc3_bad_header;
	}

	/*
	if(h.is_do_endian())
	{
		fc3_reverse_bytes(h);
	}
	*/
	return fc3_ok;
}

int fc3_s::loadv(FILE* fp, const bool doend)
{
	switch(h.format)
	{
	case 'a': return loadva(fp); break;
	case 'b': return loadvb(fp,doend); break;
	case 'c': return loadvc(fp,doend); break;
	case 'd': return loadvd(fp,doend); break;
	}
	return fc3_bad_header;
}
int fc3_s::savev(FILE* fp)
{
	// we always write our native endian
	switch(h.format)
	{
	case 'a': return saveva(fp); break;
	case 'b': return savevb(fp); break;
	case 'c': return savevc(fp); break;
	case 'd': return savevd(fp); break;
	}
	return fc3_bad_header;
}

int fc3_s::load( const char* fn )
{
	size_t n=0;
	if(NULL==fn)
	{
		return fc3_bad_arg;
	}

	FILE* fp = NULL;
	fp = fopen(fn,"rb");
	if(NULL==fp)
	{
		return fc3_fail_open;
	}
	f_s f(fp);

	// HEADER
	if(0>loadh(fp))
	{
		return fc3_bad_header;
	}

	bool doend = h.is_do_endian();
	if(doend)
	{
		fc3_reverse_bytes(h);
	}

	// HEAP
	if(false==heap())
	{
		unheap();
		return fc3_fail_heap;
	}

	// VERTS
	if(0>loadv(fp,doend))
	{
		unheap();
		return fc3_fail_heap;
	}

	// TRIS
	n = fread( &pt[0], sizeof(pt[0]), h.ntris, fp );
	if(n != h.ntris)
	{
		unheap();
		return fc3_fail_read;
	}
	if(doend)
		for(int i=0;i<h.ntris;i++)
			pt[i].reverse_bytes();

	// IMAGE
	size_t npix = h.cheight * h.cwidth;
	if(0!=npix)
	{
		n = fread( &pc[0], sizeof(pc[0]), npix, fp );
		if(n != npix)
		{
			unheap();
			return fc3_fail_read;
		}
	}
	if(doend)
		for(int i=0;i<npix;i++)
			fc3_reverse_bytes(&pc[i],sizeof(pc[0]));

	return fc3_ok;
}

int fc3_s::save( const char* fn )
{
	if(NULL==fn)
	{
		return fc3_bad_arg;
	}

	FILE* fp = fopen(fn,"wb");

	if(NULL==fp)
	{
		return fc3_fail_open;
	}
	f_s f(fp);
	size_t n=0;

	////////////
	normalize();
	////////////

	n = fwrite( &h, sizeof(h), 1, fp );
	if(1!=n)
	{
		return fc3_fail_write;
	}

	if(0>savev(fp))
	{
		return fc3_fail_write;
	}

	n = fwrite( &pt[0], sizeof(pt[0]), h.ntris, fp );
	if(n!=h.ntris)
	{
		return fc3_fail_write;
	}
	size_t npix = h.cheight*h.cwidth;

	if(h.is_image_hash_one())
		npix=1*2;

	if(NULL!=pc && 0<npix)
	{
		n = fwrite( &pc[0], sizeof(pc[0]), npix, fp );
		if(n!=npix)
		{
			return fc3_fail_write;
		}
	}

	fclose(fp);

	return fc3_ok;
}


void fc3_s::set_scaling_exponents( const signed char v, const signed char t )
{
	h.set_scaling_exponents(v,t);
}

void fc3_s::set_vertex_radius( const double r )
{
	int i; for(i= -128; i <= +126; i++) if(pow(2.0,i)>=r)break; h.vscale = (char)i;
}
void fc3_s::set_texture_radius( const double r )
{
	int i; for(i= -128; i <= +126; i++) if(pow(2.0,i)>=r)break; h.tscale = (char)i;
}

void fc3_s::find_bounds( fc3_vec_t& vhi, fc3_vec_t& vlo )
{
	vlo.zero();
	vhi.zero();
	if(0==h.nverts)return;
	vhi=vlo=pv[0].v;
	for(unsigned i=0;i<h.nverts;i++)
	{
		pv[i].v.maxv(vhi);
		pv[i].v.minv(vlo);
	}
}

void fc3_s::find_bounds( double& vhi, double& thi )
{
	vhi=0;
	thi=0;
	double a;
	for(unsigned i=0;i<h.nverts;i++)
	{
		a=pv[i].v.maxabs();if(a>vhi)vhi=a;
		a=pv[i].t.maxabs();if(a>thi)thi=a;
	}
}

void fc3_s::calculate_scaling_exponents()
{
	double vhi=0;
	double thi=0;
	find_bounds(vhi,thi);

	//double ve = ceilf(log(vhi)/log(2.0));
	//double te = ceilf(log(thi)/log(2.0));

	set_vertex_radius(vhi);
	set_texture_radius(thi);
}

void fc3_s::denormalize( const bool to_meters )
{
	if(1==h.unitlen)
	if(0==h.vscale && 0==h.tscale)
		return;

	float vmul = (float)pow(2.0,h.vscale);
	float tmul = (float)pow(2.0,h.tscale);

	if(true == to_meters)
	{
		if(0.0 != h.unitlen) // only if non-zero unitlen
		if(1.0 != h.unitlen) // only if not already meters
			vmul *= (float)h.unitlen;
		h.unitlen=1.0; // select meters
	}

	for(unsigned i=0;i<h.nverts;i++)
	{
		pv[i].v*=vmul;
		pv[i].t*=tmul;
		pv[i].n.normalize();
	}
	h.vscale=0; // 2^0 == 1
	h.tscale=0; // 2^0 == 1

}

void fc3_s::denormalize_t()
{
	if(0==h.tscale)
		return;

	float tmul = (float)pow(2.0,h.tscale);

	for(unsigned i=0;i<h.nverts;i++)
	{
		pv[i].t*=tmul;
	}
	h.tscale=0; // 2^0 == 1
}

void fc3_s::normalize()
{
	denormalize();
	calculate_scaling_exponents();
	float vmul = (float)(1.0 / pow(2.0,h.vscale));
	float tmul = (float)(1.0 / pow(2.0,h.tscale));
	for(unsigned i=0;i<h.nverts;i++)
	{
		pv[i].v*=vmul;
		pv[i].t*=tmul;
		pv[i].n.normalize();
	}
}

void fc3_s::unheap()
{
	if(NULL!=pt) delete pt;
	if(NULL!=pc) delete pc;
	if(NULL!=pv) delete pv;
	pt=NULL;
	pc=NULL;
	pv=NULL;
}

/*
void fc3_s::setv( const size_t index, const fc3_vert_b& v )
{
	pv[index]=v;
}
*/

/*
void fc3_s::sett( const size_t index, const fc3_tri_s t )
{
	pt[index]=t;
}
*/

struct axis_convert_s
{
	int what[3]; // expected 0 and 1 and 2 in any of 6 orders
	int sign[3]; // expected +1 or -1

	void reset()
	{
		what[0]=0;
		what[1]=1;
		what[2]=2;
		sign[0]=sign[1]=sign[2]= +1;
	}
	axis_convert_s()
	{
		reset();
	}

	bool is_valid() const
	{
		int bits=0;
		for(int i=0;i<3;i++)
		{
			if( what[i] < 0 || 2 < what[i] )
				return false;
			if( (+1 != sign[i]) && (-1 != sign[i]) )
				return false;
			bits |= (1<<what[i]);
		}
		if(7!=bits)
			return false;
		return true;
	}
};


/*
	This utility function takes a before F and after T axis string
	each denoting an axis orientation using the letters 
		R L U D B F right left up down back fore
	the resulting ac struct will contain a simple algorithm for axis order and sign change
	the ac is used by conver_axis below
*/
int fc3_s::calculate_convert_axis( const string& F, const string& T, axis_convert_s& ac )
{
	ac.reset();

	if(F==T)
		return fc3_ok;

	// expecting only axis characters ( R L U D F B ) right left up down fore back
	// converting to XYZ and xyz
	string f,t;
	if(false==normalize_axis_string(F,f))return -1;
	if(false==normalize_axis_string(T,t))return -2;

	// axis strings are now XYZ positive, xyz negative
	// char case is sign

	// this super-quick look up table works if and only if the above strings are valid ascii printables
	char cf[256]; memset(cf,0,sizeof(cf));
	cf[ toupper(f[0]) ] = 0;
	cf[ toupper(f[1]) ] = 1;
	cf[ toupper(f[2]) ] = 2;
	cf[ tolower(f[0]) ] = 0;
	cf[ tolower(f[1]) ] = 1;
	cf[ tolower(f[2]) ] = 2;

	ac.what[ 0 ] = cf[ t[0] ];
	ac.what[ 1 ] = cf[ t[1] ];
	ac.what[ 2 ] = cf[ t[2] ];

	// change sign when axis letters not same case
	if(t[0] != f[ ac.what[0] ]) ac.sign[0] = -1;
	if(t[1] != f[ ac.what[1] ]) ac.sign[1] = -1;
	if(t[2] != f[ ac.what[2] ]) ac.sign[2] = -1;

	return fc3_ok;
}

/*
	caller must have already correctly calculated the axis_convert_s
	function will fail if requested conversion is invalid
*/
int fc3_s::convert_axis( const axis_convert_s& ac )
{
	if(false == ac.is_valid())
		return -1;

	double f[3];
	double t[3];

	// count the number of signs that are negating
	int nneg=0;
	for(int i=0;i<3;i++)
		if (0 >ac.sign[i])
			nneg++;
	int nequ=0;
	for(int i=0;i<3;i++)
		if (i==ac.what[i])
			nequ++;
	if(2==nequ)
		return -2;

	// first the swaps
	for(unsigned i=0;i<h.nverts;i++)
	{
		f[0]=pv[i].v.x;
		f[1]=pv[i].v.y;
		f[2]=pv[i].v.z;

		for(int k=0;k<3;k++)
			t[k] = f[ac.what[k]]*ac.sign[k];

		pv[i].v.x=(float)t[0];
		pv[i].v.y=(float)t[1];
		pv[i].v.z=(float)t[2];
	}

	// then the sign changes
	for(unsigned i=0;i<h.nverts;i++)
	{
		f[0]=pv[i].n.x;
		f[1]=pv[i].n.y;
		f[2]=pv[i].n.z;

		for(int k=0;k<3;k++)
			t[k] = f[ac.what[k]]*ac.sign[k];

		pv[i].n.x=(float)t[0];
		pv[i].n.y=(float)t[1];
		pv[i].n.z=(float)t[2];
	}

	#if 1
	// then tri order changes if nneg is odd
	bool handchange=false;
	if(0==nequ || 3==nequ) handchange = (nneg&1); else
	if(1==nequ) handchange = !(nneg&1);

	if(handchange)
	for(unsigned i=0;i<h.ntris;i++)
	{
		swap( pt[i].a, pt[i].b );
	}
	#endif

	return fc3_ok;
}

int fc3_s::convert_axis_to( const string& t )
{
	string f; f+=h.xaxis; f+=h.yaxis; f+=h.zaxis;
	axis_convert_s ac;
	if(0>calculate_convert_axis( f,t, ac ))
		return -1;
	if(0>convert_axis(ac))
		return -2;
	h.xaxis=t[0];
	h.yaxis=t[1];
	h.zaxis=t[2];
	return fc3_ok;
}

int fc3_s::unit_test()
{
	fc3_s f;
	f.h.init();
	f.h.nverts=1;
	f.heap();

	f.pv[0].v.x=1;
	f.pv[0].v.y=2;
	f.pv[0].v.z=3;

	fc3_vec_t vsave = f.pv[0].v;
	fc3_vec_t nsave = f.pv[0].n;

	string F;
	F += f.h.xaxis;
	F += f.h.yaxis;
	F += f.h.zaxis;

	// make all 48 possible
	string sax[48]; int nsax=0;
	char cax[4];
	for(int i='A';i<='Z';i++)
	for(int j='A';j<='Z';j++)
	for(int k='A';k<='Z';k++)
	{
		cax[0]=i;
		cax[1]=j;
		cax[2]=k;
		cax[3]=0;
		string T = cax;
		if(false == is_valid_axis_string(T))
			continue;

		sax[nsax++]=T;
	}
	for(int i=0;i<48;i++)
	for(int j=0;j<48;j++)
	{
		f.convert_axis_to( sax[i] );
		f.convert_axis_to( sax[j] );
		f.convert_axis_to( F );
		if(f.pv[0].v != vsave)
		{
			return fc3_fail_test;
		}
		if(f.pv[0].n != nsave)
		{
			return fc3_fail_test;
		}
	}

	bool isok = fc3t::is_valid_types();
	if(false == isok)
		return fc3_fail_test;

	return fc3_ok;
}


bool fc3_s::is_valid_axis_string( const string& s )
{
	if(3 != s.length())
		return false;
	int bits=0;
	for(int i=0;i<3;i++)
	{
		char c = toupper(s[i]);
		switch(c)
		{
		case 'R': // right
		case 'L': // left
			bits |= 1; break;

		case 'U': // up
		case 'D': // down
			bits |= 2; break;

		//case 'N': // north
		//case 'S': // south
		case 'F': // fore
		case 'B': // back
			bits |= 4; break;

		default:
			return false;
		}
	}
	if(7!=bits)
		return false;
	return true;
}

bool fc3_s::normalize_axis_string( const string& F, string& f )
{
	if(false==is_valid_axis_string(F))
		return false;

	f="";

	for(int i=0;i<3;i++)
	{
		switch(toupper(F[i]))
		{
		case 'R': f+='X'; break; // right
		case 'L': f+='x'; break; // left
		case 'U': f+='Y'; break; // up
		case 'D': f+='y'; break; // down
		case 'B': f+='Z'; break; // back
		case 'F': f+='z'; break; // fore
		//case 'S': f+='Z'; break; // south
		//case 'N': f+='z'; break; // north
		}
	}
	return true;
}


bool fc3_s::heap()
{
	unheap();
	pv = new fc3_vert_s[ h.nverts ];
	pt = new fc3_tri_s[ h.ntris ];
	if(h.is_image_hash_one())
		pc = new unsigned int[ 1*2 ]; // 64bit hash value
	else
		pc = new unsigned int[ h.cheight * h.cwidth ];
	return ( true
		&& NULL!=pv
		&& NULL!=pt
		&& NULL!=pc
		);
}
bool fc3_s::heap( const size_t ntris, const size_t cw, const size_t ch )
{
	unheap();
	h.ntris=(unsigned int)ntris;
	h.nverts=(unsigned int)ntris*3;
	h.cheight=(unsigned short)ch;
	h.cwidth=(unsigned short)cw;
	return heap();
}

bool fc3_s::set_image( const unsigned long int* pbits, const size_t cw, const size_t ch )
{
	if(NULL!=pc)
		delete pc;
	h.cheight=(unsigned short)ch;
	h.cwidth=(unsigned short)cw;
	size_t npix = cw*ch;
	pc=new unsigned int[ npix ];
	if(NULL==pc)
		return false;
	size_t nbytes = npix * sizeof(pc[0]);
	memcpy( pc,pbits,nbytes );
	return true;
}

bool fc3_s::set_image( const unsigned int* pbits, const size_t cw, const size_t ch )
{
	return set_image( (const unsigned long int*)pbits, cw,ch );
}

void fc3_s::get_dimensions( double& width, double& height, double& length ) const
{
	width=height=length=0;
	double x=0;
	double X=0;
	double y=0;
	double Y=0;
	double z=0;
	double Z=0;
	for(unsigned i=0;i<h.nverts;i++)
	{
		const double xx = pv[i].v.x;
		const double yy = pv[i].v.y;
		const double zz = pv[i].v.z;
		if(xx<x)x=xx; else if(xx>X)X=xx;
		if(yy<y)y=yy; else if(yy>Y)Y=yy;
		if(zz<z)z=zz; else if(zz>Z)Z=zz;
	}
	width	= X-x;
	height	= Y-y;
	length	= Z-z;
}

double fc3_s::get_extent() const
{
	double ex=0;
	for(unsigned i=0;i<h.nverts;i++)
	{
		double v;
		v = fabs(pv[i].v.x); if(v>ex)ex=v;
		v = fabs(pv[i].v.y); if(v>ex)ex=v;
		v = fabs(pv[i].v.z); if(v>ex)ex=v;
	}
	return ex;
}

double fc3_s::get_radius() const
{
	double ex=0;
	for(unsigned i=0;i<h.nverts;i++)
	{
		double v = pv[i].v.length();
		if(v>ex)ex=v;
	}
	return ex;
}

/*
	All lossy functions expect normalized data
	where all values are within [-1.0,+1.0].
*/

using namespace fc3t;
double get_lossy_a( const double V) { const double one= (intd)(              0x7f); double v=V*one; inta i=(inta)v; v=(double)i; v/=one; return fabs(V-v); }
double get_lossy_b( const double V) { const double one= (intd)(            0x7fff); double v=V*one; intb i=(intb)v; v=(double)i; v/=one; return fabs(V-v); }
double get_lossy_c( const double V) { const double one= (intd)(        0x7fffffff); double v=V*one; intc i=(intc)v; v=(double)i; v/=one; return fabs(V-v); }
double get_lossy_d( const double V) { const double one= (intd)(0x7fffffffffffffff); double v=V*one; intd i=(intd)v; v=(double)i; v/=one; return fabs(V-v); }

void fc3_s::calc_lossy_v( double& a, double& b, double& c, double& d )
{
	a=b=c=d=0;
	double dv=0;

	for(unsigned i=0;i<h.nverts;i++)
	{
		dv=get_lossy_a(pv[i].v.x); if(dv>a)a=dv;
		dv=get_lossy_b(pv[i].v.x); if(dv>b)b=dv;
		dv=get_lossy_c(pv[i].v.x); if(dv>c)c=dv;
		dv=get_lossy_d(pv[i].v.x); if(dv>d)d=dv;

		dv=get_lossy_a(pv[i].v.y); if(dv>a)a=dv;
		dv=get_lossy_b(pv[i].v.y); if(dv>b)b=dv;
		dv=get_lossy_c(pv[i].v.y); if(dv>c)c=dv;
		dv=get_lossy_d(pv[i].v.y); if(dv>d)d=dv;

		dv=get_lossy_a(pv[i].v.z); if(dv>a)a=dv;
		dv=get_lossy_b(pv[i].v.z); if(dv>b)b=dv;
		dv=get_lossy_c(pv[i].v.z); if(dv>c)c=dv;
		dv=get_lossy_d(pv[i].v.z); if(dv>d)d=dv;
	}
}

#if 0
void fc3_s::calc_lossy_n( double& a, double& b, double& c, double& d )
{
	a=b=c=d=0;

	for(unsigned i=0;i<h.nverts;i++)
	{
		if(is_lossy_a(pv[i].n.x)) a++;
		if(is_lossy_b(pv[i].n.x)) b++;
		if(is_lossy_c(pv[i].n.x)) c++;
		if(is_lossy_d(pv[i].n.x)) d++;

		if(is_lossy_a(pv[i].n.y)) a++;
		if(is_lossy_b(pv[i].n.y)) b++;
		if(is_lossy_c(pv[i].n.y)) c++;
		if(is_lossy_d(pv[i].n.y)) d++;

		if(is_lossy_a(pv[i].n.z)) a++;
		if(is_lossy_b(pv[i].n.z)) b++;
		if(is_lossy_c(pv[i].n.z)) c++;
		if(is_lossy_d(pv[i].n.z)) d++;
	}
}

void fc3_s::calc_lossy_t( int& a, int& b, int& c, int& d )
{
	a=b=c=d=0;

	for(unsigned i=0;i<h.nverts;i++)
	{
		if(is_lossy_a(pv[i].t.x)) a++;
		if(is_lossy_b(pv[i].t.x)) b++;
		if(is_lossy_c(pv[i].t.x)) c++;
		if(is_lossy_d(pv[i].t.x)) d++;

		if(is_lossy_a(pv[i].t.y)) a++;
		if(is_lossy_b(pv[i].t.y)) b++;
		if(is_lossy_c(pv[i].t.y)) c++;
		if(is_lossy_d(pv[i].t.y)) d++;
	}
}
#endif

#if 0

/*
void debug_dups( const fc3_s& fc3 )
{
	int cnt=0;
	for(unsigned i=0;i<fc3.h.nverts-1;i++)
	{
		cnt=0;
		for(unsigned k=i+1;k<fc3.h.nverts-0;k++)
		{
			fc3_vec_t v = fc3.pv[i].v - fc3.pv[k].v;
			double dist = v.length();
			const double cusp = 0;
			if(cusp >= dist)
				cnt++;
		}
		if(cnt)
			printf("dup count = %d\n",cnt);
	}
}
*/

/*
This function expects shared vertices
*/
void fc3_s::calc_normals()
{
	debug_dups(*this);
	// now the normals
	fc3_vec_t a,b,c,xf,xt,xp;
	int nv = h.nverts;
	int nt = h.ntris;
	for(int i=0;i<nv;i++)
	{
		fc3_vec_t sum; sum.zero();
		// for each vertex
		for(int j=0;j<nt;j++)
		if(pt[j].isv(i)) // triangle j references vertex i
		{
			fc3_tri_s t = pt[j];
			a = pv[t.a].v;
			b = pv[t.b].v;
			c = pv[t.c].v;
			fc3_vec_t norm = t.calc_norm(i,a,b,c);
			sum+=norm;
		}
		sum.normalize();
		pv[i].n=sum;
	}
}
#endif


static int _cdecl fc3_vert_comp( const void* aa, const void* bb )
{
	fc3_vert_s& a= *(fc3_vert_s*)aa;
	fc3_vert_s& b= *(fc3_vert_s*)bb;
	fc3_vec_t& va = a.v;
	fc3_vec_t& vb = b.v;
	if(va.x < vb.x) return -1; else
	if(va.x > vb.x) return +1; else
	if(va.y < vb.y) return -1; else
	if(va.y > vb.y) return +1; else
	if(va.z < vb.z) return -1; else
	if(va.z > vb.z) return +1; else
	{
		fc3_vec_t& ta = a.t;
		fc3_vec_t& tb = b.t;
		if(ta.x < tb.x) return -1; else
		if(ta.x > tb.x) return +1; else
		if(ta.y < tb.y) return -1; else
		if(ta.y > tb.y) return +1; else
		if(ta.z < tb.z) return -1; else
		if(ta.z > tb.z) return +1;
	}
	return 0;
}

int fc3_s::sort_verts()
{
	int nv = h.nverts;
	int* pi = new int[nv];
	if(NULL==pi)
		return -1;
	for(int i=0;i<nv;i++)
		pi[i] = -1;

	// co-opting the z text coord as convenient storage
	for(int i=0;i<nv;i++)
		pv[i].t.z=i;
	qsort( pv, nv, sizeof(pv[0]), fc3_vert_comp );

	for(int i=0;i<nv;i++)
		pi[ (int)(pv[i].t.z) ] = i;

	for(int i=0;i<nv;i++)
	if(0>pi[i])
	{
		return -2;
	}

	int nt = h.ntris;
	for(int i=0;i<nt;i++)
	{
		pt[i].a = pi[ pt[i].a ];
		pt[i].b = pi[ pt[i].b ];
		pt[i].c = pi[ pt[i].c ];
	}

	for(int i=0;i<nv;i++)
		pv[i].t.z=0;

	delete pi;

	return 0;
}

int fc3_s::undup()
{
	sort_verts();

	int* pli = new int[h.nverts+1];
	if(NULL==pli)
		return -1;

	for(int i=0;i<h.nverts;i++)
		pli[i]=i;

	int ndup=0;
	int base=0;
	int tail=0;
	int head=0;
	for(tail=0; tail<h.nverts ;tail=head)
	{
		ndup=1; // self is dup of self so solo dup will get repointed below
		for(head=tail+1;head<h.nverts && (pv[base].same_dup( pv[head]));head++,ndup++);
		// head is on new vert
		for(int i=0;i<ndup;i++) // solo dup will get repointed
			pli[tail+i]=base;

		pv[++base] = pv[head]; // copy diff vert to next base slot
	}

	for(int i=0;i<h.ntris;i++)
	{
		//  new            old
		pt[i].a = pli[ pt[i].a ]; // translate
		pt[i].b = pli[ pt[i].b ]; // translate
		pt[i].c = pli[ pt[i].c ]; // translate
	}

	delete pli;

	// ok now base is the number of new verts
	return (h.nverts=base);
}

void fc3_s::calc_normals()
{
	for(int i=0;i<h.nverts;i++)
		pv[i].n.zero();

	fc3_vec_t a,b,c;
	for(int i=0;i<h.ntris;i++)
	{
		fc3_tri_s t = pt[i];
		a = pv[t.a].v;
		b = pv[t.b].v;
		c = pv[t.c].v;
		pv[t.a].n += t.calc_norm(t.a,a,b,c);
		pv[t.b].n += t.calc_norm(t.b,a,b,c);
		pv[t.c].n += t.calc_norm(t.c,a,b,c);
	}

	for(int i=0;i<h.nverts;i++)
		pv[i].n.normalize();

	sort_verts();
	int head=0;
	int tail=0;
	for(tail=0;tail<h.nverts;tail++)
	{
		fc3_vec_t sum;
		for(head=tail; pv[head].v == pv[tail].v; head++)
			sum += pv[head].n;
		sum.normalize();
		for(int k=head-1;k>=tail;k--)
			pv[k].n=sum;
	}

}

void fc3_s::negate_normals()
{
	// now the normals
	fc3_vec_t a,b,c,xf,xt,xp;
	int nv = h.nverts;
	for(int i=0;i<nv;i++)
	{
		pv[i].n *= -1;
	}
}

/*

*/
fc3_error_t fc3_s::make_tetrahedron()
{
	reset();
	h.init();
	h.nverts=4;
	h.ntris=4;
	if(false==heap())
		return fc3_fail_heap;

	const double q = 1.0 / sqrt(2.0);

	//0 RIGHT
	pv[0].v.x = +1;
	pv[0].v.y = +q;
	pv[0].v.z = +0;

	//1 LEFT
	pv[1].v.x = -1;
	pv[1].v.y = +q;
	pv[1].v.z = +0;

	//2 FRONT
	pv[2].v.x = +0;
	pv[2].v.y = -q;
	pv[2].v.z = -1;

	//3 BACK
	pv[3].v.x = +0;
	pv[3].v.y = -q;
	pv[3].v.z = +1;

	for(int i=0;i<h.nverts;i++)
	{
		pv[i].t.x=0.5;
		pv[i].t.y=0.5;
		pv[i].t.z=0.0;

		pv[i].n = pv[i].v;
		pv[i].n.normalize();
	}

	// FRONT
	pt[0].a = 2; // front
	pt[0].b = 1; // left
	pt[0].c = 0; // right

	// BACK
	pt[1].a = 3; // back
	pt[1].b = 0; // right
	pt[1].c = 1; // left

	// LEFT
	pt[2].a = 1; // left
	pt[2].b = 2; // front
	pt[2].c = 3; // back

	// RIGHT
	pt[3].a = 0; // right
	pt[3].b = 3; // back
	pt[3].c = 2; // front

	return fc3_ok;
}

fc3_vert_s calc_average( const fc3_vert_s& a, const fc3_vert_s& b )
{
	fc3_vert_s c;
	c.v = (a.v + b.v) * 0.5;
	c.t = (a.t + b.t) * 0.5;
	c.n = (a.n + b.n) * 0.5;
	return c;
}

void fc3_s::add_face( const fc3_vert_s& va, const fc3_vert_s& vb, const fc3_vert_s& vc )
{
	int a = h.nverts+0;
	int b = h.nverts+1;
	int c = h.nverts+2;
	h.nverts+=3;

	pv[a]=va;
	pv[b]=vb;
	pv[c]=vc;

	int n = h.ntris;
	pt[n].a = a;
	pt[n].b = b;
	pt[n].c = c;
	h.ntris+=1;

}

void fc3_s::add_quad( const fc3_vert_s& va, const fc3_vert_s& vb, const fc3_vert_s& vc, const fc3_vert_s& vd )
{
	fc3_vert_s vo;
	vo += va;
	vo += vb;
	vo += vc;
	vo += vd;

	vo.v *= 0.25;
	vo.n.normalize();
	vo.t *= 0.25;

	add_face(vo,va,vb);
	add_face(vo,vb,vc);
	add_face(vo,vc,vd);
	add_face(vo,vd,va);
}


// this will increase the number of triangle faces by a factor of 6
// nverts will be current faces * 4 * 3
fc3_error_t fc3_s::split_faces()
{
	// save a copy
	fc3_s fc3=(*this);

	reset();

	h.init();
	int NT = h.ntris = fc3.h.ntris*4;
	int NV = h.nverts = fc3.h.ntris*4*3;
	if(false==heap())
		return fc3_fail_heap;

	// these are now index counters
	h.nverts = h.ntris = 0;

	for(int i=0;i<fc3.h.ntris;i++)
	{
		fc3_tri_s t = fc3.pt[i];

		fc3_vert_s a = fc3.pv[t.a];
		fc3_vert_s b = fc3.pv[t.b];
		fc3_vert_s c = fc3.pv[t.c];

		fc3_vert_s aa = calc_average(a,b);
		fc3_vert_s bb = calc_average(b,c);
		fc3_vert_s cc = calc_average(c,a);

		add_face(a,aa,cc); // top
		add_face(b,bb,aa); // left
		add_face(c,cc,bb); // right
		add_face(aa,bb,cc); // middle

	}

	if(h.nverts != NV)
	{
		assert(false);
	}
	if(h.ntris != NT)
	{
		assert(false);
	}

	fc3.reset();

	return fc3_ok;
}

fc3_error_t fc3_s::make_sphere( const int nsplits )
{
	fc3_error_t e;
	e = make_tetrahedron();
	if(fc3_ok != e)
		return e;

	for(int i=0;i<nsplits;i++)
	{
		e = split_faces();
		if(fc3_ok != e)
			return e;
	}

	for(int i=0;i<h.nverts;i++)
	{
		pv[i].v.normalize();
		pv[i].n=pv[i].v;
	}

	// coordinates
	fc3_vec_t c; c.x=0.5; c.y=0.5; c.z=0;
	for(int i=0;i<h.nverts;i++)
	{
		fc3_vert_s& v = pv[i];
		v.t.x = 0.5f + v.v.x/2;
		v.t.y = 0.5f + v.v.y/2;
		v.t.z = 0;
		if(0>v.v.z) // front
		{
			v.t.y /= 2;
			v.t.y += 0.0;
		}
		else // back
		{
			v.t.y /= 2;
			v.t.y += 0.5;
		}
	}

	return fc3_ok;
}

fc3_error_t fc3_s::make_cube( const int nsplits )
{
	int diam = (1<<nsplits); // nquads
	int facequads = diam*diam;
	int facetris = facequads*4;
	int faceverts = facetris*3;

	reset();
	h.init();
	h.ntris = facetris*6;
	h.nverts = faceverts*6;

	if(false==heap())
		return fc3_fail_heap;

	int ntris=h.ntris;
	int nverts=h.nverts;
	h.ntris=h.nverts=0;

	fc3_vert_s a,b,c,d;

	float fstep = 1.0f / (float)diam;

	// standard algebraic cartesian layout y axis is vertical, x axis is horizontal ... origin lower left

	// NATURAL
	// negx face
	a.zero(); a.n.x = -1; d=c=b=a;
	for(int j=0;j<diam;j++) // rows Y
	for(int i=0;i<diam;i++) // cols Z
	{
		a.t.y = fstep*j; a.t.x = fstep*i;
		a.v.y = fstep*j; a.v.z = fstep*i;
		b=a; b.v.z += fstep; b.t.x += fstep;
		c=b; c.v.y += fstep; c.t.y += fstep;
		d=c; d.v.z -= fstep; d.t.x -= fstep;
		add_quad(a,b,c,d);
	}

	// NATURAL
	// negy face
	a.zero(); a.n.y = -1; d=c=b=a;
	for(int j=0;j<diam;j++) // rows Z
	for(int i=0;i<diam;i++) // cols X
	{
		a.t.y = fstep*j; a.t.x = fstep*i;
		a.v.z = fstep*j; a.v.x = fstep*i;
		b=a; b.v.x += fstep; b.t.x += fstep;
		c=b; c.v.z += fstep; c.t.y += fstep;
		d=c; d.v.x -= fstep; d.t.x -= fstep;
		add_quad(a,b,c,d);
	}

	// NATURAL
	// negz face
	a.zero(); a.n.z = -1; d=c=b=a;
	for(int j=0;j<diam;j++) // rows Y
	for(int i=0;i<diam;i++) // cols X
	{
		a.t.y = fstep*j; a.t.x = fstep*i;
		a.v.y = fstep*j; a.v.x = fstep*i;
		b=a; b.v.x += fstep; b.t.x += fstep;
		c=b; c.v.y += fstep; c.t.y += fstep;
		d=c; d.v.x -= fstep; d.t.x -= fstep;
		add_quad(d,c,b,a);
	}

	// opposite verts
	for(int i=0;i<h.nverts;i++)pv[h.nverts+i]=pv[i]; // copy verts

	for(int i=0;i<faceverts;i++){int at=h.nverts+0*faceverts+i;pv[at].v.x = +1;pv[at].n*=-1;} // posx fixup
	for(int i=0;i<faceverts;i++){int at=h.nverts+1*faceverts+i;pv[at].v.y = +1;pv[at].n*=-1;} // posy fixup
	for(int i=0;i<faceverts;i++){int at=h.nverts+2*faceverts+i;pv[at].v.z = +1;pv[at].n*=-1;} // negz fixup


	// tris
	for(int i=0;i<h.ntris;i++)
	{
		int to = h.ntris+i;
		pt[to]=pt[i];
		pt[to].a += h.nverts;
		pt[to].b += h.nverts;
		pt[to].c += h.nverts;
		swap(pt[to].b,pt[to].c);
	}

	h.nverts=nverts;
	h.ntris=ntris;

	for(int i=0;i<h.nverts;i++)
	{
		pv[i].v *= 2;
		pv[i].v.x -= 1;
		pv[i].v.y -= 1;
		pv[i].v.z -= 1;

		//pv[i].t.x = 0.5;
		//pv[i].t.y = 0.5;
		//pv[i].t.z = 0.0;
	}

	return fc3_ok;
}



/*
	must be in a normalized state
*/

/*
	impose a 4x4x4=64 array onto the model space
	use a 64bit register to record the solid bits
	where vertices are found
*/
unsigned long long int fc3_s::calc_solid()
{
	unsigned long long bits=0;
	for(int i=0;i<h.nverts;i++)
	{
		unsigned x = (pv[i].v.x+1)*2; if(3<x)x=3;
		unsigned y = (pv[i].v.y+1)*2; if(3<y)y=3;
		unsigned z = (pv[i].v.z+1)*2; if(3<z)z=3;
		unsigned ofs = (z<<4) | (y<<2) | (x<<0);
		bits |= (1llu<<ofs);
	}
	return bits;
}

void fc3_s::sphere_normals()
{
	for(int i=0;i<h.nverts;i++)
	{
		pv[i].n = pv[i].v;
		pv[i].n.normalize();
	}
}

void fc3_s::sphere_texture()
{
	sphere_normals();
	const double pi = acos(-1);
	for(int i=0;i<h.nverts;i++)
	{
		double lon = (atan2( pv[i].n.x, pv[i].n.z ) + pi);
		lon /= (2*pi);
		pv[i].t.x = lon;

		double lat = asin( pv[i].n.y ) + pi/2;
		lat /= (1*pi);

		pv[i].t.y = lat;
	}
}

void fc3_s::center()
{
	fc3_vec_t vhi,vlo,vavg;
	find_bounds(vhi,vlo);
	vavg = (vhi + vlo) * 0.5f;
	for(int i=0;i<h.nverts;i++)
		pv[i].v -= vavg;
}

void fc3_s::sink()
{
	normalize();
	fc3_vec_t vhi,vlo;
	find_bounds(vhi,vlo);
	double dy =  -1.0 - vlo.y;
	for(int i=0;i<h.nverts;i++)
	{
		pv[i].v.y += dy;
	}
}

void fc3_s::swell()
{
	// swell the whole model outward to the edge of the math
	fc3_vec_t vhi,vlo;
	find_bounds(vhi,vlo);
	double xfact = 1.0 / max(fabs(vhi.x),fabs(vlo.x));
	double yfact = 1.0 / max(fabs(vhi.y),fabs(vlo.y));
	double zfact = 1.0 / max(fabs(vhi.z),fabs(vlo.z));
	for(int i=0;i<h.nverts;i++)
	{
		pv[i].v.x *= xfact;
		pv[i].v.y *= yfact;
		pv[i].v.z *= zfact;
	}
}

void fc3_s::swellx()
{
	fc3_vec_t vhi,vlo;
	find_bounds(vhi,vlo);
	double xfact = 1.0 / max(fabs(vhi.x),fabs(vlo.x));
	for(int i=0;i<h.nverts;i++)
		pv[i].v.x *= xfact;
}
void fc3_s::swelly()
{
	fc3_vec_t vhi,vlo;
	find_bounds(vhi,vlo);
	double yfact = 1.0 / max(fabs(vhi.y),fabs(vlo.y));
	for(int i=0;i<h.nverts;i++)
		pv[i].v.y *= yfact;
}
void fc3_s::swellz()
{
	fc3_vec_t vhi,vlo;
	find_bounds(vhi,vlo);
	double zfact = 1.0 / max(fabs(vhi.z),fabs(vlo.z));
	for(int i=0;i<h.nverts;i++)
		pv[i].v.z *= zfact;
}


void fc3_s::growx( const float factor )
{
	if(factor)for(int i=0;i<h.nverts;i++)
		pv[i].v.x *= factor;
}
void fc3_s::growy( const float factor )
{
	if(factor)for(int i=0;i<h.nverts;i++)
		pv[i].v.y *= factor;
}
void fc3_s::growz( const float factor )
{
	if(factor)for(int i=0;i<h.nverts;i++)
		pv[i].v.z *= factor;
}
void fc3_s::grow( const float factor )
{
	growx(factor);
	growy(factor);
	growz(factor);
}

fc3_error_t fc3_s::add( const fc3_s& fa )
{
	// make copy of this
	fc3_s ft=*this;

	int ofs=0;

	reset();
	h.init();
	h.nverts = ft.h.nverts + fa.h.nverts;
	h.ntris = ft.h.ntris + fa.h.ntris;

	if(false==heap())
		return fc3_fail_heap;

	ofs=0;
	for(int i=0;i<ft.h.nverts;i++)
		pv[i+ofs] = ft.pv[i];
	ofs=ft.h.nverts;
	for(int i=0;i<fa.h.nverts;i++)
		pv[i+ofs] = fa.pv[i];

	ofs=0;
	for(int i=0;i<ft.h.ntris;i++)
		pt[i+ofs] = ft.pt[i];
	ofs=ft.h.ntris;
	for(int i=0;i<fa.h.ntris;i++)
		pt[i+ofs] = fa.pt[i];

	ofs=ft.h.ntris;
	for(int k=0;k<fa.h.ntris;k++)
	{
		int i = k+ofs;
		pt[i].a += ft.h.nverts;
		pt[i].b += ft.h.nverts;
		pt[i].c += ft.h.nverts;
	}

	set_image((unsigned long*)ft.pc,ft.h.cwidth,ft.h.cheight);

	return fc3_ok;
}

void fc3_s::move( const fc3_vec_t dv )
{
	for(int i=0;i<h.nverts;i++)
		pv[i].v += dv;
}

void fc3_s::movex( const float dv )
{
	for(int i=0;i<h.nverts;i++)
		pv[i].v.x += dv;
}
void fc3_s::movey( const float dv )
{
	for(int i=0;i<h.nverts;i++)
		pv[i].v.y += dv;
}
void fc3_s::movez( const float dv )
{
	for(int i=0;i<h.nverts;i++)
		pv[i].v.z += dv;
}

fc3_error_t fc3_s::monochrome( const unsigned int pixel )
{
	if(false==set_image( (unsigned long*)(&pixel), 1,1 ))
		return fc3_fail_heap;
	return fc3_ok;
}

unsigned long long fc3_s::calc_image_hash_one
(
	 unsigned int* pbits
	,const size_t H
	,const size_t W
)
{
	if(NULL==pbits)
		return 0;
	size_t npix = H*W;
	if(0==npix)
		return 0;
	if(1==npix)
		return pbits[0];

	unsigned long long hsiz=0; hsiz=H; hsiz<<=16; hsiz|=W;
	unsigned long long hadd=hsiz;
	unsigned long long hxor=hsiz;

	for(int h=0;h<H;h++)
	for(int w=0;w<W;w++)
	{
		// encode the pixel at xy
		// hi order 32bits is the rgba value
		// lo order 32bits is the pixel coord
		unsigned long long hpix = pbits[W*h+w];
		hpix <<= 16;
		hpix |= h;
		hpix <<= 16;
		hpix |= w;

		// the add/xor operators are not order sensitive
		hadd += hpix; // accumulate
		hxor ^= hpix; // accumulate

		// the left/right bit shift operators are order sensitive
		int na = hadd&1;
		int nx = hxor&1;
		if(nx) hadd>>=1; else hadd<<=1;
		if(na) hxor>>=1; else hxor<<=1;

	}
	unsigned long long hash = hadd + hxor;
	return hash;
}
unsigned long long fc3_s::calc_image_hash_one() const
{
	return calc_image_hash_one( pc,h.cheight,h.cwidth );
}

fc3_error_t fc3_s::set_image_hash_one( const unsigned long long int hash )
{
	if(false==set_image( (unsigned long*)(&hash), 2,1 ))
		return fc3_fail_heap;

	h.cheight=0	; // force height of 0
	h.cwidth=1	; // force width of 1

	return fc3_ok;
}

unsigned long long int fc3_s::get_image_hash_one() const
{
	if(NULL!=pc)
	if(true==h.is_image_hash_one())
		return *((unsigned long long int*)(pc));

	return 0;
}

double fc3_s::yfloor()
{
	fc3_vec_t vhi,vlo;
	find_bounds( vhi,vlo );
	double dy = 0-vlo.y;
	movey( dy ); // vertical move that will put the lowest y value on the zero plane
	return dy;
}

fc3_error_t fc3_s::add_verts( const fc3_vert_s* va, const int nv )
{
	fc3_s fc;
	fc.reset();
	fc.h.nverts = h.nverts+nv;
	fc.pv = new fc3_vert_s[ fc.h.nverts ];
	if(NULL==fc.pv)
		return fc3_fail_heap;

	for(int i=0;i<h.nverts;i++)
		fc.pv[i] = pv[i];
	for(int i=0;i<nv;i++)
		fc.pv[i+h.nverts] = va[i]; 

	delete pv;
	pv=NULL;
	pv = fc.pv;
	fc.pv=NULL;
	h.nverts=fc.h.nverts;

	return fc3_ok;

}

void fc3_s::radial_clamp( const double rad )
{
	for(int i=0;i<h.nverts;i++)
	{
		fc3_vec_t v = pv[i].v;
		v.y = 0;
		double len = v.length();
		if(rad < len)
		{
			v.normalize();
			v *= (float)rad;
			pv[i].v.x = v.x;
			pv[i].v.z = v.z;
		}
	}
}

/*
	This is for unit testing purposes.
	We save a single triangle fc3
	having opposite endianess to our native endianess.
*/
fc3_error_t fc3_s::save_triangle_endian(const char* fn, const bool doend)
{
	if(NULL==fn)return fc3_bad_arg;

	reset();
	h.reset();
	h.init();
	h.format='b';
	h.nverts=3;
	h.ntris=1;
	h.cheight=1;
	h.cwidth=1;
	if(!heap())
		return fc3_fail_heap;

	fc3_vertf_s v;
	v.zero();
	v.n.z = -1; // normal is facing opengl -Z
	v.t.x = 0.5f;
	v.t.y = 0.5f;
	pv[0]=v;
	pv[1]=v; pv[1].v.y=1;
	pv[2]=v; pv[2].v.x=1;
	pt[0].a = 0;
	pt[0].b = 1;
	pt[0].c = 2;
	pc[0]=0xffbb7733; // aarrggbb

	FILE *fp=fopen(fn,"wb");
	if(NULL==fp)
		return fc3_fail_open;

	// header
	if(doend)
		fc3_reverse_bytes(h);
	if(1!=fwrite(&h,sizeof(h),1,fp))
	{
		fclose(fp);
		return fc3_fail_write;
	}
	// NO h. refs beyond here because it is no longer native

	// verts
	for(int i=0;i<3;i++)
	{
		fc3_vert_b vv;
		vv.set( pv[i].v, pv[i].n, pv[i].t );
		if(doend)
			vv.reverse_bytes();
		if(!fwrite(&vv,sizeof(vv),1,fp))
		{
			fclose(fp);
			return fc3_fail_write;
		}
	}

	// tris
	if(doend)
		pt[0].reverse_bytes();
	if(1!=fwrite(pt,sizeof(pt[0]),1,fp))
	{
		fclose(fp);
		return fc3_fail_write;
	}

	// image
	if(doend)
		fc3_reverse_bytes(&pc[0],sizeof(pc[0]));
	if(1!=fwrite(pc,sizeof(pc[0]),1,fp))
	{
		fclose(fp);
		return fc3_fail_write;
	}

	fclose(fp);
	return fc3_ok;
}

bool fc3_s::unit_test_endian( const std::string& path )
{
	string fnnat = path + "natend_tri.fc3"; // native endian
	string fnrev = path + "revend_tri.fc3"; // reverse endian

	fc3_s fc3;
	if(fc3_ok != fc3.save_triangle_endian(fnnat.c_str(),false))
		return false;

	if(fc3_ok != fc3.save_triangle_endian(fnrev.c_str(),true))
		return false;

	fc3_s a,b;
	size_t z;

	if(fc3_ok != a.load(fnnat.c_str()))
		return false;
	if(fc3_ok != b.load(fnrev.c_str()))
		return false;

	// expecting binary match between a and b
	z = sizeof(a.h);
	if(memcmp(&a.h,&b.h,z))
		return false;
	z = sizeof(a.pv[0])*a.h.nverts;
	if(memcmp(a.pv,b.pv,z))
		return false;
	z = sizeof(a.pt[0])*a.h.ntris;
	if(memcmp(a.pt,b.pt,z))
		return false;
	z = sizeof(pc[0])*a.h.cheight*a.h.cwidth;
	if(memcmp(a.pc,b.pc,z))
		return false;

	return true;
}

#if 0
/*
	imagine a circle centered at a distance of 'radius' above(+) / below(-) the model
	visit each vertex and slide it along the radial line toward the center
	such that the distance from center is what the vertical flat distance was.
*/
void fc3_s::spherey(const double radius)
{
	fc3_vecf_t C(0,radius,0);
	for(int i=0;i<h.nverts;i++)
	{
		fc3_vecf_t delta = pv[i].v - C;
		double r = fabs(delta.y);
		delta.normalize();
		delta*=r;
		pv[i].v = delta;
	}
	calc_normals();
}
#endif

/*
	imagine an axis aligned with the fore/back Z axis
	the axis is situated a distance of 'radius' above(+) below(-) from the model center
	we want to curl all the points around the axis
*/
void fc3_s::curly(const double rady)
{
	double sign = (0<rady) ? -1 : +1;
	const double pi = acos(-1);
	fc3_vecf_t C(0,rady,0);
	for(int i=0;i<h.nverts;i++)
	{
		fc3_vecf_t p = pv[i].v;
		C.z = p.z;
		fc3_vecf_t V = p - C;
		double r = fabs(V.y);
		double radians = p.x / r;
		V.x = sin(radians)*r;
		V.y = cos(radians)*r*sign;
		pv[i].v = (C+V);
	}
	calc_normals();
}

// expect image
int fc3_s::make_from_image()
{
	if(pv) {delete pv; pv=NULL;}
	if(pt) {delete pt; pt=NULL;}

	int dx=h.cwidth;
	int dy=h.cheight;

	h.nverts = dx*dy;
	h.ntris = (dx-1)*(dy-1)*2;

	pv = new fc3_vert_s[ h.nverts ];
	pt = new fc3_tri_s[ h.ntris ];
	if(!pv || !pt)
	{
		return fc3_fail_heap;
	}

	double zpix2 = 0.5 / (double)dy;
	double xpix2 = 0.5 / (double)dx;

	int k=0;
	for(int z=0;z<dy;z++)
	for(int x=0;x<dx;x++)
	{
		unsigned int pix = pc[k];
		unsigned char* puc = (unsigned char*)&pix;
		float avg = (puc[0]+puc[1]+puc[2]); avg/=3;
		pv[k].v.z=z;
		pv[k].v.x=x;
		pv[k].v.y= avg / 32;
		pv[k].t.x = xpix2 + (float)x/dx;
		pv[k].t.y = zpix2 + (float)z/dy;
		pv[k].t.z = 0;
		pv[k].n.x=0;
		pv[k].n.y=1;
		pv[k].n.z=0;
		k++;
	}

	k=0;
	for(int z=0;z<dy-1;z++)
	for(int x=0;x<dx-1;x++)
	{
		pt[k].a = z*dx+x;
		pt[k].b = (z+1)*dx+(x+1);
		pt[k].c = (z+0)*dx+(x+1);
		k++;
		pt[k].a = z*dx+x;
		pt[k].b = (z+1)*dx+(x+0);
		pt[k].c = (z+1)*dx+(x+1);
		k++;
	}

	return fc3_ok;
}

int fc3_s::swap_cbytes( const unsigned a, const unsigned b )
{
	if(NULL==pc)
		return fc3_bad_arg; // no image
	if(3<a || 3<b)
		return fc3_bad_arg;

	unsigned char* puc = (unsigned char*)pc;
	int npix = h.cheight*h.cwidth;
	for(int i=0;i<npix;i++)
	{
		std::swap(puc[a],puc[b]);
		puc += 4;
	}
	return fc3_ok;
}

int fc3_s::invert_color()
{
	if(NULL==pc)
		return fc3_bad_arg; // no image

	int npix = h.cheight*h.cwidth;
	for(int i=0;i<npix;i++)
		pc[i] ^= 0x00ffffff;
	return fc3_ok;
}

int fc3_s::get_header( const char* fn, fc3_header_s& h )
{
	if(NULL==fn)
		return fc3_bad_arg;
	h.reset();
	FILE* fp = fopen(fn,"rb");
	size_t result = fread(&h, sizeof(fc3_header_s), 1, fp );
	fclose(fp);
	if(1!=result)
		return fc3_fail_read;
	return fc3_ok;
}

