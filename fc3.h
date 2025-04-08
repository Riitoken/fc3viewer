//fc3.h

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

#ifndef _fc3_h_
#define _fc3_h_

#include <string>
#include <map>
#include <math.h>
#include <vector>
#include <string>
#include "abcdef.h"

typedef floc fc3f;

typedef enum {

	 fc3_ok				= +0
	,fc3_bad_arg		= -1
	,fc3_bad_header		= -2
	//...
	,fc3_fail_heap		= -10
	,fc3_fail_open		= -11
	,fc3_fail_close		= -12
	,fc3_fail_read		= -13
	,fc3_fail_write		= -14
	,fc3_fail_test		= -15

	,fc3_fatal_error	= -20

} fc3_error_t;

void fc3_reverse_bytes( void* v, const int nbytes );

/*
	demand single byte alignment
	prevent the compiler from padding any fields in the header
	header is explicitly padded below
*/
#pragma pack(1)

struct fc3_header_flags_s
{
	/////////////////////////////////////////////////////////////////
	// LSB
	/////////////////////////////////////////////////////////////////
	untb	alpha		: 1; // model wants alpha blended
	untb	collision	: 1; // collision enabled
	untb	lofill		: 5; // future use
	untb	lsb0		: 1; // must always be 0 for lsb byte
	/////////////////////////////////////////////////////////////////
	// MSB
	/////////////////////////////////////////////////////////////////
	untb	hifill		: 6; // future use
	untb	exact		: 1; // true for exact, false for lossy
	untb	msb1		: 1; // must always be 1 for msb byte
	/////////////////////////////////////////////////////////////////
};

struct fc3_header_s
{
	// 4 byte file signature
	signed char sigF;	// 'F'
	signed char sigC;	// 'C'
	signed char sig3;	// '3'
	signed char sigv;	// 'a'

	// 4 bytes orientation and data size
	//		3 ascii bytes - axis orientation string
	//		1 ascii bytes - internal data format
	signed char	xaxis; // ascii default is 'R' for rightward
	signed char	yaxis; // ascii default is 'U' for upward
	signed char	zaxis; // ascii default is 'B' for backward
	signed char	format; // ascii case-insensitive indicator for vertex data size: a=1byte, b=2byte, c=4byte, d=8byte, e=16byte, etc., etc.

	// 4 bytes
	unsigned short int endian; // default: is msb='e' (ascii 101), lsb='E' (ascii 69), hex value 0x6545
	signed char vscale; // 2^vscale is maximum original value of any vertex coord
	signed char tscale; // 2^tscale is maximum original value of any texture coord

	// 4 bytes texture size
	unsigned short int cwidth	; // canvas texture width in pixels.  1 pixel = 4 bytes.
	unsigned short int cheight	; // canvas texture height in pixels.  1 pixel = 4 bytes.

	// 8 bytes counts
	unsigned long int nverts	; // total number of vertices in the file.  first vertex is 0, final vertex is nverts-1
	unsigned long int ntris	; // total number of triangle faces.  first triangle is 0, final triangle is ntris-1

	// 8 bytes
	double  unitlen; // the metric length of a one unit  i.e. (0.0254 meters == 1 inch)

	//////////////////////
	// 32 BYTES EXACTLY //
	//////////////////////

	fc3_header_s();
	void reset();
	bool is_valid_file_signature()const;
	void init();
	void set_scaling_exponents( const char v, const char t );
	bool is_endian_valid() const;
	bool is_endian_correct() const;
	bool is_do_endian() const;
	bool is_eE() const; // legacy endian values

	size_t get_npix() const;

	static double get_unit_meters( const std::string& name );
	static const char* get_unit_name( const double meters );

	bool is_image_hash_one() const;
	bool is_image_hash() const;

	size_t calc_file_size() const;

	// FLAGS clever overloaded use of the endian bits
	void enable_flags();
	bool has_flags() const;
	void set_flag_bit( const int b, const bool v );
	bool get_flag_bit( const int b ) const;
	int get_alpha() const;
	int get_collision() const;
	int get_exact() const;

	void set_alpha( const bool v );
	void set_collision( const bool v );
	void set_exact( const bool v );
	void get_endian( unta& msb, unta& lsb ) const;
	void set_endian( unta msb, unta lsb );

}; // end class fc3_header

#pragma pack()

template<typename T>
struct fc3_vector_template
{
	T x;
	T y;
	T z;

	// explicit conversion
	explicit operator T*() const { return (&x); }

	fc3_vector_template() : x(0),y(0),z(0) {}

	fc3_vector_template( const fc3_vector_template<float>& a ) : x((T)a.x), y((T)a.y), z((T)a.z) {}
	fc3_vector_template( const fc3_vector_template<double>& a ) : x((T)a.x), y((T)a.y), z((T)a.z) {}

	fc3_vector_template( const float a[3] )  : x((T)a[0]), y((T)a[1]), z((T)a[2]) {}
	fc3_vector_template( const double a[3] )  : x((T)a[0]), y((T)a[1]), z((T)a[2]) {}

	fc3_vector_template( const float fx, const float fy, const float fz )  : x((T)fx), y((T)fy), z((T)fz) {}
	//fc3_vector_template( const double fx, const double fy, const double fz )  : x((T)fx), y((T)fy), z((T)fz) {}

	void from( const float* f ) { x=(T)(f[0]); y=(T)(f[1]); z=(T)(f[2]); }
	void from( const double* f ) { x=(T)(f[0]); y=(T)(f[1]); z=(T)(f[2]); }

	void set( const float a[3] ) { from(a); }
	void set( const double a[3] ) { from(a); }

	void get( float* a ) const { a[0]=(T)x; a[1]=(T)y; a[2]=(T)z; }
	void get( double* a ) const { a[0]=(T)x; a[1]=(T)y; a[2]=(T)z; }

	void zero() 
	{
		x=y=z=0;
	}

	/*
	bool operator < ( const T& a ) const
	{
		if(x < a.x)return true;
		if(y < a.y)return true;
		if(z < a.z)return true;
		return false;
	}
	*/

	T maxabs() const
	{
		T a=(T)fabs(x);
		T b=(T)fabs(y);
		T c=(T)fabs(z);
		T d=(T)((a>b)?(a):(b));
		d=(d>c)?(d):(c);
		return d;
	}

	void maxv( fc3_vector_template& m )
	{
		if(x>m.x)m.x=x;
		if(y>m.y)m.y=y;
		if(z>m.z)m.z=z;
	}
	void minv( fc3_vector_template& m )
	{
		if(x<m.x)m.x=x;
		if(y<m.y)m.y=y;
		if(z<m.z)m.z=z;
	}

	fc3_vector_template operator - ( const fc3_vector_template &other ) const
	{
		fc3_vector_template vResult;

		vResult.x = x - other.x;
		vResult.y = y - other.y;
		vResult.z = z - other.z;

		return vResult;
	}
	fc3_vector_template operator + ( const fc3_vector_template &other ) const
	{
		fc3_vector_template vResult;

		vResult.x = x + other.x;
		vResult.y = y + other.y;
		vResult.z = z + other.z;

		return vResult;
	}

	fc3_vector_template& operator += ( const fc3_vector_template &other )
	{
		x += other.x;
		y += other.y;
		z += other.z;

		return *this;
	}

	fc3_vector_template& operator -= ( const fc3_vector_template &other )
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;

		return *this;
	}

	fc3_vector_template operator * ( const T scalar ) const
	{
		fc3_vector_template vResult;

		vResult.x = x * scalar;
		vResult.y = y * scalar;
		vResult.z = z * scalar;

		return vResult;
	}

	fc3_vector_template& operator *= ( const T scalar )
	{
		x *= scalar;
		y *= scalar;
		z *= scalar;
		return *this;
	}

	static T dot( const fc3_vector_template &v1,  const fc3_vector_template &v2 )
	{
		return( v1.x * v2.x + v1.y * v2.y + v1.z * v2.z  );
	}

	// cross product
	static fc3_vector_template cross( const fc3_vector_template &v1,  const fc3_vector_template &v2 )
	{
		fc3_vector_template vCrossProduct;
		vCrossProduct.x =  v1.y * v2.z - v1.z * v2.y;
		vCrossProduct.y =  v1.z * v2.x - v1.x * v2.z;
		vCrossProduct.z =  v1.x * v2.y - v1.y * v2.x;
		return vCrossProduct;
	}

	T length() const
	{
		return( (T)sqrt( x*x + y*y + z*z ) );
	}

	void normalize()
	{
		T fLength = length();

		if(0!=fLength)
		{
			x /= fLength;
			y /= fLength;
			z /= fLength;
		}
	}

	bool operator==( const fc3_vector_template& v ) const
	{
		return (x==v.x && y==v.y && z==v.z);
	}
	bool operator!=( const fc3_vector_template& v ) const
	{
		return (x!=v.x || y!=v.y || z!=v.z);
	}

};

typedef fc3_vector_template<double> fc3_vecd_t;
typedef fc3_vector_template<float> fc3_vecf_t;

template<typename T>
struct fc3_vertex_template
{
	fc3_vector_template<T> v;
	fc3_vector_template<T> n;
	fc3_vector_template<T> t;

	fc3_vertex_template ( const fc3_vertex_template<float>& a ) : v(a.v), n(a.n), t(a.t) {}
	fc3_vertex_template ( const fc3_vertex_template<double>& a ) : v(a.v), n(a.n), t(a.t) {}

	fc3_vertex_template () {}

	fc3_vertex_template operator + ( const fc3_vertex_template &a ) const
	{
		fc3_vertex_template<T> v = *this;
		v.v += a.v;
		v.n += a.n;
		v.t += a.t;
		return v;
	}
	fc3_vertex_template& operator += ( const fc3_vertex_template &a )
	{
		v += a.v;
		n += a.n;
		t += a.t;
		return *this;
	}

	bool operator == ( const fc3_vertex_template &a ) const
	{
		return(true
			&& v==a.v
			&& n==a.n
			&& t==a.t
			);
	}

	// for deciding duplicates
	static bool same_vt( const fc3_vertex_template &a,  const fc3_vertex_template &b  )
	{
		return(true
			&& a.v==b.v
			&& a.t==b.t
			);
	}
	bool same_dup( const fc3_vertex_template &a ) const
	{
		return same_vt( *this, a );
	}

	void zero()
	{
		v.zero();
		n.zero();
		t.zero();
	}
};

typedef fc3_vertex_template<float> fc3_vertf_s;
typedef fc3_vertex_template<double> fc3_vertd_s;

struct axis_convert_s;

struct fc3_tri_s
{
	unsigned long int a;
	unsigned long int b;
	unsigned long int c;

	bool isv( const unsigned long int t )
	{
		return (t==a || t==b || t==c);
	}

	fc3_vecf_t calc_norm( const int i, const fc3_vecf_t& va,  const fc3_vecf_t& vb,  const fc3_vecf_t& vc )
	{
		fc3_vecf_t xf,xt,norm;
		if((int)a==i){xf=vb-va;xt=vc-va;} else
		if((int)b==i){xf=vc-vb;xt=va-vb;} else
		if((int)c==i){xf=va-vc;xt=vb-vc;} 
		norm = fc3_vecf_t::cross(xf,xt);
		//norm.normalize();
		return norm;
	}

	void reverse_bytes()
	{
		fc3_reverse_bytes( &a, sizeof(a) );
		fc3_reverse_bytes( &b, sizeof(b) );
		fc3_reverse_bytes( &c, sizeof(c) );
	}

};

/*
	32bit floats have a 24bit mantissa
	2^24 = 16777216  (sixteen million)
	2 meters / 16777216 = about 1/10th of a micron
	For the majority of models, floats are more than enough resolution

	To use double precision internally:
	force the define of FC3_DOUBLES below
	and recompile the fc3 library and native applications
*/
#undef FC3_DOUBLES
#ifdef FC3_DOUBLES
// INTERNAL DOUBLES
typedef fc3_vecd_t fc3_vec_t;
typedef fc3_vertd_s fc3_vert_s;
#else
// INTERNAL FLOATS
typedef fc3_vecf_t fc3_vec_t;
typedef fc3_vertf_s fc3_vert_s;
#endif

struct fc3_s
{
	fc3_header_s h;
	fc3_vert_s* pv;
	fc3_tri_s* pt;
	unsigned int* pc;

	std::vector<std::string> vso; // the string name of the object id
	std::vector<unta> vto; // the object id for each triangle.

	int loadva(FILE* fp);
	int loadvb(FILE* fp,const bool doend);
	int loadvc(FILE* fp,const bool doend);
	int loadvd(FILE* fp,const bool doend);
	int saveva(FILE* fp);
	int savevb(FILE* fp);
	int savevc(FILE* fp);
	int savevd(FILE* fp);

	/**
	* Put the object into clean state, deallocate memory, reset the header.
	**/
	void reset();
	fc3_s();
	fc3_s( const fc3_s& rval );
	~fc3_s();

	static bool is_valid_axis_string( const std::string& s );

	static bool normalize_axis_string( const std::string& F, std::string& f );

	static int unit_test();
	static int calculate_convert_axis( const std::string& f, const std::string& t, axis_convert_s& ac );
	int convert_axis( const axis_convert_s& ac );

	/**
	* Convert the model orientation to the setting given by the parameter\n
	* @param t - a 3 character string that identifies the target orientation.\n
	*	These 3 characters corrspond to the order inherently present in the coordinate value XYZ.\n
	*	These directions are from the natural POV of a human head model\n
	*	'R'	= Right (ear)\n
	*	'L'	= Left (ear)\n
	*	'U'	= Up (crown)\n
	*	'D'	= Down (torso)\n
	*	'F'	= Fore (nose)\n
	*	'B'	= Back (hind)\n
	* @return -1 : the existing axis is invalid or the axis given is invalid\n
	* @return -2 : internal error\n
	* @return fc3_ok : success\n
	**/
	int convert_axis_to( const std::string& t );

	/**
	* Convert the model orientation to the OpenGL default Right Up Back (RUB)
	**/
	int convert_axis_to_opengl() { return convert_axis_to("RUB"); }

	void set_scaling_exponents( const signed char v, const signed char t );

	void set_vertex_radius( const double r );
	void set_texture_radius( const double r );
	void find_bounds( double& vhi, double& thi );
	void find_bounds( fc3_vec_t& vhi, fc3_vec_t& vlo );
	void calculate_scaling_exponents();

	double get_vertex_radius() const;

	/**
	* Optionally called after calling load().  Apply the vscale exponent to the vertex coordinates. Apply the tscale exponent to the texture coordinates.\n
	* Not called by load().\n
	* @param to_meters	-	convert to meters.\n
	**/
	void denormalize( const bool to_meters = true );
	void denormalize_t(); // just the texture coords

	/**
	* Calculate the correct vscale exponent from the vertext coordinates.\n
	* Calculate the correct tscale exponent from the texture coordinates.\n
	* Convert all vertex coordinates to the interval [-1.0,+1.0].\n
	* Convert all texture coordinates to the interval [-1.0,+1.0].\n
	* Automatically called by the save() function.\n
	**/
	void normalize();

	void unheap();
	bool heap();
	bool heap( const size_t ntris, const size_t cw, const size_t ch );
	bool set_image( const unsigned long int* pbits, const size_t cw, const size_t ch );
	bool set_image( const unsigned int* pbits, const size_t cw, const size_t ch );

	void take_image( unsigned int* &pbits, int& w, int& h );
	void give_image( unsigned int* &pbits, int& w, int& h );

	void get_dimensions( double& width, double& height, double& length ) const;
	double get_extent() const;
	double get_radius() const;
	void calc_lossy_v( double& a, double& b, double& c, double& d );
	//void calc_lossy_n( double& a, double& b, double& c, double& d );
	//void calc_lossy_t( double& a, double& b, double& c, double& d );

	/**
	* For each vertex, calculate a new normal as the summed average of the triangle face normals in which the vertex participates
	**/
	void calc_normals();
	void renorm();
	void negate_normals();
	int sort_verts();

	struct f_s{FILE* fp;f_s(FILE* a):fp(a){} ~f_s(){fclose(fp);}};

	int loadh( FILE* fp);
	int loadv(FILE* fp,const bool doend);
	int savev(FILE* fp);

	/**
	* Load an fc3 model from a file\n
	* @param fn					- the name of the fc3 file
	* @return fc3_bad_arg		- fn is null
	* @return fc3_fail_open		- could not open the file
	* @return fc3_bad_header	- could not read header or bad header format
	* @return fc3_fail_heap		- could not allocate memory
	* @return fc3_fail_read		- could not read the file
	* @return fc3_ok			- success
	**/
	int load( const char* fn );

	/**
	* Save an fc3 model to a file\n
	* @param fn					- the name of the fc3 file
	* @return fc3_bad_arg		- fn is null
	* @return fc3_fail_open		- could not open the file
	* @return fc3_fail_write	- could not write the file
	* @return fc3_ok			- success
	**/
	int save( const char* fn );

	int export_obj( const std::string& folder, const std::string& name );

	int export_obj_mtl( const std::string& fn, const std::string& mtl="" );
	int import_obj( const std::string& fn );

	//static int split_obj( string& path, const string& objnam );

	fc3_error_t make_tetrahedron();
	fc3_error_t make_sphere( const int nsplits );
	fc3_error_t split_faces();
	fc3_error_t make_cube( const int nsplits );

	fc3_error_t make_globe( const int nlayers );

	fc3_error_t make_circle( const int nsides );

	void add_face( const fc3_vert_s& a, const fc3_vert_s& b, const fc3_vert_s& c );
	void add_quad( const fc3_vert_s& va, const fc3_vert_s& vb, const fc3_vert_s& vc, const fc3_vert_s& vd );

	unsigned long long int calc_solid();

	void sphere_normals();
	void sphere_texture();
	void horizontal_texture();
	void vertical_texture();

	void center();
	void swell();
	void swellx();
	void swelly();
	void swellz();

	void rise();
	void sink();
	void grow( const float f );
	void growx( const float f );
	void growy( const float f );
	void growz( const float f );

	void movex( const float v );
	void movey( const float v );
	void movez( const float v );

	void move( const fc3_vec_t dv );

	fc3_error_t add( const fc3_s& fa );

	/**
	* Combine vertices having same spatial coord and same texture coord.
	* @return - the new number of vertices.
	**/
	int undup();

	fc3_error_t monochrome( const unsigned int pixel );
	fc3_error_t middle_gray() { return monochrome( 0xff808080 ); }

	static unsigned long long calc_image_hash_one
	(
		 unsigned int* pbits
		,const size_t H
		,const size_t W
	);

	unsigned long long calc_image_hash_one() const;

	unsigned long long int get_image_hash_one() const;
	fc3_error_t set_image_hash_one( const unsigned long long int hash );

	double yfloor();
	double yceil();

	fc3_error_t add_verts( const fc3_vert_s* va, const int nv );

	void radial_clamp( const double rad );
	fc3_error_t save_triangle_endian(const char* fn, const bool doend);

	static bool unit_test_endian( const std::string& path );

	void curly(const double radius);

	int make_from_image();

	int swap_cbytes( const unsigned a, const unsigned b );
	int invert_color();

	static int get_header( const char* fn, fc3_header_s& h );

	static int copy_image( const fc3_s& fr, fc3_s& to );

	bool set_alpha( const unsigned int* pbits, const size_t cw, const size_t ch );

	int count_unused_verts() const;

	bool is_inside( const fc3_vec_t& p );

	bool is_intersected( const fc3_tri_s& tri, const fc3_vec_t& v, const fc3_vec_t& dir );

	int make_parts();

	fc3_error_t add( const std::string& fn );

	size_t get_image_size() const;

	struct mtl_s
	{
		std::map< std::string, std::map<std::string,std::string> > m;
		std::map< std::string, std::map<std::string,std::string> >::iterator it;
	};

	static int load_materials(  const std::string& fn, mtl_s& mtl );

	static double get_metric_radius( const std::string& fn );

	void add_channel( const int ch, const int delta );

	int set_objects( const std::vector<std::string>& ps, const std::vector<unta>& pv );

	void set_exact( const bool v );
	bool get_exact() const;

	struct fc3img_s
	{
		untc* pbits;
		untc  width;
		untc  height;
		fc3img_s() : pbits(NULL),width(0),height(0){}
	};

	int get_triangle_bounds( const untc tndx,  fc3_vec_t& vhi, fc3_vec_t& vlo );

	void delete_vert();
	void delete_tri();
	void delete_img();

	fc3_error_t make_platform( const int dimquads, const int curvequads );

}; // end struct fc3_s


#if 1
template<typename T>
struct array2d
{
	T* pa;
	int nrow;
	int ncol;
	bool heap;

	int get_nitem() const
	{
		return nrow*ncol;
	}
	int get_index( const int r, const int c ) const
	{
		return ncol*r+c;
	}

	void reset()
	{
		if(NULL!=pa && heap)
			delete pa;
		pa=NULL;
		heap=false;
		nrow=ncol=0;
	}

	void init( const int nr, const int nc )
	{
		reset();
		pa=new T[(nrow=nr)*(ncol=nc)];
		heap=true;
	}

	array2d() : pa(NULL),nrow(0),ncol(0),heap(false) {}

	array2d( const int nr, const int nc ) : pa(NULL), nrow(nr), ncol(nc),heap(false)
	{
		init(nr,nc);
	}

	array2d( T* p, const int nr, const int nc ) : pa(p), nrow(nr), ncol(nc),heap(false)
	{
	}

	array2d( const T& fr ) : array2d(fr.nrow,fr.ncol)
	{
		memcpy(pa,fr.pa,nrow*ncol*sizeof(pa[0]));
	}
	~array2d()
	{
		reset();
	}
	bool isoob( const int r, const int c ) const
	{
		bool oob = (false
			|| 0>r
			|| 0>c
			|| nrow <= r
			|| ncol <= c
		);
		if(oob)
			return true;
		else
			return false;
	}

	const T& getac( const int row, const int col ) const
	{
		if(pa && !isoob(row,col))
		{
			return pa[get_index(row,col)];
		}
		else
		{
			return pa[0];
		}
	}

	T& geta( const int row, const int col )
	{
		if(pa && !isoob(row,col))
		{
			return pa[get_index(row,col)];
		}
		else
		{
			return pa[0];
		}
	}

	T get( const int row, const int col )
	{
		if(pa && !isoob(row,col))
		{
			return pa[get_index(row,col)];
		}
		else
		{
			return pa[0];
		}
	}
	void set( const int row, const int col, const T& val )
	{
		if(pa && !isoob(row,col))
		{
			pa[get_index(row,col)]=val;
		}
		else
		{
			;
		}
	}
};
#endif

#endif
