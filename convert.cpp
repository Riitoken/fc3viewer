// convert.cpp

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
#include <string>
#include <vector>

using namespace std;

#pragma warning( disable : 4018 ) // '<': signed/unsigned mismatch
#pragma warning( disable : 4244 ) // conversion from <big> to <small>, possible loss of data
#pragma warning( disable : 4267 ) // conversion from 'size_t' to <smaller>, possible loss of data

int write_mtl( const string& fn, const string& mat )
{
	if(mat<=" ")
		return 0;

	FILE* fp = fopen( fn.c_str(),"w");
	if(NULL==fp)
		return -1;

	fprintf(fp,"# fc3 export\n");
	fprintf(fp,"# %s.mtl\n",mat.c_str());
	fprintf(fp,"newmtl %s\n",mat.c_str());
	fprintf(fp,"illum 1\n");
	fprintf(fp,"Kd 1 1 1\n"); // diffuse full
	fprintf(fp,"Ka 1 1 1\n"); // ambient full
//	fprintf(fp,"Ks 0.500 0.500 0.500\n");
//	fprintf(fp,"Ns 96.078431\n");
	fprintf(fp,"map_Kd %s.png\n",mat.c_str()); // diffuse texture

	fclose(fp);
	return 0;
}

string remove_suffix( const string& s, const string& delim )
{
	std::size_t found = s.find_last_of(delim);
	string pre = s.substr(0,found);
	return pre;
}

int fc3_s::export_obj_mtl( const string& fn, const string& mtl )
{
	FILE* fp = fopen( fn.c_str(),"w");
	if(NULL==fp)
		return -1;

	if(mtl > " ")
	{
		string mfn = mtl+".mtl";
		fprintf(fp,"# fc3 export\n");
		fprintf(fp,"mtllib %s\n",mfn.c_str());
		fprintf(fp,"usemtl %s\n",mtl.c_str());
		fprintf(fp,"\n");
		string path = remove_suffix(fn,"/\\");
		mfn = path + "/" + mfn;
		write_mtl(mfn,mtl);
	}


	int zlo=h.nverts+1;
	int zhi=0;
	for(int i=0;i<h.nverts;i++)
	{
		fc3_vert_s& v = pv[i];
		if(v.t.z > zhi)zhi=v.t.z;
		if(v.t.z < zlo)zlo=v.t.z;
		fprintf(fp,"v %f %f %f\n",v.v.x,v.v.y,v.v.z);
		fprintf(fp,"vn %f %f %f\n",v.n.x,v.n.y,v.n.z);
		fprintf(fp,"vt %f %f\n",v.t.x,v.t.y);
		fprintf(fp,"\n");
	}

	for(int z=zlo;z<=zhi;z++)
	{
		bool flag=false;
		for(int i=0;i<h.ntris;i++)
		{
			fc3_tri_s& t = pt[i];

			if(z != pv[ t.a ].t.z)
				continue;

			if(!flag)
				fprintf(fp,"o main.%d\n",z);
			flag=true;

			int a = t.a+1;
			int b = t.b+1;
			int c = t.c+1;
			fprintf(fp,"f %d/%d/%d %d/%d/%d %d/%d/%d\n"
				,a,a,a
				,b,b,b
				,c,c,c
			);
		}
	}
	fclose(fp);

	return 0;
}

int fc3_s::export_obj( const string& folder, const string& name )
{
	string fn = folder + "/" + name + ".obj";
	return export_obj_mtl(fn,name);
}

// string splitter - slightly slow but very generic
void split_strings( char* pc, vector<string>& vs, const char delim=0 )
{
	vs.clear();
	if(NULL==pc)
		return;

	string s;

	size_t len = strlen(pc);
	if(delim)
		{ for(int i=0;i<len;i++) if(delim==pc[i]) pc[i]=0; }
	else
		{ for(int i=0;i<len;i++) if(' '>=pc[i]) pc[i]=0; }

	int k=0;
	while(k<len)
	{
		while((k<len) && (0==pc[k])) k++;
		if(k>=len)break;
		s = (&pc[k]); k+= s.length();
		if(0<s.length())
			vs.push_back(s);
	}
}

// to parse obj "f" records e.g.
//		f 1/2/3 4/5/6 7/8/9
void split_triple( const char* s, int& vv, int& vt, int& vn )
{
	char f[3][32]; f[0][0]=0; f[1][0]=0; f[2][0]=0;
	int k=0;
	for(int i=0;i<3;i++)
	{
		int j;
		for(j=0;j<31 && s[k] && s[k]!='/';j++)
		{
			f[i][j] = s[k++];
		}
		f[i][j]=0;
		k++;
	}
	vv = atoi(f[0]);
	vt = atoi(f[1]);
	vn = atoi(f[2]);
}

#define vector3f fc3_vec_t

int fc3_s::import_obj( const string& fn )
{
	FILE* fp = fopen(fn.c_str(),"r");
	if(NULL==fp)
	{
		return -1;
	}

	vector3f zero;
	vector<vector3f> vv; vv.push_back(zero);
	vector<vector3f> vn; vn.push_back(zero);
	vector<vector3f> vt; vt.push_back(zero);

	struct face_s { int a,b,c; face_s() : a(0),b(0),c(0) {} };
	typedef vector<face_s> vf_t;
	vector<vf_t> vf;

	vector<string> vs;

	char buf[4096]={0};
	vector3f v;
	while(!feof(fp))
	{
		buf[0]=0;
		if(NULL==fgets(buf,sizeof(buf),fp))
			break;
		split_strings( buf, vs );
		if(0>=vs.size())
			continue;

		if(vs[0] == "v")
		{
			if(4 > vs.size())
			{
				return -2;
			}
			v.x = atof( vs[1].c_str() );
			v.y = atof( vs[2].c_str() );
			v.z = atof( vs[3].c_str() );
			vv.push_back( v );
		}
		else
		if(vs[0] == "vt")
		{
			if(3 > vs.size())
			{
				return -3;
			}
			v.x = atof( vs[1].c_str() );
			v.y = atof( vs[2].c_str() );
			v.z = 0;
			vt.push_back( v );
		}
		else
		if(vs[0] == "vn")
		{
			if(4 > vs.size())
			{
				return -4;
			}
			v.x = atof( vs[1].c_str() );
			v.y = atof( vs[2].c_str() );
			v.z = atof( vs[3].c_str() );
			vn.push_back( v );
		}
		else
		if(vs[0] == "f")
		{
			if(4 > vs.size())
			{
				return -5;
			}
			vf_t tt; tt.resize( vs.size() -1 );
			vf.push_back(tt);
			int ndx = vf.size()-1;
			for(int i=1;i<vs.size();i++)
			{
				face_s f;
				split_triple( vs[i].c_str(), f.a, f.b, f.c );
				vf[ndx][i-1] = f;
			}
		}
	}
	fclose(fp);

	vector<fc3_vertf_s> vfc3;
	vector< vector<int> > vface;
	vector< fc3_tri_s > vtri;
	vector<int> vi;

	fc3_vertf_s fc3v;
	for(int i=0;i<vf.size();i++)
	{
		vi.resize(vf[i].size());
		vface.push_back(vi);
		for(int k=0;k<vf[i].size();k++)
		{
			face_s f = vf[i][k];
			vector3f a = vv[f.a];
			vector3f b = vt[f.b];
			vector3f c = vn[f.c];
			fc3v.v.from(&a.x);
			fc3v.t.from(&b.x);
			fc3v.n.from(&c.x);
			vface[i][k]=vfc3.size();
			vfc3.push_back(fc3v);
		}
	}

	for(int i=0;i<vface.size();i++)
	{
		fc3_tri_s tri;
		tri.a = vface[i][0];
		for(int k=2;k<vface[i].size();k++)
		{
			tri.b = vface[i][k-1];
			tri.c = vface[i][k-0];
			vtri.push_back(tri);
		}
	}

	fc3_s fc3;
	fc3.h.init();
	fc3.h.nverts = vfc3.size();
	fc3.h.ntris = vtri.size();
	fc3.heap();
	for(int i=0;i<fc3.h.nverts;i++)
		fc3.pv[i] = vfc3[i];
	for(int i=0;i<fc3.h.ntris;i++)
		fc3.pt[i] = vtri[i];

	return 0;
}

#if 0
void path_suffix( string& path )
{
	if(path <= " ")
	{
		path = "./";
	}
	else
	{
		char c = path[path.length()-1];
		if( '\\' != c && '/' != c)
			path += "/";
	}
}

int fc3_s::split_obj( string& path, const string& objnam )
{
	path_suffix(path);
	string fn = path + objnam;
	FILE* fp = fopen(fn.c_str(),"r");
	if(NULL==fp)
		return fc3_fail_open;



	fclose(fp);

}
#endif


#undef vector3f
