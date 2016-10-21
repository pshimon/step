/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "surf.h"
#include "geom.h"
int main(int argc,char * argv[]) {
    int ret,i;
    Flt * dist;
    TSurf * s=ALLOC_MEM0(TSurf,1);
    ret=mk_tetrahedron(s);
    dist=ALLOC_MEM(Flt,s->nt);
    ret=check_normals(dist,s);
    printf("tetrahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);
    //ret=prinTSurf(s,"tetrahedron.tsa");
    ret=writeTSurf(s,"tetrahedron.tsb");

    ret=mk_hexahedron(s);
    dist=ALLOC_MEM(Flt,s->nt);
    ret=check_normals(dist,s);
    printf("hexahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);

    //ret=prinTSurf(s,"hexahedron.tsa");
    ret=writeTSurf(s,"hexahedron.tsb");

    ret=mk_octahedron(s);
   dist=ALLOC_MEM(Flt,s->nt);
    ret=check_normals(dist,s);
    printf("octahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);

    //ret=prinTSurf(s,"octahedron.tsa");
    ret=writeTSurf(s,"octahedron.tsb");

    ret=mk_dodecahedron(s);
    dist=ALLOC_MEM(Flt,s->nt);
    ret=check_normals(dist,s);
    printf("dodecahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);

    //ret=prinTSurf(s,"dodecahedron.tsa");
    ret=writeTSurf(s,"dodecahedron.tsb");

    ret=mk_icosahedron(s);

  dist=ALLOC_MEM(Flt,s->nt);
    ret=check_normals(dist,s);
    printf("icosahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);
  

    //ret=prinTSurf(s,"icosahedron.tsa");
    ret=writeTSurf(s,"icosahedron.tsb");

    return 0;
}
