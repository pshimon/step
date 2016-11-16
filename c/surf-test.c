/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "tsurf.h"

int main(int argc,char * argv[]) {
    int ret,i;
    Flt * dist;
    TSurf * s=ALLOC_MEM0(TSurf,1);
    ret=mkTetrahedron(s);
    dist=ALLOC_MEM(Flt,s->nt);
    ret=checkNormals(dist,s);
    printf("tetrahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);
    //ret=prinTSurf(s,"tetrahedron.tsa");
    ret=writeTSurf(s,"tetrahedron.tsb");

    ret=mkHexahedron(s);
    dist=ALLOC_MEM(Flt,s->nt);
    ret=checkNormals(dist,s);
    printf("hexahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);

    //ret=prinTSurf(s,"hexahedron.tsa");
    ret=writeTSurf(s,"hexahedron.tsb");

    ret=mkOctahedron(s);
   dist=ALLOC_MEM(Flt,s->nt);
    ret=checkNormals(dist,s);
    printf("octahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);

    //ret=prinTSurf(s,"octahedron.tsa");
    ret=writeTSurf(s,"octahedron.tsb");

    ret=mkDodecahedron(s);
    dist=ALLOC_MEM(Flt,s->nt);
    ret=checkNormals(dist,s);
    printf("dodecahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);

    //ret=prinTSurf(s,"dodecahedron.tsa");
    ret=writeTSurf(s,"dodecahedron.tsb");

    ret=mkIcosahedron(s);

  dist=ALLOC_MEM(Flt,s->nt);
    ret=checkNormals(dist,s);
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
