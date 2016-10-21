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
    t_surf * s=ALLOC_MEM0(t_surf,1);
    ret=mk_tetrahedron(s);
    dist=ALLOC_MEM(Flt,s->nt);
    ret=check_normals(dist,s);
    printf("tetrahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);
    //ret=print_surf(s,"tetrahedron.tsa");
    ret=write_surf(s,"tetrahedron.tsb");

    ret=mk_hexahedron(s);
    dist=ALLOC_MEM(Flt,s->nt);
    ret=check_normals(dist,s);
    printf("hexahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);

    //ret=print_surf(s,"hexahedron.tsa");
    ret=write_surf(s,"hexahedron.tsb");

    ret=mk_octahedron(s);
   dist=ALLOC_MEM(Flt,s->nt);
    ret=check_normals(dist,s);
    printf("octahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);

    //ret=print_surf(s,"octahedron.tsa");
    ret=write_surf(s,"octahedron.tsb");

    ret=mk_dodecahedron(s);
    dist=ALLOC_MEM(Flt,s->nt);
    ret=check_normals(dist,s);
    printf("dodecahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);

    //ret=print_surf(s,"dodecahedron.tsa");
    ret=write_surf(s,"dodecahedron.tsb");

    ret=mk_icosahedron(s);

  dist=ALLOC_MEM(Flt,s->nt);
    ret=check_normals(dist,s);
    printf("icosahedron\n");
    if(ret) 
	fprintf(stdout,"bad trg %d\n",ret);
    else 
	for(i=0;i<s->nt;i++) printf("\t%d\t%e\n",i,dist[i]);
    FREE_MEM(dist);
  

    //ret=print_surf(s,"icosahedron.tsa");
    ret=write_surf(s,"icosahedron.tsb");

    return 0;
}
