#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include "geometry.h"
#include "TRI2STL.h"
#include "smooth.h"
#include "Util.h"
#include "trimesh.h"
#include "Spacing_Field.h"
#include "PList.h"
#include "Octree_Storage.h"

//global variables
extern FILE *in_f, *jou_f, *out_f, *debug_f; /*input output journal files global*/

void tri2stl_lib(geometry *geom, int in)
{
  int bdim = 132;
  char *sname = new char[bdim];
  const int bdim2 = 400;
  char buff[bdim2];
  int i, j, k, l, m, n;
  int n0, n1, n2;
  Vector v0, v1, norm;
  FILE *ucd_f = NULL;
  //declared here for use w or w/out restart
  POLYMESH *finalmesh = new POLYMESH();
  
  finalmesh->nelems = 0;
  finalmesh->nn = geom->n_gnodes;
  finalmesh->nb = geom->ngb;
  
  fprintf(out_f,"\nAllocating for geom boundary.\n");
  fflush(out_f);
  
  finalmesh->boundary = new BOUNDARY_MESH[geom->ngb];
  
  //set all names and elists
  for (i = 0; i < geom->ngb; i++)
    {
    finalmesh->boundary[i].name = new char[400];
    sprintf(finalmesh->boundary[i].name,"%s",geom->g_bname[i+1]);
    finalmesh->boundary[i].elist = new List();
    }

  fprintf(out_f,"\nAllocating for geom boundary elements.\n");
  fflush(out_f);
  
  if ((finalmesh->nelems+geom->n_gfacets) > finalmesh->element_dim)
    my_mem(finalmesh->element_dim, finalmesh->nelems, &(finalmesh->element), (finalmesh->nelems+geom->n_gfacets), ELEMENT_CHUNK);
  
  //since all elements are one face, set up nf and f_n now
  for (k = finalmesh->nelems; k < (finalmesh->nelems+geom->n_gfacets); k++)
    {
    finalmesh->element[k].Initialize();
    finalmesh->element[k].nf = 1;
    finalmesh->element[k].f_n = new List*[1];
    finalmesh->element[k].f_n[0] = new List();
    }
    
  fprintf(out_f,"\nStoring facets as boundary objects.\n");
  fflush(out_f);
  
  //loop through geom facets, noting all unique nodes at this point
  //set up geom as elements, then tris as elements
  for (i=1; i <= geom->ngb; i++)
    {
    for (j=geom->n_begin[i]; j < geom->n_end[i]; j++)
      {
      n0 = geom->g_facet[j].nodes[0];
      n1 = geom->g_facet[j].nodes[1];
      n2 = geom->g_facet[j].nodes[2];
      
      //take the resulting triangles and set up tri
      finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n0);
      finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n1);
      finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n2);
      
      finalmesh->boundary[i-1].elist->Add_To_List(finalmesh->nelems); 
      //inc nelems
      finalmesh->nelems++;
      }
    }
  
  //add in geom nodes
  for (i = 0; i < geom->n_gnodes; i++)
    {
    if ((finalmesh->nn + 1) > finalmesh->node_dim)
      my_mem(finalmesh->node_dim, finalmesh->nn, &(finalmesh->node), (finalmesh->nn + 1), NODE_CHUNK);
      
    finalmesh->node[i].vert = geom->g_vert[i];
    }
      
  //finally, output stl file
  sname[0] = '\0';
  
  sprintf(sname,"TRI2STL.stl");

  if ((ucd_f=fopen(sname,"w")) == NULL)
    {
    fprintf(stderr,"\nCouldn't open surface mesh output file %s",sname);
    fflush(stderr);
    exit(0);
    }
   
  //STL format
  fprintf(ucd_f, "solid for_pointwise\n");
  fflush(ucd_f);
    
  Vector areasum = Vector(0.0,0.0,0.0);
    
  //put in triangles and normals
  for (i = 0; i < finalmesh->nb; i++)
    {
    for (j = 0; j < finalmesh->boundary[i].elist->max; j++)
      {
      n0 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[0];
      n1 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[1];
      n2 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[2];
      
      //find normal vector
      v0 = Vector(finalmesh->node[n0].vert,finalmesh->node[n1].vert);
      v1 = Vector(finalmesh->node[n0].vert,finalmesh->node[n2].vert);
      norm = v0 % v1;
      norm.normalize();
        
      fprintf(ucd_f, "facet normal %16.10e %16.10e %16.10e\n",norm[0],norm[1],norm[2]);
        
      //do area summation for inner bd only..reset norm in process
      norm = (v0%v1)*0.5;
      areasum += norm;
        
      fprintf(ucd_f, "  outer loop\n");
      fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n0].vert[0],finalmesh->node[n0].vert[1],finalmesh->node[n0].vert[2]);
      fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n1].vert[0],finalmesh->node[n1].vert[1],finalmesh->node[n1].vert[2]);
      fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n2].vert[0],finalmesh->node[n2].vert[1],finalmesh->node[n2].vert[2]);
      fprintf(ucd_f, "  endloop\n");
      fprintf(ucd_f, "endfacet\n");
      }
    }
      
  fprintf(ucd_f, "endsolid for_pointwise\n");
  fflush(ucd_f);
    
  fprintf(out_f, "Area summation (x, y, z) = %16.10e %16.10e %16.10e\n",areasum[0],areasum[1],areasum[2]);
  fflush(out_f);
    
  //close ucd file
  fclose(ucd_f);
  
  delete finalmesh;
  
  fprintf(out_f,"\nFinished with output STL file.\n");
  fflush(out_f);

  //NO CAN DO...NO VOL ELEM 
  //reset file name before volume mesh output
  //sname[0] = '\0';
  
  //sprintf(sname,"TRI2STL.cgns");

  //output vol mesh
  //finalmesh->smooth_io(1,geom,sname);
  
  geom->Plotfile(in);
  
  fprintf(out_f,"\nFinished with output Fieldview file.\n");
  fflush(out_f);

  return;
}
