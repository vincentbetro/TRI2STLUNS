#include "Point.h"

int CGNS_read(char *name, int &nn, Point **p, int &nb, char ***b_name,
              int **nt, int ****tri_conn,
              int **nq, int ****quad_conn,
              int **ngon, int ****ngon_conn,
              int &ntet, int ***tet_conn,
              int &npyr, int ***pyr_conn,
              int &npri, int ***pri_conn,
              int &nhex, int ***hex_conn,
              int &nply, int ****poly_conn);

int CGNS_write(char *name, int nn, Point *p, int nb, char **b_name,
               int *nt, int ***tri_conn,
               int *nq, int ***quad_conn,
               int *ngon, int ***ngon_conn,
               int ntet, int **tet_conn,
               int npyr, int **pyr_conn,
               int npri, int **pri_conn,
               int nhex, int **hex_conn,
               int nply, int ***poly_conn);

int P_CGNS_read(int parallel, char *fname, int &nn, Point **p, int &nb, char ***b_name,
              int **nt, int ****tri_conn,
              int **nq, int ****quad_conn,
              int **ngon, int ****ngon_conn,
              int &ntet, int ***tet_conn,
              int &npyr, int ***pyr_conn,
              int &npri, int ***pri_conn,
              int &nhex, int ***hex_conn,
              int &nply, int ****poly_conn,
              int ***nmap, int ***tri_map, int ***quad_map, int ***ngon_map,
              int **tet_map, int **pyr_map, int **pri_map, int **hex_map, int **poly_map);

int P_CGNS_write(int parallel, char *fname, int nn, Point *p, int nb, char **b_name,
               int *nt, int ***tri_conn,
               int *nq, int ***quad_conn,
               int *ngon, int ***ngon_conn,
               int ntet, int **tet_conn,
               int npyr, int **pyr_conn,
               int npri, int **pri_conn,
               int nhex, int **hex_conn,
               int nply, int ***poly_conn,
               int **nmap, int **tri_map, int **quad_map, int **ngon_map,
               int *tet_map, int *pyr_map, int *pri_map, int *hex_map, int *poly_map);
