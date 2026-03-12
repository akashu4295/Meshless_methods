#include "structures.h"
#include "kdtree.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <stdio.h>

static double dist_sq(double *a1, double *a2)
{
    double d0 = a1[0] - a2[0];
    double d1 = a1[1] - a2[1];
    double d2 = a1[2] - a2[2];
    return d0*d0 + d1*d1 + d2*d2;
}

// KD-tree creation

void* create_kdtree(PointStructure* ps)
{
    void* ptree = kd_create(3);
    for (int i = 0; i < ps->num_nodes; i++) {
        if (!ps->corner_tag[i]) {
            assert(kd_insert3(ptree, ps->x[i], ps->y[i], ps->z[i],
                &ps->point_index[i]) == 0);
        }
    }
    return ptree;
}

void* create_kdtree_without_boundarynodes(PointStructure* ps)
{
    void* ptree = kd_create(3);
    for (int i = 0; i < ps->num_nodes; i++) {
        if (!ps->boundary_tag[i]) {
            assert(kd_insert3(ptree, ps->x[i], ps->y[i], ps->z[i],
                           &ps->point_index[i]) == 0);
        }
    }
    return ptree;
}

void free_kdtree(void* ptree)
{
    kd_free(ptree);
}

// Neighbour search

int* find_neighbours(double* p, void* ptree, double radius, int num_cloud_points)
{
    const double INF = DBL_MAX;
    const int MAX_EXPANDS = 50;

    double pos[3], dist;
    struct kdres *presults;

    double* distance = malloc(num_cloud_points * sizeof(double));
    int* ind = malloc(num_cloud_points * sizeof(int));

    if (!distance || !ind) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < num_cloud_points; i++) {
        distance[i] = INF;
        ind[i] = -1;
    }

    for (int expand = 0; expand < MAX_EXPANDS; expand++) {
        presults = kd_nearest_range(ptree, p, radius);
        while (!kd_res_end(presults)) {
            int* pch = (int*)kd_res_item(presults, pos);
            dist = sqrt(dist_sq(p, pos));

            for (int j = 0; j < num_cloud_points; j++) {
                if (dist < distance[j]) {
                    for (int k = num_cloud_points - 1; k > j; k--) {
                        distance[k] = distance[k - 1];
                        ind[k] = ind[k - 1];
                    }
                    distance[j] = dist;
                    ind[j] = *pch;
                    break;
                }
            }
            kd_res_next(presults);
        }
        kd_res_free(presults);
        if (distance[num_cloud_points - 1] < INF) {
            free(distance);
            return ind;
        }

        radius *= 1.5;
    }

    fprintf(stderr, "KD-tree search failed: insufficient neighbors\n");
    free(distance);
    free(ind);
    return NULL;
}

//   Cloud construction

void* create_kdtree_no_corners(PointStructure* ps)
{
    void* ptree = kd_create(3);
    for (int i = 0; i < ps->num_nodes; i++) {
        if (!ps->corner_tag[i]) {
            kd_insert3(ptree, ps->x[i], ps->y[i], ps->z[i],
                       &ps->point_index[i]);
        }
    }
    return ptree;
}

void* create_kdtree_interior_only(PointStructure* ps)
{
    void* ptree = kd_create(3);
    for (int i = 0; i < ps->num_nodes; i++) {
        if (!ps->boundary_tag[i] && !ps->corner_tag[i]) {
            kd_insert3(ptree, ps->x[i], ps->y[i], ps->z[i],
                       &ps->point_index[i]);
        }
    }
    return ptree;
}

void find_cloud_index(PointStructure* ps)
{
    int n = ps->num_cloud_points;
    int N = ps->num_nodes;

    ps->cloud_index = malloc(N * n * sizeof(int));
    if (!ps->cloud_index) abort();

    double radius = ps->d_avg * n;
    double pt[3];

    /* -------- Interior nodes (exclude corners) -------- */
    void* ptree_all = create_kdtree_no_corners(ps);

    for (int i = 0; i < N; i++) {
        if (ps->corner_tag[i]) continue;
        if (ps->boundary_tag[i]) continue;

        pt[0] = ps->x[i];
        pt[1] = ps->y[i];
        pt[2] = ps->z[i];

        int* neigh = find_neighbours(pt, ptree_all, radius, n);

        for (int j = 0; j < n; j++)
            ps->cloud_index[i*n + j] = neigh[j];

        free(neigh);
    }
    free_kdtree(ptree_all);

    // Boundary nodes → interior-only cloud 
    void* ptree_int = create_kdtree_interior_only(ps);

    for (int i = 0; i < ps->num_nodes; i++) {
        if (ps->corner_tag[i]) continue;
        if (!ps->boundary_tag[i]) continue;

        pt[0] = ps->x[i];
        pt[1] = ps->y[i];
        pt[2] = ps->z[i];

        int* neigh = find_neighbours(pt, ptree_int, radius, n - 1);

        ps->cloud_index[i*n] = i;  /* self */
        for (int j = 1; j < n; j++)
            ps->cloud_index[i*n + j] = neigh[j - 1];

        free(neigh);
    }
    free_kdtree(ptree_int);
}

// Nearest-point mapping 

int* find_nearest_point(PointStructure* ps1,
                        PointStructure* ps2,
                        int num_cloud_points)
{
    int* neighbour = malloc(ps1->num_nodes * sizeof(int));
    if (!neighbour) abort();

    void* ptree = create_kdtree(ps2);
    double radius = ps2->d_avg * 10.0;
    double pt[3];

    for (int i = 0; i < ps1->num_nodes; i++) {
        if (!ps1->corner_tag[i]) {
            pt[0] = ps1->x[i];
            pt[1] = ps1->y[i];
            pt[2] = ps1->z[i];

            int* tmp = find_neighbours(pt, ptree, radius, 3);
            neighbour[i] = tmp[0];
            free(tmp);
        } else {
            neighbour[i] = -1;
        }
    }

    free_kdtree(ptree);
    return neighbour;
}