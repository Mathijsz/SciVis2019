// Usage: Drag with the mouse to add smoke to the fluid. This will also move a "rotor" that disturbs
//        the velocity field at the mouse location. Press the indicated keys to change options
//--------------------------------------------------------------------------------------------------

#include <rfftw.h>              //the numerical simulation FFTW library
#include <stdio.h>              //for printing the help text
#include <math.h>

#include <QOpenGLFunctions>
#include <QDebug>
#include <QMatrix4x4>
#include <QVector>
#include <QVector2D>
#include <QVector3D>
#include "fluids.h"

namespace fluids {

//--- SIMULATION PARAMETERS ------------------------------------------------------------------------
const int DIM = 50;             //size of simulation grid
int DIM_X = 50;
int DIM_Y = 50;
double dt = 0.4;                //simulation time step
float visc = 0.001;             //fluid viscosity
fftw_real *vx, *vy;             //(vx,vy)   = velocity field at the current moment
fftw_real *vx0, *vy0;           //(vx0,vy0) = velocity field at the previous moment
fftw_real *fx, *fy;             //(fx,fy)   = user-controlled simulation forces, steered with the mouse
fftw_real *rho, *rho0;          //smoke density at the current (rho) and previous (rho0) moment
rfftwnd_plan plan_rc, plan_cr;  //simulation domain discretization

// Divergence of data
fftw_real *divd;

//--- VISUALIZATION PARAMETERS ---------------------------------------------------------------------
int   winWidth, winHeight;      //size of the graphics window, in pixels
int   color_dir = 0;            //use direction color-coding or not
float vec_scale = 1000;         //scaling of hedgehogs
int   enable_smoke = 0;           //draw the smoke or not
int   draw_vecs = 1;            //draw the vector field or not

int   scalar_col = 0;           //method for scalar coloring
int   frozen = 0;               //toggles on/off the animation


bool enable_bands = false;
bool autoscale_colormaps = false;
bool enable_isolines = false;
bool enable_repeats = false;
bool enable_bounded_isolines = false;
bool enable_heightmap = false;

bool enable_shading = false;

float isoline = 0.5;
float upper_isoline = 0.75;
float lower_isoline = 0.25;
int isoline_count = 2;

int height_scale = 20;
float streamtube_scale = 0.1;

int bands = 2;
int repeat_levels = 1;
float min_col = 0.0;
float max_col = 1.0;
vis_data_type color_data_type = DENSITY_RHO;
vis_data_type vector_data_type = VELOCITY_V;
vis_data_type heightmap_data_type = DENSITY_RHO;
interpol_type interpolation = NEAREST_NEIGHBOR;
glyph_type glyph_shape = HEDGEHOGS;

QMatrix4x4 model, view, projection;
QVector2D last_mouse_pos;

QVector3D lighting;
QVector3D material;
QVector<QVector<QVector2D>> streamtubes;
// {{(x0,y0),{x1,y1)}, }

//------ SIMULATION CODE STARTS HERE -----------------------------------------------------------------

//init_simulation: Initialize simulation data structures as a function of the grid size 'n'.
//                 Although the simulation takes place on a 2D grid, we allocate all data structures as 1D arrays,
//                 for compatibility with the FFTW numerical library.
void init_simulation(int n)
{
    int i; size_t dim;

    dim     = n * 2*(n/2+1)*sizeof(fftw_real);        //Allocate data structures
    vx       = (fftw_real*) malloc(dim);
    vy       = (fftw_real*) malloc(dim);
    vx0      = (fftw_real*) malloc(dim);
    vy0      = (fftw_real*) malloc(dim);
    dim     = n * n * sizeof(fftw_real);
    fx      = (fftw_real*) malloc(dim);
    fy      = (fftw_real*) malloc(dim);
    rho     = (fftw_real*) malloc(dim);
    rho0    = (fftw_real*) malloc(dim);
    plan_rc = rfftw2d_create_plan(n, n, FFTW_REAL_TO_COMPLEX, FFTW_IN_PLACE);
    plan_cr = rfftw2d_create_plan(n, n, FFTW_COMPLEX_TO_REAL, FFTW_IN_PLACE);

    divd = (fftw_real*) malloc(dim);

    for (i = 0; i < n * n; i++)                      //Initialize data structures to 0
    { vx[i] = vy[i] = vx0[i] = vy0[i] = fx[i] = fy[i] = rho[i] = rho0[i] = 0.0f; }
}

void destroy_simulation()
{
    free(vx);
    free(vy);
    free(vx0);
    free(vy0);
    free(fx);
    free(fy);
    free(rho);
    free(rho0);
    free(divd);
    rfftwnd_destroy_plan(plan_rc);
    rfftwnd_destroy_plan(plan_cr);
}


//FFT: Execute the Fast Fourier Transform on the dataset 'vx'.
//     'dirfection' indicates if we do the direct (1) or inverse (-1) Fourier Transform
void FFT(int direction,void* vx)
{
    if(direction==1) rfftwnd_one_real_to_complex(plan_rc,(fftw_real*)vx,(fftw_complex*)vx);
    else             rfftwnd_one_complex_to_real(plan_cr,(fftw_complex*)vx,(fftw_real*)vx);
}

int clamp(float x)
{ return ((x)>=0.0?((int)(x)):(-((int)(1-(x))))); }

float max(float x, float y)
{ return x > y ? x : y; }

float min(float x, float y)
{ return x <= y ? x : y; }

//solve: Solve (compute) one step of the fluid flow simulation
void solve(int n, fftw_real* vx, fftw_real* vy, fftw_real* vx0, fftw_real* vy0, fftw_real visc, fftw_real dt)
{
    fftw_real x, y, x0, y0, f, r, U[2], V[2], s, t;
    int i, j, i0, j0, i1, j1;

    //    Testing density fixed point
//        rho0[25+25*n] = 20.0;
//        rho[25+25*n] = 20.0;

    for (i=0;i<n*n;i++)
    { vx[i] += dt*vx0[i]; vx0[i] = vx[i]; vy[i] += dt*vy0[i]; vy0[i] = vy[i]; }

    for ( x=0.5f/n,i=0 ; i<n ; i++,x+=1.0f/n )
       for ( y=0.5f/n,j=0 ; j<n ; j++,y+=1.0f/n )
       {
          x0 = n*(x-dt*vx0[i+n*j])-0.5f;
          y0 = n*(y-dt*vy0[i+n*j])-0.5f;
          i0 = clamp(x0); s = x0-i0;
          i0 = (n+(i0%n))%n;
          i1 = (i0+1)%n;
          j0 = clamp(y0); t = y0-j0;
          j0 = (n+(j0%n))%n;
          j1 = (j0+1)%n;
          vx[i+n*j] = (1-s)*((1-t)*vx0[i0+n*j0]+t*vx0[i0+n*j1])+s*((1-t)*vx0[i1+n*j0]+t*vx0[i1+n*j1]);
          vy[i+n*j] = (1-s)*((1-t)*vy0[i0+n*j0]+t*vy0[i0+n*j1])+s*((1-t)*vy0[i1+n*j0]+t*vy0[i1+n*j1]);
       }

    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            {  vx0[i+(n+2)*j] = vx[i+n*j]; vy0[i+(n+2)*j] = vy[i+n*j]; }

    FFT(1,vx0);
    FFT(1,vy0);

    for (i=0;i<=n;i+=2)
    {
        x = 0.5f*i;
        for (j=0;j<n;j++)
        {
            y = j<=n/2 ? (fftw_real)j : (fftw_real)j-n;
            r = x*x+y*y;
            if ( r==0.0f ) continue;
            f = (fftw_real)exp(-r*dt*visc);
            U[0] = vx0[i  +(n+2)*j]; V[0] = vy0[i  +(n+2)*j];
            U[1] = vx0[i+1+(n+2)*j]; V[1] = vy0[i+1+(n+2)*j];

            vx0[i  +(n+2)*j] = f*((1-x*x/r)*U[0]     -x*y/r *V[0]);
            vx0[i+1+(n+2)*j] = f*((1-x*x/r)*U[1]     -x*y/r *V[1]);
            vy0[i+  (n+2)*j] = f*(  -y*x/r *U[0] + (1-y*y/r)*V[0]);
            vy0[i+1+(n+2)*j] = f*(  -y*x/r *U[1] + (1-y*y/r)*V[1]);
        }
    }

    FFT(-1,vx0);
    FFT(-1,vy0);


    f = 1.0/(n*n);
    for (i=0;i<n;i++)
       for (j=0;j<n;j++)
       { vx[i+n*j] = f*vx0[i+(n+2)*j]; vy[i+n*j] = f*vy0[i+(n+2)*j]; }

    if (color_data_type == DIVERGENCE) {
        fftw_real *dx = vx;
        fftw_real *dy = vy;

        switch (vector_data_type) {
            case FORCE_FIELD_F:
                dx = fx;
                dy = fy;
                break;
            default:
                break;
        }

        // No boundary
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                divd[i+n*j] = dx[((i+1) % n)+n*j] - dx[(i-1+n) % n + n*j] + dy[i+n*((j+1) % n)] - dy[i+n*((j-1+n)%n)];
            }
     }

    for (int t = 0; t < streamtubes.size(); t++) {
        auto prev = streamtubes[t].last();

        streamtubes[t] += {(float)(prev.y() + streamtube_scale*get_data_interpol(&get_vx_idx, prev.y(), prev.x(), DIM_X, DIM_Y)),
                           (float)(prev.x() + streamtube_scale*get_data_interpol(&get_vy_idx, prev.y(), prev.x(), DIM_X, DIM_Y))};
    }

}


// diffuse_matter: This function diffuses matter that has been placed in the velocity field. It's almost identical to the
// velocity diffusion step in the function above. The input matter densities are in rho0 and the result is written into rho.
void diffuse_matter(int n, fftw_real *vx, fftw_real *vy, fftw_real *rho, fftw_real *rho0, fftw_real dt)
{
    fftw_real x, y, x0, y0, s, t;
    int i, j, i0, j0, i1, j1;

    for ( x=0.5f/n,i=0 ; i<n ; i++,x+=1.0f/n )
        for ( y=0.5f/n,j=0 ; j<n ; j++,y+=1.0f/n )
        {
            x0 = n*(x-dt*vx[i+n*j])-0.5f;
            y0 = n*(y-dt*vy[i+n*j])-0.5f;
            i0 = clamp(x0);
            s = x0-i0;
            i0 = (n+(i0%n))%n;
            i1 = (i0+1)%n;
            j0 = clamp(y0);
            t = y0-j0;
            j0 = (n+(j0%n))%n;
            j1 = (j0+1)%n;
            rho[i+n*j] = (1-s)*((1-t)*rho0[i0+n*j0]+t*rho0[i0+n*j1])+s*((1-t)*rho0[i1+n*j0]+t*rho0[i1+n*j1]);
        }
}

//set_forces: copy user-controlled forces to the force vectors that are sent to the solver.
//            Also dampen forces and matter density to get a stable simulation.
void set_forces(void)
{
    int i;
    for (i = 0; i < DIM * DIM; i++)
    {
        rho0[i]  = 0.995 * rho[i];
        fx[i] *= 0.85;
        fy[i] *= 0.85;
        vx0[i]    = fx[i];
        vy0[i]    = fy[i];
    }
}


//do_one_simulation_step: Do one complete cycle of the simulation:
//      - set_forces:
//      - solve:            read forces from the user
//      - diffuse_matter:   compute a new set of velocities
//      - gluPostRedisplay: draw a new visualization frame
void do_one_simulation_step(void)
{
    if (!frozen)
    {
        set_forces();
        solve(DIM, vx, vy, vx0, vy0, visc, dt);
        diffuse_matter(DIM, vx, vy, rho, rho0, dt);
#ifdef USE_GLUT
        glutPostRedisplay();
#endif
    }
}


//------ VISUALIZATION CODE STARTS HERE -----------------------------------------------------------------

float rescale(float value)
{
    if (min_col >= max_col) {
        min_col = max_col - 0.0001;
    }
    if (value < min_col)
        value=min_col;
    if (value > max_col)
        value=max_col;
    return (value - min_col) / (max_col-min_col);
}

//rainbow: Implements a color palette, mapping the scalar 'value' to a rainbow color RGB
void rainbow(float value,float* R,float* G,float* B)
{
    const float dx=0.8;
    value = (6-2*dx)*value+dx;
    *R = max(0.0,(3-fabs(value-4)-fabs(value-5))/2);
    *G = max(0.0,(4-fabs(value-2)-fabs(value-4))/2);
    *B = max(0.0,(3-fabs(value-1)-fabs(value-2))/2);
}

void with_banding(color_func f, float value, float *R, float *G, float *B, int banding_levels)
{
    value = rescale(value);
    if (enable_bands && banding_levels > 0) {
        banding_levels--;
        value *= banding_levels;
        value = (int)value;
        value /= banding_levels;
    }

    if (enable_repeats) {
        value = round(fmod(value * repeat_levels, 1));
    }

    (*f)(value, R, G, B);
}

void red_to_white(float value, float *R, float *G, float *B)
{
    *R = 1;
    *G = value;
    *B = value;
}

void blue_to_red_via_white(float value, float *R, float *G, float *B)
{
    *R = min(1.0, 2*value);
    *G = 1-fabs(value*2-1);
    *B = min(1.0, 2 - 2*value);
}

void white_to_black(float value, float *R, float *G, float *B)
{
    *R = *G = *B = value;
}

color_func get_color_func(colormap col)
{
    switch (col) {
        case COLOR_RAINBOW:
             return &rainbow;
        case COLOR_RED_TO_WHITE:
             return &red_to_white;
        case COLOR_BLUE_TO_RED_VIA_WHITE:
             return &blue_to_red_via_white;
        case COLOR_BLACKWHITE:
        default:
             return &white_to_black;
    }
}

//set_colormap: Sets three different types of colormaps
void set_colormap(float vy)
{
    float R,G,B;
    color_func f = get_color_func((colormap)scalar_col);
    with_banding(f, vy, &R, &G, &B, bands);
    glColor3f(R,G,B);
}

QVector3D get_color(float vy)
{
    float R,G,B;
    color_func f = get_color_func((colormap)scalar_col);
    with_banding(f, vy, &R, &G, &B, bands);
    return {R, G, B};
}

fftw_real get_vis_data(int idx)
{
    switch (color_data_type) {
        case VELOCITY_V:
            return sqrt(vx[idx]*vx[idx] + vy[idx]*vy[idx]);
        case FORCE_FIELD_F:
            return sqrt(fx[idx]*fx[idx] + fy[idx]*fy[idx]);
        case DIVERGENCE:
            return divd[idx];
        case DENSITY_RHO:
        default:
            return rho[idx];
    }
}

fftw_real get_height_data(int idx)
{
    switch (heightmap_data_type) {
        case VELOCITY_V:
            return sqrt(vx[idx]*vx[idx] + vy[idx]*vy[idx]);
        case FORCE_FIELD_F:
            return sqrt(fx[idx]*fx[idx] + fy[idx]*fy[idx]);
        case DIVERGENCE:
            return divd[idx];
        case DENSITY_RHO:
        default:
            return rho[idx];
    }
}

fftw_real get_vx_idx(int idx)
{
    return vx[idx];
}

fftw_real get_vy_idx(int idx)
{
    return vy[idx];
}

fftw_real get_vec_data_x(int idx)
{
    switch (vector_data_type) {
        case VELOCITY_V:
            return vx[idx];
        case FORCE_FIELD_F:
        default:
            return fx[idx];
    }
}

fftw_real get_vec_data_y(int idx)
{
    switch (vector_data_type) {
        case VELOCITY_V:
            return vy[idx];
        case FORCE_FIELD_F:
        default:
            return fy[idx];
    }
}

float get_interpolated_value(float q1, float q2, float dx)
{
    return (q2 - q1)*dx + q1;
}

fftw_real get_data_interpol(fftw_real (*f)(int), float y, float x, int dimx, int dimy)
{
    x = (x / dimx) * DIM;
    y = (y / dimy) * DIM;
    switch (interpolation) {
        case BILINEAR: {
            int x1 = floorf(x);
            int y1 = floorf(y);
            int x2 = ceilf(x);
            int y2 = ceilf(y);
            float r1 = get_interpolated_value(f(x1 + DIM * y1), f(x1 + DIM * y2), fmod(y, 1));
            float r2 = get_interpolated_value(f(x2 + DIM * y1), f(x2 + DIM * y2), fmod(y, 1));
            return get_interpolated_value(r1, r2, fmod(x, 1));
        }
        case NEAREST_NEIGHBOR:
        default:
            return f(roundf(y) * DIM + roundf(x));
        }
}

//direction_to_color: Set the current color by mapping a direction vector (x,y), using
//                    the color mapping method 'method'. If method==1, map the vector direction
//                    using a rainbow colormap. If method==0, simply use the white color
void direction_to_color(int y, int x, int method, int dimx, int dimy)
{
    float r,g,b,f;
    if (method)
    {
        float val_x = get_data_interpol(&get_vec_data_x, y, x, dimx, dimy);
        float val_y = get_data_interpol(&get_vec_data_y, y, x, dimx, dimy);

        f = atan2(val_y,val_x) / 3.1415927 + 1;
        r = f;
        if(r > 1) r = 2 - r;
        g = f + .66667;
        if(g > 2) g -= 2;
        if(g > 1) g = 2 - g;
        b = f + 2 * .66667;
        if(b > 2) b -= 2;
        if(b > 1) b = 2 - b;
    }
    else
    {
        color_func func = get_color_func((colormap)scalar_col);
        with_banding(func, get_data_interpol(&get_vis_data, y, x, dimx, dimy), &r, &g, &b, bands);
    }
    glColor3f(r,g,b);
}

// Given 2 values of points and 1 value in between, return the relative coord of the value
float value_to_coord(float a, float b, float c)
{
    if (b-a == 0) return 0;
    return fabs(c - a) / fabs(b - a);
}

void draw_isolines(fftw_real hn, fftw_real wn, float isoline_value, int dimx, int dimy)
{
    int i, j, k;
    float R, G, B;
    glBegin(GL_LINES);
    (get_color_func((colormap)scalar_col))(rescale(isoline_value), &R, &G, &B);
    glColor3f(R, G, B);
    fftw_real points[4];

    for (i = 0; i < dimx - 1; i++) {

        points[0] = get_data_interpol(&get_vis_data, 0, i, dimx, dimy);
        points[3] = get_data_interpol(&get_vis_data, 0, i+1, dimx, dimy);

        for (j = 0; j < dimy - 1; j++)
        {
            unsigned char code = 0;
            points[1] = get_data_interpol(&get_vis_data, j+1, i, dimx, dimy);
            points[2] = get_data_interpol(&get_vis_data, j+1, i+1, dimx, dimy);
            for (k = 0; k < 4; k++) {
                if (points[k] > isoline_value)
                    code += pow(2, k);
            }
            std::vector<float> vertex_x;
            std::vector<float> vertex_y;

            switch (code) {
            case 1:
            case 14:
                vertex_x.push_back(j);
                vertex_y.push_back(i + value_to_coord(points[0], points[3], isoline_value));
                vertex_x.push_back(j + value_to_coord(points[0], points[1], isoline_value));
                vertex_y.push_back(i);
                break;
            case 2:
            case 13:
                vertex_x.push_back(j + value_to_coord(points[0], points[1], isoline_value));
                vertex_y.push_back(i);
                vertex_x.push_back(j+1);
                vertex_y.push_back(i + value_to_coord(points[1], points[2], isoline_value));
                break;
            case 3:
            case 12:
                vertex_x.push_back(j);
                vertex_y.push_back(i + value_to_coord(points[0], points[3], isoline_value));
                vertex_x.push_back(j+1);
                vertex_y.push_back(i + value_to_coord(points[1], points[2], isoline_value));
                break;
            case 4:
            case 11:
                vertex_x.push_back(j + 1);
                vertex_y.push_back(i + value_to_coord(points[1], points[2], isoline_value));
                vertex_x.push_back(j + value_to_coord(points[3], points[2], isoline_value));
                vertex_y.push_back(i + 1);
                break;
            case 5:
            case 10:
                vertex_x.push_back(j + value_to_coord(points[0], points[1], isoline_value));
                vertex_y.push_back(i);
                vertex_x.push_back(j);
                vertex_y.push_back(i + value_to_coord(points[0], points[3], isoline_value));

                vertex_x.push_back(j + 1);
                vertex_y.push_back(i + value_to_coord(points[1], points[2], isoline_value));
                vertex_x.push_back(j + value_to_coord(points[3], points[2], isoline_value));
                vertex_y.push_back(i + 1);
                break;
            case 6:
            case 9:
                vertex_x.push_back(j + value_to_coord(points[0], points[1], isoline_value));
                vertex_y.push_back(i);
                vertex_x.push_back(j + value_to_coord(points[3], points[2], isoline_value));
                vertex_y.push_back(i + 1);
                break;
            case 7:
            case 8:
                vertex_x.push_back(j);
                vertex_y.push_back(i + value_to_coord(points[0], points[3], isoline_value));
                vertex_x.push_back(j + value_to_coord(points[3], points[2], isoline_value));
                vertex_y.push_back(i + 1);
                break;
            }

            for (int l = 0; l < vertex_x.size(); l++) {
                float py = wn + (fftw_real)vertex_y[l]*wn;
                float px = hn + (fftw_real)vertex_x[l]*hn;
                float height = enable_heightmap ? height_scale * get_data_interpol(&get_height_data, vertex_x[l], vertex_y[l], dimx, dimy) : 0;
                glVertex3f(py, px, height);
            }
            points[0] = points[1];
            points[3] = points[2];
        }
    }
    glEnd();
}

void draw_hedgehogs(fftw_real hn, fftw_real wn, int dimx, int dimy, float hedge_scale = 1)
{
    int i, j;
    glBegin(GL_LINES);              //draw velocities
    for (i = 0; i < dimx; i++)
        for (j = 0; j < dimy; j++)
        {

            float height = enable_heightmap ? height_scale * get_data_interpol(&get_height_data, j, i, dimx, dimy) : 0;
            direction_to_color(j, i, color_dir, dimx, dimy);
            glVertex3f(wn + (fftw_real)i * wn, hn + (fftw_real)j * hn, height);
            glVertex3f((wn + (fftw_real)i * wn) + vec_scale * hedge_scale * get_data_interpol(&get_vec_data_x, j, i, dimx, dimy),
                       (hn + (fftw_real)j * hn) + vec_scale * hedge_scale * get_data_interpol(&get_vec_data_y, j, i, dimx, dimy), height);
        }
    glEnd();
}


void draw_cones(fftw_real hn, fftw_real wn, int dimx, int dimy, float offset = 0)
{
    int i, j;
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBegin(GL_TRIANGLES);
    for (i = 0; i < dimx; i++)
        for (j = 0; j < dimy; j++)
        {
            direction_to_color(j, i, color_dir, dimx, dimy);
            fftw_real x = get_data_interpol(&get_vec_data_x, j, i, dimx, dimy);
            fftw_real y = get_data_interpol(&get_vec_data_y, j, i, dimx, dimy);
            float height = enable_heightmap ? height_scale * get_data_interpol(&get_height_data, j, i, dimx, dimy) : 0;
            glVertex3f((wn + (fftw_real)i * wn) + vec_scale * (x*0.2 + offset * x),
                       (hn + (fftw_real)j * hn) + vec_scale * (y*0.2 + offset * y), height);
            glVertex3f((wn + (fftw_real)i * wn) + vec_scale * (y*0.05 + offset * x),
                       (hn + (fftw_real)j * hn) + vec_scale * (-x*0.05 + offset * y), height);
            glVertex3f((wn + (fftw_real)i * wn) + vec_scale * (-y*0.05 + offset * x),
                       (hn + (fftw_real)j * hn) + vec_scale * (x*0.05 + offset * y), height);
        }
    glEnd();
}

void draw_arrows(fftw_real hn, fftw_real wn, int dimx, int dimy)
{
    draw_hedgehogs(hn, wn, dimx, dimy,  0.8f);
    draw_cones(hn, wn, dimx, dimy, 0.8f);
}

void draw_smoke(fftw_real hn, fftw_real wn, int dimx, int dimy)
{
    int i, j;
    float val0, val1, val2, val3;
    double px0, py0, px1, py1, px2, py2, px3, py3;
    QVector<QVector3D> vertices, surface_normals, vertex_normals, colors;

    int nvertices = (dimx-1)*(dimy-1)*6;
    vertices.reserve(nvertices);
    surface_normals.reserve(nvertices/3);
    vertex_normals.resize(nvertices);
    colors.reserve(nvertices);


    for (j = 0; j < dimy - 1; j++)            //draw smoke
    {
        for (i = 0; i < dimx - 1; i++)
        {
            px0 = wn + (fftw_real)i * wn;
            py0 = hn + (fftw_real)j * hn;
            val0 = get_data_interpol(&get_vis_data, j, i, dimx, dimy);


            px1 = wn + (fftw_real)i * wn;
            py1 = hn + (fftw_real)(j + 1) * hn;
            val1 = get_data_interpol(&get_vis_data, (j + 1), i, dimx, dimy);


            px2 = wn + (fftw_real)(i + 1) * wn;
            py2 = hn + (fftw_real)(j + 1) * hn;
            val2 = get_data_interpol(&get_vis_data, j + 1, i + 1, dimx, dimy);


            px3 = wn + (fftw_real)(i + 1) * wn;
            py3 = hn + (fftw_real)j * hn;
            val3 = get_data_interpol(&get_vis_data, j, i + 1, dimx, dimy);

            float h0, h1, h2, h3;
            h0 = h1 = h2 = h3 = 0;

            if (enable_heightmap) {
                h0 = height_scale * get_data_interpol(&get_height_data, j, i, dimx, dimy);
                h1 = height_scale * get_data_interpol(&get_height_data, j+1, i, dimx, dimy);
                h2 = height_scale * get_data_interpol(&get_height_data, j+1, i+1, dimx, dimy);
                h3 = height_scale * get_data_interpol(&get_height_data, j, i+1, dimx, dimy);

            }
            vertices.push_back({(float)px0, (float)py0, (float)h0});
            vertices.push_back({(float)px1, (float)py1, (float)h1});
            vertices.push_back({(float)px2, (float)py2, (float)h2});

            surface_normals +=
                QVector3D::crossProduct(
                    vertices.end()[-2] - vertices.end()[-3],
                    vertices.end()[-1] - vertices.end()[-2]
                ).normalized() * -1;

            vertices.push_back({(float)px0, (float)py0, (float)h0});
            vertices.push_back({(float)px2, (float)py2, (float)h2});
            vertices.push_back({(float)px3, (float)py3, (float)h3});

            surface_normals +=
                QVector3D::crossProduct(
                    vertices.end()[-2] - vertices.end()[-3],
                    vertices.end()[-1] - vertices.end()[-2]
                ).normalized() * -1;

            colors.push_back(get_color(val0));
            colors.push_back(get_color(val1));
            colors.push_back(get_color(val2));

            colors.push_back(get_color(val0));
            colors.push_back(get_color(val2));
            colors.push_back(get_color(val3));

        }
    }

    // Iterate over every lower-left vertex point, then:
    // - find 6 adjacent face normals
    // - average them to obtain vertex normal
    // - put vertex normal at indices that share the same normal.
    // Indexing logic somewhat convoluted

    int dx = 2*(dimx-1);
    int dy = dx*(dimy-1);

    for (i = 0; i < surface_normals.size(); i += 2) {
        QVector3D average = surface_normals[i] + surface_normals[i+1] + surface_normals[i % dx != 0 ? i-1 : i-1+dx];
        int wrap_idx = i - dx + (i < dx ? dy : 0);
        average += surface_normals[wrap_idx] + surface_normals[wrap_idx - 1] + surface_normals[wrap_idx - 2];
        average /= 6;
        average.normalize();

        int xwrap = i % dx ? i*3 - 1 : i*3 - 1 + 3*dx;
        vertex_normals[i*3] = average;
        vertex_normals[i*3 + 3] = average;
        vertex_normals[xwrap] = average;

        int ywrap = i*3 - 3*dx + 1 + (i < dx ? 3*dy : 0);
        ywrap = (ywrap + 3*dy) % (3*dy);

        int xwrap_ywrap = ywrap + (i % dx == 0 ? 3*dx : 0);

        vertex_normals[ywrap] = average;
        vertex_normals[xwrap_ywrap - 3] = average;
        vertex_normals[xwrap_ywrap - 5] = average;

    }

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBegin(GL_TRIANGLES);
    for (i = 0; i < nvertices; i++) {
        glNormal3f(vertex_normals[i].x(), vertex_normals[i].y(), vertex_normals[i].z());
        glColor3f(colors[i].x(), colors[i].y(), colors[i].z());
        glVertex3f(vertices[i].x(), vertices[i].y(), vertices[i].z());
    }
    glEnd();

}

void scale_colormap()
{
    fftw_real n, current_min = INFINITY, current_max = -INFINITY;
    for (int i = 0; i < DIM*DIM; i++) {
        n = get_vis_data(i);
        if (n < current_min) current_min = n;
        else if (n > current_max) current_max = n;
    }
    min_col = current_min;
    max_col = current_max;

}

void draw_mouse(int mx, int my, int dimx, int dimy)
{
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBegin(GL_TRIANGLES);
    glColor3f(1,(float)20/255,(float)147/255);
    float height = 10 + (enable_heightmap ? height_scale * get_data_interpol(&get_height_data, my, mx, dimx, dimy) : 0.0);

    float x = (float)mx;
    float y = (float)my;
    QVector<QVector3D> points = {
        {x, y, height},
        {x+5, y+5, height+7.5f},
        {x-5, y+5, height+7.5f},

        {x, y, height},
        {x+5, y-5, height+7.5f},
        {x+5, y+5, height+7.5f},

        {x-5, y-5, height+7.5f},
        {x+5, y-5, height+7.5f},
        {x, y, height},

        {x-5, y+5, height+7.5f},
        {x-5, y-5, height+7.5f},
        {x, y, height},

        {x+5, y-5, height+7.5f},
        {x-5, y+5, height+7.5f},
        {x+5, y+5, height+7.5f},

        {x-5, y-5, height+7.5f},
        {x-5, y+5, height+7.5f},
        {x+5, y-5, height+7.5f}
    };

    QVector<QVector3D> normals;
    for (int i = 0; i < points.size(); i += 3) {
        normals += -1 * QVector3D::crossProduct(points[i+1] - points[i], points[i+2] - points[i+1]).normalized();
    }

    for (int i = 0; i < points.size(); i++) {
        glNormal3f(normals[i/3].x(), normals[i/3].y(), normals[i/3].z());
        glVertex3f(points[i].x(), points[i].y(), points[i].z());
    }

    glEnd();
}

void draw_streamtubes()
{
    int points = 20;
    float radius = 10;
    for (int t = 0; t < streamtubes.size(); t++) {
        glBegin(GL_TRIANGLE_STRIP);
        streamtubes[t];
        for (int p = 0; p < points; p++) {

        }
        glEnd();
    }

}

//visualize: This is the main visualization function
void visualize(int dimx, int dimy)
{
    fftw_real  wn = (fftw_real)winWidth / (fftw_real)(dimx + 1);   // Grid cell width
    fftw_real  hn = (fftw_real)winHeight / (fftw_real)(dimy + 1);  // Grid cell heigh

    if (autoscale_colormaps)
        scale_colormap();

    if (enable_smoke)
        draw_smoke(hn, wn, dimx, dimy);

    if (draw_vecs) {
        switch (glyph_shape) {
            case CONES:
                draw_cones(hn, wn, dimx, dimy);
                break;
            case ARROWS:
                draw_arrows(hn, wn, dimx, dimy);
                break;
            case HEDGEHOGS:
            default:
                draw_hedgehogs(hn, wn, dimx, dimy);
        }
    }

    if (enable_isolines)
        draw_isolines(hn, wn, isoline, dimx, dimy);

    if (enable_bounded_isolines && lower_isoline < upper_isoline)
        for (float l = lower_isoline; l <= upper_isoline; l += (upper_isoline - lower_isoline)/(isoline_count-1))
            draw_isolines(hn, wn, l, dimx, dimy);

    draw_mouse(last_mouse_pos.x(), winHeight- last_mouse_pos.y(), dimx, dimy);

}


//------ INTERACTION CODE STARTS HERE -----------------------------------------------------------------

//display: Handle window redrawing events. Simply delegates to visualize().
void display(int dimx, int dimy)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMultMatrixf((view * model).data());

    visualize(dimx, dimy);
    glFlush();
}

//reshape: Handle window resizing (reshaping) events
void reshape(int w, int h)
{
    winWidth = w; winHeight = h;
    initialize_env();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glMultMatrixf((view * model).data());

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    projection.perspective(60, (float)winWidth/winHeight, -10, 10);
    glMultMatrixf(projection.data());
//    glOrtho(0.0, (GLdouble)w, 0.0, (GLdouble)h, -10, 50);
    glFrustum(0.0, (GLdouble)w, 0.0, (GLdouble)h, -1, 1);

    glViewport(0.0f, 0.0f, (GLfloat)w, (GLfloat)h);
}

void reset_viewpoint()
{
    model.setToIdentity();
    view.setToIdentity();
    projection.setToIdentity();
    view.lookAt({(float)winWidth/2, (float)winHeight/2, (float)((winHeight/2)/(sqrt(3)/3)) }, {(float)winWidth/2, (float)winHeight/2, 0}, {0,1,0});
    lighting = {0,0,300};
}

void initialize_env()
{
    reset_viewpoint();
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glDepthMask(GL_TRUE);
    glDepthRange(1, -1);
    glClearDepth(1.0);

    toggle_shading();
}

//keyboard: Handle key presses
void keyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 't': dt -= 0.001; break;
        case 'T': dt += 0.001; break;
        case 'c': color_dir = 1 - color_dir; break;
        case 'S': vec_scale *= 1.2; break;
        case 's': vec_scale *= 0.8; break;
        case 'V': visc *= 5; break;
        case 'v': visc *= 0.2; break;
        case 'x': enable_smoke = 1 - enable_smoke;
//			if (draw_smoke==0) draw_vecs = 1; break;
        case 'y': draw_vecs = 1 - draw_vecs;
//			if (draw_vecs==0) draw_smoke = 1; break;
        case 'm': scalar_col++; if (scalar_col>=NUMCOLS) scalar_col=COLOR_BLACKWHITE; break;
        case 'a': frozen = 1-frozen; break;
        case 'r': reset_viewpoint(); break;
        case 'q': exit(0);
    }
}



// drag: When the user drags with the mouse, add a force that corresponds to the direction of the mouse
//       cursor movement. Also inject some new matter into the field at the mouse location.
void drag(int mx, int my)
{
    int xi,yi,X,Y; double  dx, dy, len;
    static int lmx=0,lmy=0;             //remembers last mouse location

    // Compute the array index that corresponds to the cursor location
    xi = (int)clamp((double)(DIM + 1) * ((double)mx / (double)winWidth));
    yi = (int)clamp((double)(DIM + 1) * ((double)(winHeight - my) / (double)winHeight));

    X = xi; Y = yi;

    if (X > (DIM - 1))  X = DIM - 1; if (Y > (DIM - 1))  Y = DIM - 1;
    if (X < 0) X = 0; if (Y < 0) Y = 0;

    // Add force at the cursor location
    my = winHeight - my;
    dx = mx - lmx; dy = my - lmy;
    len = sqrt(dx * dx + dy * dy);
    if (len != 0.0) {  dx *= 0.1 / len; dy *= 0.1 / len; }
    fx[Y * DIM + X] += dx;
    fy[Y * DIM + X] += dy;
    rho[Y * DIM + X] = 10.0f;
    lmx = mx; lmy = my;
}

#ifdef USE_GLUT
//main: The main program
int main(int argc, char **argv)
{
    printf("Fluid Flow Simulation and Visualization\n");
    printf("=======================================\n");
    printf("Click and drag the mouse to steer the flow!\n");
    printf("T/t:   increase/decrease simulation timestep\n");
    printf("S/s:   increase/decrease hedgehog scaling\n");
    printf("c:     toggle direction coloring on/off\n");
    printf("V/v:   increase decrease fluid viscosity\n");
    printf("x:     toggle drawing matter on/off\n");
    printf("y:     toggle drawing hedgehogs on/off\n");
    printf("m:     toggle thru scalar coloring\n");
    printf("a:     toggle the animation on/off\n");
    printf("q:     quit\n\n");

    fflush(stdout);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(500,500);
    glutCreateWindow("Real-time smoke simulation and visualization");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(do_one_simulation_step);
    glutKeyboardFunc(keyboard);
    glutMotionFunc(drag);

    init_simulation(DIM);   //initialize the simulation data structures
    glutMainLoop();         //calls do_one_simulation_step, keyboard, display, drag, reshape
    return 0;
}
#endif

void toggle_shading()
{
    if (glIsEnabled(GL_LIGHTING))
        return;

    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    float light_ambient[] = {0.1, 0.1, 0.1, 1.0};
    float light_diffuse[] = {1.0, 1.0, 1.0, 1.0};
    float light_specular[] = {1.0, 1.0, 1.0, 1.0};

    float ambient[] = {1.0f, 1.0f, 1.0f, 1.0f};
    float diffuse[] = {0.1f, 0.5f, 0.7f, 1.0f};
    float specular[] = {0.8f, 0.6f, 0.7f, 1.0f};
    float data[4] = {lighting[0], lighting[1], lighting[2], 1.0};
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, data);

    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.0f);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);

    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);
}

void rotate(int dx, int dy)
{
    view.translate(winWidth/2, winHeight/2, 0);
    view.rotate(QQuaternion::fromAxisAndAngle(QVector3D(1, 0, 0), dy));
    view.rotate(QQuaternion::fromAxisAndAngle(QVector3D(0, 0, 1), dx));
    view.translate(-winWidth/2,- winHeight/2, 0);
}

void scale(float scale)
{
    view.translate(winWidth/2, winHeight/2, 0);
    view.scale(scale);
    view.translate(-winWidth/2,- winHeight/2, 0);
}

void add_seed_point(int mx, int my)
{
    streamtubes += {{(float)mx,(float)my}};
}

void reset_seed_points()
{
    streamtubes.clear();
}

// End of namespace 'fluids'
}
