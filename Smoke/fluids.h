#ifndef FLUIDS_H
#define FLUIDS_H

typedef void (*color_func)(float, float*, float*, float*);

//different types of color mapping: black-and-white, rainbow, banded
typedef enum colormap {
    COLOR_BLACKWHITE,
    COLOR_RAINBOW,
    COLOR_RED_TO_WHITE,
    COLOR_BLUE_TO_YELLOW,
    COLOR_BLUE_TO_RED_VIA_WHITE,
    NUMCOLS
} colormap;

namespace fluids {

    extern int color_dir;
    extern int scalar_col;
    extern int draw_vecs;
    extern int draw_smoke;
    extern int bands;

    void init_simulation(int n);
    void display(void);
    void reshape(int w, int h);
    void do_one_simulation_step(void);
    void drag(int mx, int my);
    void keyboard(unsigned char key, int x, int y);

    void with_banding(color_func f, float value, float* R,float* G,float* B, int levels);
    color_func get_color_func(colormap col);
}

#endif // FLUIDS_H
