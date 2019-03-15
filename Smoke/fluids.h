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

typedef enum vis_data_type {
    DENSITY_RHO,
    VELOCITY_V,
    FORCE_FIELD_F,
    NUMDATATYPES
} vis_data_type;

typedef enum interpol_type {
    NEAREST_NEIGHBOR,
    BILINEAR,
    NUMINTERPOLTYPES
} interpol_type;

typedef enum glyph_type {
    HEDGEHOGS,
    CONES,
    ARROWS,
    NUMGLYPHTYPES
} glyph_type;


namespace fluids {

    extern int color_dir;
    extern int scalar_col;
    extern int draw_vecs;
    extern int enable_smoke;
    extern int bands;
    extern bool enable_bands;
    extern float min_col;
    extern float max_col;
    extern vis_data_type color_data_type;
    extern vis_data_type vector_data_type;
    extern interpol_type interpolation;
    extern glyph_type glyph_shape;
    extern float vec_scale;

    extern int DIM_X;
    extern int DIM_Y;

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
