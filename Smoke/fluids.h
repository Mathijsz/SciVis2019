#ifndef FLUIDS_H
#define FLUIDS_H

namespace fluids {

    extern int color_dir;
    extern int scalar_col;
    extern int draw_vecs;
    extern int draw_smoke;

    void init_simulation(int n);
    void display(void);
    void reshape(int w, int h);
    void do_one_simulation_step(void);
    void drag(int mx, int my);
    void keyboard(unsigned char key, int x, int y);
}

//different types of color mapping: black-and-white, rainbow, banded
enum colormap {
    COLOR_BLACKWHITE,
    COLOR_RAINBOW,
    COLOR_BANDS,
    COLOR_RED_TO_WHITE,
    COLOR_BLUE_TO_YELLOW,
    COLOR_BLUE_TO_RED_VIA_WHITE,
    NUMCOLS
};

#endif // FLUIDS_H
