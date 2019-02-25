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

#endif // FLUIDS_H
