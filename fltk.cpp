/*#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Box.H>
#include <iostream>

int main(int argc, char **argv)
{
    std::cout << "main() invoked" << std::endl;
    Fl_Window *window = new Fl_Window(340, 180);
    Fl_Box *box = new Fl_Box(20, 40, 300, 100, "Hello, World!");
    box->box(FL_UP_BOX);
    box->labelfont(FL_BOLD + FL_ITALIC);
    box->labelsize(36);
    box->labeltype(FL_SHADOW_LABEL);
    window->end();
    window->show(argc, argv);
    return Fl::run();
}*/
#include <FL/Fl.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Window.H>

int main()
{
    Fl_Window window(200, 200, "Window title");
    Fl_Box box(0, 0, 200, 200, "Hey, I mean, Hello, World!");
    window.show();
    return Fl::run();
}