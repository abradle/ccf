#include "mainwindow.h"
#include <QApplication>

/**
 * @brief main
 *
 * Starts the OOMMPPAA server as an external process and displays a QtWebKit Window
 * with the OOMMPPAA application within.  When the QtWebKit window is closed the
 * OOMMPPAA server is also closed.
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]){
    QApplication a(argc, argv);

    MainWindow w; // Create window object

    w.showMaximized(); // Make it full screen

    w.show(); // Make it visible

    return a.exec(); // Wait for application to close
}
