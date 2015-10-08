#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVBoxLayout>
#include <QWebSettings>
#include <QGLWidget>
#include <QStatusBar>
#include <QGraphicsView>
#include <QtWebKitWidgets/QGraphicsWebView>
#include <QtWebKitWidgets/QWebView>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QProcess>
#include <QMessageBox>

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    QWebView* view; // QtWebKit component
    QNetworkAccessManager *nam; // Used to see if OOMMPPAA URL is ready
    QProcess *process; // Stores the OOMMPPAA server process
    QProcess *killAll;
    QUrl url; // Stores the OOMMPPAA server URL

    int slept; // Stores the amount of time we have been waiting for the OOMMPPAA server to be ready
    const static int sleepLimit = 45; // Amount of time we are willing to wait for the OOMMPPAA server

    void killall();

protected:
    void closeEvent(QCloseEvent *event); // Called when the window is closed

private slots:
    void finished(QNetworkReply *reply); // Called when a URL request sent to nam has received a response
    void loadUrl(); // Called directly and by a QTimer as a slot to load the OOMMPPAA URL (results in finished being called)
    void showProcessStderr();
    void showProcessStdout();
};

#endif // MAINWINDOW_H
