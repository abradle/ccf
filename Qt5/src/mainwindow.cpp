#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QWebSettings>
#include <QGLWidget>
#include <QStatusBar>
#include <QProcess>
#include <QFileInfo>
#include <QNetworkRequest>
#include <QThread>
#include <QMessageBox>
#include <QDialogButtonBox>
#include <QPushButton>
#include <QDir>
#include <QTimer>
#include <QProcess>

/**
 * @brief MainWindow::MainWindow
 *
 * Constructor called when the window is created.
 *
 * @param parent
 */
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent) {

    killall();

    slept = 0; // Set the amount of time we have waited for OOMMPPAA to load

    // Create QtWebKit component and enable plugins (required for active icm) and developer extensions
    view = new QWebView(this);
    view->settings()->setAttribute(QWebSettings::PluginsEnabled,  true);
    view->settings()->setAttribute(QWebSettings::DeveloperExtrasEnabled,true);

    // Set overflow:hidden on the div containing the object to prevent the plugin
    // scrolling off screen as the window is resized
    view->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);

    // Set an initial loading screen
    view->setHtml("<html><body><h1>Loading...</h1></body></html>");

    // Install QtWebKit component in layout
    setCentralWidget(view);

    // Set working directory to location of executable
    QDir appPath = QFileInfo( QCoreApplication::applicationFilePath() ).absoluteDir();
    QDir::setCurrent(appPath.absolutePath());

    // Create new process for OOMMPPAA server
    process = new QProcess(this);

    // Implicit based on QDir::setCurrent call but we make it explicit anyway
    process->setWorkingDirectory(appPath.absolutePath());

    QString winServer = appPath.filePath("winserver.exe");

    // Start OOMMPPAA server
    process->start(winServer);

    if(process->waitForStarted(2000)){
        // We get here if the OOMMPPAA server has started (i.e. the executable actually exists)

        // This is the URL we expect OOMMPPAA to be served from (default port)
        url = QUrl("http://localhost:9020/OOMMPPAA");

        // Create QNetworkAccessManager which we will use to test if the above URL is ready
        nam = new QNetworkAccessManager(this);

        // Everytime a URL request is made using "nam" the response is sent to our method/slot finished
        connect(nam,SIGNAL(finished(QNetworkReply*)),this,SLOT(finished(QNetworkReply*)));

        // Try loading the OOMMPPAA URL (results in a call to the method/slot finished)
        loadUrl();
    }else{
        // We get here if the OOMMPPAA server has started (i.e. the executable name it wrong etc.)
        QMessageBox *dialog = new QMessageBox(this);

        QPushButton *quitBtn = new QPushButton();
        quitBtn->setText("OK");

        // Close application when the OK button is clicked
        connect(quitBtn,SIGNAL(clicked()),qApp,SLOT(quit()));

        dialog->addButton(quitBtn,QMessageBox::AcceptRole);
        dialog->setText("OMMPPAA server has failed to start\n"+winServer);
        dialog->setWindowTitle("Server failure");
        dialog->show();
    }
}

/**
 * @brief MainWindow::~MainWindow
 *
 * Called when the window is destroyed
 */
MainWindow::~MainWindow()
{

}

/**
 * @brief MainWindow::finished
 *
 * Called when a response is received from the request to load the OOMMPPAA URL
 *
 * @param reply
 */
void MainWindow::finished(QNetworkReply *reply){
    if(reply->error() == QNetworkReply::NoError){
        // No errors so load OOMMPPAA in internal QtWebKit window
        view->load(url);
    }else{
        // We will get here if the resource is not available or another error has occurred
        if(slept > sleepLimit){
            // We get here if we have waited longer than the timeout (defaults to 45 seconds)
            QMessageBox * dialog = new QMessageBox(this);

            QPushButton *quitBtn = new QPushButton();
            quitBtn->setText("OK");

            // Close when OK is clicked
            connect(quitBtn, SIGNAL(clicked()), qApp, SLOT(quit()));

            dialog->addButton(quitBtn, QMessageBox::AcceptRole);
            dialog->setText("OOMMPPAA server has failed to start in 45 seconds");
            dialog->setWindowTitle("Server timeout");
            dialog->show();
        }else if(process->state() == QProcess::NotRunning){
            // We get here if an error has occurred and the OOMMPPAA server instance we
            // started is no longer running
            QMessageBox * dialog = new QMessageBox(this);

            QPushButton *quitBtn = new QPushButton();
            quitBtn->setText("OK");

            // Close when OK is clicked
            connect(quitBtn, SIGNAL(clicked()), qApp, SLOT(quit()));

            dialog->addButton(quitBtn, QMessageBox::AcceptRole);
            dialog->setText("OOMMPPAA server has quit");
            dialog->setWindowTitle("Server error");
            dialog->show();
        }else{
            // We get here if an error has occurred but the OOMMPPAA server is still running
            // and we have not waited more than 45 seconds
            slept += 5;

            // Wait 5 seconds and try loading the URL again
            QTimer::singleShot(5000, this, SLOT(loadUrl()));
        }
    }
}

/**
 * @brief MainWindow::loadUrl
 *
 * Called to attempt to load a URL from the OOMMPPAA server. Note that MainWindow::finished
 * will be called at the end of the request via the SIGNAL finished on the object nam
 */
void MainWindow::loadUrl(){
    nam->get(QNetworkRequest(url));
}

/**
 * @brief MainWindow::closeEvent
 *
 * Called when the window is closed.  Note that if the OOMMPPAA server process object is not
 * null the OOMMPPAA server is asked to close
 *
 * @param event
 */
void MainWindow::closeEvent(QCloseEvent *event) {
    if(process != NULL){
        // Stop the process we directly started
        process->close();

        killall();
    }

    event->accept();
}

void MainWindow::killall(){
#ifdef Q_OS_WIN
    // winserver.exe spawns another winserver.exe that isn't shutdown correctly so we must killall
    killAll = new QProcess(this);
    killAll->start("taskkill", QStringList() << "/F" << "/T" << "/IM" << "winserver.exe");

    if(!killAll->waitForFinished(10000)){
        QMessageBox * dialog = new QMessageBox(this);

        QPushButton *quitBtn = new QPushButton();
        quitBtn->setText("OK");

        // Close when OK is clicked
        connect(quitBtn, SIGNAL(clicked()), qApp, SLOT(quit()));

        dialog->addButton(quitBtn, QMessageBox::AcceptRole);
        dialog->setText("Failed to killall spawned processed\nYou must manually stop winserver.exe");
        dialog->setWindowTitle("Server error");
        dialog->show();
    }
#endif
}

void MainWindow::showProcessStderr(){
    //QMessageBox::information(this,"",killAll->readAllStandardError());
}

void MainWindow::showProcessStdout(){
    //QMessageBox::information(this, "",killAll->readAllStandardOutput());
}
