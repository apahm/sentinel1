#include <iostream>
#include <fstream>
#include <string>
#include <qspinbox.h>
#include <qaction.h>
#include "RLI.h"
#include "windows.h"
#include <QPixmap>
#include <QGraphicsPixmapItem>
#include "sentinel1_packet_decode.h"
#include "sentinel1.h"
#include <QImageWriter>

RLI::RLI(QWidget *parent): QMainWindow(parent) {
    ui.setupUi(this);

    Sentinel sentinel;
    QVector<double> x;
    QVector<double> y0;
    QVector<double> y1;

    x.resize(sentinel.refFunc.size());
    y0.resize(sentinel.refFunc.size());
    y1.resize(sentinel.refFunc.size());

    for (size_t i = 0; i < sentinel.refFunc.size(); i++) {
        x[i] = i;
        y0[i] = sentinel.refFunc[i].real();
        y1[i] = sentinel.refFunc[i].imag();
    }

    // add two new graphs and set their look:
    ui.plot1->addGraph();
    ui.plot1->graph(0)->setPen(QPen(Qt::blue)); // line color blue for first graph
    ui.plot1->addGraph();
    ui.plot1->graph(1)->setPen(QPen(Qt::red)); // line color red for second graph
    // generate some points of data (y0 for first, y1 for second graph):
    // configure right and top axis to show ticks but no labels:
    // (see QCPAxisRect::setupFullAxesBox for a quicker method to do this)
    ui.plot1->xAxis2->setVisible(true);
    ui.plot1->xAxis2->setTickLabels(false);
    ui.plot1->yAxis2->setVisible(true);
    ui.plot1->yAxis2->setTickLabels(false);
    // make left and bottom axes always transfer their ranges to right and top axes:
    connect(ui.plot1->xAxis, SIGNAL(rangeChanged(QCPRange)), ui.plot1->xAxis2, SLOT(setRange(QCPRange)));
    connect(ui.plot1->yAxis, SIGNAL(rangeChanged(QCPRange)), ui.plot1->yAxis2, SLOT(setRange(QCPRange)));
    // pass data points to graphs:
    ui.plot1->graph(0)->setData(x, y0);
    ui.plot1->graph(1)->setData(x, y1);
    // let the ranges scale themselves so graph 0 fits perfectly in the visible area:
    ui.plot1->graph(0)->rescaleAxes();
    // same thing for graph 1, but only enlarge ranges (in case graph 1 is smaller than graph 0):
    ui.plot1->graph(1)->rescaleAxes(true);
    // Note: we could have also just called customPlot->rescaleAxes(); instead
    // Allow user to drag axis ranges with mouse, zoom with mouse wheel and select graphs by clicking:
    ui.plot1->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    uint32_t azimuth = sentinel.sentinel1PacketDecode.out.size();
    uint32_t range = sentinel.sentinel1PacketDecode.out[0].size();

    _scene = new QGraphicsScene();
    ui.graphicsView->setScene(_scene);

    QImage image(azimuth, range, QImage::Format_Grayscale8);
    
    float max = 0.0;
    
    for (int i = 0; i < azimuth; ++i) {
        for (int j = 0; j < range; ++j) {
            if (std::abs(sentinel.sentinel1PacketDecode.out[i][j]) > max)
                max = std::abs(sentinel.sentinel1PacketDecode.out[i][j]);
        }
    }

    for (int i = 0; i < azimuth; ++i) {
        for (int j = 0; j < range; ++j) {
            image.setPixel(i, j, qRgb(std::abs(sentinel.sentinel1PacketDecode.out[i][j]) * 255/max,
                std::abs(sentinel.sentinel1PacketDecode.out[i][j]) * 255 / max, 
                std::abs(sentinel.sentinel1PacketDecode.out[i][j]) * 255 / max));
        }
    }

    _item = new QGraphicsPixmapItem (QPixmap::fromImage(image));
    _scene->addItem(_item);

    z = new Graphics_view_zoom(ui.graphicsView);
    z->set_modifiers(Qt::NoModifier);
    ui.graphicsView->show();
}

