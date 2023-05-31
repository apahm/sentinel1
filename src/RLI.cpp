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
#include "ipp.h"
#include <QVector>

RLI::RLI(QWidget *parent): QMainWindow(parent) {
    ui.setupUi(this);

    Sentinel sentinel;

    // add two new graphs and set their look:
    ui.plot->addGraph();
    ui.plot->graph(0)->setPen(QPen(Qt::blue)); // line color blue for first graph
    //ui.plot->graph(0)->setBrush(QBrush(QColor(0, 0, 255, 20))); // first graph will be filled with translucent blue
    ui.plot->addGraph();
    ui.plot->graph(1)->setPen(QPen(Qt::red)); // line color red for second graph
    // configure right and top axis to show ticks but no labels:
    // (see QCPAxisRect::setupFullAxesBox for a quicker method to do this)
    ui.plot->xAxis2->setVisible(true);
    ui.plot->xAxis2->setTickLabels(false);
    ui.plot->yAxis2->setVisible(true);
    ui.plot->yAxis2->setTickLabels(false);
    // make left and bottom axes always transfer their ranges to right and top axes:
    connect(ui.plot->xAxis, SIGNAL(rangeChanged(QCPRange)), ui.plot->xAxis2, SLOT(setRange(QCPRange)));
    connect(ui.plot->yAxis, SIGNAL(rangeChanged(QCPRange)), ui.plot->yAxis2, SLOT(setRange(QCPRange)));
    // pass data points to graphs:
    ui.plot->graph(0)->setData(QVector<double>::fromStdVector(sentinel.sentinel1PacketDecode.time), QVector<double>::fromStdVector(sentinel.sentinel1PacketDecode.normOfPosition));
    ui.plot->graph(1)->setData(QVector<double>::fromStdVector(sentinel.timeS), QVector<double>::fromStdVector(sentinel.position));
    // let the ranges scale themselves so graph 0 fits perfectly in the visible area:
    ui.plot->graph(0)->rescaleAxes();
    // same thing for graph 1, but only enlarge ranges (in case graph 1 is smaller than graph 0):
    ui.plot->graph(1)->rescaleAxes(true);
    ui.plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    
    uint32_t azimuth = sentinel.sentinel1PacketDecode.out.size();
    uint32_t range = sentinel.sentinel1PacketDecode.out[0].size();

    _scene = new QGraphicsScene();
    ui.graphicsView->setScene(_scene);

    QImage image(azimuth, range, QImage::Format_Grayscale8);
    
    float max = 0.0;
    
    for (int i = 0; i < azimuth; ++i) {
        for (int j = 0; j < range; ++j) {
            if (std::sqrt(std::pow(sentinel.sentinel1PacketDecode.out[i][j].re, 2) + std::pow(sentinel.sentinel1PacketDecode.out[i][j].im, 2)) > max)
                max = std::sqrt(std::pow(sentinel.sentinel1PacketDecode.out[i][j].re, 2) + std::pow(sentinel.sentinel1PacketDecode.out[i][j].im, 2));
        }
    }
    
    for (int i = 0; i < azimuth; ++i) {
        for (int j = 0; j < range; ++j) {
            image.setPixel(i, j, qRgb(std::sqrt(std::pow(sentinel.sentinel1PacketDecode.out[i][j].re, 2) + std::pow(sentinel.sentinel1PacketDecode.out[i][j].im, 2)) * 255/max,
                std::sqrt(std::pow(sentinel.sentinel1PacketDecode.out[i][j].re, 2) + std::pow(sentinel.sentinel1PacketDecode.out[i][j].im, 2)) * 255 / max,
                std::sqrt(std::pow(sentinel.sentinel1PacketDecode.out[i][j].re, 2) + std::pow(sentinel.sentinel1PacketDecode.out[i][j].im, 2)) * 255 / max));
        }
    }

    _item = new QGraphicsPixmapItem (QPixmap::fromImage(image));
    _scene->addItem(_item);
    
    z = new Graphics_view_zoom(ui.graphicsView);
    z->set_modifiers(Qt::NoModifier);
    ui.graphicsView->show();
    ui.graphicsView->rotate(90);
    
}

