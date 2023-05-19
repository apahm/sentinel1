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

RLI::RLI(QWidget *parent): QMainWindow(parent) {
    ui.setupUi(this);

    Sentinel sentinel;

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
}

