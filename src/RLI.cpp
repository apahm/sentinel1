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

    Sentinel1PacketDecode sentinel1PacketDecode;
    sentinel1PacketDecode.ReadSARParam("C:/S1A_S3_RAW__0SDH_20220710T213600_20220710T213625_044043_0541DB_56CE/S1A_S3_RAW__0SDH_20220710T213600_20220710T213625_044043_0541DB_56CE.SAFE/s1a-s3-raw-s-hh-20220710t213600-20220710t213625-044043-0541db.dat");
    
    uint32_t azimuth = sentinel1PacketDecode.out.size();
    uint32_t range = sentinel1PacketDecode.out[0].size();

    _scene = new QGraphicsScene();
    ui.graphicsView->setScene(_scene);

    QImage image(azimuth, range, QImage::Format_Grayscale8);
    
    float max = 0.0;

    for (int i = 0; i < azimuth; ++i) {
        for (int j = 0; j < range; ++j) {
            if (std::abs(sentinel1PacketDecode.out[i][j]) > max)
                max = std::abs(sentinel1PacketDecode.out[i][j]);
        }
    }

    for (int i = 0; i < azimuth; ++i) {
        for (int j = 0; j < range; ++j) {
            image.setPixel(i, j, qRgb(std::abs(sentinel1PacketDecode.out[i][j]) * 255/max, 
                std::abs(sentinel1PacketDecode.out[i][j]) * 255 / max, 
                std::abs(sentinel1PacketDecode.out[i][j]) * 255 / max));
        }
    }

    _item = new QGraphicsPixmapItem (QPixmap::fromImage(image));
    _scene->addItem(_item);

    z = new Graphics_view_zoom(ui.graphicsView);
    z->set_modifiers(Qt::NoModifier);
    ui.graphicsView->show();
}

