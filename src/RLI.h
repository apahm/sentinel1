#pragma once
#include <QtWidgets/QMainWindow>
#include "ui_RLI.h"
#include <QGraphicsScene>
#include <QImage>
#include "zoom.h"
class RLI : public QMainWindow
{
    Q_OBJECT
private:
    double _samleRate = 1'000'000'000;
    uint64_t _countOfPulse = 200;
    
    bool _digitizerInit = false;
    bool _transponderInit = false;
public:
    RLI(QWidget *parent = nullptr);
    ~RLI() {
        delete _item;
    }
private:
    QGraphicsScene* _scene = nullptr;
    QGraphicsPixmapItem* _item = nullptr;
    Graphics_view_zoom* z = nullptr;
    Ui::RLIClass ui;
};
