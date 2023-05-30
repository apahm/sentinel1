
#include <QtWidgets/QApplication>
#include "RLI.h"
#include <fstream>



int main(int argc, char* argv[]) {

	
	QApplication a(argc, argv);
    RLI w;
    w.show();


    
    return a.exec();
}

