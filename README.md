helvete och numeriska metoder
=============================

requires:

1. qwt >= 6.0.0 (http://qwt.sourceforge.net)
2. qt4 >= 4.3 (http://qt-project.org/ / http://www.qtcentre.org/forums/5-Installation-and-Deployment)

to build:

    $> qmake honm.pro
    $> make

if it cannot find your qwt:

    $> qmake -set QMAKEFEATURES /path/to/your/qwt.prf honm.pro
    $> make
