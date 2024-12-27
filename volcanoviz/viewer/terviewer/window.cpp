/*******************************************************************************
 *
 * EcoSynth - Data-driven Authoring of Large-Scale Ecosystems (Undergrowth simulator)
 * Copyright (C) 2020  J.E. Gain  (jgain@cs.uct.ac.za)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 ********************************************************************************/


/****************************************************************************
**
** Copyright (C) 2012 Digia Plc and/or its subsidiary(-ies).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Digia Plc and its Subsidiary(-ies) nor the names
**     of its contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "glwidget.h"
#include "window.h"
#include "vecpnt.h"
#include "common/str.h"

#include <cmath>
#include <string>

using namespace std;

QSize Window::sizeHint() const
{
    return QSize(1100, 1000);
}

void Window::addSlider(QVBoxLayout *layout, const QString &label,
                          float startValue, int scale, float low, float high)
{
    const float defaultValue = startValue;
    const int defaultScaled = int(std::round(defaultValue * (float) scale));

    innerlabel = new QLabel(label);
    layout->addWidget(innerlabel);

    innerslider = new QSlider(Qt::Horizontal);
    innerslider->setMinimum(int(std::ceil(low * (float) scale)));
    innerslider->setMaximum(int(std::floor(high * (float) scale)));
    innerslider->setPageStep(200);
    innerslider->setSingleStep(1);
    innerslider->setTracking(true);
    innerslider->setValue(defaultScaled);

    const float invScale = 1.0f / scale;
    sliderscale = invScale;

    connect(innerslider, &QSlider::valueChanged, [this, invScale] (int newValue)
            {
                // perspectiveView->setRadius(newValue * invScale);
                repaintAllGL();
            });
    layout->addWidget(innerslider);
}

float Window::getSliderValue()
{
    return innerslider->value() * sliderscale;
}


Window::Window(string datadir)
{
    QWidget *mainWidget = new QWidget;
    QGridLayout *mainLayout = new QGridLayout;
    int dx, dy;
    float sx, sy;

    // default rendering parameters, set using text entry
    // mirrors TRenderer settings
    // grid params
    gridIntensity = 0.8f; // 80% of base colour
    gridSepX = 2500.0f; // separation of grid lines, depends on how input data is scaled
    gridSepZ = 2500.0f; //
    gridWidth = 1.5f; // in pixels?
    
    // contour params
    contourSep = 25.f; // separation (Y direction) depends on how input data is normalized
    numContours = 1.0f / contourSep;
    contourWidth = 1.0f; // in pixels ?
    contourIntensity = 1.2f; // 130% of base colour

    // radiance scaling parameters
    radianceTransition = 0.2f;
    radianceEnhance = 2.2f;

    
    setCentralWidget(mainWidget);
    mainLayout->setColumnStretch(0, 0);
    mainLayout->setColumnStretch(1, 0);
    mainLayout->setColumnStretch(2, 1);
    
    // render panel
    renderPanel = new QWidget;
    QVBoxLayout *renderLayout = new QVBoxLayout;

    // Grid Line Widgets
    QGroupBox *gridGroup = new QGroupBox(tr("Grid Lines"));
    QCheckBox * checkGridLines = new QCheckBox(tr("Show Grid Lines"));
    checkGridLines->setChecked(false);
    QLabel *gridSepXLabel = new QLabel(tr("Grid Sep X:"));
    gridSepXEdit = new QLineEdit;
    gridSepXEdit->setFixedWidth(60);
    // gridSepXEdit->setValidator(new QDoubleValidator(0.0, 500000.0, 2, gridSepXEdit));
    gridSepXEdit->setInputMask("0000.0");
    QLabel *gridSepZLabel = new QLabel(tr("Grid Sep Z:"));
    gridSepZEdit = new QLineEdit;
    gridSepZEdit->setFixedWidth(60);
    // gridSepZEdit->setValidator(new QDoubleValidator(0.0, 500000.0, 2, gridSepZEdit));
    gridSepZEdit->setInputMask("0000.0");
    QLabel *gridWidthLabel = new QLabel(tr("Grid Line Width:"));
    gridWidthEdit = new QLineEdit;
    gridWidthEdit->setFixedWidth(60);
    // gridWidthEdit->setValidator(new QDoubleValidator(0.0, 10.0, 2, gridWidthEdit));
    gridWidthEdit->setInputMask("0.0");
    QLabel *gridIntensityLabel = new QLabel(tr("Grid Intensity:"));
    gridIntensityEdit = new QLineEdit;
    gridIntensityEdit->setFixedWidth(60);
    // gridIntensityEdit->setValidator(new QDoubleValidator(0.0, 2.0, 2, gridIntensityEdit));
    gridIntensityEdit->setInputMask("0.0");
    
    // set initial grid values
    gridSepXEdit->setText(QString::number(gridSepX, 'g', 2));
    gridSepZEdit->setText(QString::number(gridSepZ, 'g', 2));
    gridWidthEdit->setText(QString::number(gridWidth, 'g', 2));
    gridIntensityEdit->setText(QString::number(gridIntensity, 'g', 2));

    QGridLayout *gridLayout = new QGridLayout;
    gridLayout->addWidget(checkGridLines, 0, 0);
    gridLayout->addWidget(gridSepXLabel, 1, 0);
    gridLayout->addWidget(gridSepXEdit, 1, 1);
    gridLayout->addWidget(gridSepZLabel, 2, 0);
    gridLayout->addWidget(gridSepZEdit, 2, 1);
    gridLayout->addWidget(gridWidthLabel, 3, 0);
    gridLayout->addWidget(gridWidthEdit, 3, 1);
    gridLayout->addWidget(gridIntensityLabel, 4, 0);
    gridLayout->addWidget(gridIntensityEdit, 4, 1);
    gridGroup->setLayout(gridLayout);

    // Contour Widgets
    QGroupBox *contourGroup = new QGroupBox(tr("Contours"));
    QCheckBox * checkContours = new QCheckBox(tr("Show Contours"));
    checkContours->setChecked(false);
    QLabel *contourSepLabel = new QLabel(tr("Contour Sep:"));
    contourSepEdit = new QLineEdit;
    contourSepEdit->setFixedWidth(60);
    //contourSepEdit->setValidator(new QDoubleValidator(0.0, 10000.0, 2, contourSepEdit));
    contourSepEdit->setInputMask("000.0");
    QLabel *contourWidthLabel = new QLabel(tr("Contour Line Width:"));
    contourWidthEdit = new QLineEdit;
    // contourWidthEdit->setValidator(new QDoubleValidator(0.0, 10.0, 2, contourWidthEdit));
    contourWidthEdit->setInputMask("0.0");
    contourWidthEdit->setFixedWidth(60);
    QLabel *contourIntensityLabel = new QLabel(tr("Contour Intensity:"));
    contourIntensityEdit = new QLineEdit;
    contourIntensityEdit->setFixedWidth(60);
    contourIntensityEdit->setInputMask("0.0");

    // set initial contour values
    contourSepEdit->setText(QString::number(contourSep, 'g', 2));
    contourWidthEdit->setText(QString::number(contourWidth, 'g', 2));
    contourIntensityEdit->setText(QString::number(contourIntensity, 'g', 2));

    QGridLayout *contourLayout = new QGridLayout;
    contourLayout->addWidget(checkContours, 0, 0);
    contourLayout->addWidget(contourSepLabel, 1, 0);
    contourLayout->addWidget(contourSepEdit, 1, 1);
    contourLayout->addWidget(contourWidthLabel, 2, 0);
    contourLayout->addWidget(contourWidthEdit, 2, 1);
    contourLayout->addWidget(contourIntensityLabel, 3, 0);
    contourLayout->addWidget(contourIntensityEdit, 3, 1);
    contourGroup->setLayout(contourLayout);

    // Radiance
    QGroupBox *radianceGroup = new QGroupBox(tr("Radiance"));
    QLabel *radianceEnhanceLabel = new QLabel(tr("Radiance Enhancement:"));
    radianceEnhanceEdit = new QLineEdit;
    radianceEnhanceEdit->setFixedWidth(60);
    radianceEnhanceEdit->setInputMask("0.0");

    // set initial radiance values
    radianceEnhanceEdit->setText(QString::number(radianceEnhance, 'g', 2));

    QGridLayout *radianceLayout = new QGridLayout;
    radianceLayout->addWidget(radianceEnhanceLabel, 0, 0);
    radianceLayout->addWidget(radianceEnhanceEdit, 0, 1);
    radianceGroup->setLayout(radianceLayout);

    renderLayout->addWidget(gridGroup);
    renderLayout->addWidget(contourGroup);
    renderLayout->addWidget(radianceGroup);
    
   /*
    // right-hand panel for pallete
    palPanel = new QWidget;
    QVBoxLayout *palLayout = new QVBoxLayout;
   */

    // OpenGL widget
    // Specify an OpenGL 3.2 format.

    QGLFormat glFormat;
    glFormat.setVersion( 4, 1 );
    glFormat.setProfile( QGLFormat::CoreProfile );
    glFormat.setSampleBuffers( false );

    perspectiveView = new GLWidget(glFormat, datadir, this);
    
    getView().setForcedFocus(getTerrain().getFocus());
    getView().setViewScale(getTerrain().longEdgeDist());

    getTerrain().getGridDim(dx, dy);
    getTerrain().getTerrainDim(sx, sy);

    numGridX = 1.0f / gridSepX;
    numGridZ = 1.0f / gridSepZ;

    // Palette Widget
    /*
    palLayout->addWidget(perspectiveView->getPalette());
    addSlider(palLayout, tr("Active Brush Size"), perspectiveView->getRadius(), 200, 100.0f, 2500.0f);
    */

    // signal to slot connections
    connect(perspectiveView, SIGNAL(signalRepaintAllGL()), this, SLOT(repaintAllGL()));
    connect(gridSepXEdit, SIGNAL(editingFinished()), this, SLOT(lineEditChange()));
    connect(gridSepZEdit, SIGNAL(editingFinished()), this, SLOT(lineEditChange()));
    connect(gridWidthEdit, SIGNAL(editingFinished()), this, SLOT(lineEditChange()));
    connect(gridIntensityEdit, SIGNAL(editingFinished()), this, SLOT(lineEditChange()));
    connect(contourSepEdit, SIGNAL(editingFinished()), this, SLOT(lineEditChange()));
    connect(contourWidthEdit, SIGNAL(editingFinished()), this, SLOT(lineEditChange()));
    connect(contourIntensityEdit, SIGNAL(editingFinished()), this, SLOT(lineEditChange()));
    connect(radianceEnhanceEdit, SIGNAL(editingFinished()), this, SLOT(lineEditChange()));
    connect(radianceEnhanceEdit, SIGNAL(returnPressed()), this, SLOT(lineEditChange()));

    // display switches
    connect(checkContours, SIGNAL(stateChanged(int)), this, SLOT(showContours(int)));
    connect(checkGridLines, SIGNAL(stateChanged(int)), this, SLOT(showGridLines(int)));

    renderPanel->setLayout(renderLayout);
    // palPanel->setLayout(palLayout);

    mainLayout->addWidget(renderPanel, 0, 0, Qt::AlignTop);
    mainLayout->addWidget(perspectiveView, 0, 2);
    // mainLayout->addWidget(palPanel, 0, 3, Qt::AlignTop);

    createActions();
    createMenus();
    
    mainWidget->setLayout(mainLayout);
    setWindowTitle(tr("TerrainViewer"));
    mainWidget->setMouseTracking(true);
    setMouseTracking(true);
    
    renderPanel->hide();
    // palPanel->hide();

    perspectiveView->getRenderer()->setGridParams(numGridX, numGridZ, gridWidth, gridIntensity);
    perspectiveView->getRenderer()->setContourParams(numContours, contourWidth, contourIntensity);
    perspectiveView->getRenderer()->setRadianceScalingParams(radianceEnhance);
    
}

void Window::run_viewer(int numframes, int start)
{
    perspectiveView->loadScene(numframes, start);
    repaintAllGL();
}

void Window::scaleRenderParams(float scale)
{
    gridSepX = scale / 5.0f; // separation of grid lines, depends on how input data is scaled
    gridSepZ = scale / 5.0f;
    numGridX = 1.0f / gridSepX;
    numGridZ = 1.0f / gridSepZ;
    gridSepXEdit->setText(QString::number(gridSepX, 'g', 2));
    gridSepZEdit->setText(QString::number(gridSepZ, 'g', 2));

    contourSep = scale / 100.f; // separation (Y direction) depends on how input data is normalized
    numContours = 1.0f / contourSep;
    contourSepEdit->setText(QString::number(contourSep, 'g', 2));

    perspectiveView->getRenderer()->setGridParams(numGridX, numGridZ, gridWidth, gridIntensity);
    perspectiveView->getRenderer()->setContourParams(numContours, contourWidth, contourIntensity);
    perspectiveView->getRenderer()->setRadianceScalingParams(radianceEnhance);
}

void Window::keyPressEvent(QKeyEvent *e)
{
    // pass to render window
    perspectiveView->keyPressEvent(e);
}

void Window::mouseMoveEvent(QMouseEvent *event)
{
    QWidget *child=childAt(event->pos());
    QGLWidget *glwidget = qobject_cast<QGLWidget *>(child);
    if(glwidget) {
        QMouseEvent *glevent=new QMouseEvent(event->type(),glwidget->mapFromGlobal(event->globalPos()),event->button(),event->buttons(),event->modifiers());
        QCoreApplication::postEvent(glwidget,glevent);
    }
}

void Window::repaintAllGL()
{
    perspectiveView->repaint();
}

void Window::showRenderOptions()
{
    renderPanel->setVisible(showRenderAct->isChecked());
}

void Window::showPalette(bool show)
{
    palPanel->setVisible(show);
}

void Window::showContours(int show)
{
    perspectiveView->getRenderer()->drawContours(show == Qt::Checked);
    repaintAllGL();
}

void Window::showGridLines(int show)
{
    perspectiveView->getRenderer()->drawGridlines(show == Qt::Checked);
    repaintAllGL();
}

void Window::lineEditChange()
{
    bool ok;
    float val;
    float tx, ty, hr;

    tx = 1.0f; ty = 1.0f; // to fix when scale added
    hr = 1.0f;

    if(sender() == gridSepXEdit)
    {
        val = gridSepXEdit->text().toFloat(&ok);
        if(ok)
        {
            gridSepX = val;
            numGridX = tx / gridSepX; // convert separation to num grid lines
        }
    }
    if(sender() == gridSepZEdit)
    {
        val = gridSepZEdit->text().toFloat(&ok);
        if(ok)
        {
            gridSepZ = val;
            numGridZ = ty / gridSepZ;
        }
    }
    if(sender() == gridWidthEdit)
    {
        val = gridWidthEdit->text().toFloat(&ok);
        if(ok)
        {
            gridWidth = val;
        }
    }
    if(sender() == gridIntensityEdit)
    {
        val = gridIntensityEdit->text().toFloat(&ok);
        if(ok)
        {
            gridIntensity = val;
        }
    }
    if(sender() == contourSepEdit)
    {
        val = contourSepEdit->text().toFloat(&ok);
        if(ok)
        {
            contourSep = val;
            numContours = hr / contourSep;
        }
    }
    if(sender() == contourWidthEdit)
    {
        val = contourWidthEdit->text().toFloat(&ok);
        if(ok)
        {
            contourWidth = val;
        }
    }
    if(sender() == contourIntensityEdit)
    {
        val = contourIntensityEdit->text().toFloat(&ok);
        if(ok)
        {
            contourIntensity = val;
        }
    }
    if(sender() == radianceEnhanceEdit)
    {
        val = radianceEnhanceEdit->text().toFloat(&ok);
        if(ok)
        {
            radianceEnhance = val;
        }
    }

    // cerr << "val entered " << val << endl;

    // without this the renderer defaults back to factory settings at certain stages - very wierd bug
    perspectiveView->getRenderer()->setGridParams(numGridX, numGridZ, gridWidth, gridIntensity);
    perspectiveView->getRenderer()->setContourParams(numContours, contourWidth, contourIntensity);
    perspectiveView->getRenderer()->setRadianceScalingParams(radianceEnhance);
    repaintAllGL();
}

void Window::saveMaps()
{
    QFileDialog::Options options;
    QString selectedFilter;
    QString paintfilename = QFileDialog::getSaveFileName(this,
                                    tr("Save Type Maps"),
                                    "~/",
                                    tr("Text Files (*.txt)"),
                                    &selectedFilter,
                                    options);
    if (!paintfilename.isEmpty())
    {
        std::string outfile = paintfilename.toStdString();

        if(endsWith(outfile, ".txt"))
        {
            QDir dir = QFileInfo(paintfilename).absoluteDir();
            QString rootname = dir.absolutePath();

            // remove the file extension
            QFileInfo file(paintfilename);
            rootname.append("/");
            rootname.append(file.completeBaseName());
            rootname.append(".txt");
            perspectiveView->writePaintMaps(rootname.toStdString());
        }
    }
}

void Window::loadMaps()
{
    QString fileName = QFileDialog::getOpenFileName(this,
                                                     tr("Open Type Maps"),
                                                     "~/",
                                                     tr("Text File (*.txt)"));
     if (!fileName.isEmpty())
     {
         std::string infile = fileName.toUtf8().constData();

         // use file extension to determine action
         if(endsWith(infile, ".txt"))
         {
             QDir dir = QFileInfo(fileName).absoluteDir();
             QString rootname = dir.absolutePath();

             // remove the file extension
             QFileInfo file(fileName);
             rootname.append("/");
             rootname.append(file.completeBaseName());
             rootname.append(".txt");
             perspectiveView->readPaintMaps(rootname.toStdString());
             repaintAllGL();
         }
         else
         {
             cerr << "Error Window::loadMaps: attempt to open unrecognized file format" << endl;
         }
     }

}

void Window::createActions()
{
    showRenderAct = new QAction(tr("Show Terrain Options"), this);
    showRenderAct->setCheckable(true);
    showRenderAct->setChecked(false);
    showRenderAct->setStatusTip(tr("Hide/Show Rendering Options"));
    connect(showRenderAct, SIGNAL(triggered()), this, SLOT(showRenderOptions()));
    saveAct = new QAction(tr("Save Type Maps"), this);
    connect(saveAct, SIGNAL(triggered()), this, SLOT(saveMaps()));
    openAct = new QAction(tr("Open Type Maps"), this);
    connect(openAct, SIGNAL(triggered()), this, SLOT(loadMaps()));
}

void Window::createMenus()
{
    viewMenu = menuBar()->addMenu(tr("&File"));
    viewMenu->addAction(openAct);
    viewMenu->addAction(saveAct);
    viewMenu = menuBar()->addMenu(tr("&View"));
    viewMenu->addAction(showRenderAct);
}
