/*
Copyright (C) 2011 David Doria, daviddoria@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

// Custom
#include "DisplayPoints.h"
#include "Layer.h"
class RegionSelectionWidget;

// Qt
#include "ui_MainWindow.h"

// VTK
#include <vtkSmartPointer.h>
#include <vtkInteractorStyleImage.h>
#include <vtkInteractorStyleTrackballCamera.h>

// Forward declarations
class vtkPolyDataMapper;
class vtkActor;
class vtkRenderer;
class vtkVertexGlyphFilter;
class vtkPoints;
class vtkPolyData;
class vtkUnsignedCharArray;

// ITK
#include "itkImage.h"
#include "itkCovariantVector.h"

class MainWindow : public QMainWindow, private Ui::MainWindow
{
Q_OBJECT
public:
  MainWindow(QWidget *parent = 0);

  typedef itk::Image< itk::CovariantVector<unsigned char, 3>, 2 > ImageType;
  
public slots:

  void on_chkShowFullSpace_clicked();

  void on_radRGB_clicked();
  void on_radHSV_clicked();
  void on_radCIELab_clicked();
  
  void on_actionOpen_activated();

protected:
  
  void CreatePointsFromRegion(const itk::ImageRegion<2>& region);
  
  void EndSelectionSlot();
  
  void CreateColors();

  void SetupFromGUI();
  
  void SetupRGBCube();
  void SetupHSVCylinder();
  void SetupCIELab();
  
  DisplayPoints SelectedPoints;
  
  vtkSmartPointer<vtkUnsignedCharArray> Colors;
  DisplayPoints RGBPoints;
  DisplayPoints HSVPoints;
  DisplayPoints CIELabPoints;

  ImageType::Pointer Image;
  Layer ImageLayer;
  
  vtkSmartPointer<vtkRenderer> ImageRenderer;
  vtkSmartPointer<vtkRenderer> ColorSpaceRenderer;

  vtkSmartPointer<vtkInteractorStyleImage> InteractorStyleImage;
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> InteractorStyleTrackballCamera;
  
  vtkSmartPointer<RegionSelectionWidget> RegionSelector;
  
  itk::ImageRegion<2> SelectedRegion;
};

void OutputBounds(const std::string& name, double bounds[6]);
void ITKImageToVTKRGBImage(const MainWindow::ImageType::Pointer image, vtkImageData* outputImage);

#endif
