/*
Copyright (C) 2010 David Doria, daviddoria@gmail.com

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

#include "MainWindow.h"
#include "Conversions.h"
#include "RegionSelectionWidget.h"

// VTK
#include <vtkActor.h>
#include <vtkBorderRepresentation.h>
#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkImageSlice.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkWidgetRepresentation.h>
#include <vtkUnsignedCharArray.h>

// STL
#include <iostream>

// Qt
#include <QButtonGroup>
#include <QFileDialog>

// ITK
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"

MainWindow::MainWindow(QWidget *parent)
{
  // Setup the GUI and connect all of the signals and slots
  setupUi(this);

  this->InteractorStyleImage = vtkSmartPointer<vtkInteractorStyleImage>::New();
  this->InteractorStyleTrackballCamera = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  //this->InteractorStyle->TrackballStyle->AddObserver(CustomTrackballStyle::PatchesMovedEvent, this, &MainWindow::UserPatchMoved);
  
  this->Image = NULL;
  
  QButtonGroup* groupFrom = new QButtonGroup;
  groupFrom->addButton(radRGB);
  groupFrom->addButton(radHSV);
  groupFrom->addButton(radCIELab);
  radRGB->setChecked(true);

  this->ImageRenderer = vtkSmartPointer<vtkRenderer>::New();
  this->ImageRenderer->AddViewProp(this->ImageLayer.ImageSlice);
  
  this->ColorSpaceRenderer = vtkSmartPointer<vtkRenderer>::New();
  
  this->ColorSpaceRenderer->AddViewProp(this->SelectedPoints.Actor);
  
  this->InteractorStyleImage->SetCurrentRenderer(this->ImageRenderer);
  this->qvtkWidgetImage->GetRenderWindow()->GetInteractor()->SetInteractorStyle(this->InteractorStyleImage);
  
  this->InteractorStyleTrackballCamera->SetCurrentRenderer(this->ColorSpaceRenderer);
  this->qvtkWidgetColorSpace->GetRenderWindow()->GetInteractor()->SetInteractorStyle(this->InteractorStyleTrackballCamera);
  
  this->Colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  
  this->Colors->SetNumberOfComponents(3);
  this->Colors->SetName ("Colors");

  this->RGBPoints.PolyData->GetPointData()->SetScalars(this->Colors);
  this->HSVPoints.PolyData->GetPointData()->SetScalars(this->Colors);
  this->CIELabPoints.PolyData->GetPointData()->SetScalars(this->Colors);

  this->qvtkWidgetImage->GetRenderWindow()->AddRenderer(this->ImageRenderer);
  this->qvtkWidgetColorSpace->GetRenderWindow()->AddRenderer(this->ColorSpaceRenderer);

  CreateColors();
  SetupRGBCube();
  SetupHSVCylinder();
  SetupCIELab();

  SetupFromGUI();
  
  this->RegionSelector = vtkSmartPointer<RegionSelectionWidget>::New();
  this->RegionSelector->SetInteractor(this->qvtkWidgetImage->GetRenderWindow()->GetInteractor());
  this->RegionSelector->CreateDefaultRepresentation();
  this->RegionSelector->SelectableOff();
  this->RegionSelector->On();
  
  this->RegionSelector->Renderer = this->ImageRenderer;

  this->RegionSelector->AddObserver(this->RegionSelector->ChangedEvent, this, &MainWindow::EndSelectionSlot);

}

void MainWindow::EndSelectionSlot()
{
  double* lowerLeft = static_cast<vtkBorderRepresentation*>(this->RegionSelector->GetRepresentation())->GetPositionCoordinate()->GetComputedWorldValue (this->ImageRenderer);
  //std::cout << "Lower left: " << lowerLeft[0] << " " << lowerLeft[1] << std::endl;

  double* upperRight = static_cast<vtkBorderRepresentation*>(this->RegionSelector->GetRepresentation())->GetPosition2Coordinate()->GetComputedWorldValue (this->ImageRenderer);
  //std::cout << "Upper right: " << upperRight[0] << " " << upperRight[1] << std::endl;
  
  itk::Index<2> corner;
  corner[0] = std::max(round(lowerLeft[0]), 0.0d);
  corner[1] = std::max(round(lowerLeft[1]), 0.0d);
  
  itk::Size<2> size;
  size[0] = round(upperRight[0] - corner[0]);
  size[1] = round(upperRight[1] - corner[1]);
  
  SelectedRegion.SetIndex(corner);
  SelectedRegion.SetSize(size);
  
  CreatePointsFromRegion(SelectedRegion);
}

void MainWindow::CreatePointsFromRegion(const itk::ImageRegion<2>& region)
{
  if(!this->Image)
    {
    return;
    }
  itk::ImageRegionConstIterator<ImageType> imageIterator(this->Image, region);
 
  this->SelectedPoints.Points->Reset();
  this->SelectedPoints.Points->Squeeze();
 
  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName ("Colors");
  
  while(!imageIterator.IsAtEnd())
    {
    ImageType::PixelType pixel = imageIterator.Get();
    unsigned char color[3] = {static_cast<unsigned char>(pixel[0]), static_cast<unsigned char>(pixel[1]), static_cast<unsigned char>(pixel[2])};
    float color_float[3] = {static_cast<float>(pixel[0])/255.0f, static_cast<float>(pixel[1])/255.0f, static_cast<float>(pixel[2])/255.0f};
  
    float coordinate[3];
    if(this->radRGB->isChecked())
      {
      coordinate[0] = pixel[0];
      coordinate[1] = pixel[1];
      coordinate[2] = pixel[2];
      }
    else if(this->radHSV->isChecked())
      {
      RGBToShapedHSV(color_float, coordinate);
      }
    else if(this->radCIELab->isChecked())
      {
      RGBtoCIELab(color, coordinate);
      }
    this->SelectedPoints.Points->InsertNextPoint(coordinate[0], coordinate[1], coordinate[2]);
    
    colors->InsertNextTupleValue(color);
    ++imageIterator;
    }
  this->SelectedPoints.PolyData->GetPointData()->SetScalars(colors);
  this->SelectedPoints.Points->Modified();
  this->ColorSpaceRenderer->ResetCamera();
  this->qvtkWidgetColorSpace->GetRenderWindow()->Render();
}

void MainWindow::CreateColors()
{
  unsigned int spacing = 10;
  for(unsigned int r = 0; r < 256; r += spacing)
    {
    for(unsigned int g = 0; g < 256; g += spacing)
      {
      for(unsigned int b = 0; b < 256; b += spacing)
        {
        unsigned char color[3] = {r,g,b};
        this->Colors->InsertNextTupleValue(color);
        }
      }
    }
}

void MainWindow::SetupRGBCube()
{
  std::cout << "SetupRGBCube" << std::endl;
  this->RGBPoints.Points->Reset();
  this->RGBPoints.Points->Squeeze();
  for(unsigned int i = 0 ;i < this->Colors->GetNumberOfTuples(); ++i)
    {
    unsigned char color[3];
    this->Colors->GetTupleValue(i, color);
    this->RGBPoints.Points->InsertNextPoint(color[0], color[1], color[2]);
    }
}

void MainWindow::SetupCIELab()
{
  std::cout << "SetupCIELab" << std::endl;
  this->CIELabPoints.Points->Reset();
  this->CIELabPoints.Points->Squeeze();
  for(unsigned int i = 0 ;i < this->Colors->GetNumberOfTuples(); ++i)
    {
    unsigned char color[3];
    this->Colors->GetTupleValue(i, color);
    
    float cielab[3];
    RGBtoCIELab(color, cielab);

    //double color_double[3] = {static_cast<float>(color[0])/255.0f, static_cast<float>(color[1])/255.0f, static_cast<float>(color[2])/255.0f};
    //double cielab[3];
    //vtkMath::RGBToLab(color_double, cielab);

    float L = cielab[0];
    float a = cielab[1];
    float b = cielab[2];

    float xyz[3] = {L,a,b};
    this->CIELabPoints.Points->InsertNextPoint(xyz);
    }

  bool fitInRGBCube = false;
  if(fitInRGBCube)
    {
    // Translate and scale the points to they are in the same position and magnitude as the RGB cube
    double cielabBounds[6];
    this->CIELabPoints.Points->GetBounds(cielabBounds);

    double rgbBounds[6];
    this->RGBPoints.Points->GetBounds(rgbBounds);

    OutputBounds("rgbBounds", rgbBounds);

    float scale[3];
    for(unsigned int i = 0 ; i < 3; ++i)
      {
      scale[i] = (rgbBounds[2*i + 1] - rgbBounds[2*i])/(cielabBounds[2*i + 1] - cielabBounds[2*i]);
      }
    OutputBounds("cielabBounds", cielabBounds);

    std::cout << "Scale: " << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
    
    vtkSmartPointer<vtkTransform> scaleTransform = vtkSmartPointer<vtkTransform>::New();
    scaleTransform->Scale(scale);
    
    vtkSmartPointer<vtkTransformPolyDataFilter> scaleTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    scaleTransformFilter->SetInputConnection(this->CIELabPoints.PolyData->GetProducerPort());
    scaleTransformFilter->SetTransform(scaleTransform);
    scaleTransformFilter->Update();

    double scaledBounds[6];
    scaleTransformFilter->GetOutput()->GetBounds(scaledBounds);

    OutputBounds("scaledBounds", scaledBounds);
    
    float translation[3];
    for(unsigned int i = 0 ; i < 3; ++i)
      {
      translation[i] = rgbBounds[2*i] - scaledBounds[2*i];
      }
      
    vtkSmartPointer<vtkTransform> translateTransform = vtkSmartPointer<vtkTransform>::New();
    translateTransform->Translate(translation);

    vtkSmartPointer<vtkTransformPolyDataFilter> translateTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    translateTransformFilter->SetInputConnection(scaleTransformFilter->GetOutputPort());
    translateTransformFilter->SetTransform(translateTransform);
    translateTransformFilter->Update();
    
    std::cout << "translation: " << translation[0] << " " << translation[1] << " " << translation[2] << std::endl;

    this->CIELabPoints.Points->DeepCopy(translateTransformFilter->GetOutput()->GetPoints());

    double newBounds[6];
    this->CIELabPoints.Points->GetBounds(newBounds);
    OutputBounds("newBounds", newBounds);

    this->CIELabPoints.Points->Modified();
    this->CIELabPoints.PolyData->Modified();
    }
}

void MainWindow::SetupHSVCylinder()
{
  std::cout << "SetupHSVCylinder" << std::endl;
  this->HSVPoints.Points->Reset();
  this->HSVPoints.Points->Squeeze();
  for(unsigned int i = 0 ;i < this->Colors->GetNumberOfTuples(); ++i)
    {
    unsigned char color[3];
    this->Colors->GetTupleValue(i, color);
    float floatRGB[3] = {static_cast<float>(color[0])/255.0f, static_cast<float>(color[1])/255.0f, static_cast<float>(color[2])/255.0f};

    float hsv[3];
    RGBToShapedHSV(floatRGB, hsv);

    this->HSVPoints.Points->InsertNextPoint(hsv);
    }

  bool fitInRGBCube = false;
  if(fitInRGBCube)
    {
    // Translate and scale the points to they are in the same position and magnitude as the RGB cube
    double hsvBounds[6];
    this->HSVPoints.Points->GetBounds(hsvBounds);

    double rgbBounds[6];
    this->RGBPoints.Points->GetBounds(rgbBounds);

    float translation[3];
    float scale[3];
    for(unsigned int i = 0 ; i < 3; ++i)
      {
      //scale[i] = (rgbBounds[2*i] - rgbBounds[2*i + 1])/(hsvBounds[2*i] - hsvBounds[2*i + 1]);
      scale[i] = (rgbBounds[2*i + 1] - rgbBounds[2*i])/(hsvBounds[2*i + 1] - hsvBounds[2*i]);
      translation[i] = rgbBounds[2*i] - hsvBounds[2*i];
      }

    std::cout << "Scale: " << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
    std::cout << "translation: " << translation[0] << " " << translation[1] << " " << translation[2] << std::endl;
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Scale(scale);
    transform->Translate(translation);

    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilter->SetInputConnection(this->HSVPoints.PolyData->GetProducerPort());
    transformFilter->SetTransform(transform);
    transformFilter->Update();

    this->HSVPoints.Points->DeepCopy(transformFilter->GetOutput()->GetPoints());
    
    this->HSVPoints.Points->Modified();
    this->HSVPoints.PolyData->Modified();
    }
}

void MainWindow::on_chkShowFullSpace_clicked()
{
  this->RGBPoints.Actor->SetVisibility(chkShowFullSpace->isChecked());
  this->HSVPoints.Actor->SetVisibility(chkShowFullSpace->isChecked());
  this->CIELabPoints.Actor->SetVisibility(chkShowFullSpace->isChecked());
  SetupFromGUI();
  this->ColorSpaceRenderer->ResetCamera();
  this->qvtkWidgetColorSpace->GetRenderWindow()->Render();
}

void MainWindow::SetupFromGUI()
{
  if(chkShowFullSpace->isChecked())
    {
    if(radRGB->isChecked())
      {
      on_radRGB_clicked();
      }
    else if(radHSV->isChecked())
      {
      on_radHSV_clicked();
      }
    else if(radCIELab->isChecked())
      {
      on_radCIELab_clicked();
      }
    }
}

void MainWindow::on_radRGB_clicked()
{
  CreatePointsFromRegion(SelectedRegion);
  if(chkShowFullSpace->isChecked())
    {
    this->ColorSpaceRenderer->RemoveAllViewProps();
    this->ColorSpaceRenderer->AddViewProp(this->SelectedPoints.Actor);
    this->ColorSpaceRenderer->AddViewProp(this->RGBPoints.Actor);
    this->ColorSpaceRenderer->ResetCamera();
    this->qvtkWidgetColorSpace->GetRenderWindow()->Render();
    }
}

void MainWindow::on_radHSV_clicked()
{
  CreatePointsFromRegion(SelectedRegion);
  if(chkShowFullSpace->isChecked())
    {
    this->ColorSpaceRenderer->RemoveAllViewProps();
    this->ColorSpaceRenderer->AddViewProp(this->SelectedPoints.Actor);
    this->ColorSpaceRenderer->AddViewProp(this->HSVPoints.Actor);
    this->ColorSpaceRenderer->ResetCamera();
    this->qvtkWidgetColorSpace->GetRenderWindow()->Render();
    }
}

void MainWindow::on_radCIELab_clicked()
{
  CreatePointsFromRegion(SelectedRegion);
  if(chkShowFullSpace->isChecked())
    {
    this->ColorSpaceRenderer->RemoveAllViewProps();
    this->ColorSpaceRenderer->AddViewProp(this->SelectedPoints.Actor);
    this->ColorSpaceRenderer->AddViewProp(this->CIELabPoints.Actor);
    this->ColorSpaceRenderer->ResetCamera();
    this->qvtkWidgetColorSpace->GetRenderWindow()->Render();
    }
}

void OutputBounds(const std::string& name, double bounds[6])
{
  std::cout << name << " xmin: " << bounds[0] << " "
            << name << " xmax: " << bounds[1] << std::endl
            << name << " ymin: " << bounds[2] << " "
            << name << " ymax: " << bounds[3] << std::endl
            << name << " zmin: " << bounds[4] << " "
            << name << " zmax: " << bounds[5] << std::endl;
}

void MainWindow::on_actionOpen_activated()
{
  QString fileName = QFileDialog::getOpenFileName(this,
                    "OpenFile", ".", "All Files (*.*)");

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fileName.toStdString());
  reader->Update();
  
  this->Image = reader->GetOutput();
  
  ITKImageToVTKRGBImage(reader->GetOutput(), this->ImageLayer.ImageData);
  this->ImageLayer.ImageData->Modified();
  this->ImageRenderer->ResetCamera();
  
  
  vtkBorderRepresentation* representation = static_cast<vtkBorderRepresentation*>(this->RegionSelector->GetRepresentation());
//   representation->GetPositionCoordinate()->SetCoordinateSystemToWorld();
//   representation->GetPosition2Coordinate()->SetCoordinateSystemToWorld();
//   representation->SetPosition(0, 0);
//   representation->SetPosition2(20, 20);
}


// Convert a vector ITK image to a VTK image for display
void ITKImageToVTKRGBImage(const MainWindow::ImageType::Pointer image, vtkImageData* outputImage)
{
  // This function assumes an ND (with N>3) image has the first 3 channels as RGB and extra information in the remaining channels.
  
  //std::cout << "ITKImagetoVTKRGBImage()" << std::endl;
  if(image->GetNumberOfComponentsPerPixel() < 3)
    {
    std::cerr << "The input image has " << image->GetNumberOfComponentsPerPixel() << " components, but at least 3 are required." << std::endl;
    return;
    }

  // Setup and allocate the image data
  outputImage->SetNumberOfScalarComponents(3);
  outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  outputImage->AllocateScalars();

  // Copy all of the input image pixels to the output image
  itk::ImageRegionConstIteratorWithIndex<MainWindow::ImageType> imageIterator(image,image->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    for(unsigned int component = 0; component < 3; component++)
      {
      pixel[component] = static_cast<unsigned char>(imageIterator.Get()[component]);
      }

    ++imageIterator;
    }
    
  outputImage->Modified();
}
