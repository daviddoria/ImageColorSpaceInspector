#ifndef REGIONSELECTIONWIDGET_H
#define REGIONSELECTIONWIDGET_H

#include <vtkBorderWidget.h>

class RegionSelectionWidget : public vtkBorderWidget
{
public:
  static RegionSelectionWidget *New();
  vtkTypeMacro(RegionSelectionWidget, vtkBorderWidget);

  int SubclassEndSelectAction();
  
  int ChangedEvent;
  
  RegionSelectionWidget();
  vtkRenderer* Renderer;

};

#endif
