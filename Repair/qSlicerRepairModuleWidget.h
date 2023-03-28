/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

#ifndef __qSlicerRepairModuleWidget_h
#define __qSlicerRepairModuleWidget_h

// Slicer includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerRepairModuleExport.h"

class qSlicerRepairModuleWidgetPrivate;
class vtkMRMLNode;
class QDialogButtonBox;
class qSlicerRepairTabWidget;
class QAbstractButton;
/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_REPAIR_EXPORT qSlicerRepairModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerRepairModuleWidget(QWidget *parent=0);
  virtual ~qSlicerRepairModuleWidget();
  
  void enter() override;
  void exit() override;
public slots:


protected:
  QScopedPointer<qSlicerRepairModuleWidgetPrivate> d_ptr;

  virtual void setup();


private:
  Q_DECLARE_PRIVATE(qSlicerRepairModuleWidget);
  Q_DISABLE_COPY(qSlicerRepairModuleWidget);
};

#endif
