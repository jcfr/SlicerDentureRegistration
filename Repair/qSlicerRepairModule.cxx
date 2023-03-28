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

// Repair Logic includes
#include <vtkSlicerRepairLogic.h>

// Repair includes
#include "qSlicerRepairModule.h"
#include "qSlicerRepairModuleWidget.h"
#include "qSlicerRepairModule.h"
#include "qSlicerRepairModuleWidget.h"
#include "qSlicerApplication.h"
//#include "qSlicerRepairTabWidget.h"
#include "qmainwindow.h"
#include "QDialogButtonBox.h"
#include <QList>
#include "qstatusbar.h"
#include <QAbstractButton>
//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerRepairModulePrivate
{
public:
  qSlicerRepairModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerRepairModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerRepairModulePrivate::qSlicerRepairModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerRepairModule methods

//-----------------------------------------------------------------------------
qSlicerRepairModule::qSlicerRepairModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerRepairModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerRepairModule::~qSlicerRepairModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerRepairModule::helpText() const
{
  return "This is a loadable module that can be bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerRepairModule::acknowledgementText() const
{
  return "This work was partially funded by NIH grant NXNNXXNNNNNN-NNXN";
}

//-----------------------------------------------------------------------------
QStringList qSlicerRepairModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("John Doe (AnyWare Corp.)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerRepairModule::icon() const
{
  return QIcon(":/Icons/Repair.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerRepairModule::categories() const
{
  return QStringList() << "DentureRegistration";
}

//-----------------------------------------------------------------------------
QStringList qSlicerRepairModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerRepairModule::setup()
{
  this->Superclass::setup();



}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerRepairModule
::createWidgetRepresentation()
{
  return new qSlicerRepairModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerRepairModule::createLogic()
{
  return vtkSlicerRepairLogic::New();
}
