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

// Qt includes
#include <QDebug>
#include <QDialogButtonBox>
#include <qstatusbar.h>
#include <qpushbutton.h>
#include <QDockWidget>
#include "qmainwindow.h"
#include <QJsonParseError>
#include <QJsonDocument>
#include <QJsonObject>
// Slicer includes
#include "qSlicerRepairModuleWidget.h"
#include "ui_qSlicerRepairModuleWidget.h"
#include "vtkSlicerRepairLogic.h"

#include "qSlicerRepairModule.h"
#include "qSlicerRepairModuleWidget.h"
#include "qSlicerRepairModule.h"
#include "qSlicerApplication.h"
//#include "qSlicerRepairTabWidget.h"
#include "qSlicerModuleManager.h"
#include <qSlicerApplication.h>
#include <vtkMRMLScene.h>
#include <vtkMRMLNode.h>
#include <vtkMRMLModelNode.h>
#include <vtkMRMLStorageNode.h>
#include "vtkMRMLLinearTransformNode.h"
#include "qSlicerIOManager.h"
#include <vtkMRMLMarkupsFiducialNode.h>
//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerRepairModuleWidgetPrivate: public Ui_qSlicerRepairModuleWidget
{
public:
  qSlicerRepairModuleWidgetPrivate();
};

//-----------------------------------------------------------------------------
// qSlicerRepairModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerRepairModuleWidgetPrivate::qSlicerRepairModuleWidgetPrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerRepairModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerRepairModuleWidget::qSlicerRepairModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerRepairModuleWidgetPrivate )
{
}

//-----------------------------------------------------------------------------
qSlicerRepairModuleWidget::~qSlicerRepairModuleWidget()
{

}

void qSlicerRepairModuleWidget::enter()
{
	vtkMRMLMarkupsFiducialNode* dentureInteractionFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("dentureInteractionFid"));
	if (dentureInteractionFiducialNode)
	{
		dentureInteractionFiducialNode->SetDisplayVisibility(true);
	}
	vtkMRMLMarkupsFiducialNode* imageFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->GetNodeByID("vtkMRMLMarkupsFiducialNodeImage"));
	if (imageFiducialNode)
	{
		imageFiducialNode->SetDisplayVisibility(true);
	}
	vtkMRMLMarkupsFiducialNode* modelFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->GetNodeByID("vtkMRMLMarkupsFiducialNodeModel"));

	if (modelFiducialNode)
	{
		modelFiducialNode->SetDisplayVisibility(true);
	}
}

void qSlicerRepairModuleWidget::exit()
{
	vtkMRMLMarkupsFiducialNode* dentureInteractionFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("dentureInteractionFid"));
	if (dentureInteractionFiducialNode)
	{
		dentureInteractionFiducialNode->SetDisplayVisibility(false);
		
	}
	vtkMRMLMarkupsFiducialNode* imageFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->GetNodeByID("vtkMRMLMarkupsFiducialNodeImage"));
	if (imageFiducialNode)
	{
		imageFiducialNode->SetDisplayVisibility(false);
	}
	vtkMRMLMarkupsFiducialNode* modelFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->GetNodeByID("vtkMRMLMarkupsFiducialNodeModel"));

	if (modelFiducialNode)
	{
		modelFiducialNode->SetDisplayVisibility(false);
	}
}


//-----------------------------------------------------------------------------
void qSlicerRepairModuleWidget::setup()
{
  Q_D(qSlicerRepairModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();

}


