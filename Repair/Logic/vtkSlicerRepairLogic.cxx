#include "vtkSlicerRepairLogic.h"
#include "vtkSlicerRepairLogic.h"
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
#include "vtkSlicerRepairLogic.h"

// MRML includes
#include <vtkMRMLScene.h>

// VTK includes
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>

// STD includes
#include <cassert>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerRepairLogic);

//----------------------------------------------------------------------------
vtkSlicerRepairLogic::vtkSlicerRepairLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerRepairLogic::~vtkSlicerRepairLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerRepairLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerRepairLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerRepairLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerRepairLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerRepairLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerRepairLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}
//---------------------------------------------------------------------------
void vtkSlicerRepairLogic::SetState(int _state)
{
	moduleState = _state;
}
//---------------------------------------------------------------------------
int vtkSlicerRepairLogic::GetState()
{
	return moduleState;
}
//---------------------------------------------------------------------------
void vtkSlicerRepairLogic::SaveModule(std::string dir)
{
	std::cout << "vtkSlicerRepairLogic::SaveModule()..." << dir << std::endl;
	this->modulePath = dir;
	this->InvokeEvent(vtkSlicerRepairLogic::RepairModuleSaveEvent);

}
std::string vtkSlicerRepairLogic::getModulePath()
{
	return this->modulePath;
}
//---------------------------------------------------------------------------
void vtkSlicerRepairLogic::LoadModule(std::string dir)
{
	std::cout << "vtkSlicerRepairLogic::LoadModule()..." << std::endl;
	this->modulePath = dir;
	this->InvokeEvent(vtkSlicerRepairLogic::RepairModuleLoadEvent);
}
