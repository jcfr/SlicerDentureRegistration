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

// .NAME vtkSlicerRepairLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerRepairLogic_h
#define __vtkSlicerRepairLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// MRML includes

// STD includes
#include <cstdlib>

#include "vtkSlicerRepairModuleLogicExport.h"


/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_REPAIR_MODULE_LOGIC_EXPORT vtkSlicerRepairLogic :
  public vtkSlicerModuleLogic
{
public:
	enum {
		RepairModuleLoadEvent=15500,
		RepairModuleSaveEvent,
  };
  static vtkSlicerRepairLogic *New();
  vtkTypeMacro(vtkSlicerRepairLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);
  enum RepairState
  {
	  DentureRegistratingState,
	  IdleState
  };
  void SetState(int _state);
  int GetState();
  void LoadModule(std::string dir);
  void SaveModule(std::string dir);
  std::string getModulePath();
protected:
  vtkSlicerRepairLogic();
  virtual ~vtkSlicerRepairLogic();

  virtual void SetMRMLSceneInternal(vtkMRMLScene* newScene);
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes();
  virtual void UpdateFromMRMLScene();
  virtual void OnMRMLSceneNodeAdded(vtkMRMLNode* node);
  virtual void OnMRMLSceneNodeRemoved(vtkMRMLNode* node);
private:
  int moduleState;
  std::string modulePath;
  vtkSlicerRepairLogic(const vtkSlicerRepairLogic&); // Not implemented
  void operator=(const vtkSlicerRepairLogic&); // Not implemented
	
};

#endif
