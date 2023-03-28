/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

#ifndef __qSlicerRepairWidget_h
#define __qSlicerRepairWidget_h

// Qt includes
#include <QWidget>
#include <QTableWidget>
#include "qSlicerRepairModuleWidgetsExport.h"

#include <ctkVTKObject.h>
//Eigen includes
#include <itkeigen/Eigen/Eigen>
#include <itkeigen/Eigen/Dense>
#include <itkeigen/Eigen/Core>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include "vtkMRMLSegmentationNode.h"
#include "vtkMRMLScalarVolumeNode.h"
#include "vtkMRMLMarkupsROINode.h"
#include "vtkMRMLVolumeRenderingDisplayNode.h"
#include "vtkSlicerVolumeRenderingLogic.h"
#include "vtkMRMLVolumePropertyNode.h"
#include "vtkMRMLLayoutNode.h"
//#include <iostream>
#include<stdio.h>
#include <vtkMatrix4x4.h>
#include <QList>
#include <vtkVector.h>
#include <math.h>

#include "itkImageToVTKImageFilter.h"
#include "itkJPEGImageIOFactory.h"
#include "itkRGBPixel.h"
#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkVTKImageToImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkPNGImageIOFactory.h"
#include "itksys/SystemTools.hxx"
#include "itkImageRegionIterator.h"
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>  
#include <itkPNGImageIOFactory.h>  
#include <itkConnectedThresholdImageFilter.h>  
#include <itkImageSeriesReader.h>  
#include <itkGDCMImageIO.h>  
#include <itkImageToVTKImageFilter.h>
#include<algorithm>
#include<vector>
#include<iostream>
#include "QList.h"
#include "vtkVector.h"
#include<algorithm>
#include<stdio.h>
#include<vector>
#include<iostream>
#include <math.h>
#include <cmath>
#include <numeric>
#include "QList.h"
#include "vtkVector.h"
using namespace std;


class qSlicerRepairWidgetPrivate;
class QTableWidgetItem;
class QString;
class QJsonObject;

static const int kDimension = 3;
typedef itk::Image<short, kDimension > ImageType;
typedef ImageType::IndexType PixelIndexType;

/// \ingroup Slicer_QtModules_Repair
class Q_SLICER_MODULE_REPAIR_WIDGETS_EXPORT qSlicerRepairWidget
  : public QWidget
{
  Q_OBJECT
  QVTK_OBJECT
public:
  typedef QWidget Superclass;
  qSlicerRepairWidget(QWidget *parent=0);
  virtual ~qSlicerRepairWidget();
  void CreateDentureModel(double worldSeedPoint[3]);
  QList<vtkVector3d> ExtractFiducialCentroidFromCT(QString mode);

public slots:
	void setDicomProperty();
	void onImportDentureDICOM();
	void onImportPatientDICOM();

	void onEditFiducialPoints(QString mode);

	void onEditFiducialPoints(QString mode, QList<vtkVector3d> FiducialList);

	void onPatientVolumeNodeAdded(vtkObject* scene, vtkObject* node);
	void onDentureVolumeNodeAdded(vtkObject* scene, vtkObject* node);


protected slots:

protected:
  QScopedPointer<qSlicerRepairWidgetPrivate> d_ptr;

private:

  Q_DECLARE_PRIVATE(qSlicerRepairWidget);
  Q_DISABLE_COPY(qSlicerRepairWidget);
  void fiducialRegistration();
  void volumeRender();
  ImageType::Pointer  ConstructITKImage(vtkImageData* originalImageData);
  vtkImageData* RegionGrowingByConnectedThredholdFilter(vtkImageData* pImage, double worldSeedPoint[3]);

private slots:
	void onDentureVisibilityClicked(bool visibility);
	void onDentureRegistration();
	void AutoDentureRegistration();


signals:
	void moduleStateChanged();
};

class FiducialTableWidget :public QTableWidget 
{
	Q_OBJECT
	QVTK_OBJECT
public:
	FiducialTableWidget(QWidget* parent = 0);
	~FiducialTableWidget();
	void setActiveMarkupNodeID(std::string newID);

private:
	QStringList columnLabels;
	std::string activeMarkupsNodeID;
private slots:
	void onActiveMarkupTableCurrentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn);
	void onRightClickActiveMarkupTableWidget(QPoint pos);
	void onDeleteMarkupPushButtonClicked(bool confirm=true);
	void onJumpSlicesActionTriggered();
	void updateRow(int m);
	void onUpdateFiducialTable();
	void onActiveMarkupsNodePointModifiedEvent(vtkObject* caller, void* callData);
	void onActiveMarkupsNodePointAddedEvent();
	void onActiveMarkupsNodePointRemovedEvent(vtkObject* caller, void* callData);
	void onActiveMarkupTableCellClicked(QTableWidgetItem* item);
	void onActiveMarkupTableCellChanged(int row, int column);


};



class PointMatcher
{
public:
	PointMatcher();
	~PointMatcher();

	void Point_Registration(QList<vtkVector3d> DentureFiducialList_p,
		QList<vtkVector3d> DentureFiducialList_t,
		QList<vtkVector3d>& DentureFiducialList_1,
		QList<vtkVector3d>& DentureFiducialList_2);

	std::vector<double> patient_x;
	std::vector<double> patient_y;
	std::vector<double> patient_z;
	std::vector<double> tooth_x;
	std::vector<double> tooth_y;
	std::vector<double> tooth_z;
	std::vector<double> error_group;
	std::vector<int> wrong_teeth_num;
	std::vector<std::vector<double>> tooth_dis;
	std::vector<std::vector<double>> ontooth_src;
	std::vector<double> tooth_dis2; 
	std::vector<double> onpatient_dis;
	std::vector<std::vector<double>> patient_dis;
	std::vector<double> patient_point;
	std::vector<double> tar;
	std::vector<std::vector<double>> optional_number;
	std::vector<std::vector<double>> patient_closest_point;
	QList<vtkVector4d> list_ans;
	int teeth;
	int patient;
	double teeth_cen[3];
	double patient_cen[3];
	int* p;
	double min_value;
	bool on_target;

	void point_input(QList<vtkVector3d> DentureFiducialList_p,
		QList<vtkVector3d> DentureFiducialList_t, int p_size, int t_size);
	void onTooth_dis();
	void onPatientgroup_dis();
	void onPatient_dis();
	void select_point(int pos, int process, int totla, int key, std::vector<double>& a, std::vector<int> b, std::vector<int>& p);
	void select_point_2(int pos, int process, int totla, std::vector<std::vector<double>> aa, std::vector<double>& tar, std::vector<double>& ans);
	void Point_Screening(QList<vtkVector3d>& DentureFiducialList_1,
		QList<vtkVector3d>& DentureFiducialList_2);
	void Processing_special_point(QList<vtkVector3d>& DentureFiducialList_1, QList<vtkVector3d>& DentureFiducialList_2);
	void Point_Sorting(QList<vtkVector3d>& DentureFiducialList_1, QList<vtkVector3d>& DentureFiducialList_2);
	double Find_Closest(int mn, int mx, double a, std::vector<double> p);
};

#endif

