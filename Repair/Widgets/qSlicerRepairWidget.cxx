#include "qSlicerRepairWidget.h"
#include "qSlicerRepairWidget.h"
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

#include "qSlicerRepairWidget.h"
#include "ui_qSlicerRepairWidget.h"

//Slicer includes
#include "qSlicerCLIModule.h"
#include "qSlicerApplication.h"
#include "qSlicerCoreIOManager.h"
#include "vtkSlicerModelsLogic.h"
#include "vtkMRMLModelNode.h"
#include "vtkMRMLScene.h"
#include "qSlicerLayoutManager.h"
#include "vtkMRMLDisplayNode.h"
#include "qMRMLThreeDWidget.h"
#include "vtkMRMLViewNode.h"
#include "qMRMLThreeDView.h"
#include "vtkMRMLLinearTransformNode.h"
#include "qSlicerModuleManager.h"
#include "qSlicerAbstractCoreModule.h"
#include "vtkSlicerCLIModuleLogic.h"
#include "vtkMRMLCommandLineModuleNode.h"
#include "qSlicerAbstractModuleRepresentation.h"
#include<vtkMatrix4x4.h>
#include <vtkMRMLMarkupsNode.h>
#include <vtkMRMLSelectionNode.h>
#include <vtkMRMLInteractionNode.h>
#include <vtkMRMLMarkupsFiducialNode.h>
#include "vtkSlicerMarkupsLogic.h"
#include "qSlicerModulePanel.h"
#include <ctkVTKObject.h>
#include "qSlicerAbstractModuleWidget.h"
#include "qMRMLSliceWidget.h"
#include "vtkMRMLSliceNode.h"
#include "vtkSlicerRepairLogic.h"
#include "vtkMRMLAnnotationPointDisplayNode.h"
#include "vtkMRMLScalarVolumeNode.h"
#include "vtkMRMLSegmentationNode.h"
#include"vtkMRMLSegmentEditorNode.h"
#include "qMRMLSegmentEditorWidget.h"
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkIterativeClosestPointTransform.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPoints.h>
#include <vtkLandmarkTransform.h>
#include "vtkSlicerSegmentationsModuleLogic.h"
#include <QlineEdit>
#include <itkeigen/Eigen/Eigen>
#include <itkeigen/Eigen/Dense>
#include <itkeigen/Eigen/Core>
#include <QFileDialog>
#include <qmainwindow.h>
#include <qSlicerModuleSelectorToolBar.h>
#include <qmenu.h>
#include <ctkMessageBox.h>
#include <qtablewidget.h>
#include <qDebug>
#include <QHeaderView.h>
#include <QStatusBar.h>
#include <QDockWidget>
#include <qpushbutton.h>
#include <QJsonObject>
#include <ctkCollapsibleButton.h>
#include"vtkMRMLLabelMapVolumeNode.h"
#include <vtkImageThreshold.h>
#include <QTimer>
#include<stdio.h>
#include <vtkMultiObjectMassProperties.h>
#include <vtkDataArray.h>
#include <vtkFieldData.h>
#include <ostream>
#include <vtkOrientedImageData.h>
#include <PythonQt.h>
#include <QList>
#include <vtkImageCast.h>
#include <vtkITKIslandMath.h>
#include <vtkOrientedImageData.h>

#include "qSlicerSegmentEditorAbstractEffect.h"
#include "qSlicerSegmentEditorAbstractLabelEffect.h"
#include "qSlicerSegmentEditorEffectFactory.h"


//-----------------------------------------------------------------------------
class qSlicerRepairWidgetPrivate
  : public Ui_qSlicerRepairWidget
{
  Q_DECLARE_PUBLIC(qSlicerRepairWidget);
protected:
  qSlicerRepairWidget* const q_ptr;

public:
  qSlicerRepairWidgetPrivate(
    qSlicerRepairWidget& object);
  virtual void setupUi(qSlicerRepairWidget*);
};

// --------------------------------------------------------------------------
qSlicerRepairWidgetPrivate
::qSlicerRepairWidgetPrivate(
  qSlicerRepairWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerRepairWidgetPrivate
::setupUi(qSlicerRepairWidget* widget)
{
  this->Ui_qSlicerRepairWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
qSlicerRepairWidget
::qSlicerRepairWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerRepairWidgetPrivate(*this) )
{
  Q_D(qSlicerRepairWidget);
  d->setupUi(this);


  qSlicerApplication::application()->layoutManager()->setLayout(vtkMRMLLayoutNode::SlicerLayoutDual3DView);

  QObject::connect(d->importPatientDICOMButton, &QPushButton::clicked, this, &qSlicerRepairWidget::onImportPatientDICOM);
  QObject::connect(d->importDentureDICOMButton, &QPushButton::clicked, this, &qSlicerRepairWidget::onImportDentureDICOM);

  QObject::connect(d->imageAddButton, &QPushButton::clicked,
	  [this]() {onEditFiducialPoints(QString::fromStdString("Image")); });
  QObject::connect(d->modelAddButton, &QPushButton::clicked,
	  [this]() {onEditFiducialPoints(QString::fromStdString("Model")); });
  QObject::connect(d->dentureVisibilityCheckBox, &QPushButton::clicked, this, &qSlicerRepairWidget::onDentureVisibilityClicked);
  d->dentureVisibilityCheckBox->setEnabled(false);

}

//-----------------------------------------------------------------------------
qSlicerRepairWidget
::~qSlicerRepairWidget()
{
}

void qSlicerRepairWidget::onImportPatientDICOM()
{
	Q_D(qSlicerRepairWidget);
	qSlicerApplication::application()->layoutManager()->selectModule("DICOM");
	this->qvtkConnect(qSlicerApplication::application()->mrmlScene(), vtkMRMLScene::NodeAddedEvent, this, SLOT(onPatientVolumeNodeAdded(vtkObject*, vtkObject*)));
}

void qSlicerRepairWidget::onPatientVolumeNodeAdded(vtkObject* scene, vtkObject* node)
{
	Q_D(qSlicerRepairWidget);
	vtkMRMLVolumeNode* volumeNode = vtkMRMLVolumeNode::SafeDownCast(node);
	if (!volumeNode)
	{
		return;
	}
	QTimer::singleShot(0, this, SLOT(setDicomProperty()));
	this->qvtkDisconnect(qSlicerApplication::application()->mrmlScene(), vtkMRMLScene::NodeAddedEvent, this, SLOT(onPatientVolumeNodeAdded(vtkObject*, vtkObject*)));
	qSlicerApplication::application()->layoutManager()->selectModule("Repair");

}
void qSlicerRepairWidget::setDicomProperty()
{
	vtkMRMLVolumeNode* volumeNode = vtkMRMLVolumeNode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->GetFirstNodeByClass("vtkMRMLScalarVolumeNode"));
	volumeNode->SetName("OriginalCT");
	this->volumeRender();
}

void qSlicerRepairWidget::volumeRender()
{
	vtkMRMLMarkupsROINode* roiNode =
		vtkMRMLMarkupsROINode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->AddNewNodeByClassWithID("vtkMRMLMarkupsROINode", "GlobalROI", "vtkMRMLMarkupsROINode1"));
	vtkMRMLScalarVolumeNode* volumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("OriginalCT"));

	double bounds[6];
	volumeNode->GetRASBounds(bounds);

	if (roiNode == nullptr)
	{
		roiNode = vtkMRMLMarkupsROINode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->GetNodeByID("vtkMRMLMarkupsROINode1"));
	}
	else
	{
		roiNode->SetXYZ((bounds[1] + bounds[0]) / 2,
			(bounds[3] + bounds[2]) / 2,
			(bounds[5] + bounds[4]) / 2);
		roiNode->SetRadiusXYZ(
			(bounds[1] - bounds[0]) / 2 * 0.75,
			(bounds[3] - bounds[2]) / 2 * 0.75,
			(bounds[5] - bounds[4]) / 2);
		roiNode->SetDisplayVisibility(false);
	}

	qSlicerAbstractCoreModule* volumeRenderingModule =
		qSlicerCoreApplication::application()->moduleManager()->module("VolumeRendering");
	vtkSlicerVolumeRenderingLogic* volumeRenderingLogic =
		volumeRenderingModule ? vtkSlicerVolumeRenderingLogic::SafeDownCast(volumeRenderingModule->logic()) : 0;
	if (volumeRenderingLogic)
	{
		vtkSmartPointer<vtkMRMLVolumeRenderingDisplayNode> displayNode = volumeRenderingLogic->GetFirstVolumeRenderingDisplayNode(volumeNode);
		if (displayNode == nullptr)
		{
			displayNode = vtkSmartPointer<vtkMRMLVolumeRenderingDisplayNode>::Take(volumeRenderingLogic->CreateVolumeRenderingDisplayNode());
		}
		qSlicerApplication::application()->mrmlScene()->AddNode(displayNode);
		volumeNode->AddAndObserveDisplayNodeID(displayNode->GetID());
		volumeRenderingLogic->UpdateDisplayNodeFromVolumeNode(displayNode, volumeNode);
		displayNode->CroppingEnabledOn();
		displayNode->SetAndObserveROINodeID(roiNode->GetID());
		displayNode->SetVisibility(1);
		displayNode->GetVolumePropertyNode()->Copy(volumeRenderingLogic->GetPresetByName("CT-AAA2"));

		std::vector<std::string> viewNodeIDs{ qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("View1")->GetID()};
		displayNode->SetViewNodeIDs(viewNodeIDs);
	}
	qSlicerApplication::application()->layoutManager()->threeDWidget("View1")->threeDView()->resetFocalPoint();
	
}

void qSlicerRepairWidget::onImportDentureDICOM()
{
	Q_D(qSlicerRepairWidget);
	qSlicerApplication::application()->layoutManager()->selectModule("DICOM");
	this->qvtkConnect(qSlicerApplication::application()->mrmlScene(), vtkMRMLScene::NodeAddedEvent, this, SLOT(onDentureVolumeNodeAdded(vtkObject*, vtkObject *)));
}

void qSlicerRepairWidget::onDentureVolumeNodeAdded(vtkObject* scene, vtkObject* node)
{
	Q_D(qSlicerRepairWidget);

	vtkMRMLVolumeNode* volumeNode = vtkMRMLVolumeNode::SafeDownCast(node);
	if (!volumeNode)
	{
		return;
	}
	QTimer::singleShot(0, this, SLOT(AutoDentureRegistration()));
	this->qvtkDisconnect(qSlicerApplication::application()->mrmlScene(), vtkMRMLScene::NodeAddedEvent, this, SLOT(onDentureVolumeNodeAdded(vtkObject*, vtkObject*)));
	qSlicerApplication::application()->layoutManager()->selectModule("Repair");

}

void qSlicerRepairWidget::AutoDentureRegistration()
{
	Q_D(qSlicerRepairWidget);

		int NumofVolume = qSlicerApplication::application()->mrmlScene()->GetNumberOfNodesByClass("vtkMRMLScalarVolumeNode");
		vtkMRMLScalarVolumeNode* DentureVolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->GetNthNodeByClass(NumofVolume - 1, "vtkMRMLScalarVolumeNode"));
		DentureVolumeNode->SetName("DentureCT");
		vtkMRMLScalarVolumeNode* PatientVolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("OriginalCT"));

		if (!PatientVolumeNode)
		{
			int ret = ctkMessageBox::warning(this, tr("Dental"),
				tr("Please import Dicom images of a new patient first. \n"),
				ctkMessageBox::Ok);

			return;
		}

		QList<vtkVector3d> DentureOriginalList = ExtractFiducialCentroidFromCT("Denture");
		CreateDentureModel(DentureOriginalList[0].GetData());

		QList<vtkVector3d> PatientOriginalList = ExtractFiducialCentroidFromCT("Patient");

		QList<vtkVector3d> PatientFiducialList;
		QList<vtkVector3d> DentureFiducialList;

		PointMatcher matcher;
		matcher.Point_Registration(PatientOriginalList, DentureOriginalList, PatientFiducialList, DentureFiducialList);


		onEditFiducialPoints("Model", DentureFiducialList);
		onEditFiducialPoints("Image", PatientFiducialList);

		onDentureRegistration();

}

void qSlicerRepairWidget::CreateDentureModel(double worldSeedPoint[3])
{
	Q_D(qSlicerRepairWidget);

	vtkMRMLScalarVolumeNode* DentureVolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("DentureCT"));
	vtkMRMLSegmentationNode* segmentationNode = vtkMRMLSegmentationNode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->AddNewNodeByClass("vtkMRMLSegmentationNode"));
	segmentationNode->SetName("DentureModelSegmentationNode");
	segmentationNode->CreateDefaultDisplayNodes();
	segmentationNode->SetReferenceImageGeometryParameterFromVolumeNode(DentureVolumeNode);
	segmentationNode->GetSegmentation()->AddEmptySegment("DentureWhole");
	qMRMLSegmentEditorWidget* segmentEditorWidget=new qMRMLSegmentEditorWidget;
	segmentEditorWidget->setMRMLScene(qSlicerApplication::application()->mrmlScene());
	vtkMRMLSegmentEditorNode * segmentEditorNode = vtkMRMLSegmentEditorNode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->AddNewNodeByClass("vtkMRMLSegmentEditorNode"));
	segmentEditorWidget->setMRMLSegmentEditorNode(segmentEditorNode);
	segmentEditorWidget->setSegmentationNode(segmentationNode);
	segmentEditorWidget->setMasterVolumeNode(DentureVolumeNode);
	//thresh
	segmentEditorWidget->setActiveEffectByName("Threshold");
	qSlicerSegmentEditorAbstractEffect* effect = segmentEditorWidget->activeEffect();
	effect->setParameter("MinimumThreshold", "-600");
	effect->setParameter("MaximumThreshold", "5000");
	vtkOrientedImageData* masterImageData = effect->masterVolumeImageData();
	auto modifierLabelmap = effect->defaultModifierLabelmap();
	vtkSmartPointer<vtkMatrix4x4> originalImageToWorldMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
	modifierLabelmap->GetImageToWorldMatrix(originalImageToWorldMatrix);
	auto min = effect->doubleParameter("MinimumThreshold");
	auto max = effect->doubleParameter("MaximumThreshold");
	effect->saveStateForUndo();
	vtkSmartPointer<vtkImageThreshold> thresh = vtkSmartPointer<vtkImageThreshold>::New();
	thresh->SetInputData(masterImageData);
	thresh->ThresholdBetween(min, max);
	thresh->SetInValue(1);
	thresh->SetOutValue(0);
	thresh->Update();

	double seed[3] = { -worldSeedPoint[0],-worldSeedPoint[1],worldSeedPoint[2] };
	vtkImageData* RGimage=RegionGrowingByConnectedThredholdFilter(thresh->GetOutput(), seed);

	modifierLabelmap->DeepCopy(RGimage);
	effect->modifySelectedSegmentByLabelmap(modifierLabelmap, qSlicerSegmentEditorAbstractEffect::ModificationModeSet);
	effect->selectEffect("");

	segmentationNode->CreateClosedSurfaceRepresentation();
	segmentationNode->GetDisplayNode()->VisibilityOff();
	vtkSmartPointer<vtkPolyData> outputPolyData = vtkSmartPointer<vtkPolyData>::New();
	auto segmentId = segmentationNode->GetSegmentation()->GetNthSegmentID(0);
	segmentationNode->GetClosedSurfaceRepresentation(segmentId, outputPolyData);
	segmentationNode->RemoveClosedSurfaceRepresentation();

	if (qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("Denture"))
		qSlicerApplication::application()->mrmlScene()->RemoveNode(qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("Denture"));
	
	vtkSlicerModelsLogic* modelsLogic = vtkSlicerModelsLogic::SafeDownCast(qSlicerApplication::application()->moduleLogic("Models"));
	vtkMRMLModelNode* dentureModelNode = modelsLogic->AddModel(outputPolyData);
	dentureModelNode->SetName("Denture");
	dentureModelNode->GetDisplayNode()->SetColor(1, 1, 0);
	dentureModelNode->GetDisplayNode()->SetScalarVisibility(false);
	dentureModelNode->GetDisplayNode()->SetDisplayableOnlyInView(
		qSlicerApplication::application()->layoutManager()->threeDWidget("View2")->mrmlViewNode()->GetID());
	qSlicerApplication::application()->layoutManager()->threeDWidget("View2")->threeDView()->resetFocalPoint();

}

ImageType::Pointer  qSlicerRepairWidget::ConstructITKImage(vtkImageData* originalImageData)
{
	ImageType::Pointer m_itkImage = ImageType::New();
	double imageSpacing[3] = { 0.0 };
	double imageOrigin[3] = { 0.0 };
	int imageDimension[3] = { 0 };
	originalImageData->GetSpacing(imageSpacing);
	originalImageData->GetOrigin(imageOrigin);
	originalImageData->GetDimensions(imageDimension);

	const ImageType::SizeType  size = { {imageDimension[0], imageDimension[1], imageDimension[2]} }; 
	const ImageType::IndexType start = { {0,0,0} };

	ImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(start);

	typedef itk::VTKImageToImageFilter<ImageType> VTKImageToImageType;

	//Converting to ITK Image Format
	VTKImageToImageType::Pointer vtkImageToImageFilter = VTKImageToImageType::New();
	vtkImageToImageFilter->SetInput(originalImageData);
	vtkImageToImageFilter->UpdateLargestPossibleRegion();
	//vtkImageToImageFilter->Update();

	m_itkImage->SetRegions(region);
	m_itkImage->Allocate();
	m_itkImage = const_cast<itk::Image<short, kDimension>*>(vtkImageToImageFilter->GetImporter()->GetOutput());

	ImageType::SpacingType spacing;
	spacing[0] = imageSpacing[0]; // spacing along X
	spacing[1] = imageSpacing[1]; // spacing along Y
	spacing[2] = imageSpacing[2]; // spacing along Z
	m_itkImage->SetSpacing(spacing);

	ImageType::PointType newOrigin;
	newOrigin[0] = imageOrigin[0];
	newOrigin[1] = imageOrigin[1];
	newOrigin[2] = imageOrigin[2];
	m_itkImage->SetOrigin(newOrigin);

	ImageType::DirectionType direction;
	direction.SetIdentity();
	m_itkImage->SetDirection(direction);
	m_itkImage->UpdateOutputInformation();

	return m_itkImage;
};

vtkImageData* qSlicerRepairWidget::RegionGrowingByConnectedThredholdFilter(vtkImageData* pImage, double worldSeedPoint[3])
{

	ImageType::Pointer itkimage = ConstructITKImage(pImage);

	ImageType::IndexType pixelIndex;
	typedef itk::Point<double, ImageType::ImageDimension > PointType;
	PointType point;
	point[0] = worldSeedPoint[0];    // x coordinate
	point[1] = worldSeedPoint[1];    // y coordinate
	point[2] = worldSeedPoint[2];    // z coordinate

	bool isSeedPointInsideImage = itkimage->TransformPhysicalPointToIndex(point, pixelIndex);

	if (!isSeedPointInsideImage)
	{
		exit;
	}

	typedef itk::Image< short, kDimension > InternalImageType;
	typedef itk::ConnectedThresholdImageFilter< InternalImageType, InternalImageType > ConnectedFilterType;
	ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();

	connectedThreshold->SetInput(itkimage);

	connectedThreshold->SetLower(1);
	connectedThreshold->SetUpper(1);
	connectedThreshold->SetReplaceValue(1);
	connectedThreshold->SetSeed(pixelIndex);

	typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;

	//Converting Back from ITK to VTK Image for Visualization.
	ConnectorType::Pointer connector = ConnectorType::New();
	connector->SetInput(connectedThreshold->GetOutput());
	//connector->SetInput(itkimage);//reader->GetOutput()

	try
	{
		connector->Update();
	}
	catch (itk::ExceptionObject& error)
	{
		//std::cerr << "Error: " << error << std::endl;
	}
	vtkImageData* img = vtkImageData::New();
	img->DeepCopy(connector->GetOutput());

	return img;
}

QList<vtkVector3d> qSlicerRepairWidget::ExtractFiducialCentroidFromCT(QString mode)
{
	Q_D(qSlicerRepairWidget);

	vtkMRMLScalarVolumeNode* VolumeNode;

	if (mode == "Patient")
	{
		VolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("OriginalCT"));
	}
	else if (mode == "Denture")
	{
		VolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("DentureCT"));

	}

	vtkMRMLSegmentationNode* segmentationNode = vtkMRMLSegmentationNode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->AddNewNodeByClass("vtkMRMLSegmentationNode"));
	segmentationNode->CreateDefaultDisplayNodes();
	segmentationNode->SetReferenceImageGeometryParameterFromVolumeNode(VolumeNode);
	segmentationNode->GetSegmentation()->AddEmptySegment("Segment");

	qMRMLSegmentEditorWidget* segmentEditorWidget = new qMRMLSegmentEditorWidget;
	segmentEditorWidget->setMRMLScene(qSlicerApplication::application()->mrmlScene());
	vtkMRMLSegmentEditorNode* segmentEditorNode = vtkMRMLSegmentEditorNode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->AddNewNodeByClass("vtkMRMLSegmentEditorNode"));
	segmentEditorWidget->setMRMLSegmentEditorNode(segmentEditorNode);
	segmentEditorWidget->setSegmentationNode(segmentationNode);
	segmentEditorWidget->setMasterVolumeNode(VolumeNode);
	//thresh
	segmentEditorWidget->setActiveEffectByName("Threshold");
	qSlicerSegmentEditorAbstractEffect* effect = segmentEditorWidget->activeEffect();

	QString str;
	if (mode == "Patient")
	{
		str = QString("%1").arg(d->PatientSlider->value());

	}
	else if (mode == "Denture")
	{
		str = QString("%1").arg(d->DentureSlider->value());
	}
	effect->setParameter("MinimumThreshold", str);
	effect->setParameter("MaximumThreshold", "5000");

	vtkOrientedImageData* masterImageData = effect->masterVolumeImageData();
	auto modifierLabelmap = effect->defaultModifierLabelmap();
	vtkSmartPointer<vtkMatrix4x4> originalImageToWorldMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
	modifierLabelmap->GetImageToWorldMatrix(originalImageToWorldMatrix);
	auto min = effect->doubleParameter("MinimumThreshold");
	auto max = effect->doubleParameter("MaximumThreshold");
	effect->saveStateForUndo();
	vtkSmartPointer<vtkImageThreshold> thresh = vtkSmartPointer<vtkImageThreshold>::New();
	thresh->SetInputData(masterImageData);
	thresh->ThresholdBetween(min, max);
	thresh->SetInValue(1);
	thresh->SetOutValue(0);
	thresh->Update();

	modifierLabelmap->DeepCopy(thresh->GetOutput());
	effect->modifySelectedSegmentByLabelmap(modifierLabelmap, qSlicerSegmentEditorAbstractEffect::ModificationModeSet);
	effect->selectEffect("");

	segmentationNode->CreateClosedSurfaceRepresentation();
	segmentationNode->GetDisplayNode()->VisibilityOff();

	vtkSmartPointer<vtkPolyData> outputPolyData = vtkSmartPointer<vtkPolyData>::New();
	auto segmentId = segmentationNode->GetSegmentation()->GetNthSegmentID(0);
	segmentationNode->GetClosedSurfaceRepresentation(segmentId, outputPolyData);
	segmentationNode->RemoveClosedSurfaceRepresentation();

	
	vtkSmartPointer<vtkMultiObjectMassProperties> getMass = vtkSmartPointer<vtkMultiObjectMassProperties>::New();
	if (mode == "Denture")
	{
		getMass->SetInputData(outputPolyData);

		vtkSlicerModelsLogic* modelsLogic = vtkSlicerModelsLogic::SafeDownCast(qSlicerApplication::application()->moduleLogic("Models"));
		vtkMRMLModelNode* dentureFiducialModelNode = modelsLogic->AddModel(outputPolyData);
		dentureFiducialModelNode->SetName("DentureFiducialPoint");

		dentureFiducialModelNode->GetDisplayNode()->SetDisplayableOnlyInView(
			qSlicerApplication::application()->layoutManager()->threeDWidget("View2")->mrmlViewNode()->GetID());
	}
	else if (mode == "Patient")
	{
		getMass->SetInputData(outputPolyData);

	}

	getMass->Update();
	auto numObjects = getMass->GetNumberOfObjects();
	auto centroidArray = getMass->GetOutput()->GetFieldData()->GetArray("ObjectCentroids");
	auto volArray = getMass->GetOutput()->GetFieldData()->GetArray("ObjectVolumes");
	auto areaArray = getMass->GetOutput()->GetFieldData()->GetArray("ObjectAreas");


	QList<vtkVector3d> FiducialList;
	for (vtkIdType i = 0; i < numObjects; ++i)
	{
		double volume = volArray->GetTuple1(i);
		double area = areaArray->GetTuple1(i);
		double* centroid = centroidArray->GetTuple3(i);
		if (mode == "Denture")
		{
			vtkVector3d p0{ centroid[0] ,centroid[1] ,centroid[2] };
			FiducialList.append(p0);
		}
		else if (mode == "Patient")
		{
			if (volume < 10 && volume>1e-5 && pow(area, 1.5) / volume < 60)
			{
				vtkVector3d p0{ centroid[0] ,centroid[1] ,centroid[2] };
				FiducialList.append(p0);
			}
		}
	}

	return FiducialList;

}



void qSlicerRepairWidget::fiducialRegistration()
{
	Q_D(qSlicerRepairWidget);
	vtkMRMLLinearTransformNode* transformNode = vtkMRMLLinearTransformNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->AddNewNodeByClassWithID(
			"vtkMRMLLinearTransformNode", "Repair Transform", "vtkMRMLLinearTransformNodeRepair"));

	if (transformNode == nullptr)
	{
		transformNode = vtkMRMLLinearTransformNode::SafeDownCast(qSlicerApplication::application()->mrmlScene()->GetNodeByID("vtkMRMLLinearTransformNodeRepair"));
	}

	qSlicerCLIModule* regModule = static_cast<qSlicerCLIModule*> (qSlicerCoreApplication::application()->moduleManager()->module("FiducialRegistration"));//qSlicerCoreApplication
	vtkSlicerCLIModuleLogic* moduleLogic = vtkSlicerCLIModuleLogic::SafeDownCast(regModule->logic());

	vtkMRMLCommandLineModuleNode* cliNode = moduleLogic->CreateNodeInScene();
	cliNode->SetParameterAsString("saveTransform", transformNode->GetID());
	cliNode->SetParameterAsString("movingLandmarks", qSlicerApplication::application()->mrmlScene()
		->GetFirstNodeByName("ModelFiducial")->GetID());
	cliNode->SetParameterAsString("fixedLandmarks", qSlicerApplication::application()->mrmlScene()
		->GetFirstNodeByName("ImageFiducial")->GetID());
	moduleLogic->ApplyAndWait(cliNode);

	qSlicerApplication::application()->mainWindow()->statusBar()->showMessage(
		QString::fromStdString("rms = ") + QString::fromStdString(cliNode->GetParameterAsString("rms")) + "mm");
	QString bb = QString::fromStdString("rms = ") + QString::fromStdString(cliNode->GetParameterAsString("rms")) + "mm";
    d->textEdit->insertPlainText(bb );


}

void qSlicerRepairWidget::onDentureRegistration()
{
	Q_D(qSlicerRepairWidget);

	vtkMRMLModelNode* modelNode = vtkMRMLModelNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("Denture"));
	if (modelNode)
	{
		this->fiducialRegistration();

		vtkMRMLLinearTransformNode* registrationTransformNode = vtkMRMLLinearTransformNode::SafeDownCast(
			qSlicerApplication::application()->mrmlScene()->GetNodeByID("vtkMRMLLinearTransformNodeRepair"));

		if (qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("DentureTransformed"))
			qSlicerApplication::application()->mrmlScene()->RemoveNode(
				qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("DentureTransformed"));
		vtkSlicerModelsLogic* modelsLogic = vtkSlicerModelsLogic::SafeDownCast(
			qSlicerApplication::application()->moduleLogic("Models"));
		vtkMRMLModelNode* dentureTransformedModelNode = modelsLogic->AddModel(modelNode->GetPolyData());
		dentureTransformedModelNode->SetName("DentureTransformed");
		std::vector<std::string> viewNodeIDs{
			qSlicerApplication::application()->layoutManager()->sliceWidget("Red")->mrmlSliceNode()->GetID(),
			qSlicerApplication::application()->layoutManager()->sliceWidget("Green")->mrmlSliceNode()->GetID(),
			qSlicerApplication::application()->layoutManager()->sliceWidget("Yellow")->mrmlSliceNode()->GetID(),
			qSlicerApplication::application()->layoutManager()->threeDWidget("View1")->mrmlViewNode()->GetID() };
		dentureTransformedModelNode->GetDisplayNode()->SetViewNodeIDs(viewNodeIDs);
		dentureTransformedModelNode->GetDisplayNode()->SliceIntersectionVisibilityOn();
		dentureTransformedModelNode->GetDisplayNode()->SetColor(1, 1, 0);
		dentureTransformedModelNode->GetDisplayNode()->SetScalarVisibility(false);
		dentureTransformedModelNode->SetAndObserveTransformNodeID(registrationTransformNode->GetID());

		d->dentureVisibilityCheckBox->setEnabled(true);
		vtkMRMLMarkupsFiducialNode* dentureInteractionFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
			qSlicerApplication::application()->mrmlScene()->AddNewNodeByClassWithID(
				"vtkMRMLMarkupsFiducialNode", "dentureInteractionFid", "dentureInteractionFid"));
	
	}
	else
	{
		ctkMessageBox::warning(this, tr("Dental"), tr("Please import denture before applying registration."), QMessageBox::Yes);
	}
}



void qSlicerRepairWidget::onEditFiducialPoints(QString mode)
{
	Q_D(qSlicerRepairWidget);
	vtkMRMLScene* scene = qSlicerApplication::application()->mrmlScene();

	QDockWidget* dockWidget = qSlicerApplication::application()->mainWindow()->findChild<QDockWidget*>("PanelDockWidget");
	dockWidget->show();

	if (mode=="Image")
	{
		qSlicerApplication::application()->mainWindow()->statusBar()->showMessage(
			"Left click: add image fiducials, Right click: finish");

		vtkMRMLMarkupsFiducialNode* imageFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
			scene->AddNewNodeByClassWithID("vtkMRMLMarkupsFiducialNode", "ImageFiducial", "vtkMRMLMarkupsFiducialNodeImage"));
		if (imageFiducialNode == nullptr)
			imageFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(scene->GetNodeByID("vtkMRMLMarkupsFiducialNodeImage"));
		else
		{
			std::vector<std::string> viewNodeIDs{
				qSlicerApplication::application()->layoutManager()->sliceWidget("Axial")->mrmlSliceNode()->GetID(),
				qSlicerApplication::application()->layoutManager()->sliceWidget("Coronal")->mrmlSliceNode()->GetID(),
				qSlicerApplication::application()->layoutManager()->sliceWidget("Sagittal")->mrmlSliceNode()->GetID(),
				qSlicerApplication::application()->layoutManager()->threeDWidget("View1")->mrmlViewNode()->GetID() };
			imageFiducialNode->GetMarkupsDisplayNode()->SetViewNodeIDs(viewNodeIDs);
			imageFiducialNode->GetMarkupsDisplayNode()->SetSelectedColor(1, 0, 0);
			imageFiducialNode->GetMarkupsDisplayNode()->SetUseGlyphScale(false);
			imageFiducialNode->GetMarkupsDisplayNode()->SetGlyphSize(1);
			imageFiducialNode->SetMarkupLabelFormat("I-%d");
			imageFiducialNode->GetMarkupsDisplayNode()->SetTextScale(4);
			imageFiducialNode->GetMarkupsDisplayNode()->SetGlyphType(6);
			imageFiducialNode->GetMarkupsDisplayNode()->PointLabelsVisibilityOn();
			MRMLNodeModifyBlocker blocker(imageFiducialNode);
			vtkSmartPointer<vtkTextProperty> textProperty = imageFiducialNode->GetMarkupsDisplayNode()->GetTextProperty(); 
			textProperty->SetFontFamily(2);
			textProperty->SetBold(false);
			textProperty->SetItalic(false);
			textProperty->SetShadow(false);
			textProperty->SetBackgroundOpacity(0);
			d->imageFiducialTableWidget->setActiveMarkupNodeID("vtkMRMLMarkupsFiducialNodeImage");
		}
		vtkMRMLSelectionNode* selectionNode = vtkMRMLSelectionNode::SafeDownCast(scene->GetNodeByID("vtkMRMLSelectionNodeSingleton"));
		selectionNode->SetReferenceActivePlaceNodeClassName("vtkMRMLMarkupsFiducialNode");
		selectionNode->SetReferenceActivePlaceNodeID("vtkMRMLMarkupsFiducialNodeImage");
		vtkMRMLInteractionNode* interactionNode = vtkMRMLInteractionNode::SafeDownCast(
			scene->GetNodeByID("vtkMRMLInteractionNodeSingleton"));
		int placeModePersistence = 1;
		interactionNode->SetPlaceModePersistence(placeModePersistence);
		interactionNode->SetCurrentInteractionMode(1);
		interactionNode->SwitchToPersistentPlaceMode();
	}
	else if (mode=="Model")
	{
		qSlicerApplication::application()->mainWindow()->statusBar()->showMessage(
			"Left click: add model fiducials, Right click: finish");
		vtkMRMLMarkupsFiducialNode* modelFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
			scene->AddNewNodeByClassWithID("vtkMRMLMarkupsFiducialNode", "ModelFiducial", "vtkMRMLMarkupsFiducialNodeModel"));
		if (modelFiducialNode == nullptr)
			modelFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(scene->GetNodeByID("vtkMRMLMarkupsFiducialNodeModel"));
		else
		{
			modelFiducialNode->GetMarkupsDisplayNode()->SetSelectedColor(0, 0, 1);
			modelFiducialNode->GetMarkupsDisplayNode()->SetDisplayableOnlyInView(
				qSlicerApplication::application()->layoutManager()->threeDWidget("View2")->mrmlViewNode()->GetID());
			modelFiducialNode->GetMarkupsDisplayNode()->SetUseGlyphScale(false);
			modelFiducialNode->GetMarkupsDisplayNode()->SetGlyphSize(1);
			modelFiducialNode->SetMarkupLabelFormat("M-%d");
			modelFiducialNode->GetMarkupsDisplayNode()->SetTextScale(4);
			modelFiducialNode->GetMarkupsDisplayNode()->SetGlyphType(6);
			modelFiducialNode->GetMarkupsDisplayNode()->PointLabelsVisibilityOn();
			MRMLNodeModifyBlocker blocker(modelFiducialNode);
			vtkSmartPointer<vtkTextProperty> textProperty = modelFiducialNode->GetMarkupsDisplayNode()->GetTextProperty(); 
			textProperty->SetFontFamily(2);
			textProperty->SetBold(false);
			textProperty->SetItalic(false);
			textProperty->SetShadow(false);
			textProperty->SetBackgroundOpacity(0);
			d->modelFiducialTableWidget->setActiveMarkupNodeID("vtkMRMLMarkupsFiducialNodeModel");
		}

		vtkMRMLSelectionNode* selectionNode = vtkMRMLSelectionNode::SafeDownCast(scene->GetNodeByID("vtkMRMLSelectionNodeSingleton"));
		selectionNode->SetReferenceActivePlaceNodeClassName("vtkMRMLMarkupsFiducialNode");
		selectionNode->SetReferenceActivePlaceNodeID("vtkMRMLMarkupsFiducialNodeModel");
	
		vtkMRMLInteractionNode* interactionNode = vtkMRMLInteractionNode::SafeDownCast(scene->GetNodeByID("vtkMRMLInteractionNodeSingleton"));
		int placeModePersistence = 1;
		interactionNode->SetPlaceModePersistence(placeModePersistence);
		interactionNode->SetCurrentInteractionMode(1);
		interactionNode->SwitchToPersistentPlaceMode();
	}


}
void qSlicerRepairWidget::onEditFiducialPoints(QString mode, QList<vtkVector3d> FiducialList)
{
	Q_D(qSlicerRepairWidget);
	vtkMRMLScene* scene = qSlicerApplication::application()->mrmlScene();
	int p_size = FiducialList.size();


	if (mode == "Image")
	{
		vtkMRMLMarkupsFiducialNode* imageFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
			scene->AddNewNodeByClassWithID("vtkMRMLMarkupsFiducialNode", "ImageFiducial", "vtkMRMLMarkupsFiducialNodeImage"));
		if (imageFiducialNode == nullptr)
			imageFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(scene->GetNodeByID("vtkMRMLMarkupsFiducialNodeImage"));
		else
		{
			std::vector<std::string> viewNodeIDs{
				qSlicerApplication::application()->layoutManager()->threeDWidget("View1")->mrmlViewNode()->GetID() };
			imageFiducialNode->GetMarkupsDisplayNode()->SetViewNodeIDs(viewNodeIDs);
			imageFiducialNode->GetMarkupsDisplayNode()->SetSelectedColor(1, 0, 0);
			imageFiducialNode->GetMarkupsDisplayNode()->SetUseGlyphScale(false);
			imageFiducialNode->GetMarkupsDisplayNode()->SetGlyphSize(1);
			imageFiducialNode->SetMarkupLabelFormat("I-%d");
			imageFiducialNode->GetMarkupsDisplayNode()->SetTextScale(4);
			imageFiducialNode->GetMarkupsDisplayNode()->SetGlyphType(6);
			imageFiducialNode->GetMarkupsDisplayNode()->PointLabelsVisibilityOn();
			for (int i = 0; i < p_size; i++)
			{
				vtkVector3d p0 = FiducialList.at(i);
				imageFiducialNode->AddControlPointWorld(p0);
			}
			MRMLNodeModifyBlocker blocker(imageFiducialNode);
			vtkSmartPointer<vtkTextProperty> textProperty = imageFiducialNode->GetMarkupsDisplayNode()->GetTextProperty();
			textProperty->SetFontFamily(2);
			textProperty->SetBold(false);
			textProperty->SetItalic(false);
			textProperty->SetShadow(false);
			textProperty->SetBackgroundOpacity(0);
			d->imageFiducialTableWidget->setActiveMarkupNodeID("vtkMRMLMarkupsFiducialNodeImage");
		}
	}
	else if (mode == "Model")
	{
		vtkMRMLMarkupsFiducialNode* modelFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
			scene->AddNewNodeByClassWithID("vtkMRMLMarkupsFiducialNode", "ModelFiducial", "vtkMRMLMarkupsFiducialNodeModel"));
		if (modelFiducialNode == nullptr)
			modelFiducialNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(scene->GetNodeByID("vtkMRMLMarkupsFiducialNodeModel"));
		else
		{
			modelFiducialNode->GetMarkupsDisplayNode()->SetSelectedColor(0, 0, 1);
			modelFiducialNode->GetMarkupsDisplayNode()->SetDisplayableOnlyInView(
				qSlicerApplication::application()->layoutManager()->threeDWidget("View2")->mrmlViewNode()->GetID());
			modelFiducialNode->GetMarkupsDisplayNode()->SetUseGlyphScale(false);
			modelFiducialNode->GetMarkupsDisplayNode()->SetGlyphSize(1);
			modelFiducialNode->SetMarkupLabelFormat("M-%d");
			modelFiducialNode->GetMarkupsDisplayNode()->SetTextScale(4);
			modelFiducialNode->GetMarkupsDisplayNode()->SetGlyphType(6);
			modelFiducialNode->GetMarkupsDisplayNode()->PointLabelsVisibilityOn();

			for (int i = 0; i < p_size; i++)
			{
				vtkVector3d p0 = FiducialList.at(i);
				modelFiducialNode->AddControlPointWorld(p0);

			}

			MRMLNodeModifyBlocker blocker(modelFiducialNode);
			vtkSmartPointer<vtkTextProperty> textProperty = modelFiducialNode->GetMarkupsDisplayNode()->GetTextProperty();
			textProperty->SetFontFamily(2);
			textProperty->SetBold(false);
			textProperty->SetItalic(false);
			textProperty->SetShadow(false);
			textProperty->SetBackgroundOpacity(0);
			d->modelFiducialTableWidget->setActiveMarkupNodeID("vtkMRMLMarkupsFiducialNodeModel");
			qvtkConnect(scene->GetNodeByID("vtkMRMLMarkupsFiducialNodeModel"), vtkMRMLMarkupsNode::PointPositionDefinedEvent,
				this, SLOT(onIsFiducialsEqualBetweenModelAndImage()));
			qvtkConnect(scene->GetNodeByID("vtkMRMLMarkupsFiducialNodeModel"), vtkMRMLMarkupsNode::PointRemovedEvent,
				this, SLOT(onIsFiducialsEqualBetweenModelAndImage()));
		}


	}


}



void qSlicerRepairWidget::onDentureVisibilityClicked(bool visibility)
{
	Q_D(qSlicerRepairWidget);
	vtkMRMLModelNode* dentureModelNode = vtkMRMLModelNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->GetFirstNodeByName("DentureTransformed"));
	dentureModelNode->GetDisplayNode()->SetVisibility(visibility);
	
}

//-----------------------------------------------------------------------------

FiducialTableWidget::FiducialTableWidget(QWidget* parent):
	QTableWidget(parent)
{
	QObject::connect(this, SIGNAL(cellChanged(int, int)),
		this, SLOT(onActiveMarkupTableCellChanged(int, int)));

	QObject::connect(this, SIGNAL(itemClicked(QTableWidgetItem*)),
		this, SLOT(onActiveMarkupTableCellClicked(QTableWidgetItem*)));
	QObject::connect(this, SIGNAL(currentCellChanged(int, int, int, int)),
		this, SLOT(onActiveMarkupTableCurrentCellChanged(int, int, int, int)));
	this->setContextMenuPolicy(Qt::CustomContextMenu);
	QObject::connect(this, SIGNAL(customContextMenuRequested(QPoint)),
		this, SLOT(onRightClickActiveMarkupTableWidget(QPoint)));

	this->qvtkConnect(qSlicerApplication::application()->mrmlScene(), vtkMRMLScene::EndImportEvent,
		this, SLOT(onUpdateFiducialTable()));
	this->qvtkConnect(qSlicerApplication::application()->mrmlScene(), vtkMRMLScene::EndBatchProcessEvent,
		this, SLOT(onUpdateFiducialTable()));
	this->qvtkConnect(qSlicerApplication::application()->mrmlScene(), vtkMRMLScene::EndCloseEvent,
		this, SLOT(onUpdateFiducialTable()));
	this->qvtkConnect(qSlicerApplication::application()->mrmlScene(), vtkMRMLScene::EndRestoreEvent,
		this, SLOT(onUpdateFiducialTable()));
	this->setMinimumSize(QSize(0, 250));
	this->setAlternatingRowColors(true);
	this->setSelectionBehavior(QAbstractItemView::SelectRows);
	this->setSelectionMode(QAbstractItemView::ExtendedSelection);

	columnLabels << "Visible" << "Name" << "R" << "A" << "S";

	this->setColumnCount(columnLabels.size());
	this->setHorizontalHeaderLabels(columnLabels);
	this->horizontalHeader()->setFixedHeight(32);

	this->setColumnWidth(columnLabels.indexOf("Name"), 60);
	this->setColumnWidth(columnLabels.indexOf("R"), 65);
	this->setColumnWidth(columnLabels.indexOf("A"), 65);
	this->setColumnWidth(columnLabels.indexOf("S"), 65);

	QTableWidgetItem* visibleHeader = this->horizontalHeaderItem(columnLabels.indexOf("Visible"));
	visibleHeader->setText("");
	visibleHeader->setIcon(QIcon(":/Icons/Small/SlicerVisibleInvisible.png"));
	visibleHeader->setToolTip(QString("Click in this column to show/hide markups in 2D and 3D"));
	this->setColumnWidth(columnLabels.indexOf("Visible"), 30);

	this->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
	this->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
}

FiducialTableWidget::~FiducialTableWidget()
{
}

void FiducialTableWidget::setActiveMarkupNodeID(std::string newID)
{
	if (activeMarkupsNodeID != newID)
	{
		vtkMRMLMarkupsFiducialNode* activeMarkupsNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
			qSlicerApplication::application()->mrmlScene()->GetNodeByID(activeMarkupsNodeID.c_str()));
		vtkMRMLMarkupsFiducialNode* newActiveMarkupsNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
			qSlicerApplication::application()->mrmlScene()->GetNodeByID(newID));
		qvtkReconnect(activeMarkupsNode, newActiveMarkupsNode, vtkCommand::ModifiedEvent,
			this, SLOT(onUpdateFiducialTable()));
		qvtkReconnect(activeMarkupsNode, newActiveMarkupsNode, vtkMRMLMarkupsNode::PointModifiedEvent,
			this, SLOT(onActiveMarkupsNodePointModifiedEvent(vtkObject*, void*)));
		qvtkReconnect(activeMarkupsNode, newActiveMarkupsNode, vtkMRMLMarkupsNode::PointAddedEvent,
			this, SLOT(onActiveMarkupsNodePointAddedEvent()));
		qvtkReconnect(activeMarkupsNode, newActiveMarkupsNode, vtkMRMLMarkupsNode::PointRemovedEvent,
			this, SLOT(onActiveMarkupsNodePointRemovedEvent(vtkObject*, void*)));
	}
	activeMarkupsNodeID = newID;
	this->onUpdateFiducialTable();
}

void FiducialTableWidget::onActiveMarkupTableCurrentCellChanged(
	int currentRow, int currentColumn, int previousRow, int previousColumn)
{
	Q_UNUSED(currentColumn);
	Q_UNUSED(previousRow);
	Q_UNUSED(previousColumn);
}

void FiducialTableWidget::onRightClickActiveMarkupTableWidget(QPoint pos)
{
	Q_UNUSED(pos);
	QMenu menu;

	QAction* deleteFiducialAction =
		new QAction(QString("Delete fiducials"), &menu);
	menu.addAction(deleteFiducialAction);
	QObject::connect(deleteFiducialAction, SIGNAL(triggered()),
		this, SLOT(onDeleteMarkupPushButtonClicked()));

	QAction* jumpSlicesAction =
		new QAction(QString("Jump slices"), &menu);
	menu.addAction(jumpSlicesAction);
	QObject::connect(jumpSlicesAction, SIGNAL(triggered()),
		this, SLOT(onJumpSlicesActionTriggered()));

	menu.exec(QCursor::pos());
}

void FiducialTableWidget::onDeleteMarkupPushButtonClicked(bool confirm /*=true*/)
{

	vtkMRMLMarkupsFiducialNode* MarkupsNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->GetNodeByID(activeMarkupsNodeID.c_str()));
	if (!MarkupsNode)
	{
		return;
	}

	QList<QTableWidgetItem*> selectedItems =this->selectedItems();

	if (selectedItems.isEmpty())
	{
		return;
	}

	QList<int> rows;
	for (int i = 0; i < selectedItems.size(); i += columnLabels.size())
	{
		int row = selectedItems.at(i)->row();
		rows << row;
	}
	std::sort(rows.begin(), rows.end());

	if (confirm)
	{
		ctkMessageBox deleteAllMsgBox;
		deleteAllMsgBox.setWindowTitle("Delete Markups in this list?");
		QString labelText = QString("Delete ")
			+ QString::number(rows.size())
			+ QString(" Markups from this list?");
		deleteAllMsgBox.setText(labelText);

		QPushButton* deleteButton =
			deleteAllMsgBox.addButton(tr("Delete"), QMessageBox::AcceptRole);
		deleteAllMsgBox.addButton(QMessageBox::Cancel);
		deleteAllMsgBox.setDefaultButton(deleteButton);
		deleteAllMsgBox.setIcon(QMessageBox::Question);
		deleteAllMsgBox.setDontShowAgainVisible(true);
		deleteAllMsgBox.setDontShowAgainSettingsKey("Markups/AlwaysDeleteMarkups");
		deleteAllMsgBox.exec();
		if (deleteAllMsgBox.clickedButton() != deleteButton)
		{
			return;
		}
	}

	for (int i = rows.size() - 1; i >= 0; --i)
	{
		int index = rows.at(i);
		MarkupsNode->RemoveNthControlPoint(index);
	}

	this->clearSelection();
}

void FiducialTableWidget::onJumpSlicesActionTriggered()
{
	vtkMRMLMarkupsFiducialNode* MarkupsNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->GetNodeByID(activeMarkupsNodeID.c_str()));
	if (!MarkupsNode)
	{
		return;
	}

	QList<QTableWidgetItem*> selectedItems = this->selectedItems();

	if (selectedItems.isEmpty())
	{
		return;
	}

	bool jumpCentered = false;
	jumpCentered = true;
	vtkSlicerMarkupsLogic* markupsLogic = vtkSlicerMarkupsLogic::SafeDownCast(
		qSlicerApplication::application()->moduleLogic("Markups"));

	if (markupsLogic)
	{
		markupsLogic->JumpSlicesToNthPointInMarkup(MarkupsNode->GetID(), selectedItems.at(0)->row(), jumpCentered);
	}
}

void FiducialTableWidget::updateRow(int m)
{
	this->blockSignals(true);

	vtkMRMLMarkupsFiducialNode* markupsNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->GetNodeByID(activeMarkupsNodeID.c_str()));
	if (!markupsNode)
	{
		return;
	}

	QTableWidgetItem* visibleItem = new QTableWidgetItem();
	visibleItem->setData(Qt::CheckStateRole, QVariant());
	visibleItem->setFlags(visibleItem->flags() & ~Qt::ItemIsUserCheckable);
	visibleItem->setFlags(visibleItem->flags() & ~Qt::ItemIsEditable);
	if (markupsNode->GetNthControlPointVisibility(m))
	{
		visibleItem->setData(Qt::UserRole, QVariant(true));
		visibleItem->setData(Qt::DecorationRole, QPixmap(":/Icons/Small/SlicerVisible.png"));
	}
	else
	{
		visibleItem->setData(Qt::UserRole, QVariant(false));
		visibleItem->setData(Qt::DecorationRole, QPixmap(":/Icons/Small/SlicerInvisible.png"));
	}
	int visibleIndex = columnLabels.indexOf("Visible");
	if (this->item(m, visibleIndex) == nullptr ||
		this->item(m, visibleIndex)->data(Qt::UserRole) != visibleItem->data(Qt::UserRole))
	{
		this->setItem(m, visibleIndex, visibleItem);
	}

	int nameIndex = columnLabels.indexOf("Name");
	QString markupLabel = QString(markupsNode->GetNthControlPointLabel(m).c_str());
	if (this->item(m, nameIndex) == nullptr ||
		this->item(m, nameIndex)->text() != markupLabel)
	{
		this->setItem(m, nameIndex, new QTableWidgetItem(markupLabel));
	}

	double point[3] = { 0.0, 0.0, 0.0 };
	if (true)
	{
		double worldPoint[4] = { 0.0, 0.0, 0.0, 1.0 };
		markupsNode->GetNthControlPointPositionWorld(m, worldPoint);
		for (int p = 0; p < 3; ++p)
		{
			point[p] = worldPoint[p];
		}
	}
	else
	{
		markupsNode->GetNthControlPointPosition(m, point);
	}
	int rColumnIndex = columnLabels.indexOf("R");
	for (int p = 0; p < 3; p++)
	{
		QString coordinate = QString::number(point[p], 'f', 3);
		if (this->item(m, rColumnIndex + p) == nullptr ||
			this->item(m, rColumnIndex + p)->text() != coordinate)
		{
			this->setItem(m, rColumnIndex + p, new QTableWidgetItem(coordinate));
		}
	}

	this->blockSignals(false);
}

void FiducialTableWidget::onUpdateFiducialTable()
{
	vtkMRMLMarkupsFiducialNode* markupsNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->GetNodeByID(activeMarkupsNodeID.c_str()));
	if (!markupsNode)
	{
		return;
	}
	int numberOfPoints = markupsNode->GetNumberOfControlPoints();
	if (this->rowCount() != numberOfPoints)
	{
		this->setRowCount(numberOfPoints);
	}
	for (int m = 0; m < numberOfPoints; m++)
	{
		this->updateRow(m);
	}
}

void FiducialTableWidget::onActiveMarkupsNodePointModifiedEvent(vtkObject* caller, void* callData)
{
	if (caller == nullptr)
	{
		return;
	}

	int* nPtr = reinterpret_cast<int*>(callData);
	int n = (nPtr ? *nPtr : -1);
	if (n >= 0)
	{
		this->updateRow(n);
	}
	else
	{
		this->onUpdateFiducialTable();
	}
}

void FiducialTableWidget::onActiveMarkupsNodePointAddedEvent()
{

	int newRow = this->rowCount();
	this->insertRow(newRow);

	this->updateRow(newRow);

	this->setCurrentCell(newRow, 0);

	this->onUpdateFiducialTable();
}

void FiducialTableWidget::onActiveMarkupsNodePointRemovedEvent(vtkObject* caller, void* callData)
{

	if (caller == nullptr)
	{
		return;
	}

	int* nPtr = reinterpret_cast<int*>(callData);
	int n = (nPtr ? *nPtr : -1);
	if (n >= 0)
	{
		this->removeRow(n);
	}
	else
	{
		this->onUpdateFiducialTable();
	}
}

void FiducialTableWidget::onActiveMarkupTableCellClicked(QTableWidgetItem* item)
{
	if (item == nullptr)
	{
		return;
	}
	int column = item->column();
	if (column == columnLabels.indexOf(QString("Visible")))
	{
		if (item->data(Qt::UserRole) == QVariant(false))
		{
			item->setData(Qt::UserRole, QVariant(true));
		}
		else
		{
			item->setData(Qt::UserRole, QVariant(false));
		}
	}
}

void FiducialTableWidget::onActiveMarkupTableCellChanged(int row, int column)
{
	vtkMRMLMarkupsFiducialNode* markupsNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(
		qSlicerApplication::application()->mrmlScene()->GetNodeByID(activeMarkupsNodeID.c_str()));
	if (!markupsNode)
	{
		return;
	}
	if (!markupsNode)
	{
		return;
	}

	int n = row;

	QTableWidgetItem* item = this->item(row, column);
	if (!item)
	{
		qDebug() << QString("Unable to find item in table at ") + QString::number(row) + QString(", ") + QString::number(column);
		return;
	}
	if (column == columnLabels.indexOf("Visible"))
	{
		bool flag = item->data(Qt::UserRole) == QVariant(true) ? true : false;
		if (flag)
		{
			item->setData(Qt::DecorationRole, QPixmap(":/Icons/Small/SlicerVisible.png"));
		}
		else
		{
			item->setData(Qt::DecorationRole, QPixmap(":/Icons/Small/SlicerInvisible.png"));
		}
		markupsNode->SetNthControlPointVisibility(n, flag);
	}
	else if (column == columnLabels.indexOf("Name"))
	{
		std::string name = std::string(item->text().toUtf8());
		markupsNode->SetNthControlPointLabel(n, name);
	}
	else if (column == columnLabels.indexOf("R") ||
		column == columnLabels.indexOf("A") ||
		column == columnLabels.indexOf("S"))
	{
		double newPoint[3] = { 0.0, 0.0, 0.0 };
		if (this->item(row, columnLabels.indexOf("R")) == nullptr ||
			this->item(row, columnLabels.indexOf("A")) == nullptr ||
			this->item(row, columnLabels.indexOf("S")) == nullptr)
		{
			return;
		}
		newPoint[0] = this->item(row, columnLabels.indexOf("R"))->text().toDouble();
		newPoint[1] = this->item(row, columnLabels.indexOf("A"))->text().toDouble();
		newPoint[2] = this->item(row, columnLabels.indexOf("S"))->text().toDouble();

		double point[3] = { 0.0, 0.0, 0.0 };
		if (true)
		{
			double worldPoint[4] = { 0.0, 0.0, 0.0, 1.0 };
			markupsNode->GetNthControlPointPositionWorld(n, worldPoint);
			for (int p = 0; p < 3; ++p)
			{
				point[p] = worldPoint[p];
			}
		}
		else
		{
			markupsNode->GetNthControlPointPosition(n, point);
		}

		double minChange = 0.001;
		if (fabs(newPoint[0] - point[0]) > minChange ||
			fabs(newPoint[1] - point[1]) > minChange ||
			fabs(newPoint[2] - point[2]) > minChange)
		{
			if (true)
			{
				markupsNode->SetNthControlPointPositionWorld(n, newPoint[0], newPoint[1], newPoint[2]);
			}
			else
			{
				markupsNode->SetNthControlPointPositionFromArray(n, newPoint);
			}
		}
		else
		{
		}
	}
	else
	{
		qDebug() << QString("Cell Changed: unknown column: ") + QString::number(column);
	}
}







PointMatcher::PointMatcher()
{
	on_target = false;
	min_value = 0;
	p = NULL;
}

PointMatcher::~PointMatcher()
{

}

void PointMatcher::Point_Registration(QList<vtkVector3d> DentureFiducialList_p,
	QList<vtkVector3d> DentureFiducialList_t,
	QList<vtkVector3d>& DentureFiducialList_1,
	QList<vtkVector3d>& DentureFiducialList_2)
{
	this->point_input(DentureFiducialList_p, DentureFiducialList_t, DentureFiducialList_p.size(), DentureFiducialList_t.size());
	this->onTooth_dis();
	this->onPatient_dis();
	this->Point_Screening(DentureFiducialList_1, DentureFiducialList_2);
}

void PointMatcher::point_input(QList<vtkVector3d> DentureFiducialList_p,
	QList<vtkVector3d> DentureFiducialList_t, int p_size, int t_size)
{
	for (int i = 0; i < p_size; i++)
	{
		patient_x.push_back(DentureFiducialList_p[i][0]);
		patient_y.push_back(DentureFiducialList_p[i][1]);
		patient_z.push_back(DentureFiducialList_p[i][2]);
	}
	for (int i = 0; i < t_size; i++)
	{
		tooth_x.push_back(DentureFiducialList_t[i][0]);
		tooth_y.push_back(DentureFiducialList_t[i][1]);
		tooth_z.push_back(DentureFiducialList_t[i][2]);
	}
}

void PointMatcher::onTooth_dis()
{
	teeth = tooth_x.size();
	for (int i = 0; i < teeth - 1; i++)
	{
		for (int j = 1; j < teeth - i; j++)
		{
			double dis = sqrt((double)pow((tooth_x[i] - tooth_x[i + j]), 2)
				+ (double)pow((tooth_y[i] - tooth_y[i + j]), 2)
				+ (double)pow((tooth_z[i] - tooth_z[i + j]), 2));
			tooth_dis2.push_back(dis);
		}
	}
	sort(tooth_dis2.begin(), tooth_dis2.end());
	for (int i = 0; i < teeth - 1; i++)
	{
		std::vector<double> group;
		for (int j = 0; j < i + 1; j++)
		{
			group.push_back(0);
		}

		for (int j = 1; j < teeth - i; j++)
		{
			double dis = sqrt((double)pow((tooth_x[i] - tooth_x[i + j]), 2)
				+ (double)pow((tooth_y[i] - tooth_y[i + j]), 2)
				+ (double)pow((tooth_z[i] - tooth_z[i + j]), 2));
			group.push_back(dis);
		}
		tooth_dis.push_back(group);
	}

	for (int i = 0; i < teeth; i++)
	{
		std::vector<double> group_src;
		for (int j = 0; j < teeth; j++)
		{
			if (i > j) {
				group_src.push_back(tooth_dis[j][i]);
			}
			else if (i < j) {
				group_src.push_back(tooth_dis[i][j]);
			}
		}
		ontooth_src.push_back(group_src);
	}
}

void PointMatcher::onPatient_dis()
{
	patient = patient_x.size();
	for (int i = 0; i < patient - 1; i++)
	{
		std::vector<double> group;
		for (int j = 0; j < i + 1; j++)
		{
			group.push_back(0);
		}

		for (int j = 1; j < patient - i; j++)
		{
			double dis = sqrt((double)pow((patient_x[i] - patient_x[i + j]), 2)
				+ (double)pow((patient_y[i] - patient_y[i + j]), 2)
				+ (double)pow((patient_z[i] - patient_z[i + j]), 2));
			group.push_back(dis);
		}
		patient_dis.push_back(group);
	}
}

void PointMatcher::onPatientgroup_dis()
{
	std::vector<double> onpatient;

	double error = 0;
	int min_error = 0;

	for (int i = 0; i < patient_point.size() - 1; i++)
	{
		for (int j = 1; j < patient_point.size() - i; j++)
		{
			onpatient.push_back(patient_dis[patient_point[i]][patient_point[i + j]]);
		}
	}
	sort(onpatient.begin(), onpatient.end());

	for (int i = 0; i < tooth_dis2.size(); i++)
	{
		double err = pow((onpatient[i] - tooth_dis2[i]), 2);
		error = err + error;
	}
	error_group.emplace_back(error);

	if (error <= min_value || min_value == 0)
	{
		min_value = error;
		on_target = true;
	}

}

void PointMatcher::select_point(int pos, int process, int totla, int key, std::vector<double>& a, std::vector<int> b, std::vector<int>& p)
{
	if (process > totla)
		return;
	if (pos == key)
	{
		for (int i = 0; i < key; i++)
		{
			this->patient_point.push_back(b[i]);
		}
		this->onPatientgroup_dis();
		if (on_target == true)
		{
			for (int i = 0; i < key; i++)
			{
				p[i] = patient_point[i];
			}
			on_target = false;
		}
		patient_point.clear();
		return;
	}
	else
	{
		select_point(pos, process + 1, totla, key, a, b, p);
		if (process >= this->patient)
		{
			select_point(pos + 1, process + 1, totla, key, a, b, p);
		}
		else
		{
			b[pos] = a[process];
			select_point(pos + 1, process + 1, totla, key, a, b, p);
		}
	}
}

void PointMatcher::Point_Screening(QList<vtkVector3d>& DentureFiducialList_1,
	QList<vtkVector3d>& DentureFiducialList_2)
{
	for (int i = 0; i < patient; i++)
	{
		std::vector<double> onpatient_src;
		for (int j = 0; j < patient; j++)
		{
			if (i > j) {
				onpatient_src.push_back(patient_dis[j][i]);
			}
			else if (i < j) {
				onpatient_src.push_back(patient_dis[i][j]);
			}
		}
		double sum = 0;
		double min_error = 0;
		double min_num = 0;
		sort(onpatient_src.begin(), onpatient_src.end());
		for (int k = 0; k < ontooth_src.size(); k++)
		{
			for (int m = 0; m < ontooth_src[k].size(); m++)
			{
				double num = 0;
				num = Find_Closest(0, onpatient_src.size() - 1, ontooth_src[k][m], onpatient_src);
				sum += abs(num - ontooth_src[k][m]);
			}
			if (sum < min_error || min_error == 0)
			{
				min_error = sum;
				min_num = k;
			}
			sum = 0;
		}
		std::vector<double> group;
		group.push_back(i);
		group.push_back(min_num);
		group.push_back(min_error);
		patient_closest_point.push_back(group);
	}

	int number = 0;
	bool err = false;
	for (int n = 0; n < teeth; n++)
	{
		for (int m = 0; m < patient_closest_point.size(); m++)
		{
			if (patient_closest_point[m][1] == n)
			{
				err = true;
			}
		}
		if (err == true)
		{
			number++;
		}
		else {
			wrong_teeth_num.push_back(n);
		}
		err = false;
	}
	if (number != teeth) {
		this->Processing_special_point(DentureFiducialList_1, DentureFiducialList_2);
		return;
	}

	for (int i = 0; i < teeth; i++)
	{
		std::vector<double> group;
		optional_number.push_back(group);
	}

	for (int i = 0; i < patient; i++)
	{
		for (int j = 0; j < teeth; j++)
		{
			if (patient_closest_point[i][1] == j) {
				optional_number[j].push_back(i);
			}
		}
	}
	std::vector<double> ans;
	for (int i = 0; i < teeth; i++)
	{
		tar.push_back(0);
		ans.push_back(0);
	}

	select_point_2(0, 0, optional_number.size(), optional_number, tar, ans);
	for (int i = 0; i < ans.size(); i++)
	{
		vtkVector4d p0 = { ans[i], patient_x[ans[i]], patient_y[ans[i]], patient_z[ans[i]] };
		list_ans.append(p0);
	}
	this->Point_Sorting(DentureFiducialList_1, DentureFiducialList_2);
}

void PointMatcher::select_point_2(int pos, int process, int total,
	std::vector<std::vector<double>> aa, std::vector<double>& tar, std::vector<double>& ans)
{
	if (process > total)
		return;
	if (process != total && pos >= aa[process].size())
		return;
	if (process == total)
	{
		for (int i = 0; i < process; i++)
		{
			this->patient_point.push_back(tar[i]);
		}
		sort(patient_point.begin(), patient_point.end());
		this->onPatientgroup_dis();
		if (on_target == true)
		{
			for (int i = 0; i < process; i++)
			{
				ans[i] = patient_point[i];
			}
			on_target = false;
		}
		patient_point.clear();
		return;
	}
	else
	{
		select_point_2(pos + 1, process, total, aa, tar, ans);
		if (process >= total)
		{
			select_point_2(pos, process + 1, total, aa, tar, ans);
		}
		else
		{
			tar[process] = aa[process][pos];
			select_point_2(0, process + 1, total, aa, tar, ans);
		}
	}
}

double PointMatcher::Find_Closest(int mn, int mx, double a, std::vector<double> p)
{
	if (a >= p[mx])
		return p[mx];
	else if (a <= p[mn])
		return p[mn];

	if (mx - mn == 1)
	{
		int t1 = a - p[mn], t2 = p[mx] - a;
		if (t1 > t2)
			return p[mx];
		else        //(t1<=t2)
			return p[mn];
	}
	else
	{
		int mid = (mn + mx) / 2;
		if (a > p[mid])
			return Find_Closest(mid, mx, a, p);
		else if (a < p[mid])
			return Find_Closest(mn, mid, a, p);
		else
			return p[mid];
	}
}

void PointMatcher::Point_Sorting(QList<vtkVector3d>& DentureFiducialList_1, QList<vtkVector3d>& DentureFiducialList_2)
{
	for (int i = 0; i < list_ans.size(); i++)
	{
		for (int j = 0; j < list_ans.size(); j++)
		{
			if (patient_closest_point[list_ans[j][0]][1] == i)
			{
				vtkVector3d ans0 = { list_ans[j][1], list_ans[j][2], list_ans[j][3] };
				DentureFiducialList_1.append(ans0);
			}
		}
	}

	for (double i = 0, k = 0; i < tooth_x.size() + wrong_teeth_num.size(); i++, k++)
	{
		for (int j = 0; j < wrong_teeth_num.size(); j++)
		{
			if (wrong_teeth_num[j] == i)
			{
				i++;
			}
		}
		vtkVector3d ans1 = { tooth_x[k], tooth_y[k], tooth_z[k] };
		DentureFiducialList_2.append(ans1);
	}
}

void PointMatcher::Processing_special_point(QList<vtkVector3d>& DentureFiducialList_1, QList<vtkVector3d>& DentureFiducialList_2)
{
	for (int i = 0; i < wrong_teeth_num.size(); i++)
	{
		std::vector<double>::iterator it_x = std::find(tooth_x.begin(), tooth_x.end(), tooth_x[wrong_teeth_num[i]]);
		tooth_x.erase(it_x);
		std::vector<double>::iterator it_y = std::find(tooth_y.begin(), tooth_y.end(), tooth_y[wrong_teeth_num[i]]);
		tooth_y.erase(it_y);
		std::vector<double>::iterator it_z = std::find(tooth_z.begin(), tooth_z.end(), tooth_z[wrong_teeth_num[i]]);
		tooth_z.erase(it_z);
	}
	tooth_dis2.clear();
	tooth_dis.clear();
	ontooth_src.clear();
	patient_closest_point.clear();
	this->onTooth_dis();
	this->Point_Screening(DentureFiducialList_1, DentureFiducialList_2);
}