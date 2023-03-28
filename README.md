# SlicerDentureRegistration
This is an extension module to realize the registration of patients' CT images and the corresponding three-dimensional denture models. It is available in 3D Slicer.
![](https://raw.githubusercontent.com/gongmingjun/SlicerDentureRegistration/master/DentureRegistration.jpg)

# Usage
1. Adjust the lower thresholds of patients' CT and dentures' CT for the segmentation of radiopaque markers.
2. Load patients' CT images stored in the form of DICOM by clicking the button "Import Patient DICOM". The patients need to wear customized dentures with radiopaque markers for CT scanning. The picture shows an example of a reasonable patients' CT images.
![](https://raw.githubusercontent.com/gongmingjun/SlicerDentureRegistration/master/patientCT.png)

3. Load dentures' CT images stored in the form of DICOM by clicking the button "Import Denture DICOM". The CT scanning should be taken purely for the customized dentures.The picture shows an example of a reasonable dentures' CT images.
![](https://raw.githubusercontent.com/gongmingjun/SlicerDentureRegistration/master/dentureCT.png)

4. Once the dentures' CT images are loaded, the registration will be automatically performed. The coordinates of the extracted fiducial markers (radiopaque markers) and the rms of the registration will be shown on the panel.
