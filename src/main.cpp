#include <cstdlib>
#include <iostream>
#include <string>

#include <vtkActor.h>
#include <vtkCleanPolyData.h>
#include <vtkCurvatures.h>
#include <vtkOBJReader.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkParametricTorus.h>
#include <vtkParametricFunctionSource.h>

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
    return EXIT_FAILURE;
  }

  std::string filename = argv[1];

  auto reader = vtkSmartPointer<vtkOBJReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  auto polyData = reader->GetOutput();

  auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetInputData(polyData);
  cleaner->Update();

  auto curvatureFilter = vtkSmartPointer<vtkCurvatures>::New();
  curvatureFilter->SetInputData(cleaner->GetOutput());
  // Options: ToMean, ToGaussian, ToMaximum, ToMinimum
  curvatureFilter->SetCurvatureTypeToGaussian();
  curvatureFilter->Update();

  auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(curvatureFilter->GetOutput());

  auto actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  auto renderer = vtkSmartPointer<vtkRenderer>::New();
  auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  auto interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

  renderWindow->AddRenderer(renderer);
  interactor->SetRenderWindow(renderWindow);
  renderer->AddActor(actor);
  renderer->SetBackground(0.1, 0.1, 0.1);

  renderWindow->SetSize(800, 600);
  renderWindow->Render();
  interactor->Start();

  return EXIT_SUCCESS;
}
