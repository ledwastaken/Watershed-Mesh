#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkOBJReader.h>

#include <iostream>
#include <string>
#include <cstdlib>

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

  auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(reader->GetOutputPort());

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
