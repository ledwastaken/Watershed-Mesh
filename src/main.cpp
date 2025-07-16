#include <algorithm>
#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <vtkActor.h>
#include <vtkCleanPolyData.h>
#include <vtkCurvatures.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkLookupTable.h>
#include <vtkNamedColors.h>
#include <vtkOBJReader.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>

// Special labels for the watershed algorithm
#define UNLABELED -2
#define IN_QUEUE -1
#define WATERSHED 0

// Random seed for reproducible results
#define RANDOM_SEED 123

// Priority queue element for watershed flooding
struct QueueElement
{
  double priority;
  vtkIdType pointId;

  bool operator>(const QueueElement& other) const
  {
    return priority > other.priority;
  }
};

// Find all neighboring points of a given point in the mesh
void GetPointNeighbors(vtkPolyData* polyData, vtkIdType pointId,
                       vtkIdList* neighborIds)
{
  vtkNew<vtkIdList> cellIds;
  polyData->GetPointCells(pointId, cellIds);
  neighborIds->Reset();

  std::set<vtkIdType> neighborSet;

  for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i)
  {
    vtkIdType cellId = cellIds->GetId(i);
    vtkNew<vtkIdList> pointIdsInCell;
    polyData->GetCellPoints(cellId, pointIdsInCell);

    for (vtkIdType j = 0; j < pointIdsInCell->GetNumberOfIds(); ++j)
    {
      vtkIdType currentPointId = pointIdsInCell->GetId(j);
      if (currentPointId != pointId)
      {
        neighborSet.insert(currentPointId);
      }
    }
  }

  for (const auto& id : neighborSet)
  {
    neighborIds->InsertNextId(id);
  }
}

// Main watershed segmentation class
class WatershedSegmentation
{
public:
  WatershedSegmentation(vtkPolyData* polyData, const char* scalarArrayName)
    : m_PolyData(polyData)
    , m_Scalars(polyData->GetPointData()->GetArray(scalarArrayName))
  {
    if (!m_Scalars)
    {
      throw std::runtime_error("Scalar array not found in polydata.");
    }

    m_NumPoints = m_PolyData->GetNumberOfPoints();
    m_Labels = vtkSmartPointer<vtkIdTypeArray>::New();
    m_Labels->SetName("WatershedLabels");
    m_Labels->SetNumberOfComponents(1);
    m_Labels->SetNumberOfTuples(m_NumPoints);
    m_PolyData->GetPointData()->AddArray(m_Labels);
  }

  int Run()
  {
    // Initialize all points as unlabeled
    for (vtkIdType i = 0; i < m_NumPoints; ++i)
    {
      m_Labels->SetValue(i, UNLABELED);
    }

    // Find and label seed points
    int numLabels = FindAndLabelSeeds();
    std::cout << "Found " << numLabels - 1 << " seed regions." << std::endl;

    // Main flooding loop: process points in order of increasing curvature
    while (!m_PriorityQueue.empty())
    {
      QueueElement current = m_PriorityQueue.top();
      m_PriorityQueue.pop();

      vtkIdType currentPointId = current.pointId;

      if (m_Labels->GetValue(currentPointId) != IN_QUEUE)
      {
        continue;
      }

      // Check which labeled regions this point neighbors
      vtkNew<vtkIdList> neighborIds;
      GetPointNeighbors(m_PolyData, currentPointId, neighborIds);

      std::set<vtkIdType> neighborLabels;
      for (vtkIdType i = 0; i < neighborIds->GetNumberOfIds(); ++i)
      {
        vtkIdType neighborId = neighborIds->GetId(i);
        vtkIdType label = m_Labels->GetValue(neighborId);
        if (label > 0)
        {
          neighborLabels.insert(label);
        }
      }

      // Assign label based on neighboring regions
      if (neighborLabels.size() > 1)
      {
        // This point is on the boundary between multiple regions -> watershed
        m_Labels->SetValue(currentPointId, WATERSHED);
      }
      else if (neighborLabels.size() == 1)
      {
        // This point belongs to a single region
        vtkIdType singleLabel = *neighborLabels.begin();
        m_Labels->SetValue(currentPointId, singleLabel);

        // Add unlabeled neighbors to queue with slight priority randomization
        for (vtkIdType i = 0; i < neighborIds->GetNumberOfIds(); ++i)
        {
          vtkIdType neighborId = neighborIds->GetId(i);
          if (m_Labels->GetValue(neighborId) == UNLABELED)
          {
            m_Labels->SetValue(neighborId, IN_QUEUE);
            double priority = m_Scalars->GetTuple1(neighborId);
            m_PriorityQueue.push({ priority, neighborId });
          }
        }
      }
      else
      {
        // No labeled neighbors yet, mark as watershed for now
        m_Labels->SetValue(currentPointId, WATERSHED);
      }
    }

    // Mark any remaining unlabeled points as watershed
    for (vtkIdType i = 0; i < m_NumPoints; ++i)
    {
      if (m_Labels->GetValue(i) < 0)
      {
        m_Labels->SetValue(i, WATERSHED);
      }
    }

    return numLabels;
  }

private:
  // Select seed points using low-curvature regions
  int FindAndLabelSeeds()
  {
    // Use points with low curvature
    double range[2];
    m_Scalars->GetRange(range);
    double threshold = range[0] + (range[1] - range[0]) * 0.2;

    std::vector<vtkIdType> candidates;
    for (vtkIdType i = 0; i < m_NumPoints; ++i)
    {
      if (m_Scalars->GetTuple1(i) <= threshold)
      {
        candidates.push_back(i);
      }
    }

    std::cout << "Found " << candidates.size() << " candidate seed points"
              << std::endl;

    // Randomly shuffle and select seeds for better segmentation
    std::mt19937 rng(RANDOM_SEED);
    std::shuffle(candidates.begin(), candidates.end(), rng);

    int targetNumSeeds = std::min(static_cast<int>(candidates.size()),
                                 std::max(5, static_cast<int>(m_NumPoints / 30)));

    std::vector<vtkIdType> selectedSeeds(candidates.begin(),
                                         candidates.begin() + targetNumSeeds);

    // Improve spatial distribution of seeds
    std::vector<vtkIdType> finalSeeds;
    for (const auto& seed : selectedSeeds)
    {
      bool tooClose = false;
      for (const auto& existingSeed : finalSeeds)
      {
        // Minimum distance between seeds
        int minDistance = m_NumPoints / 15;
        if (std::abs(static_cast<int>(seed - existingSeed)) < minDistance)
        {
          tooClose = true;
          break;
        }
      }
      
      if (!tooClose)
      {
        finalSeeds.push_back(seed);
      }
      
      // Limit number of seeds to a reasonable maximum
      if (finalSeeds.size() >= 20)
      {
        break;
      }
    }

    std::cout << "Selected " << finalSeeds.size() << " well-distributed seeds"
              << std::endl;

    // Assign unique labels to seeds
    vtkIdType currentLabel = 1;
    for (const auto& seedId : finalSeeds)
    {
      m_Labels->SetValue(seedId, currentLabel);
      currentLabel++;
    }

    // Initialize priority queue with neighbors of seeds
    for (const auto& seedId : finalSeeds)
    {
      vtkNew<vtkIdList> neighborIds;
      GetPointNeighbors(m_PolyData, seedId, neighborIds);
      for (vtkIdType j = 0; j < neighborIds->GetNumberOfIds(); ++j)
      {
        vtkIdType neighborId = neighborIds->GetId(j);
        if (m_Labels->GetValue(neighborId) == UNLABELED)
        {
          m_Labels->SetValue(neighborId, IN_QUEUE);
          double priority = m_Scalars->GetTuple1(neighborId);
          m_PriorityQueue.push({ priority, neighborId });
        }
      }
    }

    return currentLabel;
  }

  vtkPolyData* m_PolyData;
  vtkDataArray* m_Scalars;
  vtkIdType m_NumPoints;
  vtkSmartPointer<vtkIdTypeArray> m_Labels;
  std::priority_queue<QueueElement, std::vector<QueueElement>,
                      std::greater<QueueElement>>
      m_PriorityQueue;
};

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
    return EXIT_FAILURE;
  }

  std::string filename = argv[1];

  // Load and clean mesh data
  auto reader = vtkSmartPointer<vtkOBJReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetInputData(reader->GetOutput());
  cleaner->Update();

  // Compute Gaussian curvature
  auto curvatureFilter = vtkSmartPointer<vtkCurvatures>::New();
  curvatureFilter->SetInputConnection(cleaner->GetOutputPort());
  curvatureFilter->SetCurvatureTypeToGaussian();
  curvatureFilter->Update();

  vtkPolyData* polyData = curvatureFilter->GetOutput();

  // Print curvature statistics
  vtkDataArray* curvatureArray =
      polyData->GetPointData()->GetArray("Gauss_Curvature");
  if (curvatureArray)
  {
    double range[2];
    curvatureArray->GetRange(range);
    std::cout << "Gaussian curvature range: [" << range[0] << ", " << range[1]
              << "]" << std::endl;
  }

  // Run watershed segmentation
  std::cout << "Starting watershed segmentation..." << std::endl;
  int numLabels = 0;
  try
  {
    WatershedSegmentation watershed(polyData, "Gauss_Curvature");
    numLabels = watershed.Run();
  }
  catch (const std::exception& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Watershed segmentation finished. Found " << numLabels
            << " labels (including watershed)." << std::endl;

  // Print segmentation statistics
  vtkIdTypeArray* labels = vtkIdTypeArray::SafeDownCast(
      polyData->GetPointData()->GetArray("WatershedLabels"));
  std::vector<int> labelCounts(numLabels + 1, 0);
  for (vtkIdType i = 0; i < polyData->GetNumberOfPoints(); ++i)
  {
    int label = labels->GetValue(i);
    if (label >= 0 && label <= numLabels)
    {
      labelCounts[label]++;
    }
  }

  std::cout << "Label statistics:" << std::endl;
  std::cout << "Watershed lines (label 0): " << labelCounts[0] << " points"
            << std::endl;
  for (int i = 1; i <= numLabels; ++i)
  {
    if (labelCounts[i] > 0)
    {
      std::cout << "Region " << i << ": " << labelCounts[i] << " points"
                << std::endl;
    }
  }

  // Create RGB color array for visualization
  auto colorArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colorArray->SetName("Colors");
  colorArray->SetNumberOfComponents(3);
  colorArray->SetNumberOfTuples(polyData->GetNumberOfPoints());

  auto colors = vtkSmartPointer<vtkNamedColors>::New();
  
  // Generate random colors for each region (deterministic with fixed seed)
  std::mt19937 colorRng(RANDOM_SEED); // Use same seed as other RNG for consistency
  std::uniform_int_distribution<int> colorDist(50, 255); // Avoid very dark colors
  
  // Always keep watershed lines black
  std::vector<std::array<unsigned char, 3>> colorMap;
  colorMap.push_back({0, 0, 0}); // Label 0 (watershed) -> Black
  
  // Generate random colors for each region found
  for (int i = 1; i < numLabels; ++i)
  {
    std::array<unsigned char, 3> color = {
      static_cast<unsigned char>(colorDist(colorRng)),
      static_cast<unsigned char>(colorDist(colorRng)),
      static_cast<unsigned char>(colorDist(colorRng))
    };
    colorMap.push_back(color);
  }
  
  std::cout << "Generated " << colorMap.size() - 1 << " random colors for regions" << std::endl;

  // Assign colors to each point according to its label
  for (vtkIdType i = 0; i < polyData->GetNumberOfPoints(); ++i)
  {
    int label = labels->GetValue(i);
    if (label >= 0 && label < static_cast<int>(colorMap.size()))
    {
      colorArray->SetTuple3(i, colorMap[label][0], colorMap[label][1], colorMap[label][2]);
    }
    else
    {
      colorArray->SetTuple3(i, 128, 128, 128); // Gray for unrecognized labels
    }
  }
  
  polyData->GetPointData()->SetScalars(colorArray);

  // Setup VTK visualization pipeline
  auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(polyData);
  mapper->SetScalarModeToUsePointData();
  mapper->SetColorModeToDirectScalars();
  mapper->SetInterpolateScalarsBeforeMapping(false);

  auto actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetInterpolationToFlat();
  actor->GetProperty()->LightingOff();

  auto renderer = vtkSmartPointer<vtkRenderer>::New();
  auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  auto interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("3D Watershed Segmentation");
  interactor->SetRenderWindow(renderWindow);
  renderer->AddActor(actor);
  renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());

  renderWindow->SetSize(800, 600);
  int screenWidth = 1920;
  int screenHeight = 1080;
  int middle_x = (screenWidth - renderWindow->GetSize()[0]) / 2;
  int middle_y = (screenHeight - renderWindow->GetSize()[1]) / 2;
  renderWindow->SetPosition(middle_x, middle_y);
  renderWindow->Render();
  interactor->Start();

  return EXIT_SUCCESS;
}