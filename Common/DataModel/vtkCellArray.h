/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCellArray.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkCellArray - object to represent cell connectivity
// .SECTION Description
// vtkCellArray is a supporting object that explicitly represents cell
// connectivity. The cell array structure is a raw integer list
// of the form: (n,id1,id2,...,idn, n,id1,id2,...,idn, ...)
// where n is the number of points in the cell, and id is a zero-offset index
// into an associated point list.
//
// Advantages of this data structure are its compactness, simplicity, and
// easy interface to external data.  However, it is totally inadequate for
// random access.  This functionality (when necessary) is accomplished by
// using the vtkCellTypes and vtkCellLinks objects to extend the definition of
// the data structure.
//
// .SECTION See Also
// vtkCellTypes vtkCellLinks

#ifndef __vtkCellArray_h
#define __vtkCellArray_h

#include "vtkCommonDataModelModule.h" // For export macro
#include "vtkObject.h"

#include "vtkIdTypeArray.h" // Needed for inline methods
#include "vtkCell.h" // Needed for inline methods

class VTKCOMMONDATAMODEL_EXPORT vtkCellArray : public vtkObject
{
public:
  vtkTypeMacro(vtkCellArray,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Instantiate cell array (connectivity list).
  static vtkCellArray *New();

  // Description:
  // Allocate enough memory to fit the source object
  void Reserve(const vtkCellArray* src)
  {
    this->Points->Allocate(src->Points->GetNumberOfTuples());
    this->CellIndex->Resize(src->CellIndex->GetNumberOfTuples());
  }


  // Description:
  // Allocate memory specifying the number of points
  int Allocate(const vtkIdType sz, int extSize = 1000)
  {
    // return this->Points->Allocate(sz) &&
    //   this->CellIndex->Resize(sz/3); // assume triangles
  }

  // Description:
  // Allocate memory specifying the number of cells and
  // the number of points per cell
  int Reserve(vtkIdType numberOfCells, vtkIdType cellsPerPoint)
  {

    return this->Points->Allocate(numberOfCells * cellsPerPoint) &&
      this->CellIndex->Resize(numberOfCells+1);
  }

  // Description:
  // Allocate memory specifying the number of cells and
  // the number of points
  int ReserveAll(vtkIdType numberOfCells, vtkIdType numberOfPoints)
  {
    return this->Points->Allocate(numberOfPoints) &&
      this->CellIndex->Resize(numberOfCells+1);
  }


  // Description:
  // Free any memory and reset to an empty state.
  void Initialize();

  // Description:
  // Get the number of cells in the array.
  vtkIdType GetNumberOfCells() const
  {
    return this->CellIndex->GetNumberOfTuples() - 1;
  }

  // Description:
  // Get total number of points in the array
  vtkIdType GetNumberOfPoints() const
  {
    assert(this->Points->GetNumberOfTuples() ==
           this->CellIndex->GetValue(this->CellIndex->GetNumberOfTuples() - 1));
    return this->Points->GetNumberOfTuples();
  }

  // Description:
  // A cell traversal methods that is more efficient than vtkDataSet traversal
  // methods.  InitTraversal() initializes the traversal of the list of cells.
  void InitTraversal() {this->TraversalId=0;};

  // Description:
  // A cell traversal methods that is more efficient than vtkDataSet traversal
  // methods.  GetNextCell() gets the next cell in the list. If end of list
  // is encountered, 0 is returned. A value of 1 is returned whenever
  // npts and pts have been updated without error.
  int GetNextCell(vtkIdType& npts, vtkIdType* &pts);

  // Description:
  // A cell traversal methods that is more efficient than vtkDataSet traversal
  // methods.  GetNextCell() gets the next cell in the list. If end of list
  // is encountered, 0 is returned.
  int GetNextCell(vtkIdList *pts);

  // Description:
  // Get the size of the allocated connectivity data structure.
  vtkIdType GetSize()
  {
    return this->Points->GetSize() + this->CellIndex->GetSize();
  }

  // Description:
  // Retrieves a cell given a cellId.
  // The returned pointer can be used to access all the points of a cell
  // as all cell's points are contiguous.
  void GetCellFromId(vtkIdType cellId, vtkIdType &npts, vtkIdType* &pts) const;

  // Description:
  // Internal method used to retrieve a cell point count given a cellId
  inline vtkIdType GetCellPointCountFromId(vtkIdType cellId) const;

  // Description:
  // Internal method used to retrieve a cell given an offset into
  // the internal array.
  void GetCellFromId(vtkIdType cellId, vtkIdList* pts);

  // Description:
  // Insert a cell object. Return the cell id of the cell.
  vtkIdType InsertNextCell(vtkCell *cell);

  // Description:
  // Create a cell by specifying the number of points and an array of point
  // id's.  Return the cell id of the cell.
  template<typename T> vtkIdType InsertNextCell(vtkIdType npts, const T* pts);

  // Description:
  // Create a cell by specifying a list of point ids. Return the cell id of
  // the cell.
  vtkIdType InsertNextCell(vtkIdList *pts);

  // Description:
  // Create cells by specifying count, and then adding points one at a time
  // using method InsertCellPoint(). If you don't know the count initially,
  // use the method UpdateCellCount() to complete the cell. Return the cell
  // id of the cell.
  vtkIdType InsertNextCell(int npts);

  // Description:
  // Used in conjunction with InsertNextCell(int npts) to add another point
  // to the list of cells.
  void InsertCellPoint(vtkIdType id);

  // Description:
  // Used in conjunction with InsertNextCell(int npts) and InsertCellPoint() to
  // update the number of points defining the cell.
  void UpdateCellCount(int npts);

  // Description:
  // Get/Set the current traversal location.
  vtkIdType GetTraversalId()
    {return this->TraversalId;}
  void SetTraversalId(vtkIdType cellId)
    {this->TraversalId = cellId;}

  // Description:
  // Special method inverts ordering of current cell. Must be called
  // carefully or the cell topology may be corrupted.
  void ReverseCellFromId(vtkIdType cellId);

  // Description:
  // Replace the point ids of the cell with a different list of point ids.
  void ReplaceCellFromId(vtkIdType cellId, int npts, const vtkIdType *pts);

  // Description:
  // Returns the size of the largest cell. The size is the number of points
  // defining the cell.
  int GetMaxCellSize();

  // Description:
  // Perform a deep copy (no reference counting) of the given cell array.
  void DeepCopy(vtkCellArray *ca);
  // Description:
  // Append all cells in 'other' with their points offset with 'otherStartPoint'
  void Append(const vtkCellArray *other, vtkIdType otherStartPoint);

  // Description:
  // Copies cells from/to a format where each cell is represented by
  // the number of points followed by the point ids. For example,
  // for two cells which 3 and 4 points respectively, 'cells' is
  // (3, p0_0, p0_1, p0_2, 4, p1_0, p1_1, p1_2, p1_3) where pi_j are point ids.
  template<typename T>
  void AppendFromCountPointsFormat(vtkIdType ncells, T *cells)
  {
    for (int i = 0; i < ncells; ++i)
      {
      T npoints = *cells++;
      this->InsertNextCell(npoints);
      for (int j = 0; j < npoints; ++j)
        {
        this->InsertCellPoint(*cells++);
        }
      }
  }

  template<typename T>
  void CopyFromCountPointsFormat(vtkIdType ncells, T *cells)
  {
    this->Initialize();
    AppendFromCountPointsFormat(ncells, cells);
  }

  template<typename T>
  void CopyToCountPointsFormat(T *cells) const
  {
    vtkIdType* cellIndex = this->CellIndex->GetPointer(0);
    for (vtkIdType i = 0; i < this->GetNumberOfCells(); ++i, ++cellIndex)
      {
      vtkIdType npts = *(cellIndex + 1) - *cellIndex;
      vtkIdType* pts = this->Points->GetPointer(*cellIndex);
      *(cells++) = npts;
      for (vtkIdType j = 0; j < npts; ++j)
        {
        *(cells++) = *(pts++);
        }
      }
  }

  // Description:
  // Returns 1 if the cell array has cells with the same number of points,
  // and in that case it sets 'numberOfPoins'. Otherwise it returns 0 and
  // leaves numberOfpoints unchanged.
  int IsHomogeneous(int* numberOfPoints) const;

  // Description:
  // Reuse list. Reset to initial condition.
  void Reset();

  // Description:
  // Reclaim any extra memory.
  void Squeeze()
  {
    this->Points->Squeeze();
    this->CellIndex->Squeeze();
  }

  // Description:
  // Return the memory in kilobytes consumed by this cell array. Used to
  // support streaming and reading/writing data. The value returned is
  // guaranteed to be greater than or equal to the memory required to
  // actually represent the data represented by this object. The
  // information returned is valid only after the pipeline has
  // been updated.
  unsigned long GetActualMemorySize();

protected:
  vtkCellArray();
  ~vtkCellArray();

  // keeps track of the current point location for InsertCellPoint
  vtkIdType InsertPointLocation;
  vtkIdType TraversalId;   //keep track of traversal position
  // Description:
  // Points of cells. For instance, for two cells c0 with 3 points and
  // c1 with 4 points this array will contain the following point ids:
  // (p0_1, p0_2, p0_3, p1_1, p1_2, p1_3, p1_4)
  vtkIdTypeArray *Points;
  // Description:
  // Array that stores, for each cell, the index of the first point of the cell.
  // This array has an extra entry that points to the end of the Points array
  // which is used as a sentinel to avoid an extra test.
  // For instance for two cells c0 with 3 points and c1 with 4 points this array
  // will contain (0, 3, 7)
  vtkIdTypeArray* CellIndex;

private:
  vtkCellArray(const vtkCellArray&);  // Not implemented.
  void operator=(const vtkCellArray&);  // Not implemented.
};


//----------------------------------------------------------------------------
template<typename T>
inline vtkIdType vtkCellArray::InsertNextCell(vtkIdType npts, const T* pts)
{
  // add the location of the first point and the sentinel
  vtkIdType newCellIndex = this->Points->GetNumberOfTuples();
  vtkIdType *cellIndex = this->CellIndex->WritePointer(
    this->CellIndex->GetNumberOfTuples() - 1, 2);
  *cellIndex = newCellIndex;
  *(cellIndex + 1) = newCellIndex + npts;
  // add the points
  vtkIdType *point = this->Points->WritePointer(newCellIndex, npts);
  for (int i = 0; i < npts; i++)
    {
    *point++ = *pts++;
    }
  return this->GetNumberOfCells() - 1;
}

//----------------------------------------------------------------------------
inline vtkIdType vtkCellArray::InsertNextCell(int npts)
{
  // add the location of the first point and the sentinel
  vtkIdType newCellIndex = this->Points->GetNumberOfTuples();
  vtkIdType *cellIndex = this->CellIndex->WritePointer(
    this->CellIndex->GetNumberOfTuples() - 1, 2);
  *cellIndex = newCellIndex;
  *(cellIndex + 1) = newCellIndex + npts;
  this->InsertPointLocation = newCellIndex;
  return this->GetNumberOfCells() - 1;
}

//----------------------------------------------------------------------------
inline void vtkCellArray::InsertCellPoint(vtkIdType id)
{
  this->Points->InsertValue(this->InsertPointLocation++, id);
}

//----------------------------------------------------------------------------
inline void vtkCellArray::UpdateCellCount(int npts)
{
  vtkIdType *cellIndex = this->CellIndex->GetPointer(
    this->CellIndex->GetNumberOfTuples() - 2);
  *(cellIndex + 1) = *cellIndex + npts;
}

//----------------------------------------------------------------------------
inline vtkIdType vtkCellArray::InsertNextCell(vtkIdList *pts)
{
  return this->InsertNextCell(pts->GetNumberOfIds(), pts->GetPointer(0));
}

//----------------------------------------------------------------------------
inline vtkIdType vtkCellArray::InsertNextCell(vtkCell *cell)
{
  return this->InsertNextCell(cell->GetNumberOfPoints(),
                              cell->PointIds->GetPointer(0));
}

//----------------------------------------------------------------------------
inline void vtkCellArray::Reset()
{
  this->Points->Reset();
  this->CellIndex->Reset();
  this->CellIndex->InsertValue(0, 0);
  this->InsertPointLocation = 0;
  this->TraversalId = 0;
}

//----------------------------------------------------------------------------
inline int vtkCellArray::GetNextCell(vtkIdType& npts, vtkIdType* &pts)
{
  if (this->TraversalId < this->GetNumberOfCells())
    {
    this->GetCellFromId(this->TraversalId++, npts, pts);
    return 1;
    }
  npts = 0;
  pts = 0;
  return 0;
}


//----------------------------------------------------------------------------
inline void vtkCellArray::GetCellFromId(vtkIdType cellId, vtkIdType &npts,
                                        vtkIdType* &pts) const
{
  vtkIdType* cellIndex = this->CellIndex->GetPointer(cellId);
  npts = *(cellIndex+1) - *cellIndex;
  pts = this->Points->GetPointer(*cellIndex);
}

//----------------------------------------------------------------------------
inline void vtkCellArray::ReverseCellFromId(vtkIdType cellId)
{
  vtkIdType *cellIndex = this->CellIndex->GetPointer(cellId);
  vtkIdType npts = *(cellIndex + 1) - *cellIndex;
  vtkIdType *pts=this->Points->GetPointer(*cellIndex);
  for (int i=0; i < (npts/2); ++i)
    {
    vtkIdType tmp = pts[i];
    pts[i] = pts[npts-i-1];
    pts[npts-i-1] = tmp;
    }
}

//----------------------------------------------------------------------------
inline void vtkCellArray::ReplaceCellFromId(vtkIdType cellId, int npts,
                                            const vtkIdType *pts)
{
  vtkIdType *cellIndex = this->CellIndex->GetPointer(cellId);
  assert (*(cellIndex + 1) - *cellIndex == npts);
  vtkIdType *oldPts=this->Points->GetPointer(*cellIndex);
  for (int i=0; i < npts; ++i)
    {
    oldPts[i] = pts[i];
    }
}

//----------------------------------------------------------------------------
template<>
inline void vtkCellArray::CopyToCountPointsFormat<vtkIdTypeArray>(
  vtkIdTypeArray *cells) const
{
  CopyToCountPointsFormat(
    cells->WritePointer(0, this->GetNumberOfCells() + this->GetNumberOfPoints()));
}

//----------------------------------------------------------------------------
inline vtkIdType vtkCellArray::GetCellPointCountFromId(vtkIdType cellId) const
{
  vtkIdType* cellIndex = this->CellIndex->GetPointer(cellId);
  return *(cellIndex+1) - *cellIndex;
}

//----------------------------------------------------------------------------
inline void vtkCellArray::Append(
  const vtkCellArray *other, vtkIdType startPoint)
{
  this->Points->Resize(this->GetNumberOfPoints() +
                       other->GetNumberOfPoints());
  this->CellIndex->Resize(this->GetNumberOfCells() +
                          other->GetNumberOfCells() + 1);
  for (vtkIdType i = 0; i < other->GetNumberOfCells(); ++i)
    {
    vtkIdType nPoints = 0;
    vtkIdType* points = NULL;
    other->GetCellFromId(i, nPoints, points);
    this->InsertNextCell(nPoints);
    for (vtkIdType j = 0; j < nPoints; ++j)
      {
      this->InsertCellPoint(startPoint + points[j]);
      }
    }
}


#endif
